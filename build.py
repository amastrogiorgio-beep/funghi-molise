#!/usr/bin/env python3
"""Mappa settimanale della probabilita' (euristica) di trovare funghi in Molise.

Uso:  python3 build.py [--out DIR] [--date YYYY-MM-DD]

Fonti:
  - Confine regionale: Eurostat NUTS 2024 (ITF2) via EEA discomap
  - Bosco e tipo di foglia: Copernicus HRL Dominant Leaf Type 2018 (10 m)
    + Tree Cover Density 2018, via EEA discomap ImageServer
  - Quota (DEM 90 m) e meteo giornaliero: Open-Meteo Historical API
    (archive-api.open-meteo.com), con temperatura gia' riportata alla quota del punto

I layer statici sono messi in cache (cache/static.json); ogni settimana si
scarica solo il meteo. L'indice e' EURISTICO: non e' calibrato su ritrovamenti.
"""
import argparse, datetime as dt, io, json, math, os, sys, time
import numpy as np
import requests
import tifffile
from shapely.geometry import shape, Point, box
from shapely.prepared import prep

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
EEA = "https://image.discomap.eea.europa.eu/arcgis/rest/services"
ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"
BBOX = (13.93, 41.35, 15.17, 42.08)       # lon0, lat0, lon1, lat1
CELL = 0.02                                # cella di analisi (~1.7 x 2.2 km)
WX = 0.1                                   # griglia meteo (~8 x 11 km, come i modelli)
DAYS = 40                                  # giorni di storico meteo
LAPSE = 0.0065                             # gradi C per metro

S = requests.Session()
S.headers["User-Agent"] = "funghi-molise-map/1.0"


def get(url, params=None, kind="json", tries=6):
    for i in range(tries):
        try:
            r = S.get(url, params=params, timeout=180)
            if r.status_code == 429 or r.status_code >= 500:
                raise RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
            r.raise_for_status()
            if kind == "json":
                j = r.json()
                if isinstance(j, dict) and j.get("error"):
                    raise RuntimeError(j.get("reason"))
                return j
            return r.content
        except Exception as e:
            wait = 10 * (i + 1)
            print(f"  retry {i+1}/{tries} in {wait}s: {e}", file=sys.stderr)
            time.sleep(wait)
    raise SystemExit(f"Download fallito: {url}")


# ---------------------------------------------------------------- statico
def boundary():
    j = get(f"{EEA}/CORDA_NUTS/NUTS_2024/MapServer/0/query",
            {"where": "NUTS_ID='ITF2'", "outFields": "NUTS_ID", "outSR": 4326, "f": "geojson"})
    return j["features"][0]["geometry"]


def raster(service, w=2480, h=1460):
    b = get(f"{EEA}/GioLandPublic/{service}/ImageServer/exportImage",
            {"bbox": ",".join(map(str, BBOX)), "bboxSR": 4326, "imageSR": 4326,
             "size": f"{w},{h}", "format": "tiff", "pixelType": "U8",
             "interpolation": "RSP_NearestNeighbor", "f": "image"}, kind="bytes")
    return tifffile.imread(io.BytesIO(b))


def om_points(lats, lons, start, end, daily, batch=100, pause=3):
    """Interroga l'archivio Open-Meteo per molti punti, a blocchi."""
    out = []
    for i in range(0, len(lats), batch):
        p = {"latitude": ",".join(f"{x:.4f}" for x in lats[i:i+batch]),
             "longitude": ",".join(f"{x:.4f}" for x in lons[i:i+batch]),
             "start_date": start, "end_date": end, "daily": daily,
             "timezone": "Europe/Rome"}
        j = get(ARCHIVE, p)
        out.extend(j if isinstance(j, list) else [j])
        print(f"  open-meteo {min(i+batch, len(lats))}/{len(lats)}", file=sys.stderr)
        time.sleep(pause)
    return out


def build_static(today):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "static.json")
    if os.path.exists(path):
        return json.load(open(path))
    print("Costruisco i layer statici (una tantum)...", file=sys.stderr)
    geom = boundary()
    poly = shape(geom)
    pp = prep(poly)
    dlt = raster("HRL_DominantLeafType2018")    # 0 no alberi, 1 latifoglie, 2 conifere, 255 fuori
    tcd = raster("HRL_TreeCoverDensity_2018")   # 0-100 %, 255 fuori
    H, W = dlt.shape
    lon0, lat0, lon1, lat1 = BBOX
    cells = []
    lon = lon0 + CELL / 2
    while lon < lon1:
        lat = lat0 + CELL / 2
        while lat < lat1:
            if pp.contains(Point(lon, lat)):
                c0 = int((lon - CELL/2 - lon0) / (lon1 - lon0) * W)
                c1 = int((lon + CELL/2 - lon0) / (lon1 - lon0) * W)
                r0 = int((lat1 - (lat + CELL/2)) / (lat1 - lat0) * H)
                r1 = int((lat1 - (lat - CELL/2)) / (lat1 - lat0) * H)
                d = dlt[r0:r1, c0:c1]; t = tcd[r0:r1, c0:c1]
                valid = d != 255
                n = max(int(valid.sum()), 1)
                # "bosco" = pixel con albero e copertura chiome >= 30 %
                forest = valid & (d > 0) & (t >= 30) & (t <= 100)
                broad = forest & (d == 1); conif = forest & (d == 2)
                cells.append({"lat": round(lat, 4), "lon": round(lon, 4),
                              "forest": round(forest.sum() / n, 3),
                              "broad": round(broad.sum() / n, 3),
                              "conif": round(conif.sum() / n, 3)})
            lat += CELL
        lon += CELL
    print(f"  {len(cells)} celle in Molise", file=sys.stderr)
    # quota DEM 90 m dal servizio Open-Meteo (campo 'elevation'), con ripresa
    epath = os.path.join(CACHE, "elev.json")
    elev = json.load(open(epath)) if os.path.exists(epath) else {}
    y = (today - dt.timedelta(days=2)).isoformat()
    todo = [c for c in cells if f"{c['lat']},{c['lon']}" not in elev]
    for i in range(0, len(todo), 100):
        chunk = todo[i:i+100]
        res = om_points([c["lat"] for c in chunk], [c["lon"] for c in chunk], y, y, "temperature_2m_mean", pause=15)
        for c, r in zip(chunk, res):
            elev[f"{c['lat']},{c['lon']}"] = round(r["elevation"])
        json.dump(elev, open(epath, "w"))
    for c in cells:
        c["elev"] = elev[f"{c['lat']},{c['lon']}"]
    # griglia meteo
    wx = []
    lon = lon0 + WX / 2
    while lon < lon1:
        lat = lat0 + WX / 2
        while lat < lat1:
            if poly.buffer(WX / 2).contains(Point(lon, lat)):
                wx.append({"lat": round(lat, 4), "lon": round(lon, 4)})
            lat += WX
        lon += WX
    for c in cells:
        c["wx"] = min(range(len(wx)), key=lambda k: (wx[k]["lat"]-c["lat"])**2 + ((wx[k]["lon"]-c["lon"])*0.75)**2)
    st = {"boundary": geom, "cells": cells, "wx": wx, "cell": CELL,
          "built": today.isoformat()}
    json.dump(st, open(path, "w"))
    return st


# ---------------------------------------------------------------- indice
def ramp(x, a, b):
    """0 sotto a, 1 sopra b, lineare in mezzo."""
    return float(min(1.0, max(0.0, (x - a) / (b - a))))


def trapezoid(x, a, b, c, d):
    if x <= a or x >= d: return 0.0
    if b <= x <= c: return 1.0
    return (x - a) / (b - a) if x < b else (d - x) / (d - c)


def score_cell(c, w, wx_elev):
    """Restituisce (indice 0-100, componenti). w = serie giornaliere (piu' vecchio -> ieri)."""
    P = np.array(w["precipitation_sum"], float); P[np.isnan(P)] = 0
    dz = c["elev"] - wx_elev
    Tn = np.array(w["temperature_2m_min"], float) - LAPSE * dz
    Tx = np.array(w["temperature_2m_max"], float) - LAPSE * dz
    Tm = np.array(w["temperature_2m_mean"], float) - LAPSE * dz
    E = np.array(w["et0_fao_evapotranspiration"], float); E[np.isnan(E)] = 0
    n = len(P)
    # 1) pioggia d'innesco: massima pioggia cumulata in 5 giorni, cadutati 7-25 giorni fa
    trig = 0.0
    for end in range(n - 7, n - 26, -1):   # end = indice dell'ultimo giorno della finestra
        if end - 4 < 0: break
        trig = max(trig, P[end-4:end+1].sum())
    # 2) umidita' recente: bilancio pioggia - evapotraspirazione ultimi 14 giorni
    bal14 = P[-14:].sum() - E[-14:].sum()
    rain14 = P[-14:].sum()
    s_trig = ramp(trig, 10, 45)
    s_wet = ramp(bal14, -15, 20)
    W = 0.6 * s_trig + 0.4 * s_wet
    # 3) temperatura ultimi 10 giorni alla quota della cella
    tm10 = float(np.nanmean(Tm[-10:]))
    s_temp = trapezoid(tm10, 3, 9, 17, 24)
    frost = int((Tn[-10:] <= 0).sum()); hard = int((Tn[-10:] <= -3).sum())
    f_frost = (0.8 ** frost) * (0.6 ** hard)
    swing = float(np.nanmean(Tx[-10:] - Tn[-10:]))
    T = s_temp * f_frost
    # 4) habitat: frazione di bosco nella cella (piena a 50 %)
    Hb = ramp(c["forest"], 0.05, 0.5)
    idx = 100 * Hb * W * T
    return round(idx), {
        "trig": round(trig, 1), "rain14": round(rain14, 1), "bal14": round(bal14, 1),
        "tm10": round(tm10, 1), "tmin10": round(float(np.nanmin(Tn[-10:])), 1),
        "frost": frost, "swing": round(swing, 1),
        "W": round(W, 2), "T": round(T, 2), "H": round(Hb, 2)}


def habitat_label(c):
    if c["forest"] < 0.05: return "poco bosco"
    conif = c["conif"] > c["broad"]
    e = c["elev"]
    if conif: return "conifere (perlopiu' rimboschimenti)"
    if e >= 1000: return "latifoglie in quota (faggeta probabile)"
    if e >= 500: return "latifoglie collinari (cerrete, querceti, castagneti)"
    return "latifoglie di bassa quota (querceti, boschi ripariali)"


def species_hint(c, month):
    """Indicazione grossolana delle specie in stagione, per quota e tipo di bosco."""
    e, conif = c["elev"], c["conif"] > c["broad"]
    if conif:
        return "pinaroli e sanguinelli (Suillus, Lactarius deliciosus)"
    if month in (11, 12):
        if e >= 1100: return "stagione quasi chiusa: chiodini, trombette dei morti nei giorni miti"
        if e >= 600: return "chiodini, trombette dei morti, ultimi porcini nelle esposizioni calde"
        return "chiodini, prataioli e mazze di tamburo ai margini, porcini tardivi rari"
    if month in (9, 10):
        if e >= 900: return "porcino (B. edulis) e ovoli bianchi in faggeta"
        return "porcino nero (B. aereus), ovolo buono, mazze di tamburo"
    if month in (5, 6, 7, 8):
        return "porcino estivo (B. aestivalis), finferli, russule"
    return "prugnoli e spugnole in primavera; inverno quasi assente"


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--date", default=None, help="data di riferimento (default oggi)")
    ap.add_argument("--skip-if-recent", type=int, default=0,
                    help="esci senza fare nulla se c'e' gia' un calcolo negli ultimi N giorni")
    a = ap.parse_args()
    today = dt.date.fromisoformat(a.date) if a.date else dt.date.today()
    hist = os.path.join(HERE, "history")
    if a.skip_if_recent and os.path.isdir(hist):
        done = [f[:10] for f in os.listdir(hist) if f.endswith(".json")]
        if any(0 <= (today - dt.date.fromisoformat(d)).days < a.skip_if_recent for d in done):
            print("Calcolo gia' fatto di recente, salto.", file=sys.stderr)
            return
    st = build_static(today)
    end = today - dt.timedelta(days=1)
    start = end - dt.timedelta(days=DAYS - 1)
    print(f"Meteo {start} -> {end} su {len(st['wx'])} punti", file=sys.stderr)
    daily = "precipitation_sum,temperature_2m_min,temperature_2m_max,temperature_2m_mean,et0_fao_evapotranspiration"
    res = om_points([p["lat"] for p in st["wx"]], [p["lon"] for p in st["wx"]],
                    start.isoformat(), end.isoformat(), daily, batch=25, pause=12)
    out_cells = []
    for c in st["cells"]:
        r = res[c["wx"]]
        idx, comp = score_cell(c, r["daily"], r["elevation"])
        out_cells.append({**{k: c[k] for k in ("lat", "lon", "forest", "broad", "conif", "elev")},
                          "i": idx, **comp, "hab": habitat_label(c),
                          "sp": species_hint(c, today.month)})
    # serie regionale media per il grafico
    days = res[0]["daily"]["time"]
    P = np.nanmean([[x or 0 for x in r["daily"]["precipitation_sum"]] for r in res], axis=0)
    Tn = np.nanmean([r["daily"]["temperature_2m_min"] for r in res], axis=0)
    Tx = np.nanmean([r["daily"]["temperature_2m_max"] for r in res], axis=0)
    data = {"generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
            "ref": today.isoformat(), "start": start.isoformat(), "end": end.isoformat(),
            "cell": st["cell"], "boundary": st["boundary"], "cells": out_cells,
            "series": {"days": days, "P": [round(x, 1) for x in P],
                       "Tn": [round(x, 1) for x in Tn], "Tx": [round(x, 1) for x in Tx]}}
    # confronto con la settimana precedente
    os.makedirs(hist, exist_ok=True)
    prev = sorted(f for f in os.listdir(hist) if f.endswith(".json") and f[:10] < today.isoformat())
    if prev:
        old = json.load(open(os.path.join(hist, prev[-1])))
        om = {(o["lat"], o["lon"]): o["i"] for o in old["cells"]}
        for c in data["cells"]:
            if (c["lat"], c["lon"]) in om: c["d"] = c["i"] - om[(c["lat"], c["lon"])]
        data["prev"] = prev[-1][:10]
    slim = {"ref": data["ref"], "cells": [{"lat": c["lat"], "lon": c["lon"], "i": c["i"]} for c in data["cells"]]}
    json.dump(slim, open(os.path.join(hist, f"{today.isoformat()}.json"), "w"))
    os.makedirs(a.out, exist_ok=True)
    json.dump(data, open(os.path.join(a.out, "data.json"), "w"), separators=(",", ":"))
    tpl = open(os.path.join(HERE, "template.html"), encoding="utf-8").read()
    html = tpl.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
    open(os.path.join(a.out, "mappa-funghi-molise.html"), "w", encoding="utf-8").write(html)
    f = [c for c in data["cells"] if c["forest"] >= 0.05]
    print(f"OK: {len(data['cells'])} celle, {len(f)} boscate, indice medio boschi "
          f"{np.mean([c['i'] for c in f]):.0f}, max {max(c['i'] for c in f)}", file=sys.stderr)


if __name__ == "__main__":
    main()
