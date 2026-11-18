# Funghi in Molise

Mappa settimanale di un indice euristico (non tarato) sulla probabilita' di trovare funghi nei boschi del Molise.

- `build.py` scarica i dati e calcola l'indice; `template.html` e' la pagina.
- `.github/workflows/mappa.yml` gira ogni giovedi' mattina (con secondo tentativo il venerdi') e pubblica su GitHub Pages.
- `cache/` contiene i layer statici (confine, bosco, quota), `history/` lo storico settimanale dell'indice.

Fonti: Open-Meteo Historical API, Copernicus HRL Dominant Leaf Type e Tree Cover Density 2018, Eurostat NUTS 2024 via EEA.
Il metodo e i limiti sono spiegati in fondo alla pagina della mappa.
