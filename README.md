# May 2026 Kharg oil slick

A Data Desk research notebook tracking a suspected oil spill at Iran's main crude export terminal at Kharg Island, using the Copernicus Sentinel-2 (optical) and Sentinel-1 (radar) record from 6 May 2026 onwards.

**Published notebook:** https://research.datadesk.eco/kharg-slick/

## Layout

- `docs/index.html` — notebook source (Observable Notebook Kit 2.0)
- `docs/assets/scenes/` — rendered per-pass PNGs
- `data/scenes.json` — manifest of fetched Sentinel scenes
- `data/slicks.geojson` — automated SAR-derived slick polygons + manual S2 sheen traces
- `data/infrastructure.geojson` — Mapstand pipelines/platforms/terminals across the Persian Gulf
- `scripts/fetch_scenes.py` — fetches imagery from the Copernicus Data Space Ecosystem, segments slicks from S1 VV, renders the PNGs and writes the manifests

## Rebuilding the data

```sh
export CDSE_S3_ACCESS_KEY=...
export CDSE_S3_SECRET_KEY=...
make data    # pull new scenes, re-segment, re-render
make build   # build the static site to docs/.observable/dist/
make preview # local preview
```

Each S1 pass is segmented by log-amplitude thresholding (~4.4 dB below the scene's 70th-percentile background, despeckled and morphologically closed); scenes whose largest detection exceeds 250 km² are rejected as likely calm-wind contamination. Manual S2 sheen traces live in `data/may-*-georef.geojson`.

## Credits

Imagery: [Copernicus Data Space Ecosystem](https://dataspace.copernicus.eu/). Infrastructure: © [Mapstand](https://www.mapstand.com/). Basemap © OpenStreetMap contributors, © CARTO.
