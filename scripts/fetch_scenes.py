#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "boto3>=1.34",
#     "rasterio>=1.3",
#     "numpy>=1.26",
#     "pillow>=10.0",
#     "requests>=2.31",
#     "matplotlib>=3.8",
#     "pyproj>=3.6",
#     "scipy>=1.11",
#     "shapely>=2.0",
# ]
# ///
"""
Fetch all Sentinel-1 (SAR) and Sentinel-2 (optical) scenes over Kharg Island
between 2026-05-04 and today from the Copernicus Data Space Ecosystem, clip
each one to a fixed AOI, and render a PNG for the notebook.

Reproducibility:
  - Scenes are discovered via the public CDSE OData catalogue.
  - Pixel data is pulled directly from CDSE S3 (eodata.dataspace.copernicus.eu).
  - All processing parameters (AOI, render gain, output size) are constants below.

Required environment variables:
  CDSE_S3_ACCESS_KEY, CDSE_S3_SECRET_KEY

Outputs:
  docs/assets/scenes/<sensing_date>_<platform>_<product>.png
  data/scenes.json   (manifest consumed by the notebook)
"""
from __future__ import annotations

import json
import math
import os
import sys
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path

import boto3
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import rasterio.features
import rasterio.transform
import scipy.ndimage as ndi
import xml.etree.ElementTree as ET
from botocore.config import Config
from PIL import Image
from pyproj import Transformer
from rasterio.warp import Resampling, reproject, transform_bounds
from rasterio.windows import from_bounds
from shapely.geometry import Point, mapping
from shapely.geometry import shape as shapely_shape
from shapely.ops import transform as shapely_transform
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parent.parent
ASSETS = ROOT / "docs" / "assets" / "scenes"
MANIFEST = ROOT / "data" / "scenes.json"
INFRA = ROOT / "data" / "infrastructure.json"
INFRA_GEOJSON = ROOT / "data" / "infrastructure.geojson"
SLICKS = ROOT / "data" / "slicks.geojson"

# Wider AOI used for slick segmentation — fixed across all S1 scenes so the
# resulting polygons live in a common reference frame for the overview map.
SEGMENT_AOI = dict(west=49.85, south=28.55, east=50.55, north=29.40)
SEGMENT_RES_M = 50.0

# Output projection — UTM zone 39N covers all of the AOI.
DST_CRS = "EPSG:32639"

# Mapstand API: WMS+KML for vector overlay layers.
MAPSTAND_WMS = "https://app.mapstand.com/geoserver/ows/mps"
MAPSTAND_KEY = "63abb313-a3ce-4dd8-b4de-d093f897018a"
MAPSTAND_LAYERS = {
    "pipeline": "mps_mapping_pipeline",
    "platform": "mps_mapping_platform",
    "terminal": "mps_mapping_terminal",
    "floatingfacility": "mps_mapping_floatingfacility",
}

# Catalogue search AOI — a tight bbox around Kharg used purely to discover
# intersecting scenes; the actual render framing is computed separately and
# can be wider than this without affecting scene discovery (the S2 T39RVN
# tile and the S1 IW swaths are much larger than this bbox).
CATALOG_AOI = dict(west=50.00, south=28.75, east=50.45, north=29.40)
# Mapstand fetch spans the whole Persian Gulf so the map can show the spill
# in regional infrastructure context. Tiled fetch handles the area.
MAPSTAND_AOI = dict(west=47.5, south=23.0, east=57.0, north=30.5)
MAPSTAND_TILE_M = 200_000  # 200 km tiles in EPSG:3857

# The render AOI is derived programmatically from the union of detected slick
# polygons plus Kharg Island, padded outwards — so every scene is framed
# consistently around the actual observed extent of the spill. The pad is
# generous because the optical sheen on calm days extends well beyond the
# SAR-thresholded core polygons that define the union.
KHARG_LON, KHARG_LAT = 50.32, 29.25
AOI_PAD_DEG = 0.15  # ~15 km

START = "2026-05-04T00:00:00.000Z"
END = (datetime.now(timezone.utc).replace(hour=23, minute=59, second=59).strftime("%Y-%m-%dT%H:%M:%S.000Z"))

# Output pixel density (m/pixel) — pixel dimensions are computed per framing.
TARGET_M_PER_PX = 35.0

S3_ENDPOINT = "https://eodata.dataspace.copernicus.eu"
ODATA = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"


@dataclass
class Scene:
    date: str           # YYYY-MM-DD
    sensing_time: str   # ISO8601 with timezone
    sensor: str         # "Sentinel-2" | "Sentinel-1"
    platform: str       # e.g. "S2C", "S1A"
    product_id: str
    product_name: str
    s3_path: str
    orbit_direction: str | None
    relative_orbit: int | None
    png: str            # path relative to docs/


def odata_query(filter_clause: str) -> list[dict]:
    url = (
        ODATA
        + "?$filter="
        + urllib.parse.quote(filter_clause)
        + "&$expand=Attributes&$top=100&$orderby=ContentDate/Start"
    )
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r).get("value", [])


def catalogue_scenes() -> list[dict]:
    a = CATALOG_AOI
    aoi_poly = (
        f"POLYGON(({a['west']} {a['south']},{a['east']} {a['south']},"
        f"{a['east']} {a['north']},{a['west']} {a['north']},"
        f"{a['west']} {a['south']}))"
    )
    base = (
        f"OData.CSC.Intersects(area=geography'SRID=4326;{aoi_poly}') "
        f"and ContentDate/Start ge {START} and ContentDate/Start lt {END}"
    )
    # Kharg sits inside Sentinel-2 MGRS tile T39RVN. The adjacent T39RUN also
    # intersects the AOI but cuts off the island; ignore it. Orbit R006 is the
    # descending pass that captures the AOI in full — orbit R106 only catches
    # a sliver of the east edge and is filtered out.
    s2 = odata_query(
        f"Collection/Name eq 'SENTINEL-2' and contains(Name,'MSIL2A') "
        f"and contains(Name,'_T39RVN_') and contains(Name,'_R006_') and {base}"
    )
    s1 = odata_query(
        f"Collection/Name eq 'SENTINEL-1' and contains(Name,'GRDH') "
        f"and not endswith(Name,'_COG.SAFE') and {base}"
    )
    return s2 + s1


def s3_client():
    return boto3.client(
        "s3",
        endpoint_url=S3_ENDPOINT,
        aws_access_key_id=os.environ["CDSE_S3_ACCESS_KEY"],
        aws_secret_access_key=os.environ["CDSE_S3_SECRET_KEY"],
        region_name="default",
        config=Config(signature_version="s3v4"),
    )


def list_s3(s3, prefix: str) -> list[str]:
    keys: list[str] = []
    token = None
    while True:
        kw = dict(Bucket="eodata", Prefix=prefix, MaxKeys=1000)
        if token:
            kw["ContinuationToken"] = token
        r = s3.list_objects_v2(**kw)
        keys.extend(o["Key"] for o in r.get("Contents", []))
        token = r.get("NextContinuationToken")
        if not token:
            break
    return keys


def download(s3, key: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    s3.download_file("eodata", key, str(tmp))
    tmp.rename(dest)
    return dest


def compute_aoi(slick_features: list[dict]) -> dict:
    """Bounds of all slick geometries plus Kharg, padded by AOI_PAD_DEG."""
    geoms = [shapely_shape(f["geometry"]) for f in slick_features]
    geoms.append(Point(KHARG_LON, KHARG_LAT).buffer(0.012))  # ~1.3 km buffer round Kharg
    union = unary_union(geoms)
    xmin, ymin, xmax, ymax = union.bounds
    return dict(
        west=xmin - AOI_PAD_DEG,
        south=ymin - AOI_PAD_DEG,
        east=xmax + AOI_PAD_DEG,
        north=ymax + AOI_PAD_DEG,
    )


def render_size(aoi: dict) -> tuple[int, int]:
    """Pixel dimensions for an AOI at TARGET_M_PER_PX."""
    utm_w, utm_s, utm_e, utm_n = transform_bounds(
        "EPSG:4326", DST_CRS, aoi["west"], aoi["south"], aoi["east"], aoi["north"]
    )
    return (
        int(round((utm_e - utm_w) / TARGET_M_PER_PX)),
        int(round((utm_n - utm_s) / TARGET_M_PER_PX)),
    )


def render_s2(s3, item: dict, cache: Path, aoi: dict) -> np.ndarray | None:
    """Sentinel-2 L2A: clip the pre-computed True Color Image (TCI) to AOI."""
    safe_prefix = item["S3Path"].lstrip("/").replace("eodata/", "") + "/"
    keys = list_s3(s3, safe_prefix + "GRANULE/")
    tci = next((k for k in keys if k.endswith("TCI_10m.jp2")), None)
    if not tci:
        print(f"  ! no TCI band found for {item['Name']}")
        return None
    local = download(s3, tci, cache / Path(tci).name)
    out_w, out_h = render_size(aoi)
    with rasterio.open(local) as src:
        win_bounds = transform_bounds(
            "EPSG:4326", src.crs, aoi["west"], aoi["south"], aoi["east"], aoi["north"]
        )
        window = from_bounds(*win_bounds, transform=src.transform)
        rgb = src.read(window=window, out_shape=(3, out_h, out_w), resampling=Resampling.bilinear)
    img = rgb.astype(np.float32) / 255.0
    img = np.clip(img * 1.35, 0, 1) ** 0.85
    return (img * 255).astype(np.uint8).transpose(1, 2, 0)


def render_s1(s3, item: dict, cache: Path, aoi: dict) -> np.ndarray | None:
    """Sentinel-1 GRD: download VV measurement, warp to UTM, clip, log-stretch."""
    safe_prefix = item["S3Path"].lstrip("/").replace("eodata/", "") + "/"
    keys = list_s3(s3, safe_prefix + "measurement/")
    vv = next((k for k in keys if "-vv-" in k.lower() and k.endswith(".tiff")), None)
    if not vv:
        print(f"  ! no VV measurement for {item['Name']}")
        return None
    local = download(s3, vv, cache / Path(vv).name)
    out_w, out_h = render_size(aoi)

    with rasterio.open(local) as src:
        gcps, gcp_crs = src.gcps
        if not gcps:
            print(f"  ! {vv} has no GCPs")
            return None
        dst_w, dst_s, dst_e, dst_n = transform_bounds(
            "EPSG:4326", DST_CRS, aoi["west"], aoi["south"], aoi["east"], aoi["north"]
        )
        transform = rasterio.transform.from_bounds(dst_w, dst_s, dst_e, dst_n, out_w, out_h)
        dst = np.zeros((out_h, out_w), dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_crs=gcp_crs,
            gcps=gcps,
            dst_crs=DST_CRS,
            dst_transform=transform,
            resampling=Resampling.bilinear,
            src_nodata=0,
            dst_nodata=0,
        )

    # Convert DN to dB-ish log amplitude, then percentile-stretch for display.
    amp = dst
    valid = amp > 0
    if not valid.any():
        return None
    log = np.zeros_like(amp)
    log[valid] = np.log10(amp[valid])
    # 1st–99.5th percentile contrast stretch on the marine background.
    lo, hi = np.percentile(log[valid], [1, 99.5])
    img = np.clip((log - lo) / max(hi - lo, 1e-6), 0, 1)
    img[~valid] = 0
    img = (img * 255).astype(np.uint8)
    return np.stack([img, img, img], axis=-1)


def write_png(arr: np.ndarray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(arr).save(path, "PNG", optimize=True)


# ---------------------------------------------------------------------------
# Mapstand infrastructure overlay
# ---------------------------------------------------------------------------

KML_NS = {"k": "http://www.opengis.net/kml/2.2"}


def _parse_kml_coords(text: str) -> list[list[float]]:
    """KML coords are 'lon,lat[,alt]' space-separated."""
    pts = []
    for tok in text.strip().split():
        parts = tok.split(",")
        pts.append([float(parts[0]), float(parts[1])])
    return pts


def _parse_kml_placemark(pm) -> dict | None:
    """Return GeoJSON-style {type, coords} for the geometry in a placemark.
    Mapstand wraps the real geometry inside a MultiGeometry alongside a label
    Point — we want the real geometry, not the label."""
    mg = pm.find(".//k:MultiGeometry", KML_NS)
    container = mg if mg is not None else pm
    ls = container.findall(".//k:LineString/k:coordinates", KML_NS)
    if ls:
        lines = [_parse_kml_coords(e.text) for e in ls if e.text]
        if len(lines) == 1:
            return {"type": "LineString", "coords": lines[0]}
        return {"type": "MultiLineString", "coords": lines}
    poly = container.find(".//k:Polygon", KML_NS)
    if poly is not None:
        outer = poly.find(".//k:outerBoundaryIs/k:LinearRing/k:coordinates", KML_NS)
        if outer is not None and outer.text:
            return {"type": "Polygon", "coords": [_parse_kml_coords(outer.text)]}
    pt = container.find(".//k:Point/k:coordinates", KML_NS)
    if pt is not None and pt.text:
        c = _parse_kml_coords(pt.text)
        if c:
            return {"type": "Point", "coords": c[0]}
    return None


def _mapstand_tiles_3857() -> list[tuple[float, float, float, float]]:
    tx = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    x_min, y_min = tx.transform(MAPSTAND_AOI["west"], MAPSTAND_AOI["south"])
    x_max, y_max = tx.transform(MAPSTAND_AOI["east"], MAPSTAND_AOI["north"])
    tiles = []
    y = y_min
    while y < y_max:
        x = x_min
        while x < x_max:
            tiles.append((x, y, min(x + MAPSTAND_TILE_M, x_max), min(y + MAPSTAND_TILE_M, y_max)))
            x += MAPSTAND_TILE_M
        y += MAPSTAND_TILE_M
    return tiles


def fetch_infrastructure() -> dict:
    """Pull Mapstand vectors across the whole Gulf via tiled WMS+KML requests."""
    if INFRA.exists():
        return json.loads(INFRA.read_text())

    tiles = _mapstand_tiles_3857()
    out: dict[str, list[dict]] = {}
    for label, layer in MAPSTAND_LAYERS.items():
        feats: list[dict] = []
        seen_ids: set[str] = set()
        for i, (tx1, ty1, tx2, ty2) in enumerate(tiles, 1):
            params = {
                "SERVICE": "WMS",
                "VERSION": "1.1.1",
                "REQUEST": "GetMap",
                "FORMAT": "application/vnd.google-earth.kml+xml",
                "LAYERS": layer,
                "SRS": "EPSG:3857",
                "BBOX": f"{tx1},{ty1},{tx2},{ty2}",
                "WIDTH": "2048",
                "HEIGHT": "2048",
                "apikey": MAPSTAND_KEY,
            }
            url = MAPSTAND_WMS + "?" + urllib.parse.urlencode(params)
            req = urllib.request.Request(url, headers={"Referer": "https://app.mapstand.com/"})
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    body = r.read()
            except Exception as e:  # noqa: BLE001
                print(f"  ! Mapstand {layer} tile {i}/{len(tiles)}: {e}")
                continue
            try:
                tree = ET.fromstring(body)
            except ET.ParseError:
                continue
            for pm in tree.findall(".//k:Placemark", KML_NS):
                pid = pm.get("id") or ""
                if pid and pid in seen_ids:
                    continue
                g = _parse_kml_placemark(pm)
                if not g:
                    continue
                if pid:
                    seen_ids.add(pid)
                feats.append(g)
        out[label] = feats
        print(f"  Mapstand {label}: {len(feats)} features (across {len(tiles)} tiles)")
    INFRA.parent.mkdir(parents=True, exist_ok=True)
    INFRA.write_text(json.dumps(out, indent=2))

    # Also emit a MapLibre-friendly FeatureCollection.
    features: list[dict] = []
    for kind, items in out.items():
        for item in items:
            features.append({
                "type": "Feature",
                "properties": {"kind": kind},
                "geometry": {"type": item["type"], "coordinates": item["coords"]},
            })
    INFRA_GEOJSON.write_text(json.dumps({
        "type": "FeatureCollection",
        "features": features,
    }))
    return out


OVERLAY_COLOR = "#00ff66"  # bright green for Mapstand vectors

# Manual sheen tracing from S2 turned out to be too noisy to be useful — the
# 6 May morning S2 sheen is essentially the same footprint as the afternoon
# SAR core a few hours later, and the 11 May S2 sheen is filamentary and
# largely obscured by cloud cover. The map therefore shows only SAR-derived
# core slick polygons; the S2 observations are covered narratively in the
# day-by-day sections.
MANUAL_SHEENS: dict[tuple[str, str], list[list[list[float]]]] = {}


def segment_s1(s3, item: dict, cache: Path) -> list[dict]:
    """Extract slick polygons (WGS84 GeoJSON geometries) from one S1 GRD scene.

    Method: reproject the VV GRD to UTM at 50 m/px, despeckle with a 5×5 median
    filter, compute the 60th-percentile log-amplitude (an estimator of the
    marine background that is robust to large slicks pulling the mean down),
    threshold at ~6 dB below background, exclude land/vessel neighbourhoods,
    close small gaps, drop components smaller than 0.5 km², and vectorise.
    """
    safe_prefix = item["S3Path"].lstrip("/").replace("eodata/", "") + "/"
    keys = list_s3(s3, safe_prefix + "measurement/")
    vv = next((k for k in keys if "-vv-" in k.lower() and k.endswith(".tiff")), None)
    if not vv:
        return []
    local = download(s3, vv, cache / Path(vv).name)

    a = SEGMENT_AOI
    dst_w, dst_s, dst_e, dst_n = transform_bounds(
        "EPSG:4326", DST_CRS, a["west"], a["south"], a["east"], a["north"]
    )
    out_w = int(round((dst_e - dst_w) / SEGMENT_RES_M))
    out_h = int(round((dst_n - dst_s) / SEGMENT_RES_M))

    with rasterio.open(local) as src:
        gcps, gcp_crs = src.gcps
        if not gcps:
            return []
        transform = rasterio.transform.from_bounds(dst_w, dst_s, dst_e, dst_n, out_w, out_h)
        amp = np.zeros((out_h, out_w), dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=amp,
            src_crs=gcp_crs,
            gcps=gcps,
            dst_crs=DST_CRS,
            dst_transform=transform,
            resampling=Resampling.bilinear,
            src_nodata=0,
            dst_nodata=0,
        )

    valid = amp > 1
    if int(valid.sum()) < 1000:
        return []

    log_amp = np.zeros_like(amp)
    log_amp[valid] = np.log10(amp[valid])
    log_smooth = ndi.median_filter(log_amp, size=5)

    # 70th percentile estimates the marine background when extensive slicks
    # pull the median down. SAR amplitude → dB = 20·log10(amp), so 5 dB ≈ 0.25
    # in log10(amp).
    bg = float(np.percentile(log_smooth[valid], 70))
    threshold = bg - 0.22  # ~4.4 dB below background
    mask = (log_smooth < threshold) & valid

    # Exclude only the immediate vicinity of very bright targets (Kharg, large
    # vessels) — 2-iter dilation at 50 m/px = a ~100 m buffer, enough to
    # discard sidelobe artefacts but not so wide that it eats the slick.
    bright = (log_amp > bg + 0.40) & valid
    mask = mask & ~ndi.binary_dilation(bright, iterations=2)

    # Close small gaps within the slick; skip an opening step (it shrinks
    # narrow filaments and trailing fragments that are real slick features).
    mask = ndi.binary_closing(mask, iterations=2)

    lbl, n = ndi.label(mask)
    min_pixels = int(0.25e6 / (SEGMENT_RES_M ** 2))  # 0.25 km²
    sizes = ndi.sum(mask, lbl, range(n + 1))
    keep = sizes >= min_pixels
    keep[0] = False
    mask = keep[lbl]

    to_wgs = Transformer.from_crs(DST_CRS, "EPSG:4326", always_xy=True).transform
    polys: list[dict] = []
    for geom_dict, val in rasterio.features.shapes(
        mask.astype(np.uint8), mask=mask, transform=transform
    ):
        if val != 1:
            continue
        geom = shapely_shape(geom_dict)
        if geom.area < 0.5e6:
            continue
        geom = geom.simplify(80, preserve_topology=True)
        geom_wgs = shapely_transform(to_wgs, geom)
        # Spatial filter — the slick originates at Kharg (50.32°E, 29.25°N)
        # and drifts south with prevailing winds. Anything noticeably north of
        # the island or far west of it is open Persian Gulf and almost
        # certainly a calm-water false positive.
        c = geom_wgs.centroid
        if c.y > 29.30 or c.x < 49.97:
            continue
        polys.append(mapping(geom_wgs))
    return polys


def composite(raster: np.ndarray, aoi: dict, infra: dict, path: Path) -> None:
    """Render the raster with a 0.1° lat/lon graticule and Mapstand overlay."""
    out_h, out_w = raster.shape[:2]
    utm_w, utm_s, utm_e, utm_n = transform_bounds(
        "EPSG:4326", DST_CRS, aoi["west"], aoi["south"], aoi["east"], aoi["north"]
    )
    to_utm = Transformer.from_crs("EPSG:4326", DST_CRS, always_xy=True)

    fig = plt.figure(figsize=(out_w / 200, out_h / 200), dpi=200)
    fig.patch.set_alpha(0)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("none")
    ax.imshow(
        raster,
        extent=(utm_w, utm_e, utm_s, utm_n),
        origin="upper",
        interpolation="bilinear",
        aspect="auto",
    )
    ax.set_xlim(utm_w, utm_e)
    ax.set_ylim(utm_s, utm_n)
    ax.set_axis_off()

    # 0.1° lat/lon graticule. Labels stroked in black for legibility on
    # both light (cloud) and dark (sea) backgrounds.
    grid_kw = dict(color="#e0e0e0", linewidth=0.5, alpha=0.35, zorder=2)
    stroke = [path_effects.Stroke(linewidth=1.2, foreground="black", alpha=0.85),
              path_effects.Normal()]
    label_kw = dict(color="white", fontsize=6.5, alpha=0.95, family="sans-serif",
                    path_effects=stroke)

    def gridline(lons, lats, **kw):
        xs, ys = to_utm.transform(lons, lats)
        ax.plot(xs, ys, **kw)

    step = 0.1
    lat0 = math.floor(aoi["south"] / step) * step
    lon0 = math.floor(aoi["west"] / step) * step
    lats = [round(lat0 + i * step, 2) for i in range(int((aoi["north"] - lat0) / step) + 2)]
    lons = [round(lon0 + i * step, 2) for i in range(int((aoi["east"] - lon0) / step) + 2)]
    sample_x = np.linspace(aoi["west"], aoi["east"], 50)
    for lat in lats:
        if aoi["south"] + 0.005 <= lat <= aoi["north"] - 0.003:
            gridline(sample_x, np.full_like(sample_x, lat), **grid_kw)
            x, y = to_utm.transform(aoi["west"] + 0.005, lat + 0.003)
            ax.text(x, y, f"{lat:.1f}°N", verticalalignment="bottom",
                    horizontalalignment="left", **label_kw)
    sample_y = np.linspace(aoi["south"], aoi["north"], 50)
    for lon in lons:
        if aoi["west"] + 0.003 <= lon <= aoi["east"] - 0.005:
            gridline(np.full_like(sample_y, lon), sample_y, **grid_kw)
            x, y = to_utm.transform(lon + 0.003, aoi["south"] + 0.005)
            ax.text(x, y, f"{lon:.1f}°E", verticalalignment="bottom",
                    horizontalalignment="left", **label_kw)

    # Mapstand overlay, all in one neon green at 1px.
    def to_utm_xy(coords):
        lons = [c[0] for c in coords]
        lats = [c[1] for c in coords]
        return to_utm.transform(lons, lats)

    # matplotlib linewidth is in points (1/72in); at 200dpi a 1px line is ~0.36pt.
    line_kw = dict(color=OVERLAY_COLOR, linewidth=0.4, alpha=0.95, zorder=3)
    point_kw = dict(facecolors="none", edgecolors=OVERLAY_COLOR, linewidths=0.4,
                    alpha=0.95, zorder=4)

    for g in infra.get("pipeline", []):
        lines = g["coords"] if g["type"] == "MultiLineString" else [g["coords"]]
        for ln in lines:
            xs, ys = to_utm_xy(ln)
            ax.plot(xs, ys, **line_kw)

    for g in infra.get("platform", []):
        if g["type"] == "Point":
            x, y = to_utm.transform(g["coords"][0], g["coords"][1])
            ax.scatter([x], [y], s=14, marker="s", **point_kw)

    for g in infra.get("floatingfacility", []):
        if g["type"] == "Point":
            x, y = to_utm.transform(g["coords"][0], g["coords"][1])
            ax.scatter([x], [y], s=16, marker="o", **point_kw)

    for g in infra.get("terminal", []):
        if g["type"] == "Polygon":
            for ring in g["coords"]:
                xs, ys = to_utm_xy(ring)
                ax.plot(xs, ys, **line_kw)
        elif g["type"] == "Point":
            x, y = to_utm.transform(g["coords"][0], g["coords"][1])
            ax.scatter([x], [y], s=22, marker="D", **point_kw)

    fig.savefig(path, dpi=200, pil_kwargs={"optimize": True})
    plt.close(fig)


def main() -> int:
    cache = ROOT / "data" / "scenes_cache"
    cache.mkdir(parents=True, exist_ok=True)
    s3 = s3_client()

    print("Fetching Mapstand infrastructure overlay...")
    infra = fetch_infrastructure()

    items = catalogue_scenes()
    print(f"Found {len(items)} candidate products.")

    # Dedupe by sensor + sensing minute (S2 captures often have two granules).
    unique: list[dict] = []
    seen: set[str] = set()
    for item in sorted(items, key=lambda x: x["ContentDate"]["Start"]):
        sensor = "Sentinel-2" if item["Name"].startswith("S2") else "Sentinel-1"
        key = f"{sensor}|{item['ContentDate']['Start'][:16]}"
        if key in seen:
            continue
        seen.add(key)
        unique.append(item)

    # --- Pass 1: segment SAR scenes, collect slick polygons ---------------
    slick_features: list[dict] = []
    for item in unique:
        sensing = item["ContentDate"]["Start"]
        sensor = "Sentinel-2" if item["Name"].startswith("S2") else "Sentinel-1"
        date = sensing[:10]
        platform = item["Name"][:3]
        if sensor == "Sentinel-1":
            print(f"  segmenting {date} {platform}...")
            for geom in segment_s1(s3, item, cache):
                slick_features.append({
                    "type": "Feature",
                    "geometry": geom,
                    "properties": {
                        "date": date, "sensing_time": sensing,
                        "sensor": sensor, "platform": platform,
                        "kind": "core",
                        "source": "S1 VV threshold (~4.4 dB below scene background)",
                    },
                })
        if (date, sensor) in MANUAL_SHEENS:
            for ring in MANUAL_SHEENS[(date, sensor)]:
                slick_features.append({
                    "type": "Feature",
                    "geometry": {"type": "Polygon", "coordinates": [ring]},
                    "properties": {
                        "date": date, "sensing_time": sensing,
                        "sensor": sensor, "platform": platform,
                        "kind": "sheen",
                        "source": "manual trace from S2 true colour",
                    },
                })

    SLICKS.write_text(json.dumps({
        "type": "FeatureCollection", "features": slick_features,
    }))
    print(f"Wrote {SLICKS.relative_to(ROOT)} with {len(slick_features)} slick polygons.")

    # --- Compute common AOI from slick polygon union + Kharg --------------
    aoi = compute_aoi(slick_features)
    print(
        f"Common render AOI: {aoi['west']:.3f}–{aoi['east']:.3f}°E, "
        f"{aoi['south']:.3f}–{aoi['north']:.3f}°N"
    )

    # --- Pass 2: render every scene at the common AOI ---------------------
    scenes: list[Scene] = []
    for item in unique:
        sensing = item["ContentDate"]["Start"]
        name = item["Name"]
        platform = name[:3]
        sensor = "Sentinel-2" if platform.startswith("S2") else "Sentinel-1"
        date = sensing[:10]
        attrs = {a["Name"]: a.get("Value") for a in item.get("Attributes", [])}
        png_name = f"{date}_{platform}_{name.split('_')[5] if sensor == 'Sentinel-2' else name.split('_')[7]}.png"
        png_path = ASSETS / png_name

        print(f"[{sensing}] {platform} {sensor}")
        if png_path.exists() and png_path.stat().st_size > 0:
            print(f"  - cached {png_path.name}")
        else:
            try:
                if sensor == "Sentinel-2":
                    arr = render_s2(s3, item, cache, aoi)
                else:
                    arr = render_s1(s3, item, cache, aoi)
            except Exception as e:
                print(f"  ! failed: {e}")
                continue
            if arr is None:
                continue
            composite(arr, aoi, infra, png_path)
            print(f"  + wrote {png_path.name}")

        scenes.append(Scene(
            date=date,
            sensing_time=sensing,
            sensor=sensor,
            platform=platform,
            product_id=item["Id"],
            product_name=name,
            s3_path=item["S3Path"],
            orbit_direction=attrs.get("orbitDirection"),
            relative_orbit=attrs.get("relativeOrbitNumber"),
            png=f"assets/scenes/{png_name}",
        ))

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps([asdict(s) for s in scenes], indent=2))
    print(f"Wrote {MANIFEST.relative_to(ROOT)} with {len(scenes)} scenes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
