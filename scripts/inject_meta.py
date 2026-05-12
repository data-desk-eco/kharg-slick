"""Inject Open Graph + Twitter Card meta tags into the built notebook HTML."""
from pathlib import Path

URL = "https://research.datadesk.eco/kharg-slick/"
TITLE = "May 2026 Kharg oil slick"
DESCRIPTION = (
    "Tracking a suspected oil spill at Iran's main crude export terminal "
    "across Sentinel-1 and Sentinel-2 passes, 6–11 May 2026."
)
IMAGE = URL + "assets/scenes/2026-05-08_S1C_00F57D.png"

META = f"""<meta property="og:title" content="{TITLE}">
<meta property="og:description" content="{DESCRIPTION}">
<meta property="og:image" content="{IMAGE}">
<meta property="og:url" content="{URL}">
<meta property="og:type" content="article">
<meta name="twitter:card" content="summary_large_image">
"""

ANCHOR = '<meta name="viewport" content="width=device-width, initial-scale=1.0">'

path = Path("docs/.observable/dist/index.html")
html = path.read_text()
if "og:title" in html:
    raise SystemExit("meta tags already injected")
if ANCHOR not in html:
    raise SystemExit(f"anchor not found: {ANCHOR}")
path.write_text(html.replace(ANCHOR, ANCHOR + "\n" + META))
print(f"injected meta tags into {path}")
