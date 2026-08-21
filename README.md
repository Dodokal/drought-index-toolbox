# Drought Index Toolbox

An ArcGIS Python toolbox that turns CHIRPS rainfall rasters into SPI drought maps. I built it for eastern Ethiopia, but it works anywhere you have monthly rainfall.

The problem it solves: raw rainfall numbers don't tell you much on their own. 40 mm in a month might be perfectly normal in one place and a disaster 50 km away. SPI fixes that by scoring every pixel against its own history instead of against some fixed threshold.

## What's in here

**`DroughtIndexToolbox.pyt`** — three tools:

- **01 SPI Calculator** — monthly rainfall in, SPI rasters out (1, 3, 6 or 12-month)
- **02 Drought Classifier** — SPI into the seven standard severity classes
- **03 Shortage Hotspots** — how often each pixel has sat in drought across the whole record

**`download_chirps.py`** — pulls CHIRPS monthly rainfall and clips it to your study box, so you're not sitting on 500 Africa-wide GeoTIFFs.

## Getting the data

```
python download_chirps.py
```

Edit the folders and the bounding box at the top first. It skips files it already has, so if it dies halfway through 1994 just run it again.

One thing that will probably bite you: the CHIRPS server has served an expired certificate more than once, and Python refuses to download anything when that happens. The script tries a proper certificate bundle first and falls back to an unverified connection if that fails. It prints which route it took. That fallback is fine here because it's public read-only data from a known URL, but don't copy that pattern for anything with a login.

## Running the toolbox

In ArcGIS Pro: Catalog → Toolboxes → Add Toolbox → pick the `.pyt`. Point tool 01 at your rainfall folder, choose a timescale, run. Tools 02 and 03 read whatever 01 wrote.

Needs ArcGIS Pro 3.x with Spatial Analyst. numpy and scipy already ship with Pro, so there's nothing to install.

## The thing I wish I'd understood earlier

The hotspot map does not show you where it's driest.

SPI is relative to each pixel's own climatology. A permanently arid place isn't "in drought" for being dry, only for being dry *by its own standards*. So what you actually get is a map of where rainfall is most erratic. My hotspot values landed between 12.5% and 19.5%, clustered around 16%, which is exactly what the maths predicts (a standard normal puts 15.87% of months below −1). Narrow spread, and not a ranking of thirst.

If you want absolute water scarcity, pair this with an aridity index. I dropped mine from the project, the PET data I wanted had no clean download, sat on a different grid, and needed an extra annual-summing step. Three problems for an optional tool, and I ran out of time.

## Other things worth knowing

Don't give it fewer than about 30 years. The gamma distribution is fitted per pixel *per calendar month*, so a 2-year record gives you two samples per month, and you're fitting noise. Pixels with under five wet years for a month come back as NoData on purpose.

The whole raster stack loads into memory. Fine for a clipped study area. Not fine for a continent.

Zero-rainfall months are handled with the mixed distribution, which matters more than it sounds, dry-season months out toward Jijiga are frequently exactly zero, and a naive gamma fit falls over on those.


## Credit

Standard SPI recipe: McKee et al. (1993), gamma fitted with Thom's (1958) approximation, zero-rainfall handling from Edwards & McKee (1997). Rainfall from CHIRPS v2.0 (Funk et al. 2015).

MIT licensed. Use it, break it, tell me what broke.
