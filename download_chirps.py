"""
download_chirps.py -- fetch CHIRPS v2.0 *official monthly* rainfall and clip
it to the Dire Dawa - Harar - Jijiga corridor.

This is the recommended route: you get the real monthly product with the
original file naming, which DroughtIndexToolbox.pyt parses out of the box.

Run it inside the ArcGIS Pro Python so arcpy can do the clipping, e.g.:

  "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" download_chirps.py

(Or open the Python window in Pro and run it from there.)

If a download fails halfway, just run the script again -- it skips anything
already on disk.

Kalid Hassen -- data prep for the drought toolbox.
"""

import os
import gzip
import shutil
import urllib.request

# ---------------- settings: edit these four lines ----------------
YEARS      = range(1981, 2026)                 # full record for a solid fit
RAW_DIR    = r"C:\drought\chirps_raw"          # full-Africa tifs land here
CLIP_DIR   = r"C:\drought\chirps"              # clipped tifs -> feed tool 01
BBOX       = "40.5 8.0 43.5 10.5"              # W S E N, deg (DD-Harar-Jijiga)
# -----------------------------------------------------------------

BASE = "https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_monthly/tifs"

os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(CLIP_DIR, exist_ok=True)

try:
    import arcpy
    HAVE_ARCPY = True
except ImportError:
    HAVE_ARCPY = False
    print("arcpy not found -- files will be downloaded but NOT clipped.")
    print("Clip them afterwards, e.g. with GDAL:")
    print("  gdal_translate -projwin 40.5 10.5 43.5 8.0 in.tif out.tif")

downloaded = clipped = skipped = failed = 0

for year in YEARS:
    for month in range(1, 13):
        name = f"chirps-v2.0.{year}.{month:02d}.tif"
        clip_path = os.path.join(CLIP_DIR, name)
        if os.path.exists(clip_path):
            skipped += 1
            continue

        raw_path = os.path.join(RAW_DIR, name)
        if not os.path.exists(raw_path):
            url = f"{BASE}/{name}.gz"
            gz_path = raw_path + ".gz"
            try:
                urllib.request.urlretrieve(url, gz_path)
            except Exception as e:
                # recent months may simply not exist yet -- that's fine
                print(f"  could not fetch {name}: {e}")
                failed += 1
                continue
            with gzip.open(gz_path, "rb") as f_in, open(raw_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_path)
            downloaded += 1

        if HAVE_ARCPY:
            arcpy.management.Clip(raw_path, BBOX, clip_path,
                                  nodata_value="-9999")
            os.remove(raw_path)        # keep only the small clipped tif
            clipped += 1

    print(f"{year} done  (downloaded {downloaded}, clipped {clipped}, "
          f"skipped {skipped}, unavailable {failed})")

print("\nAll finished.")
print(f"Point tool '01 SPI Calculator' at: {CLIP_DIR}")
print("Pattern: chirps*.tif   (the default *.tif also works)")
