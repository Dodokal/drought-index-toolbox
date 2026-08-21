# -*- coding: utf-8 -*-
"""
DroughtIndexToolbox.pyt
=======================
SPI-based drought monitoring for eastern Ethiopia (or anywhere, really).

The idea: take a folder of monthly CHIRPS rainfall GeoTIFFs and turn them
into something a decision-maker can read --

    01  SPI Calculator     monthly rainfall  ->  SPI rasters (3/6/12-month)
    02  Drought Classifier SPI rasters       ->  severity classes (1-7)
    03  Shortage Hotspots  SPI rasters       ->  "% of months in drought" map

SPI follows the classic McKee et al. (1993) / Edwards & McKee (1997) recipe:
accumulate rainfall over k months, fit a gamma distribution per pixel and
per calendar month (Thom's approximation, so it vectorises cleanly), handle
zero-rainfall months with the mixed distribution, then map probabilities to
the standard normal.

Needs: ArcGIS Pro 3.x. numpy and scipy both ship with Pro, so nothing to
install. Rainfall folder is expected to hold one GeoTIFF per month with the
year and month somewhere in the name (CHIRPS names like
chirps-v2.0.1981.01.tif work out of the box).

Kalid Hassen -- Geospatial Python final project, summer 2026.
"""

import os
import re
import glob
import numpy as np
import arcpy
from scipy import special
from scipy.stats import norm

# CHIRPS ships -9999 as nodata; we also write it back out for float rasters.
NODATA = -9999.0

# Don't try to fit a gamma distribution on fewer than this many wet years --
# the parameters get silly. With CHIRPS (1981->now) you'll have 40+ anyway.
MIN_WET_YEARS = 5

# Standard SPI class boundaries (McKee). digitize() maps values to 1..7.
CLASS_BINS = [-2.0, -1.5, -1.0, 1.0, 1.5, 2.0]
CLASS_LEGEND = {
    1: "Extreme drought   (SPI <= -2.0)",
    2: "Severe drought    (-2.0 < SPI <= -1.5)",
    3: "Moderate drought  (-1.5 < SPI <= -1.0)",
    4: "Near normal       (-1.0 < SPI <  1.0)",
    5: "Moderately wet    ( 1.0 <= SPI < 1.5)",
    6: "Very wet          ( 1.5 <= SPI < 2.0)",
    7: "Extremely wet     (SPI >= 2.0)",
}


# ---------------------------------------------------------------------------
# The maths lives up here, arcpy-free, so it can be unit-tested outside Pro.
# ---------------------------------------------------------------------------

def rolling_sum(stack, k):
    """k-month rolling sum along the time axis.

    stack: (time, rows, cols). A window is only valid if all k months are
    present -- NaN propagates through np.sum, which is exactly what we want.
    """
    out = np.full(stack.shape, np.nan, dtype="float64")
    for i in range(k - 1, stack.shape[0]):
        out[i] = stack[i - k + 1: i + 1].sum(axis=0)
    return out


def gamma_params(month_values):
    """Fit gamma per pixel for one calendar month via Thom's approximation.

    month_values: (n_years, rows, cols) of k-month rainfall sums for a single
    calendar month across the record. Zeros are handled separately (mixed
    distribution), so the gamma is fitted to the wet years only.

    Returns (alpha, beta, q) where q is the probability of a zero month.
    """
    valid = ~np.isnan(month_values)
    n = valid.sum(axis=0).astype("float64")

    wet = np.where(month_values > 0, month_values, np.nan)
    n_wet = (~np.isnan(wet)).sum(axis=0).astype("float64")

    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(n > 0, 1.0 - n_wet / n, np.nan)

        mean = np.nanmean(wet, axis=0)
        mean_log = np.nanmean(np.log(wet), axis=0)

        # Thom (1958): A = ln(mean) - mean(ln x); then shape & scale.
        A = np.log(mean) - mean_log
        A = np.where(A > 0, A, np.nan)          # constant history -> no fit
        alpha = (1.0 + np.sqrt(1.0 + 4.0 * A / 3.0)) / (4.0 * A)
        beta = mean / alpha

    # Not enough wet years to trust the fit? Blank the parameters.
    too_dry = n_wet < MIN_WET_YEARS
    alpha[too_dry] = np.nan
    beta[too_dry] = np.nan
    return alpha, beta, q


def spi_values(sums, alpha, beta, q):
    """Rainfall sums -> SPI, given the fitted parameters for that month.

    sums can be (rows, cols) or (n_years, rows, cols); parameters broadcast.
    Mixed distribution: H = q + (1-q) * GammaCDF(x). At x = 0 the CDF term
    vanishes so H = q, which is the classic Edwards & McKee treatment.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        G = special.gammainc(alpha, sums / beta)
        H = q + (1.0 - q) * G
        # Clamp so ppf never returns +/-inf; +/-3.09 is the usual SPI cap.
        H = np.clip(H, 0.001, 0.999)
        spi = norm.ppf(H)
    spi[np.isnan(sums)] = np.nan
    return spi


def classify_spi(spi):
    """SPI -> integer severity classes 1..7 (0 = nodata)."""
    classes = np.digitize(spi, CLASS_BINS) + 1
    classes[np.isnan(spi)] = 0
    return classes.astype("uint8")


# ---------------------------------------------------------------------------
# Small arcpy conveniences shared by the tools.
# ---------------------------------------------------------------------------

MONTH_RE = re.compile(r"(\d{4})[.\-_](\d{2})")


def find_monthly_files(folder, pattern):
    """Find rasters and read (year, month) out of their names.

    Returns a list of (year, month, path) sorted in time. Raises if a file
    doesn't carry a parseable date -- better to fail loudly than mis-order
    a time series.
    """
    paths = sorted(glob.glob(os.path.join(folder, pattern)))
    if not paths:
        raise arcpy.ExecuteError(
            f"No files matching '{pattern}' in {folder}")
    dated = []
    for p in paths:
        m = MONTH_RE.search(os.path.basename(p))
        if not m:
            raise arcpy.ExecuteError(
                f"Can't find a YYYY.MM date in file name: {os.path.basename(p)}")
        dated.append((int(m.group(1)), int(m.group(2)), p))
    dated.sort()
    return dated


def load_stack(dated_files, messages):
    """Read the monthly rasters into one (time, rows, cols) array.

    Missing months in the sequence become all-NaN layers so that rolling
    windows stay aligned to the calendar -- silently skipping a month would
    shift every SPI value after it.
    """
    y0, m0, first_path = dated_files[0]
    y1, m1, _ = dated_files[-1]
    n_months = (y1 - y0) * 12 + (m1 - m0) + 1

    ref = arcpy.Raster(first_path)
    rows, cols = ref.height, ref.width
    georef = dict(
        lower_left=arcpy.Point(ref.extent.XMin, ref.extent.YMin),
        cell_w=ref.meanCellWidth, cell_h=ref.meanCellHeight,
        sr=arcpy.Describe(first_path).spatialReference,
    )

    stack = np.full((n_months, rows, cols), np.nan, dtype="float64")
    have = {(y, m): p for y, m, p in dated_files}

    idx_date = []
    missing = 0
    for i in range(n_months):
        y = y0 + (m0 - 1 + i) // 12
        m = (m0 - 1 + i) % 12 + 1
        idx_date.append((y, m))
        p = have.get((y, m))
        if p is None:
            missing += 1
            continue
        arr = arcpy.RasterToNumPyArray(p, nodata_to_value=NODATA).astype("float64")
        arr[arr == NODATA] = np.nan
        arr[arr < 0] = np.nan          # negative rainfall = bad pixels
        stack[i] = arr

    if missing:
        messages.addWarningMessage(
            f"{missing} month(s) missing from the sequence -- filled with NoData. "
            "SPI windows touching them will be NoData too.")
    messages.addMessage(
        f"Loaded {len(dated_files)} rasters ({y0}-{m0:02d} to {y1}-{m1:02d}), "
        f"grid {rows}x{cols}.")
    return stack, idx_date, georef


def save_array(arr, georef, out_path, nodata=NODATA, as_int=False):
    """Write an array back to disk with the stack's georeferencing."""
    if as_int:
        ras = arcpy.NumPyArrayToRaster(
            arr, georef["lower_left"], georef["cell_w"], georef["cell_h"],
            value_to_nodata=0)
    else:
        out = np.where(np.isnan(arr), nodata, arr).astype("float32")
        ras = arcpy.NumPyArrayToRaster(
            out, georef["lower_left"], georef["cell_w"], georef["cell_h"],
            value_to_nodata=nodata)
    arcpy.management.DefineProjection(ras, georef["sr"])
    ras.save(out_path)


# ---------------------------------------------------------------------------
# Toolbox
# ---------------------------------------------------------------------------

class Toolbox(object):
    def __init__(self):
        self.label = "Drought Index Toolbox"
        self.alias = "drought"
        self.description = ("SPI-based drought monitoring from monthly "
                            "rainfall rasters (CHIRPS-friendly).")
        self.tools = [SPICalculator, DroughtClassifier, ShortageHotspots]


class SPICalculator(object):
    def __init__(self):
        self.label = "01 SPI Calculator"
        self.description = (
            "Computes the Standardized Precipitation Index from a folder of "
            "monthly rainfall GeoTIFFs. Fits a gamma distribution per pixel "
            "and calendar month (Thom's approximation) with zero-rainfall "
            "handling, then writes one SPI raster per month: "
            "spi{k}_{YYYY}_{MM}.tif.")
        self.canRunInBackground = False

    def getParameterInfo(self):
        folder = arcpy.Parameter(
            displayName="Folder of monthly rainfall rasters",
            name="in_folder", datatype="DEFolder",
            parameterType="Required", direction="Input")

        pattern = arcpy.Parameter(
            displayName="File pattern", name="pattern",
            datatype="GPString", parameterType="Required", direction="Input")
        pattern.value = "*.tif"

        scale = arcpy.Parameter(
            displayName="Timescale (months)", name="scale",
            datatype="GPLong", parameterType="Required", direction="Input")
        scale.filter.type = "ValueList"
        scale.filter.list = [1, 3, 6, 12]
        scale.value = 3

        first_year = arcpy.Parameter(
            displayName="Only write outputs from this year on (optional)",
            name="first_year", datatype="GPLong",
            parameterType="Optional", direction="Input")

        out_folder = arcpy.Parameter(
            displayName="Output folder", name="out_folder",
            datatype="DEFolder", parameterType="Required", direction="Input")

        return [folder, pattern, scale, first_year, out_folder]

    def execute(self, parameters, messages):
        folder = parameters[0].valueAsText
        pattern = parameters[1].valueAsText
        k = int(parameters[2].value)
        first_year = parameters[3].value          # may be None
        out_folder = parameters[4].valueAsText

        dated = find_monthly_files(folder, pattern)
        stack, idx_date, georef = load_stack(dated, messages)

        messages.addMessage(f"Accumulating rainfall over {k} month(s)...")
        sums = rolling_sum(stack, k)

        # Fit per calendar month, then transform every year of that month.
        spi = np.full(sums.shape, np.nan, dtype="float64")
        months = np.array([m for (_, m) in idx_date])
        for month in range(1, 13):
            sel = np.where(months == month)[0]
            sel = sel[sel >= k - 1]               # windows that exist
            if sel.size == 0:
                continue
            month_sums = sums[sel]
            alpha, beta, q = gamma_params(month_sums)
            spi[sel] = spi_values(month_sums, alpha, beta, q)
            messages.addMessage(
                f"  month {month:02d}: fitted on {sel.size} years.")

        # Write out, one raster per month.
        written = 0
        for i, (y, m) in enumerate(idx_date):
            if i < k - 1:
                continue
            if first_year and y < first_year:
                continue
            if np.all(np.isnan(spi[i])):
                continue
            out_path = os.path.join(out_folder, f"spi{k}_{y}_{m:02d}.tif")
            save_array(spi[i], georef, out_path)
            written += 1

        messages.addMessage(f"Done -- wrote {written} SPI-{k} rasters to {out_folder}.")


class DroughtClassifier(object):
    def __init__(self):
        self.label = "02 Drought Classifier"
        self.description = (
            "Reclassifies SPI rasters into the seven standard severity "
            "classes (1 = extreme drought ... 7 = extremely wet). Takes a "
            "folder of SPI rasters (as written by tool 01) and mirrors it "
            "with classified uint8 rasters: cls_<name>.tif.")
        self.canRunInBackground = False

    def getParameterInfo(self):
        folder = arcpy.Parameter(
            displayName="Folder of SPI rasters", name="in_folder",
            datatype="DEFolder", parameterType="Required", direction="Input")
        pattern = arcpy.Parameter(
            displayName="File pattern", name="pattern",
            datatype="GPString", parameterType="Required", direction="Input")
        pattern.value = "spi*.tif"
        out_folder = arcpy.Parameter(
            displayName="Output folder", name="out_folder",
            datatype="DEFolder", parameterType="Required", direction="Input")
        return [folder, pattern, out_folder]

    def execute(self, parameters, messages):
        folder = parameters[0].valueAsText
        pattern = parameters[1].valueAsText
        out_folder = parameters[2].valueAsText

        paths = sorted(glob.glob(os.path.join(folder, pattern)))
        if not paths:
            raise arcpy.ExecuteError(f"No files matching '{pattern}' in {folder}")

        # georef from the first raster; they all share the same grid.
        ref = arcpy.Raster(paths[0])
        georef = dict(
            lower_left=arcpy.Point(ref.extent.XMin, ref.extent.YMin),
            cell_w=ref.meanCellWidth, cell_h=ref.meanCellHeight,
            sr=arcpy.Describe(paths[0]).spatialReference)

        messages.addMessage("Class legend:")
        for code, label in CLASS_LEGEND.items():
            messages.addMessage(f"  {code} = {label}")

        for p in paths:
            arr = arcpy.RasterToNumPyArray(p, nodata_to_value=NODATA).astype("float64")
            arr[arr == NODATA] = np.nan
            classes = classify_spi(arr)
            out_path = os.path.join(
                out_folder, "cls_" + os.path.basename(p))
            save_array(classes, georef, out_path, as_int=True)
        messages.addMessage(f"Classified {len(paths)} rasters into {out_folder}.")


class ShortageHotspots(object):
    def __init__(self):
        self.label = "03 Shortage Hotspots"
        self.description = (
            "Summarises a folder of SPI rasters into one 'shortage "
            "frequency' map: the percentage of months in which each pixel "
            "sat at or below a drought threshold (default SPI <= -1.0). "
            "High values = chronic shortage hotspots.")
        self.canRunInBackground = False

    def getParameterInfo(self):
        folder = arcpy.Parameter(
            displayName="Folder of SPI rasters", name="in_folder",
            datatype="DEFolder", parameterType="Required", direction="Input")
        pattern = arcpy.Parameter(
            displayName="File pattern", name="pattern",
            datatype="GPString", parameterType="Required", direction="Input")
        pattern.value = "spi*.tif"
        threshold = arcpy.Parameter(
            displayName="Drought threshold (SPI <=)", name="threshold",
            datatype="GPDouble", parameterType="Required", direction="Input")
        threshold.value = -1.0
        out_raster = arcpy.Parameter(
            displayName="Output frequency raster (%)", name="out_raster",
            datatype="DERasterDataset", parameterType="Required",
            direction="Output")
        return [folder, pattern, threshold, out_raster]

    def execute(self, parameters, messages):
        folder = parameters[0].valueAsText
        pattern = parameters[1].valueAsText
        thr = float(parameters[2].value)
        out_raster = parameters[3].valueAsText

        paths = sorted(glob.glob(os.path.join(folder, pattern)))
        if not paths:
            raise arcpy.ExecuteError(f"No files matching '{pattern}' in {folder}")

        ref = arcpy.Raster(paths[0])
        georef = dict(
            lower_left=arcpy.Point(ref.extent.XMin, ref.extent.YMin),
            cell_w=ref.meanCellWidth, cell_h=ref.meanCellHeight,
            sr=arcpy.Describe(paths[0]).spatialReference)

        drought = None
        valid = None
        for p in paths:
            arr = arcpy.RasterToNumPyArray(p, nodata_to_value=NODATA).astype("float64")
            arr[arr == NODATA] = np.nan
            if drought is None:
                drought = np.zeros(arr.shape, dtype="float64")
                valid = np.zeros(arr.shape, dtype="float64")
            drought += (arr <= thr) & ~np.isnan(arr)
            valid += ~np.isnan(arr)

        with np.errstate(invalid="ignore", divide="ignore"):
            freq = np.where(valid > 0, 100.0 * drought / valid, np.nan)

        save_array(freq, georef, out_raster)
        messages.addMessage(
            f"Shortage-frequency raster written from {len(paths)} months "
            f"(threshold SPI <= {thr}). Values are % of months in drought.")
