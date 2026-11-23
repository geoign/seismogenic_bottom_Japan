"""Write the node-wise fit results to a multi-band ZSTD GeoTIFF (EPSG:4326, pixel-is-point).

Usage: python 06_write_geotiff.py <fit_npz> <out_tif> [--mask-ess 30] [--mask-se 5]
"""
import argparse
import numpy as np
from osgeo import gdal, osr

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out")
ap.add_argument("--mask-ess", type=float, default=30.0)
ap.add_argument("--mask-se", type=float, default=5.0, help="mask nodes whose s.e.(z_b) exceeds this (km)")
args = ap.parse_args()
gdal.UseExceptions()

d = np.load(args.src)
lats, lons = d["lats"], d["lons"]
step = float(abs(lons[1] - lons[0]))
bad = ~np.isfinite(d["zb"]) | (d["ess"] < args.mask_ess) | (d["zb"] < 0) | (d["zb"] > 55)
if np.isfinite(args.mask_se):
    bad |= np.isfinite(d["se_zb"]) & (d["se_zb"] > args.mask_se)
nan = np.full_like(d["zb"], np.nan)
slab_cut = d["slab_cut"] if "slab_cut" in d else nan
slab_top = d["slab_top"] if "slab_top" in d else nan
frac30 = d["frac30"] if "frac30" in d else nan

Db10 = d["zb"] + 1.2816 * d["sb"]
# quality flag (bit mask): 1 = bottom within 3 km of the slab cut-off (limited by the slab mask, not by rheology),
# 2 = more than 30 % of the weight deeper than 30 km (lower-crustal or residual slab population),
# 4 = transition width s_b > 6 km (gradual decay, bottom ill-defined), 8 = ESS < 50
qflag = np.zeros_like(d["zb"], dtype=np.float32)
qflag += 1 * (np.isfinite(slab_cut) & ((slab_cut - Db10) < 3.0))
qflag += 2 * (np.nan_to_num(frac30) > 0.3)
qflag += 4 * (d["sb"] > 6.0)
qflag += 8 * (d["ess"] < 50)
Db10_clean = np.where(qflag == 0, Db10, np.nan)

bands = [
    ("D_b10", Db10, "Base of brittle seismicity: depth where the modelled brittle rate falls to 10 % (z_b + 1.28 s_b), km"),
    ("D_b10_clean", Db10_clean, "D_b10 restricted to nodes with quality_flag == 0"),
    ("quality_flag", qflag, "Bit mask: 1 = within 3 km of slab cut-off, 2 = >30 % of weight deeper than 30 km, 4 = s_b > 6 km, 8 = ESS < 50"),
    ("z_b", d["zb"], "Centre of the lower (brittle-ductile) transition of seismicity, km"),
    ("s_b", d["sb"], "Width (1 sigma) of the lower transition, km; includes depth-location scatter"),
    ("D_b90", d["zb"] - 1.2816 * d["sb"], "Top of the lower transition (z_b - 1.28 s_b), km"),
    ("se_z_b", d["se_zb"], "Standard error of z_b from the observed information, km"),
    ("z_t", d["zt"], "Centre of the shallow onset of seismicity, km"),
    ("s_t", d["st"], "Width of the shallow onset, km"),
    ("p_outlier", d["p"], "Fraction of the uniform outlier component (0-60 km)"),
    ("D90_weighted", d["d90"], "Classic D90 of the same sample with family weights, km (for comparison)"),
    ("D95_weighted", d["d95"], "Classic D95 of the same sample with family weights, km"),
    ("D10_weighted", d["d10"], "Classic D10 of the same sample with family weights, km"),
    ("ESS", d["ess"], "Kish effective sample size (number of independent families equivalent)"),
    ("n_events", d["n"], "Number of earthquakes used at the node"),
    ("radius_km", d["radius"], "Search radius used (adaptive), km"),
    ("slab_cut_km", slab_cut, "Median depth above which events were accepted (slab top minus margin, floor 15 km); "
                              "if D_b10 is within a few km of this, the result is limited by the slab mask, not by rheology"),
    ("slab_top_km", slab_top, "Median slab-top depth under the node (shallower of Slab2 and Hirose/Nakajima models), km"),
    ("frac_below_30km", frac30, "Weight fraction of events deeper than 30 km (QA: lower-crustal or residual slab seismicity)"),
]
drv = gdal.GetDriverByName("GTiff")
ds = drv.Create(args.out, len(lons), len(lats), len(bands), gdal.GDT_Float32,
                options=["COMPRESS=ZSTD", "ZSTD_LEVEL=15", "PREDICTOR=3", "TILED=YES",
                         "BLOCKXSIZE=256", "BLOCKYSIZE=256"])
ds.SetGeoTransform((float(lons[0]) - step / 2, step, 0.0, float(lats[0]) + step / 2, 0.0, -step))
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); ds.SetProjection(srs.ExportToWkt())
ds.SetMetadata({
    "AREA_OR_POINT": "Point",
    "TITLE": "Lower limit of brittle seismicity beneath the Japanese islands (JMA unified catalogue 2001-2024)",
    "METHOD": "Weighted maximum-likelihood fit of f(z) = (1-p) Phi((z-z_t)/s_t)[1-Phi((z-z_b)/s_b)]/norm + p/60 "
              "to upper-plate hypocentres within an adaptive radius; weights = 1/family size from "
              "nearest-neighbour declustering (Zaliapin & Ben-Zion 2013)",
    "SOURCE": "JMA unified hypocenter catalogue (depth-free, ordinary events, flags K/k/A), Slab2 (Hayes et al. 2018) "
              "and Hirose/Nakajima plate models for the slab mask",
    "UNITS": "km (depth positive down), except p_outlier (fraction), ESS, n_events and radius_km",
})
for i, (name, arr, desc) in enumerate(bands, start=1):
    a = np.array(arr, dtype=np.float32).copy()
    a[bad] = np.nan
    b = ds.GetRasterBand(i)
    b.SetDescription(name); b.SetNoDataValue(float("nan")); b.SetMetadataItem("DESCRIPTION", desc)
    b.WriteArray(a); b.ComputeStatistics(False)
ds.FlushCache(); ds = None
print("written", args.out, "valid nodes:", int((~bad).sum()))
