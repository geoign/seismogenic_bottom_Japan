"""Attainable resolution of the seismogenic-bottom map.

For every node of a fine grid, find the smallest search radius at which the node reaches
ESS >= ess_min independent families and n >= nmin events (the same acceptance rule as the
fit). The result is a map of the data-limited resolution. Land nodes (GSHHS L1) are
tabulated separately because only those matter for rheological interpretation.

Usage: python 08_resolution_limits.py <declustered_npz> <out_npz> [--step 0.025]
       [--radii 2,3,4,5,7,10,15,20,30,40] [--ess-min 30] [--nmin 60] [--bbox 122,148,24,46]
       [--coast GSHHS_i_L1.shp]
"""
import argparse, time
import numpy as np
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out")
ap.add_argument("--step", type=float, default=0.025)
ap.add_argument("--radii", default="2,3,4,5,7,10,15,20,30,40")
ap.add_argument("--ess-min", type=float, default=30)
ap.add_argument("--nmin", type=int, default=60)
ap.add_argument("--bbox", default="122,148,24,46")
ap.add_argument("--coast", default="GSHHS_i_L1.shp",
                help="GSHHG shoreline shapefile (intermediate resolution, level 1) for the land mask")
args = ap.parse_args()
RADII = [float(r) for r in args.radii.split(",")]
x0, x1, y0, y1 = [float(v) for v in args.bbox.split(",")]
R = 6371.0

d = np.load(args.src)
lat, lon, w = d["lat"], d["lon"], d["weight"]
eh = np.hypot(d["err_lat_km"], d["err_lon_km"]); ez = d["err_dep"]
print("events", len(lat))
print("horizontal s.e. (km) quantiles 10/50/75/90/95:", np.round(np.percentile(eh, [10, 50, 75, 90, 95]), 2))
print("depth s.e. (km) quantiles 10/50/75/90/95:", np.round(np.percentile(ez, [10, 50, 75, 90, 95]), 2))
la, lo = np.deg2rad(lat), np.deg2rad(lon)
xyz = np.c_[R*np.cos(la)*np.cos(lo), R*np.cos(la)*np.sin(lo), R*np.sin(la)]
tree = cKDTree(xyz)

lons = np.round(np.arange(x0, x1 + 1e-9, args.step), 6); lats = np.round(np.arange(y1, y0 - 1e-9, -args.step), 6)
LA, LO = np.meshgrid(np.deg2rad(lats), np.deg2rad(lons), indexing="ij")
cxyz = np.c_[(R*np.cos(LA)*np.cos(LO)).ravel(), (R*np.cos(LA)*np.sin(LO)).ravel(), (R*np.sin(LA)).ravel()]
N = len(cxyz)
rmin = np.full(N, np.nan, dtype=np.float32); ess_at = np.full(N, np.nan, dtype=np.float32); n_at = np.zeros(N, dtype=np.int32)
pending = np.arange(N)
t0 = time.time()
for r in RADII:
    for c0 in range(0, len(pending), 20000):
        ids = pending[c0:c0 + 20000]
        lists = tree.query_ball_point(cxyz[ids], r)
        lens = np.array([len(l) for l in lists], dtype=np.int64)
        s1 = np.zeros(len(ids)); s2 = np.zeros(len(ids))
        if lens.sum() > 0:
            cat = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists if len(l)])
            ww = w[cat]
            starts = np.r_[0, np.cumsum(lens)[:-1]]
            nz = lens > 0
            s1[nz] = np.add.reduceat(ww, starts[nz]); s2[nz] = np.add.reduceat(ww * ww, starts[nz])
        ess = np.where(s2 > 0, s1 * s1 / np.where(s2 > 0, s2, 1), 0.0)
        ok = (lens >= args.nmin) & (ess >= args.ess_min)
        rmin[ids[ok]] = r; ess_at[ids[ok]] = ess[ok]; n_at[ids[ok]] = lens[ok]
    pending = pending[np.isnan(rmin[pending])]
    print(f"r = {r:4.0f} km: nodes reaching ESS >= {args.ess_min:.0f} so far {int(np.isfinite(rmin).sum()):7d}, pending {len(pending):7d}  ({time.time()-t0:.0f} s)", flush=True)

# land mask from GSHHS (rasterise polygons at the node positions)
land = np.zeros(N, dtype=bool)
try:
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions()
    ds = ogr.Open(args.coast); lyr = ds.GetLayer(0)
    lyr.SetSpatialFilterRect(x0, y0, x1, y1)
    mem = gdal.GetDriverByName("MEM").Create("", len(lons), len(lats), 1, gdal.GDT_Byte)
    mem.SetGeoTransform((float(lons[0]) - args.step / 2, args.step, 0, float(lats[0]) + args.step / 2, 0, -args.step))
    srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); mem.SetProjection(srs.ExportToWkt())
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1], options=["ALL_TOUCHED=FALSE"])
    land = mem.ReadAsArray().astype(bool).ravel()
    print("land nodes:", int(land.sum()))
except Exception as e:
    print("land mask unavailable:", e)

rmin_g = rmin.reshape(len(lats), len(lons)); ess_g = ess_at.reshape(len(lats), len(lons)); n_g = n_at.reshape(len(lats), len(lons))
land_g = land.reshape(len(lats), len(lons))
np.savez(args.out, lats=lats, lons=lons, rmin=rmin_g, ess=ess_g, n=n_g, land=land_g, radii=np.array(RADII))
print("saved", args.out)

print("\nfraction of nodes whose data-limited resolution is <= r (cumulative):")
print(f"{'r (km)':>8} {'all nodes':>10} {'land nodes':>11} {'land, Honshu-Kyushu-Shikoku (130-142E, 31-41.5N)':>50}")
core = land_g & (LO >= np.deg2rad(130)) & (LO <= np.deg2rad(142)) & (LA >= np.deg2rad(31)) & (LA <= np.deg2rad(41.5))
anyv = np.isfinite(rmin_g)
for r in RADII:
    m = np.nan_to_num(rmin_g, nan=1e9) <= r
    print(f"{r:8.0f} {m[anyv].mean()*100 if anyv.any() else 0:9.1f}% {m[land_g].mean()*100 if land_g.any() else 0:10.1f}% {m[core].mean()*100 if core.any() else 0:49.1f}%")
print(f"{'none':>8} {(~anyv).mean()*100:9.1f}% {(~anyv[land_g]).mean()*100 if land_g.any() else 0:10.1f}% {(~anyv[core]).mean()*100 if core.any() else 0:49.1f}%")
