"""Quick-look maps of the fitted seismogenic bottom.

Usage: python 05_plot_maps.py <fit_npz> <out_png_prefix> [--mask-ess 30] [--coast GSHHS_i_L1.shp]
       [--compare other_fit_npz]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("prefix")
ap.add_argument("--mask-ess", type=float, default=30.0)
ap.add_argument("--coast", default="GSHHS_i_L1.shp",
                help="GSHHG shoreline shapefile (intermediate resolution, level 1); maps are drawn without coast if missing")
ap.add_argument("--compare", default=None, help="second fit (e.g. without weights) for a difference panel")
args = ap.parse_args()
BBOX = (122, 148, 24, 46)

def coast_segments(path):
    """Coastline polygons as a list of (N,2) arrays, clipped to the map box (via OGR)."""
    try:
        from osgeo import ogr
        ds = ogr.Open(path); lyr = ds.GetLayer(0)
        lyr.SetSpatialFilterRect(BBOX[0], BBOX[2], BBOX[1], BBOX[3])
        segs = []
        for f in lyr:
            g = f.GetGeometryRef()
            polys = [g] if g.GetGeometryName() == "POLYGON" else [g.GetGeometryRef(i) for i in range(g.GetGeometryCount())]
            for pg in polys:
                for r in range(pg.GetGeometryCount()):
                    ring = pg.GetGeometryRef(r)
                    pts = np.array(ring.GetPoints())[:, :2]
                    if len(pts) > 1: segs.append(pts)
        return segs
    except Exception as e:
        print("coastline unavailable:", e); return []

COAST = coast_segments(args.coast)
d = np.load(args.src)
lats, lons = d["lats"], d["lons"]
step = float(abs(lons[1] - lons[0]))
ext = [lons[0] - step / 2, lons[-1] + step / 2, lats[-1] - step / 2, lats[0] + step / 2]
bad = ~np.isfinite(d["zb"]) | (d["ess"] < args.mask_ess) | (d["zb"] < 0) | (d["zb"] > 55)
if "se_zb" in d: bad |= np.isfinite(d["se_zb"]) & (d["se_zb"] > 5)

def axes(fig, pos):
    ax = fig.add_subplot(*pos) if isinstance(pos, tuple) else fig.add_subplot(pos)
    ax.add_collection(LineCollection(COAST, colors="k", linewidths=0.35))
    ax.set_xlim(BBOX[0], BBOX[1]); ax.set_ylim(BBOX[2], BBOX[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(36)))
    ax.grid(True, linewidth=0.3, alpha=0.4)
    return ax

def panel(ax, a, title, vmin, vmax, cmap):
    a = np.array(a, dtype=float); a[bad] = np.nan
    im = ax.imshow(a, extent=ext, origin="upper", vmin=vmin, vmax=vmax, cmap=cmap, interpolation="nearest", zorder=0)
    ax.set_title(title, fontsize=11)
    plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)

Db10 = d["zb"] + 1.2816 * d["sb"]
panels = [
    (Db10, "D_b10 = z_b + 1.28 s_b: base of brittle seismicity (km)", 0, 35, "viridis_r"),
    (d["zb"], "z_b: centre of the lower transition (km)", 0, 35, "viridis_r"),
    (d["sb"], "s_b: transition width, 1 sigma (km)", 0, 8, "magma_r"),
    (d["d90"], "classic D90 of the same weighted sample (km)", 0, 35, "viridis_r"),
    (Db10 - d["d90"], "D_b10 minus classic D90 (km)", -6, 6, "RdBu_r"),
    (d["ess"], "effective sample size (independent families)", 0, 500, "cividis"),
]
if "frac30" in d:
    panels.append((d["frac30"], "weight fraction deeper than 30 km (slab / lower-crust QA)", 0, 0.5, "Reds"))
if "slab_cut" in d:
    panels.append((d["slab_cut"] - Db10, "slab cut-off depth minus D_b10 (km; small = limited by slab mask)", 0, 30, "Greens"))
if args.compare:
    c = np.load(args.compare)
    panels.append(((c["zb"] + 1.2816 * c["sb"]) - Db10, "D_b10 without family weights minus with weights (km)", -6, 6, "RdBu_r"))
ncol = 3; nrow = int(np.ceil(len(panels) / ncol))
fig = plt.figure(figsize=(6.2 * ncol, 5.8 * nrow))
for i, (a, title, vmin, vmax, cmap) in enumerate(panels, start=1):
    panel(axes(fig, (nrow, ncol, i)), a, title, vmin, vmax, cmap)
fig.suptitle(f"Lower limit of brittle seismicity based on JMA unified catalogue 2001-2024 (upper plate, family-weighted ML fit, {step:g}-degree nodes)", fontsize=13)
fig.tight_layout()
fig.savefig(args.prefix + "_overview.png", dpi=130)
print("saved", args.prefix + "_overview.png")

limited = np.zeros_like(bad)
if "slab_cut" in d:
    limited = np.isfinite(d["slab_cut"]) & ((d["slab_cut"] - Db10) < 3.0)
if "frac30" in d:
    limited |= np.nan_to_num(d["frac30"]) > 0.3
for key, arr, lab, vmax in [("Db10", Db10, "D_b10 (km): depth where the modelled brittle rate falls to 10 %", 35),
                            ("zb", d["zb"], "z_b (km): centre of the lower transition of seismicity", 35)]:
    fig = plt.figure(figsize=(11, 10))
    ax = axes(fig, 111)
    a = np.array(arr, dtype=float); a[bad] = np.nan
    im = ax.imshow(a, extent=ext, origin="upper", vmin=0, vmax=vmax, cmap="viridis_r", interpolation="nearest", zorder=0)
    g = np.where(limited & ~bad, 1.0, np.nan)
    ax.imshow(g, extent=ext, origin="upper", cmap=matplotlib.colors.ListedColormap(["#9a9a9a"]), interpolation="nearest",
              zorder=1, alpha=0.85)
    plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02, label=lab)
    ax.set_title(lab + f", {step:g}-degree nodes, 2001-2024\n"
                 "grey: bottom set by the slab mask or >30 % of events deeper than 30 km (not a rheological limit)", fontsize=11)
    fig.tight_layout(); fig.savefig(args.prefix + f"_{key}.png", dpi=150)
    print("saved", args.prefix + f"_{key}.png")
