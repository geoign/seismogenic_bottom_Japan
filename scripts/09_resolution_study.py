"""Resolution study: how far can the node spacing and the search radius be pushed?

Usage: python 09_resolution_study.py <work_dir> <out_dir> [--coast GSHHS_i_L1.shp]

Reads from <work_dir> (whatever is present):
  resolution_limits_0p025.npz                 data-limited minimum radius per node (08)
  fit_v3_0p025.npz, fit_v3_0p1.npz            fine and coarse products (04)
  fit_full_r{R}_0p1.npz, fit_half{0,1}_r{R}_0p1.npz   fixed-radius fits, all data / two family halves
  fit_half{0,1}_0p025_{region}.npz            split-half fine maps for test regions
Writes to <out_dir>: resolution_limit_rmin.png, splithalf_{region}.png, resolution_study.md
"""
import argparse, os, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
from scipy.ndimage import uniform_filter

ap = argparse.ArgumentParser()
ap.add_argument("work"); ap.add_argument("out")
ap.add_argument("--coast", default="GSHHS_i_L1.shp",
                help="GSHHG shoreline shapefile (intermediate resolution, level 1); maps are drawn without coast if missing")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
BBOX = (122, 148, 24, 46)
Z = 1.2816
report = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    report.append(s)


def coast_segments(path, bbox):
    try:
        from osgeo import ogr
        ds = ogr.Open(path); lyr = ds.GetLayer(0)
        lyr.SetSpatialFilterRect(bbox[0], bbox[2], bbox[1], bbox[3])
        segs = []
        for f in lyr:
            g = f.GetGeometryRef()
            polys = [g] if g.GetGeometryName() == "POLYGON" else [g.GetGeometryRef(i) for i in range(g.GetGeometryCount())]
            for pg in polys:
                for r in range(pg.GetGeometryCount()):
                    pts = np.array(pg.GetGeometryRef(r).GetPoints())[:, :2]
                    if len(pts) > 1:
                        segs.append(pts)
        return segs
    except Exception as e:
        print("coastline unavailable:", e)
        return []


def valid(d, ess_min):
    bad = ~np.isfinite(d["zb"]) | (d["ess"] < ess_min) | (d["zb"] < 0) | (d["zb"] > 55)
    bad |= np.isfinite(d["se_zb"]) & (d["se_zb"] > 5)
    return ~bad


def db10(d):
    return d["zb"] + Z * d["sb"]


def extent(d):
    lats, lons = d["lats"], d["lons"]
    step = float(abs(lons[1] - lons[0]))
    return [lons[0] - step / 2, lons[-1] + step / 2, lats[-1] - step / 2, lats[0] + step / 2]


def map_axes(fig, pos, bbox, segs):
    ax = fig.add_subplot(*pos)
    ax.add_collection(LineCollection(segs, colors="k", linewidths=0.4))
    ax.set_xlim(bbox[0], bbox[1]); ax.set_ylim(bbox[2], bbox[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(0.5 * (bbox[2] + bbox[3]))))
    ax.grid(True, linewidth=0.3, alpha=0.4)
    return ax


def radius_of(path):
    return float(re.search(r"_r([0-9.]+)_", os.path.basename(path)).group(1))


def anom(x, k):
    """Deviation of a map from its k-node running mean (NaN-aware): the fine-scale texture."""
    f = np.where(np.isfinite(x), x, 0.0); c = np.isfinite(x).astype(float)
    sm = uniform_filter(f, k) / np.maximum(uniform_filter(c, k), 1e-9)
    return x - sm


W = args.work
say("# Resolution study of the seismogenic-bottom map\n")

# ---------------------------------------------------------------- A. location errors and data density
p = os.path.join(W, "resolution_limits_0p025.npz")
if os.path.exists(p):
    r = np.load(p)
    rmin, land, radii = r["rmin"], r["land"], r["radii"]
    LO, LA = np.meshgrid(r["lons"], r["lats"])
    core = land & (LO >= 130) & (LO <= 142) & (LA >= 31) & (LA <= 41.5)
    anyv = np.isfinite(rmin)
    say("## 1. Data-limited resolution (smallest radius with ESS >= 30 families and n >= 60 events)\n")
    say("Cumulative fraction of 0.025-degree nodes whose smallest admissible radius is <= r:\n")
    say("| r (km) | land nodes of Honshu, Kyushu, Shikoku (130-142E, 31-41.5N) | all nodes that reach 40 km |")
    say("| --- | --- | --- |")
    for rr in radii:
        m = np.nan_to_num(rmin, nan=1e9) <= rr
        say(f"| {rr:.0f} | {100*m[core].mean():.1f} % | {100*m[anyv].mean():.1f} % |")
    say("")
    segs = coast_segments(args.coast, BBOX)
    fig = plt.figure(figsize=(11, 10))
    ax = map_axes(fig, (1, 1, 1), BBOX, segs)
    levels = list(radii) + [60]
    cmap = ListedColormap(plt.get_cmap("viridis")(np.linspace(0, 1, len(radii))))
    norm = BoundaryNorm(levels, cmap.N)
    im = ax.imshow(rmin, extent=extent(r), origin="upper", cmap=cmap, norm=norm, interpolation="nearest", zorder=0)
    cb = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02, ticks=radii)
    cb.set_label("smallest search radius reaching ESS >= 30 (km)")
    ax.set_title("Data-limited resolution of the seismogenic-bottom map (JMA 2001-2024, upper plate, declustered)\n"
                 "colour = smallest radius at which a 0.025-degree node has 30 independent families; "
                 "blank = not reached within 40 km", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "resolution_limit_rmin.png"), dpi=150); plt.close(fig)

# ---------------------------------------------------------------- B. split-half noise per radius (0.1-degree nodes)
rows = []
for f in sorted(glob.glob(os.path.join(W, "fit_full_r*_0p1.npz")), key=radius_of):
    rr = radius_of(f)
    h0p = os.path.join(W, f"fit_half0_r{rr:g}_0p1.npz")
    h1p = os.path.join(W, f"fit_half1_r{rr:g}_0p1.npz")
    if not (os.path.exists(h0p) and os.path.exists(h1p)):
        continue
    full, h0, h1 = np.load(f), np.load(h0p), np.load(h1p)
    m = valid(full, 30) & valid(h0, 15) & valid(h1, 15)
    if m.sum() < 20:
        continue
    a, b, c = db10(h0)[m], db10(h1)[m], db10(full)[m]
    diff = a - b
    rms = float(np.sqrt(np.mean(diff ** 2)))
    mad = float(1.4826 * np.median(np.abs(diff - np.median(diff))))
    noise_full = rms / 2          # var(A - B) = 2 var(half); var(full) = var(half) / 2
    noise_full_rob = mad / 2
    corr = float(np.corrcoef(a, b)[0, 1])
    se_h = float(np.nanmedian(full["se_zb"][m]))
    # fine-scale texture: deviation from a 0.5-degree (5-node) running mean
    A2 = np.where(valid(h0, 15), db10(h0), np.nan); B2 = np.where(valid(h1, 15), db10(h1), np.nan)
    F2 = np.where(valid(full, 30), db10(full), np.nan)
    aa, bb, ff = anom(A2, 5)[m], anom(B2, 5)[m], anom(F2, 5)[m]
    corr_an = float(np.corrcoef(aa, bb)[0, 1])
    rms_an = float(np.sqrt(np.mean(ff ** 2)))
    sig_an = float(np.sqrt(max(rms_an ** 2 - noise_full ** 2, 0)))
    rows.append((rr, int(m.sum()), rms, noise_full, noise_full_rob, se_h, corr, rms_an, corr_an, sig_an))
if rows:
    say("## 2. Reproducibility between two independent halves of the catalogue (fixed radius, 0.1-degree nodes)\n")
    say("Families (declustered sequences) were split by parity of their id into two halves, and the fit was repeated on "
        "each half with ESS >= 15. The noise of the full-data estimate is RMS(half A - half B) / 2. The fine-scale "
        "texture is the deviation of D_b10 from its 0.5-degree running mean; its reproducibility is the correlation of "
        "that texture between the two halves, and its signal is its RMS in the full-data map after removing the noise.\n")
    say("| radius (km) | nodes | RMS(A-B) km | noise of D_b10 (RMS) km | noise (robust, MAD) km | "
        "median s.e.(z_b) from Hessian km | corr(A, B) | RMS texture of D_b10 km | corr of texture (A, B) | "
        "texture signal / noise |")
    say("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for rr, n, rms, nf, nfr, se_h, corr, rms_an, corr_an, sig_an in rows:
        say(f"| {rr:g} | {n} | {rms:.2f} | {nf:.2f} | {nfr:.2f} | {se_h:.2f} | {corr:.3f} | {rms_an:.2f} | "
            f"{corr_an:.3f} | {sig_an / max(nf, 1e-6):.1f} |")
    say("")

# ---------------------------------------------------------------- C. fine product vs coarse product
pf, pc = os.path.join(W, "fit_v3_0p025.npz"), os.path.join(W, "fit_v3_0p1.npz")
if os.path.exists(pf) and os.path.exists(pc):
    F, C = np.load(pf), np.load(pc)
    vf = valid(F, 30); vc = valid(C, 30)
    say("## 3. The 0.025-degree product\n")
    say(f"Valid nodes: {int(vf.sum())} (0.025 deg) vs {int(vc.sum())} (0.1 deg).\n")
    rad = F["radius"]
    say("Search radius actually used by the valid 0.025-degree nodes:\n")
    say("| radius (km) | nodes | share | median s.e.(z_b) km | median ESS |")
    say("| --- | --- | --- | --- | --- |")
    for rr in np.unique(rad[vf]):
        m = vf & (rad == rr)
        say(f"| {rr:g} | {int(m.sum())} | {100 * m.sum() / vf.sum():.1f} % | {np.nanmedian(F['se_zb'][m]):.2f} | "
            f"{np.nanmedian(F['ess'][m]):.0f} |")
    say("")
    fl, fo = F["lats"], F["lons"]; cl, co = C["lats"], C["lons"]
    iy = np.array([np.argmin(np.abs(fl - y)) for y in cl]); ix = np.array([np.argmin(np.abs(fo - x)) for x in co])
    okl = np.abs(fl[iy] - cl) < 1e-6; okx = np.abs(fo[ix] - co) < 1e-6
    IY, IX = np.meshgrid(iy, ix, indexing="ij"); OK = np.outer(okl, okx)
    fsub = db10(F)[IY, IX]; vsub = vf[IY, IX] & OK; rsub = rad[IY, IX]
    both = vsub & vc
    say(f"Of the {int(both.sum())} coincident nodes valid in both products, {int((both & (rsub >= 10)).sum())} use a "
        f"radius >= 10 km in the fine product and are therefore identical to the coarse product. Where the fine product "
        f"uses a smaller radius, D_b10(0.025 deg) - D_b10(0.1 deg) is:\n")
    say("| fine radius (km) | nodes | median difference km | robust s.d. km | abs. difference > 3 km |")
    say("| --- | --- | --- | --- | --- |")
    for rr in (3, 5, 7):
        sel = both & (rsub == rr)
        if sel.sum() == 0:
            continue
        dd = (fsub - db10(C))[sel]
        say(f"| {rr} | {int(sel.sum())} | {np.median(dd):+.2f} | {1.4826 * np.median(np.abs(dd - np.median(dd))):.2f} | "
            f"{100 * np.mean(np.abs(dd) > 3):.1f} % |")
    say("")

# ---------------------------------------------------------------- D. split-half maps of the fine product in test regions
for reg, bbox in [("kyushu", (129.5, 132.5, 31, 34)), ("tohoku_kanto", (138, 142, 35, 40))]:
    h0p = os.path.join(W, f"fit_half0_0p025_{reg}.npz")
    h1p = os.path.join(W, f"fit_half1_0p025_{reg}.npz")
    if not (os.path.exists(h0p) and os.path.exists(h1p) and os.path.exists(pf)):
        continue
    h0, h1, F = np.load(h0p), np.load(h1p), np.load(pf)
    fl, fo = F["lats"], F["lons"]
    iy = np.array([np.argmin(np.abs(fl - y)) for y in h0["lats"]]); ix = np.array([np.argmin(np.abs(fo - x)) for x in h0["lons"]])
    IY, IX = np.meshgrid(iy, ix, indexing="ij")
    full_c = np.where(valid(F, 30)[IY, IX], db10(F)[IY, IX], np.nan)
    a = np.where(valid(h0, 15), db10(h0), np.nan)
    b = np.where(valid(h1, 15), db10(h1), np.nan)
    m = np.isfinite(a) & np.isfinite(b) & np.isfinite(full_c)
    diff = (a - b)[m]
    rms = float(np.sqrt(np.mean(diff ** 2)))
    mad = float(1.4826 * np.median(np.abs(diff - np.median(diff))))
    corr = float(np.corrcoef(a[m], b[m])[0, 1])

    aa, bb = anom(a, 11)[m], anom(b, 11)[m]
    corr_an = float(np.corrcoef(aa, bb)[0, 1])
    say(f"## 4. Split-half test of the 0.025-degree product: {reg} ({bbox[0]}-{bbox[1]}E, {bbox[2]}-{bbox[3]}N)\n")
    say(f"{int(m.sum())} nodes valid in both halves. RMS(A-B) = {rms:.2f} km, robust = {mad:.2f} km, so the noise of the "
        f"full-data D_b10 is about {mad / 2:.2f} km (robust). corr(A, B) = {corr:.3f}; correlation of the anomalies from a "
        f"0.275-degree (11-node) running mean = {corr_an:.3f} (how much of the fine texture is reproduced by independent data).\n")
    segs = coast_segments(args.coast, bbox); ext = extent(h0)
    h_in = 3.6 * (bbox[3] - bbox[2]) / (bbox[1] - bbox[0]) / np.cos(np.deg2rad(0.5 * (bbox[2] + bbox[3]))) + 1.3
    fig = plt.figure(figsize=(19, h_in))
    for i, (arr, title, vmin, vmax, cmap) in enumerate([
            (full_c, "all data (product)", 0, 35, "viridis_r"), (a, "half A (even family ids)", 0, 35, "viridis_r"),
            (b, "half B (odd family ids)", 0, 35, "viridis_r"), (a - b, "half A minus half B", -8, 8, "RdBu_r")], start=1):
        ax = map_axes(fig, (1, 4, i), bbox, segs)
        im = ax.imshow(arr, extent=ext, origin="upper", vmin=vmin, vmax=vmax, cmap=cmap, interpolation="nearest", zorder=0)
        ax.set_title(f"D_b10 (km): {title}", fontsize=10)
        plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    fig.suptitle(f"Split-half reproducibility of the 0.025-degree map, {reg}: robust noise {mad / 2:.1f} km, "
                 f"anomaly correlation {corr_an:.2f}", fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(args.out, f"splithalf_{reg}.png"), dpi=150); plt.close(fig)

open(os.path.join(args.out, "resolution_study.md"), "w", encoding="utf-8").write("\n".join(report) + "\n")
print("wrote", os.path.join(args.out, "resolution_study.md"))
