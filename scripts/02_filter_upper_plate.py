"""Select well-located, natural, upper-plate (crust + overlying mantle) earthquakes.

Usage: python 02_filter_upper_plate.py <in_npz> <slab2_dir> <out_npz> [max_err_dep_km] [max_err_h_km] [hirose_surfaces_npz]

Criteria
  record type J, hypocenter evaluation '1' (depth free), auxiliary info '1' (ordinary earthquake;
  drops artificial events, eruption-related events and low-frequency events), determination flag
  K / k / A, finite depth 0-60 km, depth s.e. <= max_err_dep, horizontal s.e. <= max_err_h,
  region 122-148E / 24-46N, and depth shallower than the Slab2 slab top minus a margin
  (margin = max(5 km, Slab2 depth uncertainty); a floor of 15 km is applied so that the thin
  toe region and the Izu collision zone keep their shallow crustal seismicity).
"""
import sys, os, glob
import numpy as np
import rasterio
from scipy.ndimage import map_coordinates

src, slab_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]
max_err_dep = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
max_err_h = float(sys.argv[5]) if len(sys.argv) > 5 else 5.0

d = dict(np.load(src))
n0 = len(d["lat"])
def report(name, m):
    print(f"{name:<34s} keep {int(m.sum()):>9d} / {n0}")

m_eval = d["eval"] == ord("1"); report("depth-free solution", m_eval)
m_aux = d["aux"] == ord("1"); report("ordinary (non-LFE, non-artificial)", m_aux)
m_det = np.isin(d["det"], [ord("K"), ord("k"), ord("A")]); report("flag K/k/A", m_det)
m_dep = np.isfinite(d["dep"]) & (d["dep"] >= 0) & (d["dep"] <= 60); report("0 <= depth <= 60 km", m_dep)
err_h = np.hypot(d["err_lat_km"], d["err_lon_km"])
m_err = (d["err_dep"] <= max_err_dep) & (err_h <= max_err_h); report(f"s.e. dep<={max_err_dep}, h<={max_err_h} km", m_err)
m_box = (d["lon"] >= 122) & (d["lon"] <= 148) & (d["lat"] >= 24) & (d["lat"] <= 46); report("region box", m_box)
m = m_eval & m_aux & m_det & m_dep & m_err & m_box
report("all of the above", m)

sel = {k: v[m] for k, v in d.items()}
lon, lat, dep = sel["lon"], sel["lat"], sel["dep"]
print("depth s.e. quantiles (selected) 50/90/99 %:", np.percentile(sel["err_dep"], [50, 90, 99]))

# --- Slab2: shallowest slab top under each event (positive km down) and its uncertainty
slab_top = np.full(len(lon), np.nan); slab_unc = np.full(len(lon), np.nan)
for reg in ["kur", "izu", "ryu", "phi"]:
    fdep = glob.glob(os.path.join(slab_dir, f"{reg}_slab2_dep_*.grd"))[0]
    func = glob.glob(os.path.join(slab_dir, f"{reg}_slab2_unc_*.grd"))[0]
    with rasterio.open(fdep) as r:
        a = r.read(1).astype(float); tr = r.transform
    with rasterio.open(func) as r:
        u = r.read(1).astype(float)
    a[~np.isfinite(a)] = np.nan
    # pixel coordinates of events (row, col) for bilinear sampling
    col = (lon - tr.c) / tr.a - 0.5
    row = (lat - tr.f) / tr.e - 0.5
    inside = (col >= 0) & (col <= a.shape[1] - 1) & (row >= 0) & (row <= a.shape[0] - 1)
    v = np.full(len(lon), np.nan); vu = np.full(len(lon), np.nan)
    if inside.any():
        v[inside] = map_coordinates(a, [row[inside], col[inside]], order=1, mode="nearest")
        vu[inside] = map_coordinates(u, [row[inside], col[inside]], order=1, mode="nearest")
    top = -v  # Slab2 depths are negative (km below sea level)
    better = np.isfinite(top) & (~np.isfinite(slab_top) | (top < slab_top))
    slab_top[better] = top[better]; slab_unc[better] = vu[better]
    print(f"{reg}: {int(np.isfinite(top).sum())} events under this slab model")

# --- optional extra slab-top surfaces (Hirose/Nakajima contours gridded by 02b_hirose_surfaces.py)
if len(sys.argv) > 6:
    hs = np.load(sys.argv[6])
    hl, hla = hs["lons"], hs["lats"]
    for key in ("phs", "pac"):
        a = hs[key].astype(float)
        col = (lon - hl[0]) / (hl[1] - hl[0]); row = (lat - hla[0]) / (hla[1] - hla[0])
        inside = (col >= 0) & (col <= a.shape[1] - 1) & (row >= 0) & (row <= a.shape[0] - 1)
        top = np.full(len(lon), np.nan)
        top[inside] = map_coordinates(np.nan_to_num(a, nan=1e6), [row[inside], col[inside]], order=1, mode="nearest")
        top[top > 1e5] = np.nan
        better = np.isfinite(top) & (~np.isfinite(slab_top) | (top < slab_top))
        slab_top[better] = top[better]; slab_unc[better] = 5.0
        print(f"Hirose {key.upper()}: {int(np.isfinite(top).sum())} events under this surface, "
              f"{int(better.sum())} where it is shallower than Slab2")

margin = np.fmax(5.0, np.nan_to_num(slab_unc, nan=5.0))
cutoff = np.where(np.isfinite(slab_top), np.fmax(slab_top - margin, 15.0), 60.0)
m_slab = dep < cutoff
print(f"above slab top (with margin): keep {int(m_slab.sum())} / {len(dep)}")
for k in list(sel):
    sel[k] = sel[k][m_slab]
sel["slab_top"] = slab_top[m_slab]
sel["slab_cut"] = cutoff[m_slab]   # depth above which events were accepted at this epicentre
np.savez(out, **sel)
print("saved", out, len(sel["lat"]), "events")
print("depth quantiles 50/90/99 %:", np.percentile(sel["dep"], [50, 90, 99]))
print("magnitude quantiles 10/50/90 %:", np.nanpercentile(sel["mag"], [10, 50, 90]))
