"""Grid the Hirose/Nakajima plate-configuration contours (MRI distribution plate_data.tar.gz)
into slab-top surfaces for the Philippine Sea (PHS) and Pacific (PAC) plates.

Usage: python 02b_hirose_surfaces.py <plate_data_dir> <out_npz>
Output: 0.05-degree grids (lon 122-148, lat 24-46) of slab-top depth in km (positive down),
NaN outside the convex hull of the source points.
"""
import sys, os, glob, re
import numpy as np
from scipy.interpolate import LinearNDInterpolator

src, out = sys.argv[1], sys.argv[2]

def read_contours(folders):
    pts = []
    for fo in folders:
        for fn in glob.glob(os.path.join(src, fo, "plate_d*km*.dat")):
            dep = float(re.search(r"plate_d(\d+)km", os.path.basename(fn)).group(1))
            rows = []
            for line in open(fn, encoding="utf-8", errors="replace"):
                parts = line.replace(",", " ").split()
                if len(parts) < 2 or parts[0].startswith(">") or parts[0].startswith("#"):
                    continue  # GMT segment separators and comments
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
            if not rows: continue
            a = np.array(rows)
            pts.append(np.c_[a[:, 0], a[:, 1], np.full(len(a), dep)])
            print(f"  {fo}/{os.path.basename(fn)}: {len(a)} points at {dep:.0f} km")
    return np.vstack(pts)

print("PHS contours:")
phs = read_contours(["sw", "kanto", "izu3", "nakajima", "baba"])
print("PAC data:")
pac_raw = np.loadtxt(os.path.join(src, "PAC", "plate_combine.dat"), ndmin=2)
print("  PAC/plate_combine.dat columns:", pac_raw.shape, "first rows:", pac_raw[:3])
if pac_raw.shape[1] >= 3:
    pac = np.c_[pac_raw[:, 0], pac_raw[:, 1], np.abs(pac_raw[:, 2])]
else:
    raise SystemExit("unexpected PAC format")

lons = np.arange(122, 148 + 1e-9, 0.05); lats = np.arange(46, 24 - 1e-9, -0.05)
LO, LA = np.meshgrid(lons, lats)
def grid(p):
    f = LinearNDInterpolator(p[:, :2], p[:, 2])
    return f(LO, LA).astype(np.float32)
g_phs = grid(phs); g_pac = grid(pac)
print("PHS grid valid:", int(np.isfinite(g_phs).sum()), "PAC grid valid:", int(np.isfinite(g_pac).sum()))
np.savez(out, lons=lons, lats=lats, phs=g_phs, pac=g_pac)
for name, x, y in [("Bungo", 132.0, 33.4), ("W Ehime", 132.4, 33.6), ("Saitama", 139.5, 36.0), ("Ibaraki", 140.0, 36.5),
                   ("Aso", 131.1, 32.9), ("Tokyo", 139.7, 35.7), ("Kii", 135.5, 34.0), ("Sendai", 140.9, 38.3)]:
    ix = int(round((x - 122) / 0.05)); iy = int(round((46 - y) / 0.05))
    print(f"  {name:8s} PHS {g_phs[iy, ix]:6.1f}  PAC {g_pac[iy, ix]:6.1f}")
print("saved", out)
