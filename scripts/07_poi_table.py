"""Print fitted values at reference locations (nearest node) as a Markdown table.

Usage: python 07_poi_table.py <fit_npz> [out_md]
"""
import sys
import numpy as np

POI = [
    ("Aso caldera", 131.08, 32.88), ("Kumamoto (2016 M7.3)", 130.76, 32.79), ("Unzen", 130.30, 32.76),
    ("Kirishima", 130.86, 31.93), ("Beppu-Shimabara graben", 131.50, 33.28), ("Hiroshima", 132.46, 34.40),
    ("Tottori", 134.20, 35.40), ("Kobe / Rokko", 135.20, 34.70), ("Kii Peninsula", 135.60, 34.10),
    ("Shikoku, south of MTL", 133.50, 33.60), ("Nobi / Gifu", 136.80, 35.50), ("Matsumoto", 137.95, 36.25),
    ("Noto (2024 M7.5)", 137.25, 37.35), ("Chuetsu (2004 M6.8)", 138.85, 37.30), ("Hakone", 139.02, 35.23),
    ("Izu-Tobu", 139.10, 34.92), ("Tokyo", 139.70, 35.70), ("Tsukuba", 140.10, 36.10),
    ("Fukushima Hamadori", 140.70, 37.00), ("Sendai", 140.90, 38.30), ("Iwate-Miyagi (2008 M7.2)", 140.88, 39.03),
    ("Akita", 140.10, 39.70), ("Aomori / Hakkoda", 140.90, 40.65), ("Sapporo", 141.35, 43.05),
    ("Hidaka", 142.80, 42.60), ("Tokachi-dake", 142.69, 43.42), ("Kushiro", 144.40, 43.00),
]

d = np.load(sys.argv[1])
lats, lons = d["lats"], d["lons"]
rows = ["| Site | lon | lat | z_b km | s_b km | D_b10 km | D90w km | s.e. km | p | ESS | radius km |",
        "|---|---|---|---|---|---|---|---|---|---|---|"]
for name, x, y in POI:
    ix = int(np.argmin(np.abs(lons - x))); iy = int(np.argmin(np.abs(lats - y)))
    zb, sb, d90, se, p, ess, r = (float(d[k][iy, ix]) for k in ["zb", "sb", "d90", "se_zb", "p", "ess", "radius"])
    if np.isfinite(zb):
        rows.append(f"| {name} | {x:.2f} | {y:.2f} | {zb:.1f} | {sb:.1f} | {zb + 1.2816*sb:.1f} | {d90:.1f} | {se:.1f} | {p:.2f} | {ess:.0f} | {r:.0f} |")
    else:
        rows.append(f"| {name} | {x:.2f} | {y:.2f} | masked | | | | | | | |")
txt = "\n".join(rows)
print(txt)
if len(sys.argv) > 2:
    open(sys.argv[2], "w", encoding="utf-8").write(txt + "\n")
