"""Parse JMA unified hypocenter files (96-byte fixed-width records) into a NumPy archive.

Usage: python 01_parse_jma.py <raw_dir> <out_npz> [first_year] [last_year]
Only record type 'J' (JMA hypocenters) is kept. All other filtering is left to later steps.
"""
import sys, glob, os
import numpy as np
import pandas as pd

raw_dir, out = sys.argv[1], sys.argv[2]
y0 = int(sys.argv[3]) if len(sys.argv) > 3 else 1900
y1 = int(sys.argv[4]) if len(sys.argv) > 4 else 2100
files = sorted(f for f in glob.glob(os.path.join(raw_dir, "h????")) if y0 <= int(os.path.basename(f)[1:]) <= y1)

def num(a, lo, hi, scale):
    """Fixed-width numeric field -> float array; blanks and malformed fields -> NaN."""
    s = a[:, lo:hi].copy().view(f"S{hi-lo}").ravel()
    s = np.char.strip(s)
    v = pd.to_numeric(pd.Series(s).str.decode("ascii"), errors="coerce").to_numpy(dtype=float)
    return v / scale

def mag(a, lo):
    s = a[:, lo:lo+2].copy().view("S2").ravel()
    v = np.full(len(s), np.nan)
    first = a[:, lo]
    digit_ok = (first >= 48) & (first <= 57) | (first == 45) | (first == 32)
    plain = digit_ok & (np.char.strip(s) != b"")
    v[plain] = np.char.strip(s[plain]).astype(float) / 10.0
    # negative magnitudes coded as -1..-9 (=-0.1..-0.9) and A0..C9 (=-1.0..-3.9)
    neg1 = first == 45  # '-'
    v[neg1] = -(a[neg1, lo+1].astype(float) - 48) / 10.0
    letter = (first >= 65) & (first <= 67)
    v[letter] = -((first[letter].astype(float) - 65 + 1) + (a[letter, lo+1].astype(float) - 48) / 10.0)
    return v

parts = []
for fn in files:
    b = open(fn, "rb").read()
    lines = b.split(b"\n")
    lines = [l for l in lines if len(l) == 96]
    a = np.frombuffer(b"".join(lines), dtype=np.uint8).reshape(-1, 96)
    a = a[a[:, 0] == ord("J")]
    n = len(a)
    Y = num(a, 1, 5, 1); M = num(a, 5, 7, 1); D = num(a, 7, 9, 1)
    h = num(a, 9, 11, 1); mi = num(a, 11, 13, 1); sec = num(a, 13, 17, 100)
    sec = np.nan_to_num(sec)
    t = ((Y - 1970).astype("M8[Y]") + (M - 1).astype("m8[M]") + (D - 1).astype("m8[D]")
         + h.astype("m8[h]") + mi.astype("m8[m]") + np.round(sec * 1000).astype("m8[ms]"))
    rec = dict(
        time=t.astype("M8[ms]"),
        lat=num(a, 21, 24, 1) + num(a, 24, 28, 100) / 60.0,
        lon=num(a, 32, 36, 1) + num(a, 36, 40, 100) / 60.0,
        dep=num(a, 44, 49, 100),
        err_t=num(a, 17, 21, 100),
        err_lat_km=num(a, 28, 32, 100) * 1.852,          # minutes of latitude -> km
        err_lon_min=num(a, 40, 44, 100),
        err_dep=num(a, 49, 52, 100),
        mag=mag(a, 52),
        magtype=a[:, 54].copy(),
        ttt=a[:, 58].copy(), eval=a[:, 59].copy(), aux=a[:, 60].copy(), det=a[:, 95].copy(),
        nsta=num(a, 92, 95, 1),
    )
    rec["err_lon_km"] = rec["err_lon_min"] * 1.852 * np.cos(np.deg2rad(rec["lat"]))
    del rec["err_lon_min"]
    parts.append(rec)
    print(os.path.basename(fn), n, "records", flush=True)

allrec = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
np.savez(out, **allrec)
print("saved", out, len(allrec["lat"]), "events")
