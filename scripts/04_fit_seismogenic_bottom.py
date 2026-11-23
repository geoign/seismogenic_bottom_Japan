"""Map the lower limit of brittle seismicity by weighted maximum-likelihood fits of a
soft-bottomed depth distribution, node by node.

Usage: python 04_fit_seismogenic_bottom.py <declustered_npz> <out_npz> [--step 0.1] [--workers 20]
       [--ess-min 30] [--nmin 60] [--radii 10,15,20,30,40] [--bbox 122,148,24,46] [--no-weights] [--half 0|1]

Model for the depth density on 0-60 km at one node:
    f(z) = (1 - p) * g(z) + p / 60
    g(z) = Phi((z - z_t)/s_t) * [1 - Phi((z - z_b)/s_b)] / norm
z_t, s_t   : centre and width of the shallow onset of seismicity
z_b, s_b   : centre and width of the lower (brittle -> ductile) transition  <- target quantities
p          : fraction of a uniform "outlier" component (mislocated or slab-related events)
Each earthquake enters the likelihood with its family weight (1 / family size) so that an
aftershock sequence counts as one independent sample. The radius is widened adaptively until the
Kish effective sample size (sum w)^2 / sum w^2 reaches --ess-min.
Derived depths: D_b10 = z_b + 1.2816 s_b is where the brittle survival function drops to 10 %
(a model-based counterpart of D90); D_b90 = z_b - 1.2816 s_b is the top of the transition.
"""
import sys, argparse, time
import numpy as np
from scipy.special import ndtr
from scipy.optimize import minimize
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out")
ap.add_argument("--step", type=float, default=0.1)
ap.add_argument("--workers", type=int, default=1)
ap.add_argument("--ess-min", type=float, default=30)
ap.add_argument("--nmin", type=int, default=60)
ap.add_argument("--radii", default="10,15,20,30,40")
ap.add_argument("--bbox", default="122,148,24,46")
ap.add_argument("--no-weights", action="store_true")
ap.add_argument("--subset-npz", default=None,
                help="keep only events that also occur in this (more strictly filtered) catalogue file")
ap.add_argument("--half", type=int, default=None,
                help="use only families with (family id %% 2) == half (0 or 1): split-half reproducibility test")
args = ap.parse_args()
np.seterr(over="ignore", divide="ignore", invalid="ignore")

ZMAX = 60.0; DZ = 0.05
ZZ = np.arange(DZ / 2, ZMAX, DZ)                  # fine grid for normalisation
BIN = 0.25; EDGES = np.arange(0, ZMAX + BIN, BIN); ZC = EDGES[:-1] + BIN / 2
RADII = [float(r) for r in args.radii.split(",")]

G = {}
def init(src, no_weights):
    d = np.load(src)
    keep = np.ones(len(d["lat"]), dtype=bool)
    if args.subset_npz:
        import pandas as pd
        def keys(x):
            return pd.DataFrame(dict(t=x["time"].astype("M8[ms]").astype(np.int64), la=np.round(x["lat"] * 1e5).astype(np.int64),
                                     lo=np.round(x["lon"] * 1e5).astype(np.int64), z=np.round(x["dep"] * 100).astype(np.int64)))
        s = np.load(args.subset_npz)
        merged = keys(d).merge(keys(s).drop_duplicates().assign(ok=True), how="left", on=["t", "la", "lo", "z"])
        keep = merged["ok"].fillna(False).to_numpy(dtype=bool)
    if args.half is not None and "family" in d:
        keep &= (d["family"].astype(np.int64) % 2) == args.half
    G["lat"], G["lon"], G["dep"] = d["lat"][keep], d["lon"][keep], d["dep"][keep]
    w = np.ones(len(d["lat"])) if (no_weights or "weight" not in d) else d["weight"]
    G["w"] = w[keep]
    G["slab_cut"] = d["slab_cut"][keep] if "slab_cut" in d else np.full(keep.sum(), np.nan)
    G["slab_top"] = d["slab_top"][keep] if "slab_top" in d else np.full(keep.sum(), np.nan)
    R = 6371.0; la, lo = np.deg2rad(G["lat"]), np.deg2rad(G["lon"])
    G["xyz"] = np.c_[R*np.cos(la)*np.cos(lo), R*np.cos(la)*np.sin(lo), R*np.sin(la)]
    G["tree"] = cKDTree(G["xyz"])

def shape(zt, st, zb, sb, z):
    return ndtr((z - zt) / st) * (1.0 - ndtr((z - zb) / sb))

PMAX = 0.05   # cap on the uniform outlier fraction, so that it cannot swallow a genuine deeper population

def unpack(q):
    zt = q[0]; st = np.exp(q[1]); zb = zt + 0.5 + np.exp(q[2]); sb = np.exp(q[3])
    p = PMAX / (1.0 + np.exp(-q[4]))
    return zt, st, zb, sb, p

def nll_q(q, W):
    zt, st, zb, sb, p = unpack(q)
    if not (0.3 <= st <= 15 and 0.3 <= sb <= 15 and -10 <= zt <= 40 and zb <= 58):
        return 1e12
    g = shape(zt, st, zb, sb, ZZ); norm = g.sum() * DZ
    if norm <= 1e-12: return 1e12
    f = (1 - p) * shape(zt, st, zb, sb, ZC) / norm + p / ZMAX
    return -np.sum(W * np.log(f + 1e-300))

def nll_theta(th, W):
    zt, st, zb, sb, p = th
    g = shape(zt, st, zb, sb, ZZ); norm = g.sum() * DZ
    f = (1 - p) * shape(zt, st, zb, sb, ZC) / norm + p / ZMAX
    return -np.sum(W * np.log(f + 1e-300))

def wquant(z, w, q):
    o = np.argsort(z); cw = np.cumsum(w[o]); cw /= cw[-1]
    return np.interp(q, cw, z[o])

def fit_node(z, w):
    W = np.bincount(np.minimum((z / BIN).astype(int), len(ZC) - 1), weights=w, minlength=len(ZC))
    d10, d50, d90, d95 = wquant(z, w, [0.10, 0.5, 0.90, 0.95])
    starts = [(d10, 2.0, d90, 2.0, 0.01), (d10, 3.0, d95, 4.0, 0.02), (d10, 1.0, d50 + 0.6*(d90-d50), 3.0, 0.01)]
    best = None
    for zt, st, zb, sb, p in starts:
        pr = p / PMAX
        q0 = np.array([zt, np.log(st), np.log(max(zb - zt - 0.5, 0.5)), np.log(sb), np.log(pr / (1 - pr))])
        r = minimize(nll_q, q0, args=(W,), method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-4, maxiter=4000, maxfev=6000))
        if best is None or r.fun < best.fun: best = r
    th = np.array(unpack(best.x))
    f0 = nll_theta(th, W)
    # standard error of z_b from the observed information in the (z_b, s_b) block
    # (finite-difference Hessian; the other parameters are held at their optimum)
    se_zb = np.nan
    try:
        h = np.array([0.0, 0.0, 0.25, 0.1, 0.0]); H = np.zeros((2, 2)); ids = [2, 3]
        for a_, i in enumerate(ids):
            for b_, j in enumerate(ids):
                if b_ < a_: continue
                ei = np.zeros(5); ej = np.zeros(5); ei[i] = h[i]; ej[j] = h[j]
                H[a_, b_] = H[b_, a_] = (nll_theta(th + ei + ej, W) - nll_theta(th + ei - ej, W)
                                         - nll_theta(th - ei + ej, W) + nll_theta(th - ei - ej, W)) / (4 * h[i] * h[j])
        if np.all(np.linalg.eigvalsh(H) > 0):
            se_zb = float(np.sqrt(np.linalg.inv(H)[0, 0]))
        elif H[0, 0] > 0:
            se_zb = float(1 / np.sqrt(H[0, 0]))
    except Exception:
        pass
    return dict(zt=th[0], st=th[1], zb=th[2], sb=th[3], p=th[4], se_zb=se_zb, nll=float(f0),
                d90=d90, d95=d95, d10=d10)

def work(nodes):
    dep, w, xyz, tree = G["dep"], G["w"], G["xyz"], G["tree"]
    R = 6371.0
    out = []
    for (iy, ix, nlat, nlon) in nodes:
        la, lo = np.deg2rad(nlat), np.deg2rad(nlon)
        c = np.array([R*np.cos(la)*np.cos(lo), R*np.cos(la)*np.sin(lo), R*np.sin(la)])
        idx = np.array(tree.query_ball_point(c, RADII[-1]), dtype=np.int64)
        if len(idx) < args.nmin: continue
        dist = np.linalg.norm(xyz[idx] - c, axis=1)
        chosen = None
        for r in RADII:
            m = dist <= r
            if m.sum() < args.nmin: continue
            ww = w[idx[m]]; ess = ww.sum()**2 / np.sum(ww**2)
            if ess >= args.ess_min: chosen = (r, m, ess); break
        if chosen is None: continue
        r, m, ess = chosen
        z = dep[idx[m]]; ww = w[idx[m]]
        res = fit_node(z, ww)
        res.update(iy=iy, ix=ix, radius=r, ess=ess, n=int(m.sum()), wsum=float(ww.sum()),
                   slab_cut=float(np.nanmedian(G["slab_cut"][idx[m]])), slab_top=float(np.nanmedian(G["slab_top"][idx[m]])),
                   frac30=float(ww[z > 30].sum() / ww.sum()))
        out.append(res)
    return out

if __name__ == "__main__":
    x0, x1, y0, y1 = [float(v) for v in args.bbox.split(",")]
    lons = np.round(np.arange(x0, x1 + 1e-9, args.step), 6); lats = np.round(np.arange(y1, y0 - 1e-9, -args.step), 6)
    nodes = [(iy, ix, nlat, nlon) for iy, nlat in enumerate(lats) for ix, nlon in enumerate(lons)]
    # quick pre-screen: drop nodes with too few events within the largest radius
    init(args.src, args.no_weights)
    R = 6371.0; la, lo = np.deg2rad(lats), np.deg2rad(lons)
    LA, LO = np.meshgrid(la, lo, indexing="ij")
    cxyz = np.c_[(R*np.cos(LA)*np.cos(LO)).ravel(), (R*np.cos(LA)*np.sin(LO)).ravel(), (R*np.sin(LA)).ravel()]
    cnt = G["tree"].query_ball_point(cxyz, RADII[-1], return_length=True)
    keep = cnt >= args.nmin
    nodes = [nd for nd, k in zip(nodes, keep) if k]
    print(f"{len(lats)}x{len(lons)} grid, {len(nodes)} nodes with >= {args.nmin} events within {RADII[-1]} km", flush=True)
    t0 = time.time()
    nchunk = max(1, args.workers * 8)
    chunks = [nodes[i::nchunk] for i in range(nchunk)]
    if args.workers > 1:
        import multiprocessing as mp
        with mp.Pool(args.workers, initializer=init, initargs=(args.src, args.no_weights)) as pool:
            results = [r for part in pool.imap_unordered(work, chunks) for r in part]
    else:
        results = [r for ch in chunks for r in work(ch)]
    print(f"fitted {len(results)} nodes in {time.time()-t0:.0f} s", flush=True)
    fields = ["zt", "st", "zb", "sb", "p", "se_zb", "nll", "d10", "d90", "d95", "radius", "ess", "n", "wsum",
              "slab_cut", "slab_top", "frac30"]
    grids = {f: np.full((len(lats), len(lons)), np.nan, dtype=np.float32) for f in fields}
    for r in results:
        for f in fields: grids[f][r["iy"], r["ix"]] = r[f]
    np.savez(args.out, lats=lats, lons=lons, **grids)
    print("saved", args.out)
