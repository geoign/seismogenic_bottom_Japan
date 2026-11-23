"""Nearest-neighbour declustering (Zaliapin & Ben-Zion 2013) with cluster weights.

Usage: python 03_decluster_nn.py <in_npz> <out_npz> [b=1.0] [df=1.6] [log10_eta0 or 'auto']

For every event j the nearest neighbour (parent) i among earlier events is the one minimising
    eta_ij = t_ij [yr] * r_ij[km]^df * 10^(-b * m_i)
Search strategy: all earlier events within the neighbouring 0.25-degree cells, 20 km and 1.5 years,
plus an exhaustive pass over all earlier events with M >= 5 (their aftershock zones exceed the
local window). For parents with M < 5 the local window is far wider than any pair that could fall
below the cluster threshold (eta0 ~ 1e-5), so the clustered/background decision is unaffected.
Events with eta < eta0 are linked to their parent; connected components are families. Every family
gets total weight 1 (w_j = 1 / family size), so one aftershock sequence counts as one independent
sample of the depth distribution instead of thousands.
"""
import sys, time
import numpy as np
from numba import njit, prange
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

src, out = sys.argv[1], sys.argv[2]
b = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
df = float(sys.argv[4]) if len(sys.argv) > 4 else 1.6
eta0_arg = sys.argv[5] if len(sys.argv) > 5 else "auto"
TWIN, RLOC = 1.5, 20.0

d = dict(np.load(src))
n = len(d["lat"])
order = np.argsort(d["time"], kind="stable")
for k in d:
    d[k] = d[k][order]
t = (d["time"] - np.datetime64("2001-01-01")).astype("m8[s]").astype(np.float64) / (365.25 * 86400.0)
lat, lon = d["lat"], d["lon"]
mag = np.where(np.isfinite(d["mag"]), d["mag"], -2.0)
pm = 10.0 ** (-b * mag)
R = 6371.0
la, lo = np.deg2rad(lat), np.deg2rad(lon)
X = R * np.cos(la) * np.cos(lo); Y = R * np.cos(la) * np.sin(lo); Z = R * np.sin(la)

CELL = 0.25
cx = np.floor((lon - 120.0) / CELL).astype(np.int64)
cy = np.floor((lat - 20.0) / CELL).astype(np.int64)
ncx = int(cx.max()) + 2
cell = cy * ncx + cx
ncell = int(cell.max()) + ncx + 2
cell_order = np.lexsort((t, cell))            # sorted by cell, then time
cell_sorted = cell[cell_order]
cell_start = np.searchsorted(cell_sorted, np.arange(ncell), side="left")
cell_end = np.searchsorted(cell_sorted, np.arange(ncell), side="right")
big = np.where(mag >= 5.0)[0]
print(f"{n} events, {len(big)} parents with M>=5, {ncell} cells", flush=True)

@njit(parallel=True, fastmath=True)
def nn_search(j0, j1, t, X, Y, Z, pm, cell, ncx, cell_order, cell_start, cell_end, big, df, twin, rloc, eta, parent):
    tmin = 1.0 / (365.25 * 86400.0)   # 1 second in years
    r2max = rloc * rloc
    for j in prange(j0, j1):
        tj = t[j]; best = np.inf; bi = -1
        xj = X[j]; yj = Y[j]; zj = Z[j]
        # exhaustive pass over large earlier events
        for q in range(len(big)):
            i = big[q]
            if t[i] >= tj: break
            dt = tj - t[i]
            if dt < tmin: dt = tmin
            dx = xj-X[i]; dy = yj-Y[i]; dz = zj-Z[i]
            d2 = dx*dx+dy*dy+dz*dz
            r = np.sqrt(d2) if d2 > 0.01 else 0.1
            e = dt * r**df * pm[i]
            if e < best: best = e; bi = i
        # local pass over neighbouring cells within the time window
        c = cell[j]
        for dcy in range(-1, 2):
            for dcx in range(-1, 2):
                cc = c + dcy*ncx + dcx
                if cc < 0 or cc >= len(cell_start): continue
                lo = cell_start[cc]; hi = cell_end[cc]
                if hi <= lo: continue
                a = lo; bnd = hi
                while a < bnd:
                    m = (a + bnd) // 2
                    if t[cell_order[m]] < tj - twin: a = m + 1
                    else: bnd = m
                for k in range(a, hi):
                    i = cell_order[k]
                    if t[i] >= tj: break
                    if i == j: continue
                    dx = xj-X[i]; dy = yj-Y[i]; dz = zj-Z[i]
                    d2 = dx*dx+dy*dy+dz*dz
                    if d2 > r2max: continue
                    dt = tj - t[i]
                    if dt < tmin: dt = tmin
                    r = np.sqrt(d2) if d2 > 0.01 else 0.1
                    e = dt * r**df * pm[i]
                    if e < best: best = e; bi = i
        eta[j] = best; parent[j] = bi

eta = np.full(n, np.inf); parent = np.full(n, -1, dtype=np.int64)
t0 = time.time(); CH = 200000
for j0 in range(0, n, CH):
    j1 = min(n, j0 + CH)
    nn_search(j0, j1, t, X, Y, Z, pm, cell, ncx, cell_order, cell_start, cell_end, big, df, TWIN, RLOC, eta, parent)
    print(f"  {j1}/{n} events ({str(d['time'][j1-1])[:10]}) after {time.time()-t0:.0f} s", flush=True)
print(f"nearest-neighbour search done in {time.time()-t0:.0f} s", flush=True)
leta = np.log10(eta)
fin = np.isfinite(leta)
print(f"{int((~fin).sum())} events without any earlier neighbour in the window (treated as background)")
print("log10 eta quantiles 5/25/50/75/95 %:", np.percentile(leta[fin], [5, 25, 50, 75, 95]))

if eta0_arg == "auto":
    from sklearn.mixture import GaussianMixture
    x = leta[fin].reshape(-1, 1)
    rng = np.random.default_rng(0)
    xs = x[rng.choice(len(x), size=min(400000, len(x)), replace=False)]
    gm = GaussianMixture(2, random_state=0).fit(xs)
    mu = gm.means_.ravel(); sd = np.sqrt(gm.covariances_.ravel()); wgt = gm.weights_
    lo_, hi_ = np.sort(mu)
    grid = np.linspace(lo_, hi_, 2001)
    post = gm.predict_proba(grid.reshape(-1, 1))
    icl = int(np.argmin(mu))
    thr = grid[np.argmin(np.abs(post[:, icl] - 0.5))]
    print(f"GMM: means {mu}, sd {sd}, weights {wgt} -> threshold log10 eta0 = {thr:.2f}")
    log_eta0 = float(thr)
else:
    log_eta0 = float(eta0_arg)
clustered = fin & (leta < log_eta0) & (parent >= 0)
print(f"log10 eta0 = {log_eta0:.2f}: {int(clustered.sum())} clustered events ({100*clustered.mean():.1f} %)")

rows = np.where(clustered)[0]; cols = parent[rows]
g = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
ncomp, fam = connected_components(g, directed=False)
size = np.bincount(fam, minlength=ncomp)[fam]
w = 1.0 / size
print(f"{ncomp} families; largest {size.max()} events; effective sample size {w.sum()**2/np.sum(w**2):.0f}")
big_fams = np.argsort(-np.bincount(fam, minlength=ncomp))[:8]
for f in big_fams:
    mm = fam == f; i0 = np.argmax(np.where(mm, mag, -9))
    print(f"  family {f}: {int(mm.sum())} events, largest M{mag[i0]:.1f} at {d['time'][i0]} ({lat[i0]:.2f}N {lon[i0]:.2f}E {d['dep'][i0]:.0f} km)")

d.update(dict(t_yr=t, log10_eta=leta, parent=parent, family=fam, family_size=size, weight=w,
              log10_eta0=np.array(log_eta0), b=np.array(b), df=np.array(df)))
np.savez(out, **d)
print("saved", out)
