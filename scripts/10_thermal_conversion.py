"""Convert the base of brittle seismicity (D_b10) into heat flow, geothermal gradient and
isotherm depths, calibrated against heat-flow observations.

Usage: python 10_thermal_conversion.py <seismogenic_bottom_0p1deg.tif> <IHFC_2021_GHFDB.csv>
       <strain_rate.xyz> <out_dir> [--curie curie_depth.tif] [--coast GSHHS_i_L1.shp]
       [--k 2.5] [--a0 2.0] [--hr 10] [--smooth 1.0] [--boot 1000]

Logic
1. D_b10 (quality_flag == 0 only, lightly smoothed) is taken as the depth of an isotherm T_c.
2. T_c is not a constant. From a power-law flow law, eps = A sigma^n exp(-Q / R T), with the
   differential stress set by friction (proportional to depth):
       1/T_c = 1/T_ref + (R/Q) * [ n ln(z / z_ref) - ln(eps / eps_ref) ]
   n and Q are fixed (wet anorthite dislocation creep, Rybacki & Dresen 2000); eps is the second
   invariant of the geodetic strain rate. T_ref (z_ref = 15 km, median strain rate) is the one
   free parameter.
3. A steady-state conductive geotherm with exponentially decaying heat production,
       T(z) = T_s + (q_m / k) z + (A0 hr^2 / k) (1 - exp(-z / hr)),
   is pinned at the surface (T_s) and at (D_b10, T_c). This gives the reduced heat flow q_m, the
   surface heat flow q_s = q_m + A0 hr, the gradients and the depth of any isotherm.
4. T_ref is calibrated so that the predicted q_s matches onshore heat-flow observations (node
   medians, least absolute deviation of ln(q_obs / q_pred)); its 95 % interval is bootstrapped.

Inputs that are not redistributed with this repository:
  IHFC_2021_GHFDB.csv  Fuchs, Norden & International Heat Flow Commission (2021), The Global Heat
                       Flow Database: Release 2021, GFZ Data Services, doi:10.5880/fidgeo.2021.014
  strain_rate.xyz      Ohtate et al. (2026), Earth Planets Space 78:210,
                       doi:10.1186/s40623-026-02475-y, Supplementary Material 3
                       (columns: lon lat shear areal D; nanostrain/yr; 0.05-degree grid)
  curie_depth.tif      optional independent check: Li, Lu & Wang (2017), Sci. Rep. 7:45129,
                       doi:10.1038/srep45129 (global Curie-point depth, gridded to GeoTIFF)
"""
import argparse, os, sys
import numpy as np
import pandas as pd
from osgeo import gdal, osr
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize, minimize_scalar
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import NullFormatter, ScalarFormatter

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
gdal.UseExceptions()

ap = argparse.ArgumentParser()
ap.add_argument("tif"); ap.add_argument("heatflow"); ap.add_argument("strain"); ap.add_argument("out")
ap.add_argument("--curie", default=None, help="Curie-point depth GeoTIFF (km) for an independent check")
ap.add_argument("--coast", default="GSHHS_i_L1.shp")
ap.add_argument("--k", type=float, default=2.5, help="thermal conductivity, W/m/K")
ap.add_argument("--a0", type=float, default=2.0, help="surface heat production, microW/m3")
ap.add_argument("--hr", type=float, default=10.0, help="e-folding depth of heat production, km")
ap.add_argument("--smooth", type=float, default=1.0, help="Gaussian sigma (nodes) applied to D_b10 and ln(strain rate)")
ap.add_argument("--boot", type=int, default=1000)
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

R_GAS = 8.314
N_EXP, Q_ACT = 3.0, 356e3          # wet anorthite, dislocation creep (Rybacki & Dresen 2000)
Z_REF = 15.0                       # km
Q_RANGE = (20.0, 300.0)            # mW/m2, heat-flow observations kept for the calibration
BBOX = (122, 148, 24, 46)
report = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    report.append(s)


def nansmooth(a, sigma):
    if sigma <= 0:
        return a.copy()
    f = np.where(np.isfinite(a), a, 0.0); w = np.isfinite(a).astype(float)
    out = gaussian_filter(f, sigma) / np.maximum(gaussian_filter(w, sigma), 1e-9)
    out[~np.isfinite(a)] = np.nan
    return out


def surface_temperature(lat):
    """Mean annual ground temperature at sea level, deg C (linear in latitude across Japan)."""
    return np.clip(16.0 - 0.9 * (lat - 35.0), 5.0, 24.0)


def radiogenic(z_km, k, a0, hr):
    """Temperature contribution of crustal heat production at depth z, deg C."""
    return a0 * 1e-6 * (hr * 1e3) ** 2 / k * (1.0 - np.exp(-z_km / hr))


def t_c(z_km, lneps_rel, t_ref, beta=N_EXP * R_GAS / Q_ACT, gamma=R_GAS / Q_ACT):
    """Temperature (deg C) at the base of brittle seismicity."""
    inv = 1.0 / (t_ref + 273.15) + beta * np.log(z_km / Z_REF) - gamma * lneps_rel
    return 1.0 / inv - 273.15


def q_mantle(z_km, tc, ts, k, a0, hr):
    """Reduced heat flow (W/m2) of the geotherm through (0, ts) and (z, tc)."""
    return (tc - ts - radiogenic(z_km, k, a0, hr)) * k / (z_km * 1e3)


def q_surface(z_km, tc, ts, k, a0, hr):
    return (q_mantle(z_km, tc, ts, k, a0, hr) + a0 * 1e-6 * hr * 1e3) * 1e3      # mW/m2


def isotherm_depth(t_target, qm, ts, k, a0, hr):
    """Depth (km) at which the geotherm reaches t_target; bisection on 0-150 km."""
    lo = np.zeros_like(qm); hi = np.full_like(qm, 150.0)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        t = ts + qm / k * mid * 1e3 + radiogenic(mid, k, a0, hr)
        below = t < t_target
        lo = np.where(below, mid, lo); hi = np.where(below, hi, mid)
    z = 0.5 * (lo + hi)
    z[~np.isfinite(qm) | (qm <= 0) | (z > 149.0)] = np.nan
    return z


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


# ------------------------------------------------------------------ inputs
ds = gdal.Open(args.tif)
gt = ds.GetGeoTransform(); nx, ny = ds.RasterXSize, ds.RasterYSize
step = gt[1]
lons = gt[0] + step / 2 + step * np.arange(nx); lats = gt[3] + gt[5] / 2 + gt[5] * np.arange(ny)
LO, LA = np.meshgrid(lons, lats)
B = {ds.GetRasterBand(i + 1).GetDescription(): ds.GetRasterBand(i + 1).ReadAsArray() for i in range(ds.RasterCount)}
clean = B["quality_flag"] == 0

sr = pd.read_csv(args.strain, sep=r"\s+")
sr.columns = ["lon", "lat", "shear", "areal", "D"]
e1 = sr.areal / 2 + sr.shear; e2 = sr.areal / 2 - sr.shear; e3 = -(e1 + e2)      # incompressible
sr["eII"] = np.sqrt(0.5 * (e1 ** 2 + e2 ** 2 + e3 ** 2))                           # nanostrain/yr
LNE = griddata((sr.lon, sr.lat), np.log(sr.eII), (LO, LA), method="linear")
dist, _ = cKDTree(np.c_[sr.lon, sr.lat]).query(np.c_[LO.ravel(), LA.ravel()])
LNE[dist.reshape(LO.shape) > 0.15] = np.nan                                       # onshore GNSS coverage only
lo_c, hi_c = np.nanpercentile(LNE, [1, 99]); LNE = np.clip(LNE, lo_c, hi_c)

domain = clean & np.isfinite(LNE)
TS = surface_temperature(LA)


def prepared(depth_band):
    z = nansmooth(np.where(domain, B[depth_band], np.nan), args.smooth)
    lne = nansmooth(np.where(domain, LNE, np.nan), args.smooth)
    return z, lne - np.nanmedian(lne[np.isfinite(z)]), float(np.exp(np.nanmedian(lne[np.isfinite(z)])))


Z, LNE_REL, EPS_REF = prepared("D_b10")

hf = pd.read_csv(args.heatflow, sep=";", low_memory=False, encoding="latin-1")
for c in ["lat", "lng", "q"]:
    hf[c] = pd.to_numeric(hf[c].astype(str).str.replace(",", "."), errors="coerce")
hf = hf[(hf.lng >= BBOX[0]) & (hf.lng <= BBOX[1]) & (hf.lat >= BBOX[2]) & (hf.lat <= BBOX[3]) & hf.q.notna()
        & (hf.A_dom == "continental")].copy()
hf["ix"] = np.round((hf.lng.values - lons[0]) / step).astype(int)
hf["iy"] = np.round((hf.lat.values - lats[0]) / gt[5]).astype(int)
hf = hf[(hf.ix >= 0) & (hf.ix < nx) & (hf.iy >= 0) & (hf.iy < ny)]
hf["node"] = hf.iy * nx + hf.ix


def node_table(z, lne_rel):
    t = hf.assign(z=z[hf.iy, hf.ix], lne=lne_rel[hf.iy, hf.ix]).dropna(subset=["z", "lne"])
    g = t.groupby("node").agg(q=("q", "median"), nobs=("q", "size"), z=("z", "first"), lne=("lne", "first"),
                              lat=("lat", "mean"), lon=("lng", "mean")).reset_index()
    return g[(g.q >= Q_RANGE[0]) & (g.q <= Q_RANGE[1])].reset_index(drop=True)


def calibrate(g, k, a0, hr, nboot=0, seed=0, free=False):
    """LAD fit of T_ref (and optionally beta, gamma) to ln(q_obs / q_pred)."""
    z, lne, q = g.z.values, g.lne.values, g.q.values
    ts = surface_temperature(g.lat.values)

    def resid(p, idx):
        tc = t_c(z[idx], lne[idx], *p) if free else t_c(z[idx], lne[idx], p[0])
        return np.log(q[idx] / np.clip(q_surface(z[idx], tc, ts[idx], k, a0, hr), 1.0, None))

    def fit(idx, start):
        if free:
            return minimize(lambda p: np.abs(resid(p, idx)).sum(), start, method="Nelder-Mead",
                            options=dict(xatol=1e-7, fatol=1e-6, maxiter=4000)).x
        r = minimize_scalar(lambda t: np.abs(resid([t], idx)).sum(), bounds=(150, 800), method="bounded")
        return np.array([r.x])

    allidx = np.arange(len(g))
    best = fit(allidx, [400.0, N_EXP * R_GAS / Q_ACT, R_GAS / Q_ACT])
    boots = None
    if nboot:
        rng = np.random.default_rng(seed)
        boots = np.array([fit(rng.integers(0, len(g), len(g)), best) for _ in range(nboot)])
    return best, boots, resid(best, allidx)


K, A0, HR = args.k, args.a0, args.hr
G = node_table(Z, LNE_REL)
best, boots, res = calibrate(G, K, A0, HR, nboot=args.boot)
T_REF = float(best[0]); T_LO, T_HI = [float(v) for v in np.percentile(boots[:, 0], [2.5, 97.5])]
mad = float(1.4826 * np.median(np.abs(res - np.median(res))))

say("# Thermal conversion of the base of brittle seismicity\n")
say("Made by `scripts/10_thermal_conversion.py`. D_b10 is treated as the depth of an isotherm whose temperature T_c "
    "is calibrated against onshore heat-flow observations; a conductive geotherm through the surface and (D_b10, T_c) "
    "then gives heat flow, gradients and isotherm depths. Independent work based on the JMA unified hypocenter "
    "catalogue; not a product of JMA or of the providers of the heat-flow and strain-rate data.\n")
say("## 1. Calibration\n")
say(f"* Heat-flow observations: {len(hf)} onshore entries of the IHFC Global Heat Flow Database (Release 2021) in "
    f"{BBOX[0]}-{BBOX[1]}E, {BBOX[2]}-{BBOX[3]}N; {int(G.nobs.sum())} of them fall on {len(G)} nodes with a usable D_b10 "
    f"(quality_flag 0, onshore strain-rate coverage, node median {Q_RANGE[0]:.0f}-{Q_RANGE[1]:.0f} mW/m2).")
say(f"* Spearman correlation between observed heat flow and D_b10 at these nodes: {spearmanr(G.q, G.z).statistic:+.2f}.")
say(f"* Geotherm parameters: k = {K} W/m/K, A0 = {A0} microW/m3, hr = {HR} km, T_s = 16 - 0.9 (lat - 35) deg C.")
say(f"* T_c model: 1/T_c = 1/T_ref + (R/Q) [n ln(z/{Z_REF:.0f} km) - ln(eps/eps_ref)], n = {N_EXP:g}, Q = {Q_ACT/1e3:.0f} kJ/mol, "
    f"eps_ref = {EPS_REF:.0f} nanostrain/yr (median second invariant over the domain).")
say(f"* **Calibrated T_ref = {T_REF:.0f} deg C (95 % bootstrap interval {T_LO:.0f}-{T_HI:.0f} deg C)**: the temperature at "
    f"D_b10 when D_b10 = {Z_REF:.0f} km and the strain rate is at its median.")
say(f"* Residual ln(q_obs / q_pred): median {np.median(res):+.3f}, robust s.d. {mad:.2f} (a factor of {np.exp(mad):.2f}); "
    f"Spearman correlation of predicted and observed heat flow {spearmanr(G.q, G.q * np.exp(-res)).statistic:+.2f}. "
    f"The scatter is that of single heat-flow sites against a 10-40 km average of seismicity, so it is an upper bound "
    f"on the error of the grid.\n")

# which depth statistic corresponds to which temperature
say("The same calibration applied to other depth statistics of the same product (constant strain-rate and depth terms as above):\n")
say("| Depth statistic | Median depth over the domain (km) | Calibrated T_ref (deg C) | 95 % interval |")
say("| --- | --- | --- | --- |")
for band, label in [("D_b10", "D_b10 (10 % of the plateau rate)"), ("D90_weighted", "classic D90 (family-weighted)"),
                    ("z_b", "z_b (centre of the transition)"), ("D_b90", "D_b90 (top of the transition)")]:
    zz, ll, _ = prepared(band)
    gg = node_table(zz, ll)
    b_, bo_, _ = calibrate(gg, K, A0, HR, nboot=max(200, args.boot // 5), seed=1)
    lo_, hi_ = (T_LO, T_HI) if band == "D_b10" else np.percentile(bo_[:, 0], [2.5, 97.5])
    say(f"| {label} | {np.nanmedian(zz):.1f} | {b_[0]:.0f} | {lo_:.0f}-{hi_:.0f} |")
say("")

# is the functional form supported by the data?
bf, bof, resf = calibrate(G, K, A0, HR, nboot=max(200, args.boot // 5), seed=2, free=True)
lo_f, hi_f = np.percentile(bof, [2.5, 97.5], axis=0)
say("## 2. What the heat-flow data say about the depth and strain-rate terms\n")
say("Fitting the two coefficients freely instead of fixing them from the flow law:\n")
say("| Coefficient | Flow law (wet anorthite) | Flow law (wet quartz, n = 4, Q = 135 kJ/mol) | Free fit | 95 % interval |")
say("| --- | --- | --- | --- | --- |")
say(f"| depth term n R/Q (1/K) | {N_EXP*R_GAS/Q_ACT:.2e} | {4*R_GAS/135e3:.2e} | {bf[1]:.2e} | {lo_f[1]:.2e} to {hi_f[1]:.2e} |")
say(f"| strain-rate term R/Q (1/K) | {R_GAS/Q_ACT:.2e} | {R_GAS/135e3:.2e} | {bf[2]:.2e} | {lo_f[2]:.2e} to {hi_f[2]:.2e} |")
say("")
say(f"Mean absolute residual: {np.abs(res).mean():.3f} (flow-law coefficients, this product), {np.abs(resf).mean():.3f} (free fit). "
    "The depth term is not constrained by the data. The strain-rate term has the sign the flow law predicts but is several "
    "times larger; since hot, weak crust also strains faster, and the geodetic field (2019-2023) contains elastic strain "
    "from plate coupling and post-seismic deformation of the 2011 Tohoku-oki earthquake, the free-fit value is read as a "
    "statistical association with heat flow, not as rheology, and is not used.\n")

# sensitivity to the geotherm parameters
say("## 3. Sensitivity to the geotherm parameters\n")
say("T_ref is re-calibrated for each parameter set, so the predicted heat flow is anchored to the observations; the "
    "temperature assigned to D_b10, and with it the depth of the 350 deg C isotherm, is what moves.\n")
say("| k (W/m/K) | A0 (microW/m3) | T_ref (deg C) | median heat flow (mW/m2) | median depth of 350 deg C (km) | median depth of 580 deg C (km) |")
say("| --- | --- | --- | --- | --- | --- |")
ok = np.isfinite(Z)
for k_, a_ in [(2.0, 2.0), (2.5, 1.0), (2.5, 2.0), (2.5, 3.0), (3.0, 2.0)]:
    b_, _, _ = calibrate(G, k_, a_, HR)
    tc_ = t_c(Z, LNE_REL, b_[0]); qm_ = q_mantle(Z, tc_, TS, k_, a_, HR)
    say(f"| {k_} | {a_} | {b_[0]:.0f} | {np.nanmedian(q_surface(Z, tc_, TS, k_, a_, HR)):.0f} | "
        f"{np.nanmedian(isotherm_depth(350.0, qm_, TS, k_, a_, HR)):.1f} | {np.nanmedian(isotherm_depth(580.0, qm_, TS, k_, a_, HR)):.1f} |")
say("")

# ------------------------------------------------------------------ grids
def grids(t_ref, k=None):
    k = K if k is None else k
    tc = t_c(Z, LNE_REL, t_ref)
    qm = q_mantle(Z, tc, TS, k, A0, HR)
    qs = q_surface(Z, tc, TS, k, A0, HR)
    return dict(tc=tc, qs=qs, g_surf=qs / k, g_mean=(tc - TS) / Z,
                z350=isotherm_depth(350.0, qm, TS, k, A0, HR), z580=isotherm_depth(580.0, qm, TS, k, A0, HR))


K_ALT = 3.0                                    # alternative conductivity: cooler D_b10, deeper isotherms
T_ALT = float(calibrate(G, K_ALT, A0, HR)[0][0])
M, MLO, MHI, MALT = grids(T_REF), grids(T_LO), grids(T_HI), grids(T_ALT, K_ALT)
say("## 4. The grids\n")
say(f"{int(ok.sum())} nodes (0.1 degree; quality_flag 0 and onshore strain-rate coverage). Quantiles over these nodes:\n")
say("| Quantity | 5 % | 25 % | 50 % | 75 % | 95 % |")
say("| --- | --- | --- | --- | --- | --- |")
for lab, a in [("D_b10, smoothed input (km)", Z), ("T at D_b10 (deg C)", M["tc"]), ("surface heat flow (mW/m2)", M["qs"]),
               ("surface gradient (deg C/km)", M["g_surf"]), ("mean gradient to D_b10 (deg C/km)", M["g_mean"]),
               ("depth of 350 deg C (km)", M["z350"]), ("depth of 580 deg C (km)", M["z580"]),
               ("depth of 350 deg C / D_b10", M["z350"] / Z)]:
    qv = np.nanpercentile(a, [5, 25, 50, 75, 95])
    say(f"| {lab} | " + " | ".join(f"{v:.2f}" if qv[2] < 3 else f"{v:.0f}" if qv[2] > 100 else f"{v:.1f}" for v in qv) + " |")
say("")
say(f"With T_ref at the ends of its 95 % interval the median heat flow is {np.nanmedian(MLO['qs']):.0f} and "
    f"{np.nanmedian(MHI['qs']):.0f} mW/m2 (bands 7-8). The depth of 350 deg C depends more on the conductivity than on "
    f"the calibration error: its median is {np.nanmedian(MHI['z350']):.1f} km with T_ref at the upper end of the interval "
    f"(band 9) and {np.nanmedian(MALT['z350']):.1f} km with k = {K_ALT} W/m/K, for which T_ref recalibrates to "
    f"{T_ALT:.0f} deg C (band 10; the surface gradient is then heat flow / {K_ALT}).\n")

if args.curie and os.path.exists(args.curie):
    cu = gdal.Open(args.curie); cg = cu.GetGeoTransform(); ca = cu.GetRasterBand(1).ReadAsArray()
    cx = np.clip(np.round((LO - (cg[0] + cg[1] / 2)) / cg[1]).astype(int), 0, cu.RasterXSize - 1)
    cy = np.clip(np.round((LA - (cg[3] + cg[5] / 2)) / cg[5]).astype(int), 0, cu.RasterYSize - 1)
    CU = ca[cy, cx].astype(float)
    m = ok & np.isfinite(CU) & np.isfinite(M["z580"])
    say("## 5. Independent check against a Curie-point depth model\n")
    say(f"Li et al. (2017) global Curie-point depth (magnetite, about 580 deg C) at the same {int(m.sum())} nodes: median "
        f"{np.median(CU[m]):.1f} km (25-75 %: {np.percentile(CU[m], 25):.1f}-{np.percentile(CU[m], 75):.1f} km). The 580 deg C "
        f"depth predicted here: median {np.median(M['z580'][m]):.1f} km ({np.percentile(M['z580'][m], 25):.1f}-"
        f"{np.percentile(M['z580'][m], 75):.1f} km). Spearman correlation between the two maps: "
        f"{spearmanr(CU[m], M['z580'][m]).statistic:+.2f}; between the Curie depth and D_b10: {spearmanr(CU[m], Z[m]).statistic:+.2f}. "
        f"The Curie model carries no node-scale information here, so only the regional medians can be compared. "
        f"With k = {K_ALT} W/m/K (T_ref = {T_ALT:.0f} deg C) the predicted median is {np.median(MALT['z580'][m]):.1f} km. "
        f"The heat-flow calibration with k = {K} and the Curie median therefore bracket the temperature at D_b10 between "
        f"about {T_ALT:.0f} and {T_REF:.0f} deg C, and the 350 deg C isotherm between "
        f"{np.nanmedian(M['z350'] / Z):.2f} and {np.nanmedian(MALT['z350'] / Z):.2f} times D_b10. The heat-flow grid is "
        f"the same in both cases.\n")

say("## 6. Limits\n")
say("* Steady-state, one-dimensional conduction. Not valid next to hydrothermal systems and young intrusions (where the "
    "shallow D_b10 reflects advected heat), nor in the fore-arc cooled by the slab. The 0.025-degree product should not be "
    "converted this way.")
say("* T_c is tied to observations only through one regional constant. Lithology (mafic lower crust, D_b10 deeper than "
    "about 25 km) and pore pressure are not modelled.")
say("* Onshore heat-flow sites over-represent geothermal areas; node medians and the 20-300 mW/m2 window limit, but do "
    "not remove, the resulting upward bias of T_ref.")
say("* Depths are JMA routine depths (1-D velocity model); a systematic depth error of 1-3 km maps into T_ref.")
say("* k is constant; its decrease with temperature would steepen the geotherm at depth and make the 580 deg C depth "
    "shallower than given here.")

# ------------------------------------------------------------------ GeoTIFF
strain_out = np.where(ok, np.exp(LNE_REL) * EPS_REF, np.nan)
bands = [
    ("heat_flow_mWm2", M["qs"], "Surface heat flow of the conductive geotherm through (0, T_s) and (D_b10, T_c), mW/m2"),
    ("gradient_surface_Cperkm", M["g_surf"], "Geothermal gradient at the surface, heat flow / k, deg C/km"),
    ("gradient_mean_Cperkm", M["g_mean"], "Mean gradient between the surface and D_b10, (T_c - T_s) / D_b10, deg C/km"),
    ("depth_350C_km", M["z350"], "Depth of the 350 deg C isotherm on the calibrated geotherm, km"),
    ("depth_580C_km", M["z580"], "Depth of the 580 deg C isotherm (extrapolated below D_b10), km"),
    ("T_at_D_b10_C", M["tc"], "Temperature assigned to D_b10 (depth- and strain-rate-dependent), deg C"),
    ("heat_flow_lo_mWm2", MLO["qs"], f"Heat flow with T_ref at the lower end of its 95 % interval ({T_LO:.0f} deg C)"),
    ("heat_flow_hi_mWm2", MHI["qs"], f"Heat flow with T_ref at the upper end of its 95 % interval ({T_HI:.0f} deg C)"),
    ("depth_350C_shallow_km", MHI["z350"], f"Depth of 350 deg C with T_ref at the upper end of its interval ({T_HI:.0f} deg C)"),
    ("depth_350C_deep_km", MALT["z350"], f"Depth of 350 deg C with k = {K_ALT} W/m/K and T_ref recalibrated to {T_ALT:.0f} deg C"),
    ("D_b10_smoothed_km", Z, f"Input: D_b10 where quality_flag == 0, Gaussian-smoothed (sigma {args.smooth:g} node)"),
    ("strain_rate_nstrain_per_yr", strain_out, "Input: second invariant of the geodetic strain rate (Ohtate et al. 2026), smoothed"),
]
tif_out = os.path.join(args.out, "thermal_from_seismogenic_bottom_Japan_0p1deg.tif")
drv = gdal.GetDriverByName("GTiff")
out = drv.Create(tif_out, nx, ny, len(bands), gdal.GDT_Float32,
                 options=["COMPRESS=ZSTD", "ZSTD_LEVEL=15", "PREDICTOR=3", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"])
out.SetGeoTransform(gt)
srs = osr.SpatialReference(); srs.ImportFromEPSG(4326); out.SetProjection(srs.ExportToWkt())
out.SetMetadata({
    "AREA_OR_POINT": "Point",
    "TITLE": "Heat flow, geothermal gradient and isotherm depths converted from the lower limit of brittle seismicity "
             "beneath the Japanese islands (based on JMA unified catalogue 2001-2024)",
    "METHOD": f"D_b10 taken as an isotherm T_c with 1/T_c = 1/T_ref + (R/Q)[n ln(z/15 km) - ln(eps/eps_ref)] "
              f"(n = {N_EXP:g}, Q = {Q_ACT/1e3:.0f} kJ/mol); steady-state conductive geotherm with k = {K} W/m/K, "
              f"A0 = {A0} microW/m3, hr = {HR} km; T_ref = {T_REF:.0f} deg C (95 % interval {T_LO:.0f}-{T_HI:.0f}) calibrated "
              f"against {len(G)} node medians of onshore heat flow",
    "SOURCE": "Created by F. Ikegami based on the JMA unified hypocenter catalogue "
              "(https://www.data.jma.go.jp/eqev/data/bulletin/hypo.html), the Global Heat Flow Database Release 2021 "
              "(Fuchs, Norden & IHFC 2021, doi:10.5880/fidgeo.2021.014, CC BY 4.0) and the strain-rate field of "
              "Ohtate et al. (2026, doi:10.1186/s40623-026-02475-y, CC BY 4.0). Independent work, not a product of the "
              "Japan Meteorological Agency or of the other data providers",
    "UNITS": "see band descriptions",
})
for i, (name, arr, desc) in enumerate(bands, start=1):
    a = np.where(ok, arr, np.nan).astype(np.float32)
    b = out.GetRasterBand(i)
    b.SetDescription(name); b.SetNoDataValue(float("nan")); b.SetMetadataItem("DESCRIPTION", desc)
    b.WriteArray(a); b.ComputeStatistics(False)
out.FlushCache(); out = None
print("written", tif_out)

# ------------------------------------------------------------------ figures
CREDIT = ("Created by F. Ikegami based on the JMA unified hypocenter catalogue 2001-2024 (Source: Japan Meteorological Agency website, "
          "https://www.data.jma.go.jp/eqev/data/bulletin/hypo.html); independent work, not a JMA product.\n"
          "Heat flow: Fuchs, Norden & IHFC (2021) Global Heat Flow Database Release 2021, doi:10.5880/fidgeo.2021.014 (CC BY 4.0). "
          "Strain rate: Ohtate et al. (2026) Earth Planets Space 78:210, doi:10.1186/s40623-026-02475-y (CC BY 4.0). "
          "Coastline: GSHHG (Wessel & Smith 1996).")
segs = coast_segments(args.coast, BBOX)
ext = [lons[0] - step / 2, lons[-1] + step / 2, lats[-1] - step / 2, lats[0] + step / 2]
VIEW = (128.5, 146.5, 30.0, 46.0)


def map_axes(fig, pos):
    ax = fig.add_subplot(*pos)
    ax.add_collection(LineCollection(segs, colors="k", linewidths=0.35))
    ax.set_xlim(VIEW[0], VIEW[1]); ax.set_ylim(VIEW[2], VIEW[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(38)))
    ax.grid(True, linewidth=0.3, alpha=0.4)
    return ax


fig = plt.figure(figsize=(15, 13.6))
panels = [(M["qs"], "Surface heat flow (mW/m2); circles = observed node medians", 30, 150, "inferno"),
          (M["g_surf"], "Surface geothermal gradient (deg C/km)", 10, 60, "inferno"),
          (M["z350"], "Depth of the 350 deg C isotherm (km)", 0, 30, "viridis_r"),
          (M["tc"], "Temperature assigned to D_b10 (deg C)", T_REF - 60, T_REF + 60, "coolwarm")]
for i, (a, title, vmin, vmax, cmap) in enumerate(panels, start=1):
    ax = map_axes(fig, (2, 2, i))
    im = ax.imshow(np.where(ok, a, np.nan), extent=ext, origin="upper", vmin=vmin, vmax=vmax, cmap=cmap,
                   interpolation="nearest", zorder=0)
    if i == 1:
        ax.scatter(G.lon, G.lat, c=G.q, s=22, vmin=vmin, vmax=vmax, cmap=cmap, edgecolors="w", linewidths=0.5, zorder=3)
    ax.set_title(title, fontsize=11)
    plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
fig.suptitle("Thermal conversion of the base of brittle seismicity (D_b10) based on JMA unified catalogue 2001-2024\n"
             f"conductive geotherm pinned at D_b10; T_ref = {T_REF:.0f} deg C ({T_LO:.0f}-{T_HI:.0f}) calibrated against "
             f"{len(G)} onshore heat-flow nodes; 0.1-degree nodes", fontsize=13)
fig.text(0.01, 0.005, CREDIT, fontsize=7, va="bottom", ha="left")
fig.tight_layout(rect=(0, 0.035, 1, 1))
fig.savefig(os.path.join(args.out, "thermal_maps.png"), dpi=130); plt.close(fig)

qp = G.q.values * np.exp(-res)
fig = plt.figure(figsize=(15, 5.6))
ax = fig.add_subplot(1, 3, 1)
ax.scatter(G.z, G.q, s=14, c="0.3")
zz = np.linspace(max(4.0, G.z.min()), G.z.max(), 100)
for tr, st, lab in [(T_REF, "-", f"T_ref = {T_REF:.0f} deg C"), (T_LO, "--", f"{T_LO:.0f} / {T_HI:.0f} deg C"), (T_HI, "--", None)]:
    ax.plot(zz, q_surface(zz, t_c(zz, 0.0, tr), 14.0, K, A0, HR), "r" + st, lw=1.3, label=lab)
ax.plot(zz, q_surface(zz, 350.0, 14.0, K, A0, HR), "b:", lw=1.3, label="constant 350 deg C")
ax.set_yscale("log"); ax.set_xlabel("D_b10 (km)"); ax.set_ylabel("observed heat flow, node median (mW/m2)")
ax.set_title("Heat flow against D_b10"); ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)
ax = fig.add_subplot(1, 3, 2)
ax.scatter(qp, G.q, s=14, c="0.3"); ax.plot([30, 300], [30, 300], "r-", lw=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("predicted heat flow (mW/m2)"); ax.set_ylabel("observed (mW/m2)")
ax.set_xlim(30, 300); ax.set_ylim(18, 330)
for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_formatter(ScalarFormatter()); axis.set_minor_formatter(NullFormatter())
    axis.set_ticks([30, 50, 100, 200, 300])
ax.set_title(f"Predicted against observed (n = {len(G)}, robust s.d. of ln ratio {mad:.2f})"); ax.grid(True, which="both", alpha=0.3)
ax = fig.add_subplot(1, 3, 3)
ax.hist(boots[:, 0], bins=40, color="0.5"); ax.axvline(T_REF, color="r")
ax.set_xlabel("T_ref (deg C)"); ax.set_title(f"Bootstrap of T_ref ({args.boot} resamples)"); ax.grid(True, alpha=0.3)
fig.text(0.01, 0.005, CREDIT, fontsize=6.5, va="bottom", ha="left")
fig.tight_layout(rect=(0, 0.075, 1, 1))
fig.savefig(os.path.join(args.out, "thermal_calibration.png"), dpi=130); plt.close(fig)

open(os.path.join(args.out, "thermal_conversion.md"), "w", encoding="utf-8", newline="\n").write("\n".join(report) + "\n")
print("wrote", os.path.join(args.out, "thermal_conversion.md"))
