"""FESOM PR #1056 (FCT low-order tendency in flux form): 10-year double-precision
CORE2/JRA55 pair, main 8fdb5e89 (f28dpmain) against the PR merged with it (f28dpflux).

In double precision the PR only reorders arithmetic, so the expected outcome is a
round-off seed that grows into chaotic divergence with no systematic signal. The
script asks three things:
  1. does the field-wise RMS difference grow and saturate at the level of two
     unrelated states of the same model (decorrelation), and how fast;
  2. is any global or regional integral different by more than the difference
     series' own month-to-month scatter allows (trend test on the difference);
  3. does the heat budget close the same way in both runs.

Writes data/clim/fesom_pr1056_dp_eval.npz, a text summary and two figures.
"""
import os
import numpy as np
import xarray as xr

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

BASE = "/albedo/scratch/user/jstreffi/runtime/fesom-2.8"
RUNS = {"main": "f28dpmain", "flux": "f28dpflux"}
YEARS = range(1958, 1968)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CACHE = os.path.join(REPO, "data", "clim", "fesom_pr1056_dp_eval.npz")
SUMMARY = os.path.join(REPO, "data", "clim", "fesom_pr1056_dp_eval.txt")
VCPW = 4.2e6  # src/oce_modules.F90

mesh = xr.open_dataset(f"{BASE}/f28dpmain/run_19580101-19671231/work/fesom.mesh.diag.nc")
lat = mesh.lat.values
lon = mesh.lon.values
area3 = mesh.nod_area.values[:47]                      # (nz1, nod2), 0 below bottom
zi = np.abs(mesh.nz.values)
dz = np.diff(zi)                                       # nominal layer thickness
zmid = 0.5 * (zi[1:] + zi[:-1])
vol3 = area3 * dz[:, None]
area2 = area3[0]
ocean_area = area2.sum()


def kidx(depth):
    return int(np.argmin(np.abs(zmid - depth)))


LEVELS = {"100m": kidx(100), "1000m": kidx(1000), "3000m": kidx(3000)}
BANDS = {"0-700": (0, 700), "700-2000": (700, 2000), "2000-bot": (2000, 1e5)}
band_w = {b: vol3 * ((zmid >= lo) & (zmid < hi))[:, None] for b, (lo, hi) in BANDS.items()}
nh, sh = lat > 0, lat < 0
lab = (lat > 53) & (lat < 65) & (lon > -65) & (lon < -45)
wed = (lat > -72) & (lat < -60) & (lon > -60) & (lon < 30)


def load(run, var, year):
    p = f"{BASE}/{run}/outdata/fesom/{var}.fesom.{year}.nc"
    with xr.open_dataset(p) as d:
        return d[var].values.astype(np.float64)


def wmean(x, w):
    x = np.where(np.isnan(x), 0.0, x)
    return (x * w).sum(axis=tuple(range(1, x.ndim))) / w.sum()


def build():
    out = {}
    for tag, run in RUNS.items():
        ts = {k: [] for k in ["T_vol", "S_vol", "sst", "sss", "ssh", "fh",
                              "iceA_nh", "iceA_sh", "iceV_nh", "iceV_sh",
                              "mld_lab", "mld_wed"] + [f"T_{b}" for b in BANDS]}
        for y in YEARS:
            t = load(run, "temp", y)
            s = load(run, "salt", y)
            ts["T_vol"].append(wmean(t, vol3))
            ts["S_vol"].append(wmean(s, vol3))
            for b in BANDS:
                ts[f"T_{b}"].append(wmean(t, band_w[b]))
            ts["sst"].append(wmean(load(run, "sst", y), area2))
            ts["sss"].append(wmean(load(run, "sss", y), area2))
            ts["ssh"].append(wmean(load(run, "ssh", y), area2))
            ts["fh"].append(wmean(load(run, "fh", y), area2))
            a = load(run, "a_ice", y)
            m = load(run, "m_ice", y)
            ts["iceA_nh"].append((a * area2 * nh).sum(1) / 1e12)
            ts["iceA_sh"].append((a * area2 * sh).sum(1) / 1e12)
            ts["iceV_nh"].append((m * area2 * nh).sum(1) / 1e12)
            ts["iceV_sh"].append((m * area2 * sh).sum(1) / 1e12)
            mld = np.abs(load(run, "MLD2", y))
            ts["mld_lab"].append(wmean(mld, area2 * lab))
            ts["mld_wed"].append(wmean(mld, area2 * wed))
            print(run, y, flush=True)
        for k, v in ts.items():
            out[f"{tag}_{k}"] = np.concatenate(v)

    # field differences: area/volume-weighted RMS per month, and the reference
    # spread of the main run (deseasonalised monthly anomalies, last five years)
    fields2 = ["sst", "sss", "ssh", "a_ice"]
    rms = {f: [] for f in fields2 + list(LEVELS) + ["temp3d"]}
    keep = {f: {"main": [], "flux": []} for f in fields2 + list(LEVELS)}
    ann_sst = {"main": [], "flux": []}
    ann_t3 = {"main": [], "flux": []}
    for y in YEARS:
        for f in fields2:
            a, b = load(RUNS["main"], f, y), load(RUNS["flux"], f, y)
            rms[f].append(np.sqrt(wmean((b - a) ** 2, area2)))
            keep[f]["main"].append(a)
            keep[f]["flux"].append(b)
            if f == "sst":
                ann_sst["main"].append(a.mean(0))
                ann_sst["flux"].append(b.mean(0))
        a, b = load(RUNS["main"], "temp", y), load(RUNS["flux"], "temp", y)
        rms["temp3d"].append(np.sqrt(wmean((b - a) ** 2, vol3)))
        for name, k in LEVELS.items():
            rms[name].append(np.sqrt(wmean((b[:, k] - a[:, k]) ** 2, area3[k])))
            keep[name]["main"].append(a[:, k])
            keep[name]["flux"].append(b[:, k])
        ann_t3["main"].append(np.nan_to_num(a.mean(0)))
        ann_t3["flux"].append(np.nan_to_num(b.mean(0)))
        print("diff", y, flush=True)
    for f, v in rms.items():
        out[f"rms_{f}"] = np.concatenate(v)

    # saturation level: two independent draws from the same distribution differ by
    # sqrt(2)*sigma; sigma from the main run's deseasonalised, detrended monthly
    # anomalies over 1963-1967
    for f in keep:
        w = area2 if f in fields2 else area3[LEVELS[f]]
        x = np.concatenate(keep[f]["main"])[60:]           # (60, nod2)
        x = x.reshape(5, 12, -1)
        x = x - x.mean(0, keepdims=True)
        tt = np.arange(5)[:, None, None] - 2.0
        slope = (x * tt).sum(0, keepdims=True) / (tt ** 2).sum()
        x = (x - slope * tt).reshape(60, -1)
        var = np.nanvar(x, axis=0, ddof=3)
        out[f"sat_{f}"] = np.sqrt(2 * wmean(var[None], w)[0])

    # 1963-1967 mean SST difference and a per-node t-test on annual means
    am = np.array(ann_sst["main"])[5:]
    af = np.array(ann_sst["flux"])[5:]
    out["sst_diff_map"] = af.mean(0) - am.mean(0)
    se = np.sqrt(am.var(0, ddof=1) / 5 + af.var(0, ddof=1) / 5)
    out["sst_t"] = out["sst_diff_map"] / np.where(se > 0, se, np.nan)
    t3m = np.array(ann_t3["main"])[5:]
    t3f = np.array(ann_t3["flux"])[5:]
    d3 = t3f.mean(0) - t3m.mean(0)
    se3 = np.sqrt(t3m.var(0, ddof=1) / 5 + t3f.var(0, ddof=1) / 5)
    t3 = np.where(se3 > 0, d3 / np.where(se3 > 0, se3, 1), 0)
    out["t3_sigfrac_vol"] = (vol3 * (np.abs(t3) > 2.306)).sum() / vol3.sum()
    # zonal mean of the 3-D difference in 2-degree bands
    edges = np.arange(-80, 92, 2)
    zm = np.full((47, len(edges) - 1), np.nan)
    for j in range(len(edges) - 1):
        sel = (lat >= edges[j]) & (lat < edges[j + 1])
        w = vol3[:, sel]
        ws = w.sum(1)
        zm[:, j] = np.where(ws > 0, (d3[:, sel] * w).sum(1) / np.where(ws > 0, ws, 1), np.nan)
    out["zm_dT"] = zm
    out["zm_lat"] = 0.5 * (edges[1:] + edges[:-1])
    out["zmid"] = zmid
    np.savez(CACHE, **out)
    return dict(np.load(CACHE))


def trend_test(d):
    """OLS slope of a monthly difference series per decade, with a lag-1
    autocorrelation-corrected standard error."""
    t = np.arange(d.size) / 120.0
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    r = d - A @ coef
    r1 = np.corrcoef(r[:-1], r[1:])[0, 1]
    neff = max(3.0, d.size * (1 - r1) / (1 + r1))
    se = np.sqrt((r ** 2).sum() / (neff - 2) / ((t - t.mean()) ** 2).sum())
    return coef[0], se, neff


def heat_budget(D, tag):
    """HC change between the first and last monthly means against the surface flux
    integrated over the same interval; returns the residual in W/m2 for each sign
    convention of fh."""
    hc = D[f"{tag}_T_vol"] * vol3.sum() * VCPW
    fh = D[f"{tag}_fh"]
    days = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] * 10, float)
    days[[12 * 2 + 1, 12 * 6 + 1]] += 1                  # 1960, 1964 leap Februaries
    sec = days * 86400
    # from mid-first-month to mid-last-month: half of both end months, all between
    span = 0.5 * sec[0] + sec[1:-1].sum() + 0.5 * sec[-1]
    q = (0.5 * fh[0] * sec[0] + (fh[1:-1] * sec[1:-1]).sum() + 0.5 * fh[-1] * sec[-1]) / span
    dhc = (hc[-1] - hc[0]) / span / ocean_area
    return dhc, q


def summarise(D):
    lines = []
    P = lines.append
    P("FESOM PR #1056, double precision, CORE2/JRA55 1958-1967, f28dpflux - f28dpmain")
    P("")
    P("1. RMS difference, against sqrt(2) x sigma of the main run's monthly anomalies (free variability;\n"
      "   forced runs share the forced part of it and cannot reach this level)")
    P(f"   {'field':8s} {'month 1':>10s} {'year 1':>10s} {'year 10':>10s} {'sqrt2 sig':>11s} {'yr10/ref':>9s}")
    for f in ["sst", "sss", "ssh", "a_ice", "100m", "1000m", "3000m"]:
        r = D[f"rms_{f}"]
        sat = float(D[f"sat_{f}"])
        P(f"   {f:8s} {r[0]:10.3e} {r[:12].mean():10.3e} {r[-12:].mean():10.3e} {sat:11.3e} {r[-12:].mean() / sat:9.2f}")
    r = D["rms_temp3d"]
    P(f"   temp3d   {r[0]:10.3e} {r[:12].mean():10.3e} {r[-12:].mean():10.3e}  (volume-weighted, no saturation reference)")
    P("")
    P("2. Global and regional integrals: mean difference 1963-67 and trend of the monthly difference")
    P(f"   {'quantity':10s} {'main 1958-67 drift/dec':>23s} {'diff 63-67':>11s} {'diff trend/dec':>15s} {'2 se':>10s} {'n_eff':>6s}")
    for k, unit in [("T_vol", "K"), ("T_0-700", "K"), ("T_700-2000", "K"), ("T_2000-bot", "K"),
                    ("S_vol", "psu"), ("sst", "K"), ("sss", "psu"), ("ssh", "m"), ("fh", "W/m2"),
                    ("iceA_nh", "1e6km2"), ("iceA_sh", "1e6km2"), ("iceV_nh", "1e3km3"),
                    ("iceV_sh", "1e3km3"), ("mld_lab", "m"), ("mld_wed", "m")]:
        a, b = D[f"main_{k}"], D[f"flux_{k}"]
        d = b - a
        drift, _, _ = trend_test(a)
        s, se, ne = trend_test(d)
        flag = "  *" if abs(s) > 2 * se else ""
        P(f"   {k:10s} {drift:+23.3e} {d[60:].mean():+11.3e} {s:+15.3e} {2 * se:10.3e} {ne:6.0f}  {unit}{flag}")
    P("   (* = trend of the difference exceeds 2 standard errors)")
    P("")
    P("3. Heat budget over 10 years (W/m2 over the ocean area)")
    for tag in ["main", "flux"]:
        dhc, q = heat_budget(D, tag)
        P(f"   {tag:5s} dHC/dt = {dhc:+.4f}   mean fh = {q:+.4f}   dHC/dt - fh = {dhc - q:+.4f}   dHC/dt + fh = {dhc + q:+.4f}")
    dm, qm = heat_budget(D, "main")
    df, qf = heat_budget(D, "flux")
    P(f"   difference flux - main: dHC/dt {df - dm:+.5f}, fh {qf - qm:+.5f}")
    P("   (nominal layer thicknesses; zstar volume changes are ignored, the same in both runs)")
    P("")
    # no per-node t-test: both runs see the same forcing, so they are paired, not
    # independent samples, and a two-sample test is meaningless here
    P("4. 1963-67 mean difference fields")
    P(f"   SST difference map: area-weighted mean {wmean(D['sst_diff_map'][None], area2)[0]:+.4f} K, "
      f"RMS {np.sqrt(wmean(D['sst_diff_map'][None] ** 2, area2)[0]):.4f} K")
    txt = "\n".join(lines)
    with open(SUMMARY, "w") as fh:
        fh.write(txt + "\n")
    print(txt)


def plot(D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
    C1, C2 = "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                         "grid.linewidth": 0.6, "lines.linewidth": 1.6})
    yrs = 1958 + (np.arange(120) + 0.5) / 12

    # Figure 1: divergence and integrals
    fig, ax = plt.subplots(3, 3, figsize=(13, 9.5), constrained_layout=True)
    for a, f, lab, unit in [(ax[0, 0], "sst", "SST", "K"), (ax[0, 1], "1000m", "T at 1000 m", "K"),
                            (ax[0, 2], "a_ice", "Sea-ice concentration", "")]:
        r = D[f"rms_{f}"]
        sat = float(D[f"sat_{f}"])
        a.semilogy(yrs, r, color=C1)
        a.axhline(sat, color=INK2, ls="--", lw=1)
        a.text(yrs[-1], sat * 1.15, "free variability, sqrt(2) sigma", ha="right", va="bottom", color=INK2, fontsize=8)
        a.set_title(f"RMS difference, {lab}", loc="left", color=INK)
        a.set_ylabel(unit)
    for a, k, lab, unit, sc in [(ax[1, 0], "T_vol", "Global mean temperature", "mK", 1e3),
                                (ax[1, 1], "S_vol", "Global mean salinity", "1e-6 psu", 1e6),
                                (ax[1, 2], "sst", "Global mean SST", "mK", 1e3),
                                (ax[2, 0], "iceA_nh", "Arctic ice area", "1e3 km2", 1e3),
                                (ax[2, 1], "iceA_sh", "Antarctic ice area", "1e3 km2", 1e3),
                                (ax[2, 2], "mld_lab", "Labrador Sea MLD", "m", 1)]:
        d = (D[f"flux_{k}"] - D[f"main_{k}"]) * sc
        a.plot(yrs, d, color=C1)
        a.axhline(0, color=INK2, lw=0.8)
        a.set_title(f"{lab}, PR - main", loc="left", color=INK)
        a.set_ylabel(unit)
    fig.suptitle("FESOM PR #1056 against main, double precision, CORE2/JRA55, 1958-1967",
                 x=0.01, ha="left", color=INK, fontsize=11)
    fig.savefig(os.path.join(REPO, "plots", "fesom_pr1056_dp_divergence.png"), dpi=130)
    plt.close(fig)

    # Figure 2: where the 1963-67 SST difference is, and the zonal-mean T difference
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True,
                           gridspec_kw={"width_ratios": [1.6, 1]})
    d = D["sst_diff_map"]
    lim = np.nanpercentile(np.abs(d), 99)
    sc = ax[0].scatter(lon, lat, c=d, s=0.4, cmap="RdBu_r", vmin=-lim, vmax=lim, rasterized=True)
    ax[0].set_xlim(-180, 180)
    ax[0].set_ylim(-80, 90)
    ax[0].grid(False)
    ax[0].set_title("SST, 1963-67 mean, PR - main", loc="left", color=INK)
    fig.colorbar(sc, ax=ax[0], label="K", shrink=0.85)
    zm = D["zm_dT"]
    lim = np.nanpercentile(np.abs(zm), 99)
    pc = ax[1].pcolormesh(D["zm_lat"], D["zmid"], zm, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="nearest")
    ax[1].set_ylim(5500, 0)
    ax[1].grid(False)
    ax[1].set_xlabel("latitude")
    ax[1].set_ylabel("depth (m)")
    ax[1].set_title("Zonal-mean temperature, 1963-67, PR - main", loc="left", color=INK)
    fig.colorbar(pc, ax=ax[1], label="K", shrink=0.85)
    fig.savefig(os.path.join(REPO, "plots", "fesom_pr1056_dp_maps.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    D = dict(np.load(CACHE)) if os.path.exists(CACHE) else build()
    summarise(D)
    plot(D)
