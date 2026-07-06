#!/usr/bin/env python3
"""
Spatial correlations among the DJF WITH-minus-WITHOUT LPJG difference fields,
over boreal land. Tests how tightly the causal chain links:
  forest loss (dcvh) -> cold (dT2m) -> stronger inversion (dInv) -> decoupling (dSfc)
and confirms cloud (dCloud) is NOT part of it (r ~ 0).

Also regresses the two lever slopes:  dT2m per unit dcvh,  dInv per unit dcvh.

Usage: python correlate_djf_roughness.py [Y0 Y1]   (default 1350 1359)
Out:   /work/bb1469/a270092/eval/plots/djf_roughness_scatter.png
"""
import sys, glob, warnings
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1350, 1359)
years = list(range(Y0, Y1 + 1))
WITHOUT = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
WITH    = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs"
VN = {"2t":"2t","skt":"skt","tcc":"tcc","strd":"strd","cvh":"cvh","pl_t":"t"}

def files(run, var):
    base = WITHOUT if run == "WITHOUT" else WITH
    out = []
    for y in years:
        if run == "WITHOUT":
            cands = ([f"atm_remapped_1m_pl_t_1m_pl_{y}-{y}.nc"] if var == "pl_t"
                     else [f"atm_remapped_1m_{var}_1m_{y}-{y}.nc",
                           f"atm_remapped_1d_{var}_1d_{y}-{y}.nc"])
        else:
            cands = [f"atm_remapped_1m_{var}_{y}-{y}.nc"]
        for c in cands:
            g = glob.glob(f"{base}/{c}")
            if g: out += g; break
    return sorted(out)

def djf(run, var, plev=None):
    fs = files(run, var)
    ds = xr.open_mfdataset(fs, combine="by_coords", use_cftime=True)
    da = ds[VN[var]]
    if plev is not None: da = da.sel(pressure_levels=plev)
    td = "time_counter" if "time_counter" in da.dims else "time"
    m = da[td].dt.month
    v = da.where((m == 12) | (m == 1) | (m == 2), drop=True).mean(td).load()
    if var == "strd":  # accumulated J/m2 -> W/m2 (3600 s both runs, verified)
        v = v / 3600.0
    return v

D = {r: {v: djf(r, v, plev=(92500.0 if v == "pl_t" else None))
         for v in ("2t","skt","tcc","strd","cvh","pl_t")} for r in ("WITHOUT","WITH")}

d = lambda v: D["WITH"][v] - D["WITHOUT"][v]
fields = {
    "dcvh":   d("cvh"),
    "dT2m":   d("2t"),
    "dInv":   (D["WITH"]["pl_t"] - D["WITH"]["2t"]) - (D["WITHOUT"]["pl_t"] - D["WITHOUT"]["2t"]),
    "dSfcInv":(D["WITH"]["2t"] - D["WITH"]["skt"]) - (D["WITHOUT"]["2t"] - D["WITHOUT"]["skt"]),
    "dCloud": d("tcc"),
    "dLWdn":  d("strd"),
}
lat = D["WITH"]["2t"]["lat"]; lon = D["WITH"]["2t"]["lon"]

# boreal land mask
lm = None
lmf = glob.glob(f"{WITH}/atm_remapped_1m_lsm_{Y0}-{Y0}.nc")
if lmf:
    lm = xr.open_dataset(lmf[0])["lsm"]
    if "time_counter" in lm.dims: lm = lm.isel(time_counter=0)
    lm = lm.reset_coords(drop=True)
boreal = (lat >= 50) & (lat <= 80)
mask = boreal & ((lm > 0.5) if lm is not None else True)
W = np.cos(np.deg2rad(lat)).broadcast_like(fields["dT2m"])

def flat(da):
    return da.where(mask).values.ravel()
w = flat(W);
keys = list(fields.keys())
cols = {k: flat(fields[k]) for k in keys}
good = np.isfinite(w)
for k in keys: good &= np.isfinite(cols[k])
w = w[good]; cols = {k: cols[k][good] for k in keys}
N = good.sum()

def wcorr(x, y):
    mx = np.average(x, weights=w); my = np.average(y, weights=w)
    cov = np.average((x-mx)*(y-my), weights=w)
    sx = np.sqrt(np.average((x-mx)**2, weights=w)); sy = np.sqrt(np.average((y-my)**2, weights=w))
    return cov/(sx*sy)

def wslope(x, y):  # y on x
    mx = np.average(x, weights=w); my = np.average(y, weights=w)
    return np.average((x-mx)*(y-my), weights=w) / np.average((x-mx)**2, weights=w)

print(f"=== DJF boreal-land ({Y0}-{Y1}), area-weighted spatial correlations, N={N} cells ===\n")
print("        " + " ".join(f"{k:>8s}" for k in keys))
for a in keys:
    row = " ".join(f"{wcorr(cols[a], cols[b]):8.2f}" for b in keys)
    print(f"{a:8s}{row}")

print("\n=== Lever regressions (over boreal land) ===")
s1 = wslope(cols["dcvh"], cols["dT2m"]);  print(f"  dT2m  = {s1:+.2f} K per unit dcvh   (=> restoring cvh by +0.4 warms screen-T by {+0.4*s1:+.2f} K)")
s2 = wslope(cols["dcvh"], cols["dInv"]);  print(f"  dInv  = {s2:+.2f} K per unit dcvh   (=> restoring cvh by +0.4 changes inversion by {+0.4*s2:+.2f} K, i.e. weaker)")
s3 = wslope(cols["dInv"], cols["dT2m"]);  print(f"  dT2m  = {s3:+.2f} K per K of dInv    (inversion<->screen-T coupling)")
r_cld = wcorr(cols["dcvh"], cols["dCloud"]); print(f"\n  cloud check: r(dcvh, dCloud) = {r_cld:+.2f}  (near 0 => not a cloud effect)")

# scatter figure
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
def scat(a, xk, yk, xl, yl):
    x, y = cols[xk], cols[yk]
    a.scatter(x, y, s=3, alpha=0.25, c=w, cmap="viridis")
    m = wslope(x, y); mx = np.average(x, weights=w); my = np.average(y, weights=w)
    xs = np.array([x.min(), x.max()]); a.plot(xs, my + m*(xs-mx), "r-", lw=1.5)
    a.set_xlabel(xl); a.set_ylabel(yl); a.axhline(0, c="k", lw=.3); a.axvline(0, c="k", lw=.3)
    a.set_title(f"r = {wcorr(x, y):+.2f},  slope = {m:+.2f}")
scat(ax[0], "dcvh", "dT2m",  "Δcvh (forest cover)", "ΔT2m [K]")
scat(ax[1], "dcvh", "dInv",  "Δcvh (forest cover)", "Δ inversion T925−T2m [K]")
scat(ax[2], "dInv", "dT2m",  "Δ inversion T925−T2m [K]", "ΔT2m [K]")
fig.suptitle(f"DJF boreal-land WITH−WITHOUT LPJG ({Y0}–{Y1}): forest loss ↔ cold ↔ inversion "
             f"(colour = cos-lat weight)", fontsize=12)
fig.tight_layout(rect=[0,0,1,0.94])
OUT = "/work/bb1469/a270092/eval/plots/djf_roughness_scatter.png"
fig.savefig(OUT, dpi=130, bbox_inches="tight"); print("\nwrote", OUT)
