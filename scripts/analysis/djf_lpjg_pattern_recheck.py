"""Re-check of the WITH-vs-WITHOUT-LPJG DJF maps: does the cold follow the forest loss?

The figure plots/djf_roughness_inversion_maps.png shows six DJF differences (LPJG Baseline
minus the no-LPJG CORE3 reference) over the boreal NH.  The report reads them as one chain:
forest lost -> smoother surface -> stronger inversion -> colder screen.  By eye the forest
loss is near-uniform around the whole circumboreal belt while the cooling is concentrated in
eastern Siberia, which a single-predictor regression on a near-constant field would not
reveal.  This script puts numbers on that: a regional table, the spread of the predictor, and
the two controls the pattern argument needs -- sea ice (the two runs are separately coupled,
so their Arctic ice is a free variable, not a held-fixed boundary) and LW down.

Usage:  python3 scripts/analysis/djf_lpjg_pattern_recheck.py [Y0 Y1]   (default 1370 1379)
"""
import sys, glob, warnings
import numpy as np, xarray as xr
warnings.filterwarnings("ignore")

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1370, 1379)
years = list(range(Y0, Y1 + 1))
WITHOUT = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
WITH    = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs"
VN = {"2t": "2t", "skt": "skt", "tcc": "tcc", "strd": "strd", "cvh": "cvh", "ci": "ci", "pl_t": "t"}


def files(run, var):
    base = WITHOUT if run == "WITHOUT" else WITH
    out = []
    for y in years:
        if run == "WITHOUT":
            # the no-LPJG reference wrote cvh only daily; everything else is monthly
            c = ([f"atm_remapped_1m_pl_t_1m_pl_{y}-{y}.nc"] if var == "pl_t"
                 else [f"atm_remapped_1m_{var}_1m_{y}-{y}.nc", f"atm_remapped_1d_{var}_1d_{y}-{y}.nc"])
        else:
            c = [f"atm_remapped_1m_{var}_{y}-{y}.nc"]
        for cc in c:
            g = glob.glob(f"{base}/{cc}")
            if g: out += g; break
    return sorted(out)


def djf(run, var, plev=None):
    fs = files(run, var)
    if not fs: return None
    ds = xr.open_mfdataset(fs, combine="by_coords", use_cftime=True)
    da = ds[VN[var]]
    if plev is not None: da = da.sel(pressure_levels=plev)
    td = "time_counter" if "time_counter" in da.dims else "time"
    m = da[td].dt.month
    return da.where((m == 12) | (m == 1) | (m == 2), drop=True).mean(td).load()


V = ("2t", "skt", "tcc", "strd", "cvh", "ci", "pl_t")
D = {r: {v: djf(r, v, plev=(92500.0 if v == "pl_t" else None)) for v in V} for r in ("WITHOUT", "WITH")}
dif = lambda v: (D["WITH"][v] - D["WITHOUT"][v])
lat = D["WITH"]["2t"]["lat"].values; lon = D["WITH"]["2t"]["lon"].values
lm = xr.open_dataset(glob.glob(f"{WITH}/atm_remapped_1m_lsm_{Y0}-{Y0}.nc")[0])["lsm"]
if "time_counter" in lm.dims: lm = lm.isel(time_counter=0)
lm = lm.reset_coords(drop=True).values
LAT = np.broadcast_to(lat[:, None], lm.shape); LON = np.broadcast_to(lon[None, :], lm.shape)
W = np.cos(np.deg2rad(LAT))
LND, OCN = lm > 0.5, lm <= 0.5

dcvh = dif("cvh").values; dt2 = dif("2t").values; dci = dif("ci").values
dinv = ((D["WITH"]["pl_t"] - D["WITH"]["2t"]) - (D["WITHOUT"]["pl_t"] - D["WITHOUT"]["2t"])).values
dsfc = ((D["WITH"]["2t"] - D["WITH"]["skt"]) - (D["WITHOUT"]["2t"] - D["WITHOUT"]["skt"])).values
dcld = dif("tcc").values
s = lambda a: a / 3600.0 if abs(np.nanmean(a)) > 1000 else a
dlw = (s(D["WITH"]["strd"].values) - s(D["WITHOUT"]["strd"].values))

av = lambda f, k: float(np.average(f[k], weights=W[k])) if k.sum() else np.nan
inbox = lambda la, lb, lo, hi: (LAT >= la) & (LAT <= lb) & (((LON - lo) % 360) <= ((hi - lo) % 360))

REG = [("E Siberia   90-150E", 55, 70, 90, 150), ("C Siberia   60-90E", 55, 70, 60, 90),
       ("Far E Sib  150-180E", 55, 70, 150, 180), ("W Russia/Fenno 10-60E", 55, 70, 10, 60),
       ("Canada     240-300E", 55, 70, 240, 300), ("Alaska     190-230E", 55, 70, 190, 230)]

print(__doc__.split("Usage:")[0])
print(f"DJF {Y0}-{Y1}, LPJG Baseline minus no-LPJG CORE3 reference\n")
print(f'{"region (55-70N land)":<24}{"dcvh":>8}{"dT2m":>8}{"dinv":>8}{"dsfc":>8}{"dcld":>8}{"dLW":>8}{"cells":>7}')
for name, la, lb, lo, hi in REG:
    k = LND & inbox(la, lb, lo, hi)
    print(f'{name:<24}{av(dcvh,k):8.3f}{av(dt2,k):8.2f}{av(dinv,k):8.2f}{av(dsfc,k):8.2f}'
          f'{av(dcld,k):8.3f}{av(dlw,k):8.1f}{int(k.sum()):7d}')
kb = LND & (LAT >= 55) & (LAT <= 70)
print(f'{"BOREAL BELT 55-70N":<24}{av(dcvh,kb):8.3f}{av(dt2,kb):8.2f}{av(dinv,kb):8.2f}'
      f'{av(dsfc,kb):8.2f}{av(dcld,kb):8.3f}{av(dlw,kb):8.1f}{int(kb.sum()):7d}')

print("\nHow much does the predictor actually vary over boreal land 50-75N?")
kk = LND & (LAT >= 50) & (LAT <= 75) & np.isfinite(dcvh)
q = np.percentile(dcvh[kk], [5, 25, 50, 75, 95])
print(f'  dcvh  mean {av(dcvh,kk):+.3f}  sd {float(np.std(dcvh[kk])):.3f}   '
      f'p5/25/50/75/95 = {q[0]:+.2f} {q[1]:+.2f} {q[2]:+.2f} {q[3]:+.2f} {q[4]:+.2f}')
print(f'  dT2m  mean {av(dt2,kk):+.2f} K  sd {float(np.std(dt2[kk])):.2f} K')


def wcorr(a, b, k):
    x, y, w = a[k], b[k], W[k]
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.average((x - mx) * (y - my), weights=w) /
                 np.sqrt(np.average((x - mx) ** 2, weights=w) * np.average((y - my) ** 2, weights=w)))


print("\nWeighted correlations over boreal land 50-75N (the report quotes the first three):")
for nm, a, b in [("dcvh vs dinv", dcvh, dinv), ("dcvh vs dsfc", dcvh, dsfc), ("dT2m vs dinv", dt2, dinv),
                 ("dcvh vs dT2m", dcvh, dt2), ("dcvh vs dcloud (control)", dcvh, dcld),
                 ("dLW  vs dT2m", dlw, dt2), ("dT2m vs dsfc", dt2, dsfc)]:
    print(f"  {nm:<26}{wcorr(a,b,kk):+.3f}")

print("\nThe sea ice the two runs did NOT hold fixed (they are separately coupled):")
for nm, la, lb, lo, hi in [("Arctic 70-90N", 70, 90, 0, 359.9), ("Okhotsk/Bering 45-65N,140-200E", 45, 65, 140, 200),
                           ("Barents/Kara 65-80N,20-90E", 65, 80, 20, 90), ("Labrador/Baffin 55-75N,280-320E", 55, 75, 280, 320)]:
    k = OCN & inbox(la, lb, lo, hi)
    print(f'  {nm:<34} dci {av(dci,k):+.3f}   dT2m {av(dt2,k):+.2f} K   dLW {av(dlw,k):+6.1f} W/m2')
