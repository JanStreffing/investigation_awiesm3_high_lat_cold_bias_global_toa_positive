"""Constrained refit: hold winter saturation, then get as close to spring as possible.

WHY THIS EXISTS.  The first fit (fit_depletion_curve.py) minimised weighted RMSE over
(depth, density) bins with sqrt(n) weights.  That let the thousands of saturated
midwinter surveys dominate and left the melt corner under-fitted: forest May came out
0.950 (b=0.20) or 0.969 (b=1) against an observed 0.931 -- and the b=1 variant is
FURTHER from observation than the as-released ramp's 0.959, i.e. worse than doing
nothing in the one season the scheme exists to improve.

The two seasons are not equally negotiable, so they should not be equally weighted:
  * WINTER saturation is a hard constraint.  Snow courses give DJF cover 0.9995 with
    99.74% of surveys at exactly 10/10, and failing to reproduce it is what cost the
    tanh -15 to -20 K of soil temperature against 174676 station observations.
  * SPRING is the objective.  It is what the scheme was introduced to supply, and it
    is where every candidate so far is wrong.

So: constrain the predicted Oct-Feb climatology to stay saturated, then minimise the
error over Mar-May.  Reported as a frontier, not a single answer, because the
trade-off is the decision and it belongs to the operator.

Parameters remain (d_c0, m, b) of
    SCF = min(1, (d/d_c)**b),   d_c = d_c200*(rho/200)**m
a pure function of depth and density -- the two snow variables SURFBC_CTL already
receives.  No calendar, no new prognostic, nothing that pins it to today's climate.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

OBSD = '/work/ab0246/a270092/obs/RIHMI-WDC/data'
MON = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
TYPES = {0: 'field (pole)', 1: 'forest (les)'}
WINTER = (9, 10, 11, 0, 1)        # Oct-Feb: must stay saturated
SPRING = (2, 3, 4)                # Mar-May: the objective

with xr.open_dataset(f'{OBSD}/RIHMI-WDC_snmar.nc') as ds:
    la, lo = ds['lat'].values, ds['lon'].values
    k = (la >= 55) & (la <= 75) & (lo >= 60) & (lo <= 180)
    COV = ds['fraction_of_the_snow_course_covered_by_snow'].isel(station=k).values / 10.0
    DEP = ds['snow_depth_mean'].isel(station=k).values / 100.0
    RHO = ds['snow_density'].isel(station=k).values * 1000.0
    tm = ds['time'].values
MO = tm.astype('datetime64[M]').astype(int) % 12

DC0 = np.arange(0.004, 0.40, 0.002)   # d_c at rho=200, decorrelated from m
MM = np.arange(0.0, 5.02, 0.05)
BB = np.arange(0.02, 1.52, 0.02)

print(__doc__)
print('=' * 100)

for t in (0, 1):
    c, d, r = COV[:, :, t], DEP[:, :, t], RHO[:, :, t]
    m3 = np.repeat(MO[None, :], c.shape[0], axis=0)
    ok = np.isfinite(c) & np.isfinite(d) & np.isfinite(r) & (d > 0) & (r > 30) & (r < 600)
    c, d, r, mm = c[ok], d[ok], r[ok], m3[ok]

    # observed monthly climatology and the per-month survey weights
    obs = np.array([c[mm == i].mean() if (mm == i).sum() >= 30 else np.nan for i in range(12)])
    nmo = np.array([(mm == i).sum() for i in range(12)])

    # Collapse surveys onto a fine (d, rho) grid once, so the parameter scan is a
    # matrix product rather than 20k function evaluations per candidate.
    de = np.concatenate([np.linspace(0.005, 0.3, 40), np.linspace(0.32, 1.6, 20)])
    re = np.linspace(60, 600, 46)
    di = np.clip(np.digitize(d, de) - 1, 0, len(de) - 2)
    ri = np.clip(np.digitize(r, re) - 1, 0, len(re) - 2)
    bidx = di * (len(re) - 1) + ri
    nb = (len(de) - 1) * (len(re) - 1)
    dc = 0.5 * (de[:-1] + de[1:]); rc = 0.5 * (re[:-1] + re[1:])
    DGRID = np.repeat(dc, len(rc)); RGRID = np.tile(rc, len(dc))
    CNT = np.zeros((12, nb))
    for i in range(12):
        s = mm == i
        if s.sum():
            CNT[i] = np.bincount(bidx[s], minlength=nb)
    ROWSUM = CNT.sum(axis=1)

    print(f'\n\n### TYPE {t+1} -- {TYPES[t]}, {c.size} surveys\n')

    # scan
    rows = []
    RG = (RGRID / 200.0)
    for b in BB:
        for m in MM:
            dcv = np.outer(DC0, RG ** m)                       # (nDC0, nb)
            scf = np.clip((DGRID[None, :] / dcv) ** b, 0, 1)   # (nDC0, nb)
            pred = (scf @ CNT.T) / np.maximum(ROWSUM, 1)[None, :]   # (nDC0, 12)
            wmin = np.nanmin(pred[:, list(WINTER)], axis=1)
            serr = np.sqrt(np.nanmean((pred[:, list(SPRING)] - obs[list(SPRING)]) ** 2, axis=1))
            for i in range(len(DC0)):
                rows.append((wmin[i], serr[i], DC0[i], m, b, pred[i]))

    def best(thr):
        cand = [x for x in rows if x[0] >= thr]
        return min(cand, key=lambda x: x[1]) if cand else None

    print(f'  FRONTIER -- how much spring accuracy each winter-saturation floor costs\n')
    print(f'  {"winter floor":>13s}{"d_c0":>8s}{"m":>7s}{"b":>7s}{"spring RMSE":>13s}'
          + ''.join(f'{MON[i]:>8s}' for i in (9, 11, 2, 3, 4)))
    print(f'  {"OBSERVED":>13s}{"":8s}{"":7s}{"":7s}{"":13s}'
          + ''.join(f'{obs[i]:8.3f}' for i in (9, 11, 2, 3, 4)))
    keep = {}
    for thr in (0.0, 0.980, 0.990, 0.995, 0.998, 0.999):
        x = best(thr)
        if x is None:
            print(f'  {thr:13.3f}   -- no parameter set reaches this floor --'); continue
        keep[thr] = x
        print(f'  {thr:13.3f}{x[2]:8.3f}{x[3]:7.2f}{x[4]:7.2f}{x[1]:13.4f}'
              + ''.join(f'{x[5][i]:8.3f}' for i in (9, 11, 2, 3, 4)))

    # incumbents on the same yardstick
    def clim(f):
        v = f(DGRID, RGRID)
        return (v @ CNT.T) / np.maximum(ROWSUM, 1)
    rel = clim(lambda dd, rr: np.clip(10 * dd, 0, 1))
    tanh = clim(lambda dd, rr: np.tanh(dd / np.maximum(2.5 * 0.016 * (rr / 100) ** 1.6, 1e-6)))
    for nm, p in (('as-released', rel), ('tanh (N2)', tanh)):
        e = np.sqrt(np.nanmean((p[list(SPRING)] - obs[list(SPRING)]) ** 2))
        print(f'  {nm:>13s}{"":8s}{"":7s}{"":7s}{e:13.4f}'
              + ''.join(f'{p[i]:8.3f}' for i in (9, 11, 2, 3, 4)))

    if 0.998 in keep:
        x = keep[0.998]
        print(f'\n  RECOMMENDED (winter floor 0.998): d_c = {x[2]:.3f}*(rho/200)^{x[3]:.2f}, b = {x[4]:.2f}')
        print(f'  {"month":8s}' + ''.join(f'{MON[i]:>8s}' for i in (8,9,10,11,0,1,2,3,4)))
        print(f'  {"obs":8s}' + ''.join(f'{obs[i]:8.3f}' for i in (8,9,10,11,0,1,2,3,4)))
        print(f'  {"fit":8s}' + ''.join(f'{x[5][i]:8.3f}' for i in (8,9,10,11,0,1,2,3,4)))
        print(f'  {"n":8s}' + ''.join(f'{nmo[i]:8d}' for i in (8,9,10,11,0,1,2,3,4)))

        # The metric that exposed the defect: fraction reaching EXACTLY complete
        # cover.  tanh scores 0 here by construction at every month; the whole point
        # of min(1,.) is that this must track the observed f(10/10).
        dcv = x[2] * (r / 200.0) ** x[3]
        pf = (np.clip((d / dcv) ** x[4], 0, 1) >= 0.999)
        pt = (np.tanh(d / np.maximum(2.5 * 0.016 * (r / 100) ** 1.6, 1e-6)) >= 0.999)
        pr = (np.clip(10 * d, 0, 1) >= 0.999)
        print(f'\n  f(complete cover) -- observed vs predicted')
        print(f'  {"obs 10/10":8s}' + ''.join(
            f'{(c[mm==i]>=0.999).mean() if nmo[i] else np.nan:8.3f}' for i in (8,9,10,11,0,1,2,3,4)))
        for nm, msk in (('fit', pf), ('released', pr), ('tanh', pt)):
            print(f'  {nm:8s}' + ''.join(
                f'{msk[mm==i].mean() if nmo[i] else np.nan:8.3f}' for i in (8,9,10,11,0,1,2,3,4)))

print("""

  READING IT.  Each row is the best spring the curve can do WITHOUT letting any
  Oct-Feb month fall below the stated floor.  If spring RMSE barely rises as the floor
  is tightened from 0.98 to 0.999, the two objectives are compatible and we can have
  both -- correct winter AND a real melt signal.  If it rises steeply, the form itself
  is the limit and no parameter choice gives both, which would be a reason to revisit
  the functional form rather than keep scanning it.""")
