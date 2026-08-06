"""Fit a snow-cover depletion curve to RIHMI snow-course surveys.

WHY THIS EXISTS.  tanh(d/L) is dead: its range is open at 1, so complete snow cover
is not representable at ANY (z0, rho_new, m), while Russian snow-course surveys put
Siberian DJF cover at exactly 10/10 in 14331 of 14369 cases (99.74%).  Running it
cost -15 to -20 K of winter soil temperature against 174676 station observations,
while the as-released clipping ramp is right to +0.8 K.  But the as-released ramp
min(1, 10*d) has no spring depletion at all, which is the thing the scheme was
introduced to supply.  We need a form that keeps BOTH.

The observations say the density term is mandatory, not optional:
    October   9.1 cm, rho 131  ->  cover 0.989
    May      30.3 cm, rho 285  ->  cover 0.898
Three times the depth and LOWER cover.  No function of depth alone can do that.

THE FORM FITTED HERE
    SCF = min(1, (d / d_c)**b),    d_c = d_c0 * (rho/100)**m
  * reaches EXACTLY 1 at d = d_c, so complete cover is in the range   <- fixes the defect
  * reduces to the as-released ramp at d_c0 = 0.1, m = 0, b = 1       <- nests the incumbent
  * d_c grows with density, so an old dense spring pack patches out
    at a depth where a fresh autumn pack of the same depth does not   <- supplies the spring melt
Fitted separately for course type 1 (field) and 2 (forest), which is where the
roughness dependence should come from -- by vegetation type, not by region.

TWO CAVEATS, both stated with the numbers rather than buried.
  1. SAMPLING.  Surveys are only run when snow is present (no d=0 records), so the
     low-cover tail is under-sampled and the fit is better constrained at the
     saturation end than at melt-out.  That is the right way round for our problem.
  2. SCALE.  A snow course is 1-2 km; a TCO95 grid box is ~100 km.  Sub-grid
     heterogeneity grows with scale, so at model scale the true depletion is at
     least as strong as fitted here -- d_c fitted on courses is a LOWER bound.
     The saturation constraint is scale-free though: midwinter Siberia is
     completely covered at 1 km and at 100 km alike, which is the part that broke.
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
BOX = dict(lat0=55, lat1=75, lon0=60, lon1=180)


# ------------------------------------------------------------------ the forms --
def scf_new(d, rho, p):
    """min(1, (d/d_c)^b), d_c = d_c0*(rho/100)^m -- reaches 1 exactly."""
    dc0, m, b = p
    dc = np.maximum(dc0 * (rho / 100.0) ** m, 1e-6)
    return np.clip((d / dc) ** b, 0.0, 1.0)


def scf_released(d, rho):
    """As released: min(1, 100*d*RQSNCR) with RQSNCR = 1/10."""
    return np.clip(10.0 * d, 0.0, 1.0)


def scf_tanh(d, rho, p=(0.016, 100.0, 1.6)):
    """Niu & Yang as implemented -- cannot reach 1."""
    z0, rho_new, m = p
    return np.tanh(d / np.maximum(2.5 * z0 * (rho / rho_new) ** m, 1e-6))


# ------------------------------------------------------------------- the data --
with xr.open_dataset(f'{OBSD}/RIHMI-WDC_snmar.nc') as ds:
    la, lo = ds['lat'].values, ds['lon'].values
    k = (la >= BOX['lat0']) & (la <= BOX['lat1']) & (lo >= BOX['lon0']) & (lo <= BOX['lon1'])
    cov = ds['fraction_of_the_snow_course_covered_by_snow'].isel(station=k).values / 10.0
    dep = ds['snow_depth_mean'].isel(station=k).values / 100.0        # cm -> m
    rho = ds['snow_density'].isel(station=k).values * 1000.0          # g/cm3 -> kg/m3
    tm = ds['time'].values
mo = tm.astype('datetime64[M]').astype(int) % 12

print(__doc__.split('TWO CAVEATS')[0])
print('=' * 92)

DBINS = np.array([0, .02, .05, .08, .12, .17, .23, .30, .40, .55, .80, 1.5])
RBINS = np.array([0, 120, 150, 180, 210, 250, 300, 500])

FIT = {}
for t in (0, 1):
    c, d, r = cov[:, :, t], dep[:, :, t], rho[:, :, t]
    m3 = np.repeat(mo[None, :], c.shape[0], axis=0)
    ok = np.isfinite(c) & np.isfinite(d) & np.isfinite(r) & (d > 0) & (r > 30) & (r < 600)
    c, d, r, mm = c[ok], d[ok], r[ok], m3[ok]
    print(f'\n\n### COURSE TYPE {t+1} -- {TYPES[t]}, {c.size} surveys\n')

    # ---- empirical surface -----------------------------------------------------
    print('  EMPIRICAL mean cover by depth (rows, m) and density (cols, kg/m3)\n')
    print(f'  {"depth":>12s}' + ''.join(f'{f"{RBINS[j]}-{RBINS[j+1]}":>10s}'
                                        for j in range(len(RBINS) - 1)))
    for i in range(len(DBINS) - 1):
        sd = (d >= DBINS[i]) & (d < DBINS[i + 1])
        if sd.sum() < 20:
            continue
        row = []
        for j in range(len(RBINS) - 1):
            s = sd & (r >= RBINS[j]) & (r < RBINS[j + 1])
            row.append(f'{c[s].mean():10.3f}' if s.sum() >= 20 else f'{"":>10s}')
        print(f'  {DBINS[i]:.2f}-{DBINS[i+1]:.2f}'.rjust(14) + ''.join(row))

    # ---- fit to BIN MEANS so the saturated bulk does not swamp the melt tail ----
    bd, br, bc, bn = [], [], [], []
    for i in range(len(DBINS) - 1):
        for j in range(len(RBINS) - 1):
            s = (d >= DBINS[i]) & (d < DBINS[i + 1]) & (r >= RBINS[j]) & (r < RBINS[j + 1])
            if s.sum() >= 20:
                bd.append(d[s].mean()); br.append(r[s].mean())
                bc.append(c[s].mean()); bn.append(s.sum())
    bd, br, bc, bn = map(np.asarray, (bd, br, bc, bn))
    wgt = np.sqrt(bn)

    def grid(brange):
        best, bp = np.inf, None
        for dc0 in np.arange(0.01, 0.60, 0.0025):
            for m in np.arange(0.0, 3.02, 0.02):
                for b in brange:
                    e = np.sum(wgt * (scf_new(bd, br, (dc0, m, b)) - bc) ** 2)
                    if e < best:
                        best, bp = e, (dc0, m, b)
        return bp
    bp = grid(np.arange(0.05, 3.02, 0.05))
    bp1 = grid(np.array([1.0]))          # b fixed at 1: pure ramp, density-dependent d_c
    FIT[t] = bp
    rms = lambda f: np.sqrt(np.sum(wgt * (f - bc) ** 2) / wgt.sum())
    print(f'\n  FITTED  d_c = {bp[0]:.3f} * (rho/100)^{bp[1]:.2f}   b = {bp[2]:.2f}'
          f'      [{len(bc)} bins, {bn.sum()} surveys]')
    print(f'  b=1 variant  d_c = {bp1[0]:.3f} * (rho/100)^{bp1[1]:.2f}'
          f'   RMSE {rms(scf_new(bd, br, bp1)):.4f}')
    print(f'  weighted RMSE   new form      {rms(scf_new(bd, br, bp)):.4f}')
    print(f'                  as-released   {rms(scf_released(bd, br)):.4f}')
    print(f'                  tanh (N2)     {rms(scf_tanh(bd, br)):.4f}')
    print(f'  d_c at rho=130 (Oct): {bp[0]*(1.30)**bp[1]:.3f} m'
          f'   at rho=285 (May): {bp[0]*(2.85)**bp[1]:.3f} m')

    # ---- seasonal cycle: average the FUNCTION over surveys, never the inputs ----
    # scf(mean d, mean rho) != mean scf(d, rho); the earlier version of this check
    # made the correct fit look like it saturated in May when it does not.
    print(f'\n  MONTHLY CHECK -- mean over individual surveys, not of the means\n')
    print(f'  {"month":8s}{"n":>7s}{"depth":>8s}{"rho":>7s}{"OBS":>8s}{"new":>8s}'
          f'{"new b=1":>9s}{"released":>10s}{"tanh":>8s}')
    # Snow-year order Sep->Aug.  Summer months are printed even when nearly empty:
    # a survey count of 3 is itself the finding (courses are not run without snow),
    # and hiding them would make the melt tail look better sampled than it is.
    for mi in (8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7):
        s = mm == mi
        if s.sum() == 0:
            print(f'  {MON[mi]:8s}{0:7d}{"":8s}{"":7s}{"-- no surveys --":>8s}')
            continue
        flag = ' *' if s.sum() < 30 else ''
        print(f'  {MON[mi]:8s}{s.sum():7d}{d[s].mean():8.3f}{r[s].mean():7.0f}{c[s].mean():8.3f}'
              f'{scf_new(d[s], r[s], bp).mean():8.3f}{scf_new(d[s], r[s], bp1).mean():9.3f}'
              f'{scf_released(d[s], r[s]).mean():10.3f}{scf_tanh(d[s], r[s]).mean():8.3f}{flag}')
    print('    * fewer than 30 surveys -- not a climatology, shown for completeness')

# ------------------------------------------------------------------- summary --
print('\n\n' + '=' * 92)
print('SUMMARY -- what to implement\n')
for t in (0, 1):
    if t in FIT:
        dc0, m, b = FIT[t]
        print(f'  type {t+1} {TYPES[t]:16s}  d_c = {dc0:.3f}*(rho/100)^{m:.2f}   b = {b:.2f}')
print("""
  The two course types give the roughness dependence empirically: a forest course
  should need MORE snow to reach complete cover than an open field, because the
  understorey and litter roughness is larger.  If the fitted d_c0 for forest exceeds
  that for field, that ratio is the physical basis for scaling d_c with the model's
  own vegetation roughness -- by TYPE, not by region, and not invented.

  Implementation note: this is a drop-in replacement for the ZCVS branch at
  surfbc_ctl_mod.F90:340-355.  It needs no new prognostic and no hysteresis flag --
  the density the model already carries does the autumn/spring separation, which is
  what the tanh was trying and failing to do.""")
