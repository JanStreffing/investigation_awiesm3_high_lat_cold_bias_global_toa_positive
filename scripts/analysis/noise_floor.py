"""How big is a tuning signal allowed to be before we believe it?

Round 10 produced three results that do not behave like physics:
  * ABB8 (A1b+B2+B8) gave Siberia JJA T2m -0.167 K when its parts sum to ~+0.76
  * E1 (RVLAMSK 5->2.5) reversed B8's (5) gain instead of extending it
  * C1 (RLAM 75) and C2 (RLAM 40) are non-monotonic in RLAM

All three are explained at once if the 4-year Siberian JJA mean simply is not
precise enough to resolve differences of a few tenths of a kelvin. This script
measures that precision instead of assuming it.

Method. Every run shares the same prescribed SST and the same 1870 initial state,
so a given year is partly common to all runs (forced by that year's SST) and
partly chaotic. Decompose the per-year, per-run field as

    X[run, year] = mu + a[run] + g[year] + eps[run, year]

a[run] is the tuning signal we are after, g[year] the shared forced excursion,
and eps the internal atmospheric noise that a 4-year mean fails to average away.
A two-way ANOVA on the run x year matrix separates them. The number that matters
is then the standard error of a *difference* between two runs' 4-year means,

    SE_diff = sd(eps) * sqrt(2/n_years)

because every entry in the results table is such a difference (run minus control).

Reported on the same box, land mask, months and years as eval_round10_A.py, so
the numbers are directly comparable to that table.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
Y0, Y1 = 1872, 1875
YEARS = list(range(Y0, Y1 + 1))
JJA = [5, 6, 7]                     # same 0-based convention as eval_round10_A.py
BOX = ((55, 75), (60, 180))         # Siberia land
ACC = 3600.0

RUNS = [('control', 'amip_pi_base'), ('A1a ovl=0.10', 'amip_A1_overlap01'),
        ('A1b ovl=0.35', 'amip_A1_overlap035'), ('A2 KKland=150', 'amip_A2_kknumland150'),
        ('expA rvrs=500', 'amip_expA_rvrsmin500'),
        ('A1c depth1500', 'amip_A1c_depliqdepth1500'),
        ('B1 detrpen.45', 'amip_B1_detrpen045'),
        ('B2 convi=25', 'amip_B2_clddiffconvi25'),
        ('AB ovl+convi', 'amip_AB_ovl035_convi25'),
        ('B3 clddiff', 'amip_B3_clddiff15e6'),
        ('B4 entshalp3', 'amip_B4_entshalp3'),
        ('B5 capdcycl0', 'amip_B5_capdcycl0'),
        ('B6 lcritsnow', 'amip_B6_lcritsnow1e5'),
        ('B7 rvice.22', 'amip_B7_rvice022'),
        ('B8 lamsk5', 'amip_B8_lamsk5'),
        ('ABB8 A1b+B2+B8', 'amip_ABB8'),
        ('C1 rlam75', 'amip_C1_rlam75'),
        ('C2 rlam40', 'amip_C2_rlam40'),
        ('E1 lamsk2.5', 'amip_E1_lamsk25')]

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def box_mean(a2d, lat, lon):
    """Area-weighted mean over the Siberian land box, for one month-field."""
    yi = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    l180 = ((lon + 180) % 360) - 180
    xi = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    sub = a2d[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    L = lsm[np.ix_(np.where(yi)[0], np.where(xi)[0])]
    m = L > 0.5
    w = np.broadcast_to(np.cos(np.deg2rad(lat[yi]))[:, None], sub.shape)
    return np.average(sub[m], weights=w[m])


def per_year(run, var='2t'):
    """Siberia JJA mean of `var` for each individual year; None if incomplete."""
    vals = []
    for y in YEARS:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            return None
        ds = xr.open_dataset(f)
        a = ds[var].values
        if var != '2t' and var != 'tcc':
            a = a / ACC
        lat, lon = ds[var].lat.values, ds[var].lon.values
        ds.close()
        vals.append(np.mean([box_mean(a[m], lat, lon) for m in JJA]))
    return np.array(vals)


# ---- assemble the run x year matrix ----------------------------------------
labs, X = [], []
for lab, run in RUNS:
    v = per_year(run)
    if v is None:
        print(f'  !! {lab}: incomplete, skipped'); continue
    labs.append(lab); X.append(v)
X = np.array(X) - 273.15
n_run, n_yr = X.shape

# ---- two-way decomposition --------------------------------------------------
mu = X.mean()
a_run = X.mean(axis=1) - mu          # tuning signal
g_yr = X.mean(axis=0) - mu           # shared forced year excursion
eps = X - (mu + a_run[:, None] + g_yr[None, :])
dof = (n_run - 1) * (n_yr - 1)
sd_eps = np.sqrt((eps ** 2).sum() / dof)
se_mean = sd_eps / np.sqrt(n_yr)
se_diff = sd_eps * np.sqrt(2.0 / n_yr)

print(f'Siberia JJA T2m, {n_run} runs x {n_yr} years ({Y0}-{Y1})\n')
print(f'  shared forced year signal g[year]  : {np.round(g_yr, 3)}  (range '
      f'{g_yr.max()-g_yr.min():.2f} K)')
print(f'  internal noise sd(eps)             : {sd_eps:.3f} K   (dof={dof})')
print(f'  SE of one run 4-yr mean            : {se_mean:.3f} K')
print(f'  SE of a run-minus-control diff     : {se_diff:.3f} K')
print(f'  => 95% detection threshold         : +-{1.96*se_diff:.3f} K\n')

ctl = labs.index('control')
print(f"  {'run':16s} {'per-year JJA T2m':34s} {'mean':>7s} {'vs ctl':>8s} {'t':>6s}  verdict")
for i, lab in enumerate(labs):
    d = X[i].mean() - X[ctl].mean()
    t = d / se_diff
    if i == ctl:
        verd = '(reference)'
    else:
        verd = 'SIGNIFICANT' if abs(t) > 1.96 else ('marginal' if abs(t) > 1.0 else 'noise')
    print(f'  {lab:16s} {str(np.round(X[i], 2)):34s} {X[i].mean():7.2f} '
          f'{d:+8.3f} {t:+6.2f}  {verd}')

print(f'\n  Of {n_run-1} tuning runs, '
      f'{sum(abs(X[i].mean()-X[ctl].mean())/se_diff > 1.96 for i in range(n_run) if i != ctl)} '
      f'clear the 95% threshold.')
