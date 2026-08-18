"""Does the Southern Ocean cloud tuning transfer from AMIP to coupled, or not?

THE QUESTION.  S4 closes only 27 % of the coupled SO gap, and DMS -- scored coupled from
11F's 20 years on 2026-08-18 -- delivers only 0.19-0.44x of its AMIP CRE forcing.  Two
readings, and they imply completely different next steps:

  (a) THE CLOUD TUNING DOES NOT TRANSFER.  The coupled cloud field is worse than AMIP's,
      so the levers are being applied to a different atmosphere and their AMIP calibration
      does not carry.  Then the SO problem must be re-tuned coupled, and AMIP screening is
      of limited value for it.

  (b) THE CLOUD FIELD TRANSFERS AND THE SURFACE IS WORSE.  The coupled SO SW CRE matches
      AMIP's, but the free ocean sits at the wrong temperature and ice cover, so the
      residual error is oceanic rather than radiative.  Then the levers are fine and the
      SO bias is not a cloud problem at all.

WHY AMIP-VS-COUPLED IS THE DISCRIMINATOR.  AMIP prescribes OBSERVED SST, so its cloud
field is the model's response to a correct ocean.  Coupled uses the model's own SST.  Put
the same diagnostic on both and the difference is attributable.

WHAT IS COMPARED, identically on every run: SO 45-65S TOA SW CRE (tsr - tsrc), clear-sky
SW (tsrc), and total cloud cover (tcc).  CERES toa_cre_sw_clim anchors the CRE.

AMIP arms included so the SO lever ladder is on the SAME diagnostic as the coupled runs,
which the report's tables are not -- its AMIP band numbers sit on a different CERES value
(-55.81) than a direct read of the gridded climatology gives, and that discrepancy is
recorded but unresolved.  P5 is the adopted-stack control; Y2 is RCL_OVERLAPLIQICE 0.10,
which reached SO SW -5.029 in the round-27 ladder, the strongest SO lever found and one
that has never been run coupled.

ACCUMULATION.  IFS TOA fluxes are accumulated J/m^2.  Verified empirically 2026-08-18:
the divisor is 3600 for BOTH the AMIP daily and the coupled monthly output -- global ASR
comes out 240.4 W/m^2 on AMIP 1900 with /3600 and 10.0 with /86400.  Do not assume the
divisor follows the output frequency.

PERIOD CAVEAT, stated because it is not correctable here.  AMIP runs 1870-1917 on
transient forcing (the campaign's AMIP arms do not read NCMIPFIXYR; Krakatoa is visible
1883-86).  The coupled runs are correctly pinned at 1850.  So an AMIP-minus-coupled
difference carries a forcing-era component as well as a coupling one.  That does not
affect the (a)-vs-(b) question, which turns on whether the gap is in CRE or in the
surface, but it does mean the absolute numbers are not interchangeable.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

AMIP = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CPL = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ACC = 3600.0
SO = (-65.0, -45.0)
FLUX = ('tsr', 'tsrc', 'ttr', 'ttrc')

RUNS = [
    ('AMIP P5  (control)', f'{AMIP}/amip_P5_swemin15', 'amip'),
    ('AMIP Y2  (ovl 0.10)', f'{AMIP}/amip_Y2_ovl01', 'amip'),
    ('AMIP LY2 (ovl.10 long)', f'{AMIP}/amip_LY2_long', 'amip'),
    ('AMIP Y1  (ovl 0.35)', f'{AMIP}/amip_Y1_ovl035', 'amip'),
    ('AMIP LX3 (DMS+INPPMIN)', f'{AMIP}/amip_LX3_long', 'amip'),
    ('CPL  11E (base)', f'{CPL}/Tuning_test_11E_swemin15_K1', 'cpl'),
    ('CPL  11G (+S4)', f'{CPL}/Tuning_test_11G_inppmin50k', 'cpl'),
    # The overlap lever coupled, submitted 2026-08-18.  Each is 11G plus one namelist
    # number, so 11G is the control and no new control run is needed.  They print
    # "no tsr/tsrc output" until their first leg lands, which is not an error.
    ('CPL  11L (+ovl 0.35)', f'{CPL}/11L', 'cpl'),
    ('CPL  11M (+ovl 0.10)', f'{CPL}/11M', 'cpl'),
]


def files_for(root, var, kind):
    d = f'{root}/outdata/oifs'
    pats = ([f'{d}/atm_remapped_1d_{var}_1d_*.nc'] if kind == 'amip'
            else [f'{d}/atm_remapped_1m_{var}_*.nc'])
    out = []
    for p in pats:
        out += sorted(glob.glob(p))
    return out


def band_series(root, var, kind, nmax=30):
    """Per-year SO band mean, using the LAST nmax years available."""
    fs = files_for(root, var, kind)
    if not fs:
        return None
    fs = fs[-nmax:]
    vals = []
    div = ACC if var in FLUX else 1.0
    for f in fs:
        try:
            with xr.open_dataset(f, decode_times=False) as d:
                if var not in d:
                    continue
                a = d[var].values / div
                lat = d['lat'].values
        except OSError:
            continue
        sel = (lat >= SO[0]) & (lat < SO[1])
        w = np.cos(np.deg2rad(lat[sel]))
        vals.append(float(np.average(a.mean(axis=0)[sel, :].mean(axis=1), weights=w)))
    return np.array(vals) if vals else None


print(__doc__)
print('=' * 98)

try:
    with xr.open_dataset(CERESF) as c:
        clat = c['lat'].values
        s = (clat >= SO[0]) & (clat < SO[1])
        cw = np.cos(np.deg2rad(clat[s]))
        CERES_CRE = float(np.average(
            c['toa_cre_sw_clim'].values.mean(axis=0)[s, :].mean(axis=1), weights=cw))
        CERES_CLR = float(np.average(
            c['toa_sw_clr_t_clim'].values.mean(axis=0)[s, :].mean(axis=1), weights=cw))
except Exception as exc:
    CERES_CRE = CERES_CLR = np.nan
    print(f'  (CERES unavailable: {exc})')

print(f'\n  CERES SO 45-65S:  SW CRE {CERES_CRE:.2f}   clear-sky SW(up) {CERES_CLR:.2f}\n')
print(f'  {"run":22s} {"n":>3s} {"SW CRE":>9s} {"vs CERES":>9s} {"clrsky SW":>10s} '
      f'{"tcc":>7s}')

res = {}
for tag, root, kind in RUNS:
    tsr = band_series(root, 'tsr', kind)
    tsrc = band_series(root, 'tsrc', kind)
    tcc = band_series(root, 'tcc', kind)
    if tsr is None or tsrc is None:
        print(f'  {tag:22s}  no tsr/tsrc output')
        continue
    n = min(len(tsr), len(tsrc))
    cre = (tsr[:n] - tsrc[:n]).mean()
    clr = tsrc[:n].mean()
    cc = tcc.mean() if tcc is not None else np.nan
    res[tag] = dict(cre=cre, clr=clr, tcc=cc, n=n,
                    cre_sd=(tsr[:n] - tsrc[:n]).std(ddof=1))
    print(f'  {tag:22s} {n:3d} {cre:9.3f} {cre - CERES_CRE:+9.3f} {clr:10.3f} '
          f'{cc:7.4f}')

print('\n' + '=' * 98)
print('\nSURFACE STATE -- the other half of the question (coupled only; AMIP SST is'
      ' prescribed observed)\n')
for tag, root, kind in RUNS:
    if kind != 'cpl':
        continue
    sst = band_series(root, 'sst', kind)
    ci = band_series(root, 'ci', kind)
    if sst is not None:
        print(f'  {tag:22s} SO SST {sst.mean() - 273.15:7.3f} C', end='')
        if ci is not None:
            print(f'   sea ice {ci.mean() * 100:6.3f} %')
        else:
            print()

print('\n' + '=' * 98)
print('\nREADING\n')
p5 = res.get('AMIP P5  (control)')
e11 = res.get('CPL  11E (base)')
if p5 and e11:
    d_cre = e11['cre'] - p5['cre']
    print(f'  AMIP P5 SO SW CRE   {p5["cre"]:+8.3f}   ({p5["cre"] - CERES_CRE:+.3f} vs CERES)')
    print(f'  CPL  11E SO SW CRE  {e11["cre"]:+8.3f}   ({e11["cre"] - CERES_CRE:+.3f} vs CERES)')
    print(f'  coupled minus AMIP  {d_cre:+8.3f}\n')
    frac = (p5['cre'] - CERES_CRE) / (e11['cre'] - CERES_CRE)
    print(f'  THE PART THAT IS NOT A COUPLING PROBLEM: the AMIP control already carries'
          f'\n  {p5["cre"] - CERES_CRE:+.3f} of the coupled arm\'s {e11["cre"] - CERES_CRE:+.3f},'
          f' i.e. {100 * frac:.0f} % of the coupled SO CRE\n  error is inherited from AMIP'
          ' with a PERFECT ocean. Coupling adds the rest.\n')
    if abs(d_cre) < 2.0:
        print('  -> (b) THE CLOUD FIELD LARGELY TRANSFERS. The coupled SO CRE is close to')
        print('     AMIP\'s, so the levers are acting on a comparable atmosphere and their')
        print('     AMIP calibration is not the problem. Look at the surface terms above:')
        print('     if SST and ice are off, the residual SO error is oceanic, and no amount')
        print('     of cloud tuning will close it.')
    else:
        sign = 'WORSE' if d_cre > 0 else 'BETTER'
        print(f'  -> (a) THE CLOUD FIELD DOES NOT TRANSFER: coupled is {sign} than AMIP by')
        print(f'     {abs(d_cre):.2f} W/m2 before any lever is applied. AMIP-calibrated SO')
        print('     levers are being applied to a different atmosphere, which is enough on')
        print('     its own to explain why they under-deliver coupled.')
y2 = res.get('AMIP LY2 (ovl.10 long)') or res.get('AMIP Y2  (ovl 0.10)')
if y2 and p5:
    print(f'\n  Lever check on the SAME diagnostic: Y2 (ovl 0.10) minus P5 = '
          f'{y2["cre"] - p5["cre"]:+.3f} W/m2 SO SW CRE')
    print('  Round-27 ladder gave -5.029 cumulative. This is the strongest SO lever found')
    print('  and it has never been run coupled.')
