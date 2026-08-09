"""The energy target, re-based on a period-clean footing.

WHY THIS EXISTS.  Every band-level energy target in this campaign was measured on
`amip_pi_base` -- an 1870s atmosphere -- against CERES EBAF 07/2005-06/2015.  A PI
model scored against a present-day satellite record carries an epoch offset in every
number, and the size of that offset has never been quantified per band.  It matters
because the whole next round would be designed around the Southern Ocean's +5.81
W/m2, and because the global figure moves a lot: the PI arm reads net TOA +0.64
against a target of ~0, while the SAME MODEL run over 1990-2014 reads +2.20 against
CERES's +0.97, i.e. +1.23 too positive.  Those two framings differ by a factor of two
and only one of them is a like-for-like comparison.

`amip_presentday` (1989-2015, transient historical forcing, which is CORRECT for that
period -- unlike the PI arm, whose NCMIPFIXYR never took effect) overlaps CERES
directly.  So this script recomputes the band decomposition on that arm and reports
BOTH arms side by side, making the epoch offset an explicit measured column instead of
an unquantified caveat.

NO REMAPPING.  Band means over full latitude circles are exact on each source's own
grid, so the model stays on its grid and CERES on its 1-degree grid.  Remapping would
add error and buy nothing.

THE OLR TEST, which is the second reason to run this.  The campaign has a measured
global tropospheric cold bias (-0.7 to -2.2 K, amip_presentday vs ERA5) and a measured
OLR deficit.  A cold troposphere under-emits, so those two should be quantitatively
consistent -- and on the PI-vs-CERES numbers they are not: a Planck scaling of the
observed cold bias predicts several W/m2 of OLR deficit where only 0.74 is seen.  This
script measures the OLR deficit period-cleanly and compares it against the Planck
expectation computed from the model's own emitting temperature, so the discrepancy is
either resolved by the epoch correction or promoted to a real compensating error.

CONVENTIONS.
  model:  tsr = net SW down (+),  ttr = net LW (-, so OLR = -ttr); all ACCUMULATED
          J/m2 over the 3600 s output step -- divide by ACC or every number is 3600x
          too large.  tcc is a fraction, NOT accumulated.
  CERES:  toa_sw_all_clim / toa_lw_all_clim are UPWARD fluxes; toa_net_all_clim is
          positive DOWN.  ASR is recovered as net + lw_up (no incoming-solar variable
          in the subset).  toa_cre_sw_clim is already all-minus-clear, positive down.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

ACC = 3600.0
RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'

ARMS = [('PI  amip_pi_base 1872-1915', 'amip_pi_base', 1872, 1915),
        ('PD  amip_presentday 1990-2014', 'amip_presentday', 1990, 2014)]

# Tiling bands: contributions must sum to the global mean.  The SO overlay overlaps
# them deliberately and is excluded from the sum.
BANDS = [('90S-65S', -90, -65), ('65S-45S  <- SO', -65, -45), ('45S-30S', -45, -30),
         ('tropics 30S-30N', -30, 30), ('30N-45N', 30, 45), ('45N-65N', 45, 65),
         ('65N-90N', 65, 90)]

print(__doc__)
print('=' * 100)


# ---------------------------------------------------------------- model side
def model_fields(run, y0, y1):
    """Annual-mean 2-D fields [W/m2] and cloud fraction, averaged over the years."""
    D = f'{RT}/{run}/outdata/oifs'
    out, lat, nyr = {}, None, 0
    for v in ('tsr', 'ttr', 'tsrc', 'ttrc', 'tcc'):
        acc = []
        for y in range(y0, y1 + 1):
            f = f'{D}/atm_remapped_1m_{v}_1m_{y}-{y}.nc'
            if not os.path.exists(f):
                continue
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[v].values
                if lat is None:
                    lat = d['lat'].values
            # tcc is a fraction; the radiative fluxes are accumulated over the step
            acc.append(a.mean(axis=0) / (1.0 if v == 'tcc' else ACC))
        if not acc:
            return None, None, 0
        out[v] = np.mean(acc, axis=0)
        nyr = max(nyr, len(acc))
    return out, lat, nyr


def band_mean(field, lat, a, b):
    sel = (lat >= a) & (lat < b)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(field[sel, :].mean(axis=1), weights=w))


def area_frac(lat, a, b):
    """Cos-lat area share of a band, on this grid."""
    w = np.cos(np.deg2rad(lat))
    sel = (lat >= a) & (lat < b)
    return float(w[sel].sum() / w.sum())


# ---------------------------------------------------------------- CERES side
cds = xr.open_dataset(CERESF)
clat = cds['lat'].values
cw = np.cos(np.deg2rad(clat))


def ceres(varname):
    """Annual-mean 2-D CERES climatology field."""
    return cds[varname].mean('ctime').values


C = dict(net=ceres('toa_net_all_clim'), lwup=ceres('toa_lw_all_clim'),
         swcre=ceres('toa_cre_sw_clim'), lwcre=ceres('toa_cre_lw_clim'),
         cld=ceres('cldarea_total_daynight_clim'))
C['asr'] = C['net'] + C['lwup']          # no incoming-solar variable in the subset


def cband(key, a, b):
    sel = (clat >= a) & (clat < b)
    return float(np.average(C[key][sel, :].mean(axis=1), weights=cw[sel]))


# sanity: the CERES global numbers must reproduce the published EBAF values
gnet = float(np.average(C['net'].mean(axis=1), weights=cw))
gasr = float(np.average(C['asr'].mean(axis=1), weights=cw))
golr = float(np.average(C['lwup'].mean(axis=1), weights=cw))
gsw = float(np.average(C['swcre'].mean(axis=1), weights=cw))
glw = float(np.average(C['lwcre'].mean(axis=1), weights=cw))
gcld = float(np.average(C['cld'].mean(axis=1), weights=cw))
print(f'CERES EBAF 07/2005-06/2015 global check: net {gnet:+.2f} (published ~+0.97), '
      f'ASR {gasr:.2f} (~240.5), OLR {golr:.2f} (~240.2),')
print(f'   SW CRE {gsw:.2f} (~-45.4), LW CRE {glw:.2f} (~+25.8), cloud area {gcld:.2f} %')
ok = abs(gasr - 240.5) < 1.5 and abs(golr - 240.2) < 1.5
print('   -> ' + ('reproduces the published values; the gridded file is being read correctly'
                  if ok else '*** DOES NOT MATCH PUBLISHED EBAF -- check variable names/signs ***'))
print()

# ---------------------------------------------------------------- per arm
res = {}
for tag, run, y0, y1 in ARMS:
    m, lat, nyr = model_fields(run, y0, y1)
    if m is None:
        print(f'{tag}: no output'); continue
    net = m['tsr'] + m['ttr']
    asr = m['tsr']
    olr = -m['ttr']
    swcre = m['tsr'] - m['tsrc']
    lwcre = m['ttr'] - m['ttrc']
    cld = m['tcc'] * 100.0
    res[run] = dict(lat=lat, net=net, asr=asr, olr=olr, swcre=swcre, lwcre=lwcre,
                    cld=cld, nyr=nyr)

    g = lambda f: band_mean(f, lat, -90, 90)
    print('=' * 100)
    print(f'{tag}   ({nyr} years)')
    print('=' * 100)
    print(f'  {"":22s} {"model":>9s} {"CERES":>9s} {"diff":>9s}')
    for lbl, mv, cv in (('net TOA', g(net), gnet), ('absorbed SW', g(asr), gasr),
                        ('OLR', g(olr), golr), ('SW CRE', g(swcre), gsw),
                        ('LW CRE', g(lwcre), glw), ('cloud area [%]', g(cld), gcld)):
        print(f'  {lbl:22s} {mv:9.3f} {cv:9.3f} {mv - cv:+9.3f}')

    print(f'\n  {"band":18s} {"net":>8s} {"CERES":>8s} {"diff":>8s} {"area":>7s} '
          f'{"contrib":>8s} | {"dSWCRE":>8s} {"dLWCRE":>8s} {"dCLD":>7s}')
    tot = 0.0
    for name, a, b in BANDS:
        mv, cv = band_mean(net, lat, a, b), cband('net', a, b)
        fr = area_frac(lat, a, b)
        contrib = (mv - cv) * fr
        tot += contrib
        print(f'  {name:18s} {mv:8.2f} {cv:8.2f} {mv - cv:+8.2f} {fr * 100:6.1f}% '
              f'{contrib:+8.3f} | {band_mean(swcre, lat, a, b) - cband("swcre", a, b):+8.2f} '
              f'{band_mean(lwcre, lat, a, b) - cband("lwcre", a, b):+8.2f} '
              f'{band_mean(cld, lat, a, b) - cband("cld", a, b):+7.2f}')
    print(f'  {"sum of bands":18s} {"":8s} {"":8s} {"":8s} {"":7s} {tot:+8.3f}'
          f'   <- must equal global net diff {g(net) - gnet:+.3f}')

# ---------------------------------------------------------------- epoch offset
if len(res) == 2:
    pi, pd = res['amip_pi_base'], res['amip_presentday']
    print()
    print('=' * 100)
    print('THE EPOCH OFFSET, per band: how much of each "model error" is the reference period')
    print('=' * 100)
    print(f'  {"band":18s} {"PI-vs-CERES":>12s} {"PD-vs-CERES":>12s} {"offset":>9s}  '
          f'{"verdict":<34s}')
    for name, a, b in BANDS:
        dpi = band_mean(pi['net'], pi['lat'], a, b) - cband('net', a, b)
        dpd = band_mean(pd['net'], pd['lat'], a, b) - cband('net', a, b)
        off = dpi - dpd
        v = ('epoch-dominated -- PI target unsafe' if abs(off) > abs(dpd)
             else 'robust to the epoch' if abs(off) < 0.3 * max(abs(dpd), 1e-9)
             else 'partly epoch')
        print(f'  {name:18s} {dpi:+12.2f} {dpd:+12.2f} {off:+9.2f}  {v:<34s}')
    dgpi = band_mean(pi['net'], pi['lat'], -90, 90) - gnet
    dgpd = band_mean(pd['net'], pd['lat'], -90, 90) - gnet
    print(f'  {"GLOBAL":18s} {dgpi:+12.2f} {dgpd:+12.2f} {dgpi - dgpd:+9.2f}')
    print(f'\n  Read the PD column as the model error.  The PI column is what the campaign has')
    print(f'  been tuning against, and the offset is how much of it was never model error.')

    # ------------------------------------------------------------ the OLR test
    print()
    print('=' * 100)
    print('THE OLR TEST: is the OLR deficit consistent with the tropospheric cold bias?')
    print('=' * 100)
    dolr = band_mean(pd['olr'], pd['lat'], -90, 90) - golr
    dasr = band_mean(pd['asr'], pd['lat'], -90, 90) - gasr
    Te = (golr / 5.670374419e-8) ** 0.25
    # Global tropospheric cold bias vs ERA5, measured in vertical_bias_column.py /
    # report.tex sub:vprof on THIS SAME RUN, so it is period-clean and needs no offset.
    PROF = {1000: -0.73, 925: -0.65, 850: -0.99, 700: -1.15, 500: -1.49, 300: -2.22}
    # Weight by where the emission comes from: crude but honest -- the bulk of clear-sky
    # OLR originates between 700 and 300 hPa, so average those levels.
    dT = np.mean([PROF[p] for p in (700, 500, 300)])
    planck = 4 * 5.670374419e-8 * Te ** 3
    print(f'  measured, period-clean:   dOLR {dolr:+.3f}   dASR {dasr:+.3f}   '
          f'-> dNET {dasr - dolr:+.3f} W/m2')
    print(f'  effective emitting T from CERES OLR: {Te:.2f} K, Planck slope '
          f'{planck:.3f} W/m2/K')
    print(f'  mean model cold bias 700-300 hPa vs ERA5: {dT:+.2f} K')
    print(f'  => a purely radiative response to that cold bias predicts '
          f'dOLR = {planck * dT:+.2f} W/m2')
    gap = dolr - planck * dT
    print(f'  => measured minus predicted = {gap:+.2f} W/m2')
    print()
    if abs(gap) > 1.0:
        print(f'  *** {abs(gap):.1f} W/m2 UNACCOUNTED.  The atmosphere is {abs(dT):.1f} K cold')
        print(f'      through the emitting layer yet emits nearly the right amount, so')
        print(f'      something is offsetting it: too little cloud LW trapping, too little')
        print(f'      water vapour, or cloud tops too low.  LW CRE is the first place to look')
        print(f'      -- model {band_mean(pd["lwcre"], pd["lat"], -90, 90):.2f} vs CERES '
              f'{glw:.2f} ({band_mean(pd["lwcre"], pd["lat"], -90, 90) - glw:+.2f}).')
        print(f'      A compensating pair like that is exactly what makes a tuned global')
        print(f'      budget fragile: both errors are real and only their sum is small.')
    else:
        print('  The two are consistent within 1 W/m2: the OLR deficit IS the cold troposphere,')
        print('  and warming the troposphere to observations would remove most of the imbalance.')

cds.close()
