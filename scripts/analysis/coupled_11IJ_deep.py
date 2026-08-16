"""11I / 11J in the campaign lineage, not in isolation.

WHY THE FIRST PASS WAS TOO SHALLOW.  coupled_11IJ_eval.py scored two pair differences
against unpaired thresholds and reported "nothing resolves".  That is true but it is
not the analysis, for three reasons:

  1. NO LINEAGE.  These arms are a chain, not two isolated pairs:
         11E  base                      (swemin15_K1)
         11G  = 11E + S4                (RCL_INPPMIN 70000 -> 50000)
         11I  = 11G + land-surface repair
         11J  = 11I + Raupach roughness
     What matters for adoption is where the STACK sits, not only each increment.  A
     lever that is individually ns can still matter cumulatively, and three ns steps
     that all lean the same way are different from three that cancel.

  2. UNPAIRED THRESHOLDS ON BRANCHED RUNS.  Every arm branches from the SAME 1350
     state, so the early years are near-identical by construction and the difference
     series grows as the runs decorrelate.  Averaging a growing difference over all 40
     years dilutes it toward zero, and scoring it against 1.96*sd*sqrt(2/n) from the
     control's absolute interannual scatter throws away the pairing entirely.  The
     difference series has its OWN, much smaller, scatter.  This script reports both,
     plus the decade breakdown and a trend, so a diluted signal cannot hide.

  3. NO TARGET.  Differences between arms say nothing about whether the model is right.
     The campaign exists because of a boreal summer cold bias severe enough to lose the
     forest, so the arms are also scored against ERA5 Siberian JJA T2m and CERES SO SW
     CRE, and against the thing that actually has to be true: is the forest there?
     FORESTFPC / TREEFPC from fpc.out, which LPJ-GUESS writes per rank continuously and
     which therefore SURVIVED the XIOS finalisation crash that destroyed 11J's leg 4.

CAVEAT CARRIED FORWARD.  11I - 11G is not one variable: the _v2 soil map, the
L096 -> TCO95-land gridlist and the cold-start binary path all change together, and 2
and 3 are inseparable (the new gridlist is what adds the 307 cells the cold-start path
lets run).  11J - 11I IS one variable: same binary md5 ff9b5507, gated on ifraupachz0.
"""
import os, glob
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
ERA5_SIB = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'era5_siberia_jja.txt')
ACC = 3600.0

CHAIN = [('11E', f'{R092}/Tuning_test_11E_swemin15_K1', 'base'),
         ('11G', f'{R092}/Tuning_test_11G_inppmin50k', '+ S4 RCL_INPPMIN 50k'),
         ('11I', f'{R092}/Tuning_test_11I_v2soil', '+ land-surface repair'),
         ('11J', f'{R092}/Tuning_test_11J_v2soil_raupach', '+ Raupach roughness')]

SIB = (55.0, 75.0, 60.0, 180.0)
SO = (-65.0, -45.0)
DJF, JJA = [12, 1, 2], [6, 7, 8]
FLUX = ('tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str')
KEY = ('net TOA [W/m2]', 'Siberia JJA T2m [C]', 'Siberia DJF soil [C]',
       'Siberia JJA soil [C]', 'SO SW CRE [W/m2]', 'global T2m [C]')


def years_avail(root, var='2t'):
    D = f'{root}/outdata/oifs'
    if not os.path.isdir(D):
        return []
    ys = []
    for f in os.listdir(D):
        if f.startswith(f'atm_remapped_1m_{var}_') and f.endswith('.nc'):
            try:
                ys.append(int(f.split('_')[-1].split('-')[0]))
            except ValueError:
                pass
    return sorted(ys)


def load_year(root, var, y):
    D = f'{root}/outdata/oifs'
    div = ACC if var in FLUX else 1.0
    for f in (f'{D}/atm_remapped_1m_{var}_{y}-{y}.nc',
              f'{D}/atm_remapped_1m_{var}_1m_{y}-{y}.nc'):
        if os.path.exists(f):
            with xr.open_dataset(f, decode_times=False) as d:
                a = d[var].values / div
                lat, lon = d['lat'].values, d['lon'].values
            return (a, lat, lon) if a.shape[0] == 12 else (None, None, None)
    return None, None, None


def gmean(f2d, lat):
    return float(np.average(f2d.mean(axis=1), weights=np.cos(np.deg2rad(lat))))


def boxmean(f2d, lat, lon, box, lsm=None):
    la0, la1, lo0, lo1 = box
    ys = (lat >= la0) & (lat <= la1)
    xs = ((lon % 360) >= lo0) & ((lon % 360) <= lo1)
    sub = f2d[np.ix_(ys, xs)]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape).copy()
    if lsm is not None:
        w = np.where(lsm[np.ix_(ys, xs)] >= 0.5, w, 0.0)
    return float(np.average(sub, weights=w)) if w.sum() else np.nan


def zband(f2d, lat, a, b):
    sel = (lat >= a) & (lat < b)
    return float(np.average(f2d[sel, :].mean(axis=1),
                            weights=np.cos(np.deg2rad(lat[sel]))))


def metrics_for_year(root, y, lsm, lat, lon):
    t2m, _, _ = load_year(root, '2t', y)
    tsr, _, _ = load_year(root, 'tsr', y)
    ttr, _, _ = load_year(root, 'ttr', y)
    if t2m is None or tsr is None or ttr is None:
        return None
    tsrc, _, _ = load_year(root, 'tsrc', y)
    stl2, _, _ = load_year(root, 'stl2', y)
    d = [m - 1 for m in DJF]
    j = [m - 1 for m in JJA]
    out = {
        'net TOA [W/m2]': gmean((tsr + ttr).mean(axis=0), lat),
        'SO SW CRE [W/m2]': (zband((tsr - tsrc).mean(axis=0), lat, *SO)
                             if tsrc is not None else np.nan),
        'global T2m [C]': gmean(t2m.mean(axis=0), lat) - 273.15,
        'Siberia JJA T2m [C]': boxmean(t2m[j].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
        'Siberia DJF T2m [C]': boxmean(t2m[d].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
    }
    if stl2 is not None:
        out['Siberia DJF soil [C]'] = boxmean(stl2[d].mean(axis=0), lat, lon, SIB, lsm) - 273.15
        out['Siberia JJA soil [C]'] = boxmean(stl2[j].mean(axis=0), lat, lon, SIB, lsm) - 273.15
    return out


def series_for(root, years):
    lsm, lat, lon = load_year(root, 'lsm', years[0])
    if lsm is not None and lsm.ndim == 3:
        lsm = lsm[0]
    out = {}
    for y in years:
        m = metrics_for_year(root, y, lsm, lat, lon)
        if m:
            out[y] = m
    return out


def fpc_boxmean(root, cols):
    """Box-mean of the named fpc.out columns in the last available decade.

    LPJ-GUESS writes fpc.out per rank as the run proceeds, outside XIOS, so this
    survives a finalisation abort that destroys the netCDF output.
    """
    legs = sorted(glob.glob(f'{root}/outdata/lpj_guess/*'))
    for leg in reversed(legs):
        files = sorted(glob.glob(f'{leg}/run*/fpc.out'))
        if not files:
            continue
        acc, hdr = {c: [] for c in cols}, None
        for fn in files:
            try:
                with open(fn) as fh:
                    hdr = fh.readline().split()
                    idx = {c: hdr.index(c) for c in cols if c in hdr}
                    if not idx:
                        continue
                    last = {}
                    for line in fh:
                        p = line.split()
                        if len(p) < len(hdr):
                            continue
                        lon, la, yr = float(p[0]), float(p[1]), int(p[2])
                        if not (SIB[0] <= la <= SIB[1] and SIB[2] <= (lon % 360) <= SIB[3]):
                            continue
                        last[(lon, la)] = (yr, p)
                    for (lon, la), (yr, p) in last.items():
                        for c, i in idx.items():
                            acc[c].append(float(p[i]))
            except Exception:
                continue
        if any(acc[c] for c in cols):
            return os.path.basename(leg), {c: (float(np.mean(acc[c])) if acc[c] else np.nan)
                                           for c in cols}, len(acc[cols[0]])
    return None, {c: np.nan for c in cols}, 0


print(__doc__)
print('=' * 100)

avail = {t: years_avail(r) for t, r, _ in CHAIN}
print('\nARMS\n')
for t, r, note in CHAIN:
    ys = avail[t]
    print(f'  {t}  {len(ys):3d} yr  {ys[0]}-{ys[-1] if ys else "-"}   {note}')

# 1370-79 is excluded wherever 11G is involved: the repaired soil map slipped into its
# leg 3 alone (organic cells 233 vs 547 in its other four legs).
CLEAN40 = list(range(1350, 1370)) + list(range(1380, 1400))
COMMON = [y for y in CLEAN40 if all(y in avail[t] for t in ('11E', '11G', '11I'))]
COMMON_J = [y for y in range(1350, 1380) if all(y in avail[t] for t in
                                                ('11E', '11G', '11I', '11J'))]
COMMON_J = [y for y in COMMON_J if y < 1370]   # keep 11G clean here too

ser = {t: series_for(r, sorted(set(COMMON) | set(COMMON_J))) for t, r, _ in CHAIN}

print('\n' + '=' * 100)
print(f'\n1. THE STACK, absolute means over the 4-arm common window '
      f'({len(COMMON_J)} yr, {COMMON_J[0]}-{COMMON_J[-1]})\n')
print(f'  {"metric":24s}' + ''.join(f'{t:>11s}' for t, _, _ in CHAIN))
for k in KEY:
    row = f'  {k:24s}'
    for t, _, _ in CHAIN:
        v = [ser[t][y][k] for y in COMMON_J if y in ser[t] and k in ser[t][y]]
        row += f'{np.mean(v):11.3f}' if v else f'{"-":>11s}'
    print(row)

print('\n' + '=' * 100)
print('\n2. INCREMENTS, paired.  Each arm branches from the same 1350 state, so the')
print('   difference series is PAIRED: report its own mean and scatter, the decade')
print('   breakdown, and whether it is still growing.\n')

for ctl, trt, window in (('11E', '11G', COMMON), ('11G', '11I', COMMON),
                         ('11I', '11J', COMMON_J)):
    yrs = [y for y in window if y in ser[ctl] and y in ser[trt]]
    print(f'  --- {trt} minus {ctl},  {len(yrs)} yr {yrs[0]}-{yrs[-1]} ---')
    print(f'  {"metric":24s} {"mean d":>9s} {"sd(d)":>8s} {"t":>7s} {"unpaired thr":>13s}'
          f' {"1st half":>9s} {"2nd half":>9s} {"trend/dec":>10s}')
    for k in KEY:
        a = np.array([ser[ctl][y][k] for y in yrs if k in ser[ctl][y]])
        b = np.array([ser[trt][y][k] for y in yrs if k in ser[trt][y]])
        if a.size != b.size or a.size < 4:
            continue
        d = b - a
        n = d.size
        sd = d.std(ddof=1)
        t = d.mean() / (sd / np.sqrt(n)) if sd > 0 else np.nan
        thr_unp = 1.96 * a.std(ddof=1) * np.sqrt(2.0 / n)
        h = n // 2
        yy = np.array(yrs, dtype=float)
        slope = np.polyfit(yy, d, 1)[0] * 10.0
        mark = '*' if abs(t) > 2.0 else ' '
        print(f'  {k:24s} {d.mean():+9.3f}{mark} {sd:8.3f} {t:7.2f} {thr_unp:13.3f}'
              f' {d[:h].mean():+9.3f} {d[h:].mean():+9.3f} {slope:+10.3f}')
    print()

print('=' * 100)
print('\n3. AGAINST THE TARGETS, not against each other\n')
try:
    era = np.loadtxt(ERA5_SIB)
    era_jja = float(era[:, 1].mean()) - 273.15
    print(f'  ERA5 Siberian JJA T2m, {int(era[0,0])}-{int(era[-1,0])} mean: {era_jja:.3f} C')
    for t, _, _ in CHAIN:
        v = [ser[t][y]['Siberia JJA T2m [C]'] for y in COMMON_J if y in ser[t]]
        if v:
            print(f'    {t}  {np.mean(v):7.3f}   bias {np.mean(v) - era_jja:+7.3f} K')
    print('    (ERA5 is a modern-era reference and these are piControl arms, so the')
    print('     absolute offset carries a forcing difference; use it for RANKING.)')
except Exception as exc:
    print(f'  (ERA5 reference unavailable: {exc})')

try:
    with xr.open_dataset(CERESF) as cds:
        clat = cds['lat'].values
        csel = (clat >= SO[0]) & (clat < SO[1])
        so_obs = float(np.average(
            cds['toa_cre_sw_clim'].values.mean(axis=0)[csel, :].mean(axis=1),
            weights=np.cos(np.deg2rad(clat[csel]))))
    print(f'\n  CERES SO SW CRE: {so_obs:.2f} W/m2')
    for t, _, _ in CHAIN:
        v = [ser[t][y]['SO SW CRE [W/m2]'] for y in COMMON_J
             if y in ser[t] and 'SO SW CRE [W/m2]' in ser[t][y]]
        if v:
            print(f'    {t}  {np.mean(v):8.3f}   gap {np.mean(v) - so_obs:+7.3f}')
except Exception as exc:
    print(f'  (CERES unavailable: {exc})')

print('\n' + '=' * 100)
print('\n4. IS THE FOREST THERE?  Siberian box mean of fpc.out, last decade on disk.')
print('   This is the campaign\'s actual objective; temperature is the mechanism.\n')
COLS = ['BNE', 'BINE', 'BNS', 'IBS', 'C3G', 'TREEFPC', 'FORESTFPC', 'GRASSFPC', 'Total']
print(f'  {"arm":5s} {"leg":22s} {"cells":>6s} ' + ''.join(f'{c:>9s}' for c in COLS))
for t, r, _ in CHAIN:
    leg, vals, ncell = fpc_boxmean(r, COLS)
    if leg is None:
        print(f'  {t:5s} {"no fpc.out found":22s}')
        continue
    print(f'  {t:5s} {leg:22s} {ncell:6d} ' + ''.join(f'{vals[c]:9.3f}' for c in COLS))
print('\n  BNS = boreal needleleaf summergreen (larch), the east-Siberian dominant.')
print('  Standing finding: BNS is NOT climate-gate limited, it is outcompeted by C3G.')
print()
