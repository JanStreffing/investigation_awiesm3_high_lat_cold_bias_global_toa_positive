"""Round 33, recomputed over OCEAN ONLY.

WHY.  The round-33 numbers were area-averaged over whole latitude bands, which folds the
Antarctic ICE SHEET into the 90-60S band.  The ice sheet is brilliant white, radiatively
quite different from the sea-ice zone, and it cancels part of the signal: over the whole
cap the surface net SW bias reads -0.57 W/m2, while over ocean alone it is +7.99.  The
sea-ice question is about the ocean, so the ocean mask is the right domain and the
committed numbers understate the effect.

Recomputes the three tables that the ice conclusion rests on:
  (a) net TOA bias by band and its share of the global mean
  (b) SW CRE bias per campaign arm, polar and sub-Antarctic
  (c) the opacity split (CRE per unit cover) and the surface-side check
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')
from runs import OBS

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
Y0, Y1 = 1380, 1389
ARMS = [('11E base', f'{R}/Tuning_test_11E_swemin15_K1'),
        ('11G +S4', f'{R}/Tuning_test_11G_inppmin50k'),
        ('11N LX4', f'{R}/11N'), ('11P LX4 1990', f'{R}/11P'),
        ('11Q +RSBLB', f'{R}/11Q'), ('11R +RSBLB 1990', f'{R}/11R')]
DJF, JJA, ANN = [11, 0, 1], [5, 6, 7], list(range(12))

f = sorted(glob.glob(f'{R}/11P/outdata/oifs/atm_remapped_1m_lsm_*.nc'))[0]
with xr.open_dataset(f, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values)
    m = m[0] if m.ndim == 3 else m
    lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
ocean = m <= 0.5
RE = 6.371e6
A = np.broadcast_to((RE**2 * np.cos(np.deg2rad(lat)) * np.deg2rad(abs(lat[1] - lat[0]))
                     * 2 * np.pi / m.shape[1])[:, None], m.shape)
GLOB = A.sum()


def load(p, var):
    out = []
    for y in range(Y0, Y1 + 1):
        fp = f'{p}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(fp):
            return None
        with xr.open_dataset(fp, decode_times=False) as d:
            out.append(np.squeeze(d[[k for k in d.data_vars
                                     if 'bnds' not in k and 'bounds' not in k][0]].values))
    return np.stack(out).mean(0)


o = xr.open_dataset(OBS)


def ceres(v, mo):
    a = o[v]
    a = a.isel({a.dims[0]: mo}).mean(a.dims[0])
    la = [x for x in a.dims if 'lat' in x][0]
    lo = [x for x in a.dims if 'lon' in x][0]
    return np.asarray(a.interp({la: ('lat', lat), lo: ('lon', lon % 360)}).values, float)


def am(x, s):
    k = s & np.isfinite(x)
    return float(np.average(x[k], weights=A[k]))


print('=' * 84)
print('(a) net TOA bias vs CERES, 11R, OCEAN ONLY.  clr_t is the model-comparable')
print('    clear sky (total-region); clr_c is cloud-free-scene only.\n')
tsr, tsrc, ttr = (load(f'{R}/11R', v) for v in ('tsr', 'ttr', 'ttr'))
tsr, tsrc, ttr = load(f'{R}/11R', 'tsr'), load(f'{R}/11R', 'tsrc'), load(f'{R}/11R', 'ttr')
BANDS = [('90-60S', -90, -60), ('60-45S', -60, -45), ('45-30S', -45, -30),
         ('30S-30N', -30, 30), ('30-45N', 30, 45), ('45-60N', 45, 60), ('60-90N', 60, 90)]
print(f'{"band":>9} {"ocean%":>7} {"net bias":>9} {"contrib":>9} {"SWCRE bias":>11}')
tot = 0.0
for nm, lo, hi in BANDS:
    s = ocean & np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], m.shape)
    on = ceres('toa_net_all_clim', ANN)
    k = s & np.isfinite(on)
    fr = A[k].sum() / GLOB
    b = am((tsr + ttr).mean(0) / 3600. - on, k)
    c = am((tsr - tsrc).mean(0) / 3600.
           - (ceres('toa_sw_clr_t_clim', ANN) - ceres('toa_sw_all_clim', ANN)), k)
    tot += b * fr
    print(f'{nm:>9} {100*fr:7.1f} {b:+9.2f} {b*fr:+9.3f} {c:+11.2f}')
print(f'{"ocean tot":>9} {"":>7} {"":>9} {tot:+9.3f}')

print('\n' + '=' * 84)
print('(b) SW CRE bias vs CERES clr_t, OCEAN ONLY, by arm\n')
print(f'{"arm":>17} | {"90-60S DJF":>11} {"90-60S ANN":>11} | {"60-45S DJF":>11} {"60-45S ANN":>11}')
for lab, p in ARMS:
    a1, a2 = load(p, 'tsr'), load(p, 'tsrc')
    if a1 is None:
        print(f'{lab:>17} | incomplete'); continue
    line = f'{lab:>17} |'
    for lo, hi in ((-90, -60), (-60, -45)):
        s = ocean & np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], m.shape)
        for mo in (DJF, ANN):
            oc = ceres('toa_sw_clr_t_clim', mo) - ceres('toa_sw_all_clim', mo)
            k = s & np.isfinite(oc)
            line += f' {am((a1[mo].mean(0)-a2[mo].mean(0))/3600. - oc, k):+11.2f}'
        line += ' |'
    print(line)

print('\n' + '=' * 84)
print('(c) opacity split and the surface check, OCEAN ONLY, DJF/JJA summer\n')
print(f'{"band":>14} {"seas":>5} | {"CRE":>8} {"cover":>7} {"CRE/cov":>8} | '
      f'{"sfc netSW":>10} {"CERES":>8} {"bias":>8}')
for nm, lo, hi, mo, sn in (('60-90S ocean', -90, -60, DJF, 'DJF'),
                           ('50-70S ocean', -70, -50, DJF, 'DJF'),
                           ('60-90N ocean', 60, 90, JJA, 'JJA')):
    s = ocean & np.broadcast_to(((lat >= lo) & (lat < hi))[:, None], m.shape)
    tsr_, tsrc_, tcc_, ssr_ = (load(f'{R}/11R', v) for v in ('tsr', 'tsrc', 'tcc', 'ssr'))
    ocre = ceres('toa_sw_clr_t_clim', mo) - ceres('toa_sw_all_clim', mo)
    ocov = ceres('cldarea_total_daynight_clim', mo) / 100.
    osfc = ceres('sfc_net_sw_all_clim', mo)
    k = s & np.isfinite(ocre) & np.isfinite(ocov) & (ocov > 0.05)
    mc = am((tsr_[mo].mean(0) - tsrc_[mo].mean(0)) / 3600., k)
    mv = am(tcc_[mo].mean(0), k)
    oc, ov = am(ocre, k), am(ocov, k)
    print(f'{nm:>14} {sn:>5} | {mc:8.2f} {mv:7.3f} {mc/mv:8.1f} | '
          f'{am(ssr_[mo].mean(0)/3600., k):10.1f} {am(osfc, k):8.1f} '
          f'{am(ssr_[mo].mean(0)/3600. - osfc, k):+8.2f}')
    print(f'{"":>14} {"CERES":>5} | {oc:8.2f} {ov:7.3f} {oc/ov:8.1f} |')
