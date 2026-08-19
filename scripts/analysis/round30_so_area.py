"""Round 30 scored on its pre-registered terms: SO cloud AREA, and the tropical cost.

THE TWO QUESTIONS, both of which the campaign had assumed rather than measured.

  N1/N2  RCL_INPSEA 0.2 -> 0.10 / 0.05.  Has the only spatially selective SO lever any
         headroom?  It sat at 0.2 in all 46 runscripts that set it and was never scanned.
  W1     RCLCRIT_SEA 2.5e-4 -> 6.0e-4.  This IS round 23's S3, designed and then CANCELLED
         on an inference: by analogy with B3/RCLDIFF, a GLOBAL cloud-erosion term, it
         "would cool the tropics".  RCLCRIT_SEA acts on stratiform warm-rain autoconversion
         over sea only, which is the SO's regime while tropical rain is largely convective.

SCORED ON AREA, NOT CRE.  Round 23's explicit instruction.  CRE conflates the two thirds
of the SO error -- round 22 split it AMOUNT +4.87 (65.5 %) / OPACITY +2.57 (34.5 %) -- and
the INP branch already looked good on CRE while moving area only +0.63 pp (D2a) and
+0.43 pp (D2b) of a 6.43 pp deficit, inside the noise floor.  A lever that improves CRE
without improving area has bought opacity again and is not what is needed.

REFERENCE HYGIENE, round 23's standing rule: score area against CERES and MODIS, NOT
against ERA5, which is not an independent cloud reference for this model -- ERA5 says
81.6 % against the model's 83.1 % while CERES and MODIS say 89.5 and 89.3.

DISQUALIFIER, pre-registered: tropical SW CRE or tropical net beyond +-0.5, the tolerance
that bound DMS.  Siberia is scored too but should be ~0 for all three by construction --
all three levers branch on PLSM and cannot act over land.  A non-zero Siberian term would
mean the selectivity argument is wrong, so it is a check on the reasoning, not on the arm.

Accumulated fluxes: /3600, verified empirically -- global ASR comes out 240.4 W/m2.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

RT = '/work/bb1469/a270092/runtime/oifsamip-cy48'
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'
MODISF = '/work/ab0246/a270092/obs/MODIS/clt_MODIS_yearmean.nc'
ACC = 3600.0
Y0, Y1 = 1872, 1915
SO = (-65.0, -45.0)
TROP = (-20.0, 20.0)
SIB = (55.0, 75.0, 60.0, 180.0)
DJF, JJA = [12, 1, 2], [6, 7, 8]
FLUX = ('tsr', 'ttr', 'tsrc', 'ttrc')

ARMS = [('S4 control', 'amip_S4_inppmin50000'),
        ('N1 inpsea0.10', 'N1'),
        ('N2 inpsea0.05', 'N2'),
        ('W1 clcrit6e-4', 'W1')]


def load(root, var, y):
    # AMIP monthly output repeats the frequency before the year --
    # atm_remapped_1m_tsr_1m_1872-1872.nc -- while the coupled tree does not.
    # Getting this wrong silently yields zero years, which is how it first ran.
    for pat in (f'{RT}/{root}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc',
                f'{RT}/{root}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc',
                f'{RT}/{root}/outdata/oifs/atm_remapped_1d_{var}_1d_{y}-{y}.nc'):
        if glob.glob(pat):
            try:
                with xr.open_dataset(pat, decode_times=False) as d:
                    a = np.asarray(d[var].values, dtype=float)
                    if var in FLUX:
                        a = a / ACC
                    return a, d['lat'].values, d['lon'].values
            except Exception:
                return None, None, None
    return None, None, None


def zband(f2d, lat, a, b):
    s = (lat >= a) & (lat < b)
    w = np.cos(np.deg2rad(lat[s]))
    return float(np.average(f2d[s, :].mean(axis=1), weights=w))


def boxmean(f2d, lat, lon, box):
    la0, la1, lo0, lo1 = box
    ys = (lat >= la0) & (lat <= la1)
    xs = ((lon % 360) >= lo0) & ((lon % 360) <= lo1)
    sub = f2d[np.ix_(ys, xs)]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape)
    return float(np.average(sub, weights=w))


def series(root):
    out = {k: [] for k in ('SO area', 'SO SW CRE', 'trop SW CRE', 'trop net',
                           'net TOA', 'Sib JJA', 'Sib DJF')}
    for y in range(Y0, Y1 + 1):
        tsr, lat, lon = load(root, 'tsr', y)
        tsrc, _, _ = load(root, 'tsrc', y)
        ttr, _, _ = load(root, 'ttr', y)
        tcc, _, _ = load(root, 'tcc', y)
        t2m, _, _ = load(root, '2t', y)
        if tsr is None or tsrc is None or ttr is None:
            continue
        cre = (tsr - tsrc).mean(axis=0)
        net = (tsr + ttr).mean(axis=0)
        out['SO SW CRE'].append(zband(cre, lat, *SO))
        out['trop SW CRE'].append(zband(cre, lat, *TROP))
        out['trop net'].append(zband(net, lat, *TROP))
        out['net TOA'].append(float(np.average(net.mean(axis=1),
                                               weights=np.cos(np.deg2rad(lat)))))
        if tcc is not None:
            out['SO area'].append(zband(tcc.mean(axis=0), lat, *SO) * 100.0)
        if t2m is not None and t2m.shape[0] == 12:
            d = [m - 1 for m in DJF]
            j = [m - 1 for m in JJA]
            out['Sib JJA'].append(boxmean(t2m[j].mean(axis=0), lat, lon, SIB) - 273.15)
            out['Sib DJF'].append(boxmean(t2m[d].mean(axis=0), lat, lon, SIB) - 273.15)
    return {k: np.array(v) for k, v in out.items() if v}


def main():
    print(__doc__)
    print('=' * 96)
    S = {}
    for tag, root in ARMS:
        S[tag] = series(root)
        n = len(S[tag].get('SO SW CRE', []))
        print(f'  {tag:16s} {n:3d} years')
    if not S[ARMS[0][0]].get('SO SW CRE', np.array([])).size:
        print('\n  ABORT: zero years loaded for the control. That is a path bug, not a\n'
              '  result -- do not read an empty table as agreement between arms.')
        return
    base = ARMS[0][0]
    keys = ['SO area', 'SO SW CRE', 'trop SW CRE', 'trop net', 'net TOA',
            'Sib JJA', 'Sib DJF']

    print(f'\nDETECTION THRESHOLDS from {base} interannual scatter, 1.96*sd*sqrt(2/n)\n')
    thr = {}
    for k in keys:
        v = S[base].get(k)
        if v is None or not len(v):
            continue
        thr[k] = 1.96 * v.std(ddof=1) * np.sqrt(2.0 / len(v))
        print(f'  {k:14s} sd {v.std(ddof=1):8.4f}   threshold +-{thr[k]:.4f}')

    # observed anchors
    try:
        with xr.open_dataset(CERESF) as c:
            cl = c['lat'].values
            s = (cl >= SO[0]) & (cl < SO[1])
            w = np.cos(np.deg2rad(cl[s]))
            ceres_cre = float(np.average(
                c['toa_cre_sw_clim'].values.mean(axis=0)[s, :].mean(axis=1), weights=w))
    except Exception:
        ceres_cre = np.nan

    print(f'\n{"metric":14s} ' + ' '.join(f'{a[0]:>15s}' for a in ARMS))
    for k in keys:
        if k not in thr:
            continue
        row = f'  {k:14s}'
        for tag, _ in ARMS:
            v = S[tag].get(k)
            if v is None or not len(v):
                row += f'{"--":>16s}'
                continue
            m = v.mean()
            if tag == base:
                row += f'{m:16.3f}'
            else:
                d = m - S[base][k].mean()
                row += f'{m:11.3f}{"*" if abs(d) > thr[k] else " "}({d:+.2f})'[:16].rjust(16)
        print(row)

    print(f'\n  CERES SO SW CRE {ceres_cre:.2f};  CERES area 89.5 %, MODIS 89.3 %, '
          f'model ~83.1 %  (ERA5 81.6 %, NOT a valid reference here)')
    print('\nVERDICT PER ARM\n')
    for tag, _ in ARMS[1:]:
        area = S[tag].get('SO area')
        d_area = area.mean() - S[base]['SO area'].mean() if area is not None and len(area) else np.nan
        d_cre = S[tag]['SO SW CRE'].mean() - S[base]['SO SW CRE'].mean()
        d_tsw = S[tag]['trop SW CRE'].mean() - S[base]['trop SW CRE'].mean()
        d_tn = S[tag]['trop net'].mean() - S[base]['trop net'].mean()
        d_sj = (S[tag]['Sib JJA'].mean() - S[base]['Sib JJA'].mean()
                if 'Sib JJA' in S[tag] else np.nan)
        trop_fail = abs(d_tsw) > 0.5 or abs(d_tn) > 0.5
        print(f'  {tag}:  d(area) {d_area:+.3f} pp of the 6.4 pp deficit,  '
              f'd(SO CRE) {d_cre:+.3f}')
        print(f'      tropics SW {d_tsw:+.3f}, net {d_tn:+.3f}  -> '
              f'{"DISQUALIFIED on the tropics" if trop_fail else "tropics within +-0.5"}')
        print(f'      Siberia JJA {d_sj:+.3f} K  (expected ~0: the lever branches on PLSM '
              f'and cannot act over land)')
        if not np.isnan(d_area) and abs(d_area) < 0.7 and d_cre < -0.5:
            print('      -> buys OPACITY, not AREA: the same trap the INP branch fell into.')
        print()


if __name__ == '__main__':
    main()
