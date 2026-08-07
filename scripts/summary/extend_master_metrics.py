"""Extend the master metric table and scatter with the 09 and 10 series.

The table and Figure 2 covered only the coupled Tuning_test_06* family.  Rounds 09
and 10 branch from 06V and are the runs the campaign now actually turns on, so they
belong on the same axes.

METHOD, matched to the existing rows so the numbers are comparable:
  * 2m-T RMSD/bias -- part8 recipe exactly: cdo remapcon to r360x180, time mean over
    the run's own 30-yr window, then sqrt(cos-lat) weighting, model minus the SAME
    cached era5_clim.nc the 06 rows used (plot_t2m_bias.py wrmsd/wmean).
  * net TOA -- 30-yr mean of (tsr+ttr)/accumulation, the recipe already validated
    against the cached campaign CSV to 0.001 W/m2 in plot_net_toa_with_10series.py.

WINDOWS.  The 06 rows use 1354-1379.  09 runs 1350-1379 and 10 runs 1350-1399, so the
last 26 years (1354-1379) is the common window and is used for every new row; 10A/10B
also get their own 1370-1399 value reported separately, because the extra 20 years
change the answer and hiding that would overstate comparability.

CMPI is NOT computed here -- it needs part4_cmpi.py, which is a separate tool with its
own reference set.  Those cells stay blank rather than being filled with a lookalike.
"""
import os, csv, math, subprocess, tempfile
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

BASE = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive'
TOOL = '/work/ab0246/a270092/software/release_evaluation_tool2/output/Tuning_test_06_overview'
ERA = f'{TOOL}/t2m/era5_clim.nc'
STEP = 3600.0
R270 = '/work/bb1469/a270270/runtime/awiesm3-v3.4'
R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'

NEW = [
    ('080a', f'{R270}/Tuning_test_080a_lpjguess_Baseline_coupled_fromCRUNCEP', (1354, 1379)),
    ('09A',  f'{R270}/Tuning_test_09A_lpjguess_Baseline_coupled_fromCRUNCEP_newSeaIce', (1354, 1379)),
    ('09B',  f'{R270}/Tuning_test_09B_06T_1hCPL_MOSPP_KPPLOW_CRUNCEPinit_newSeaIce', (1354, 1379)),
    ('09C',  f'{R092}/Tuning_test_09C_06V_CRUNCEPinit_newSeaIce', (1354, 1379)),
    ('10A',  f'{R270}/Tuning_test_10A_06V_G4_CRUNCEP_plus_CERES_init_newSeaIce', (1354, 1379)),
    ('10B',  f'{R270}/Tuning_test_10B_06V_G4_snowDepletion_CRUNCEP_plus_CERES_init_newSeaIce', (1354, 1379)),
    ('10A_late', f'{R270}/Tuning_test_10A_06V_G4_CRUNCEP_plus_CERES_init_newSeaIce', (1370, 1399)),
    ('10B_late', f'{R270}/Tuning_test_10B_06V_G4_snowDepletion_CRUNCEP_plus_CERES_init_newSeaIce', (1370, 1399)),
]

era = xr.open_dataset(ERA)
evar = [v for v in era.data_vars if 't2m' in v.lower()][0]
ERA5 = np.squeeze(era[evar].values)
LAT = era['lat'].values
W = np.broadcast_to(np.sqrt(np.cos(np.deg2rad(LAT)))[:, None], ERA5.shape)


def wrmsd(d):
    return math.sqrt(np.sum(W * d ** 2) / np.sum(W))


def wmean(d):
    return float(np.sum(W * d) / np.sum(W))


def files(root, var, y0, y1):
    d = f'{root}/outdata/oifs'
    out = []
    for y in range(y0, y1 + 1):
        for pat in (f'atm_remapped_1m_{var}_1m_{y}-{y}.nc', f'atm_remapped_1m_{var}_{y}-{y}.nc'):
            if os.path.exists(f'{d}/{pat}'):
                out.append(f'{d}/{pat}'); break
    return out


def t2m_metrics(root, y0, y1):
    fs = files(root, '2t', y0, y1)
    if not fs:
        return None, None, 0
    with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tf:
        tmp = tf.name
    cmd = f'cdo -s -O timmean -remapcon,r360x180 -cat {" ".join(fs)} {tmp}'
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print('   cdo failed:', r.stderr.strip()[:200]); return None, None, 0
    with xr.open_dataset(tmp) as d:
        v = [x for x in d.data_vars if x in ('2t', 't2m')][0]
        m = np.squeeze(d[v].values)
    os.unlink(tmp)
    dif = m - ERA5
    return wrmsd(dif), wmean(dif), len(fs)


def toa_mean(root, y0, y1):
    ts, n = [], 0
    for y in range(y0, y1 + 1):
        f1 = files(root, 'tsr', y, y); f2 = files(root, 'ttr', y, y)
        if not f1 or not f2:
            continue
        with xr.open_dataset(f1[0], decode_times=False) as a:
            tsr = a['tsr'].values; lat = a['lat'].values
        with xr.open_dataset(f2[0], decode_times=False) as b:
            ttr = b['ttr'].values
        net = (tsr + ttr) / STEP
        w = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], net.shape[1:])
        ts.append(float(np.average(net.mean(axis=0), weights=w))); n += 1
    return (float(np.mean(ts)), n) if ts else (None, 0)


print(__doc__)
print('=' * 92)
print(f'{"run":10s}{"window":>12s}{"nyr":>5s}{"netTOA":>9s}{"RMSD":>8s}{"bias":>8s}')
rows = []
for lab, root, (y0, y1) in NEW:
    toa, nt = toa_mean(root, y0, y1)
    rm, bi, n2 = t2m_metrics(root, y0, y1)
    if toa is None or rm is None:
        print(f'{lab:10s}{f"{y0}-{y1}":>12s}   incomplete (TOA {nt} yr, 2t {n2} yr)'); continue
    rows.append((lab, f'{y0}-{y1}', nt, toa, rm, bi))
    print(f'{lab:10s}{f"{y0}-{y1}":>12s}{nt:5d}{toa:+9.3f}{rm:8.3f}{bi:+8.3f}')

out = f'{BASE}/data/MASTER_metrics_0910.csv'
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['run', 'window', 'nyr', 'TOA_imbalance_Wm2', 't2m_RMSD_K', 't2m_bias_K'])
    for r in rows:
        w.writerow([r[0], r[1], r[2], round(r[3], 3), round(r[4], 3), round(r[5], 3)])
print('\nSaved:', out)
print('\nBaseline for reference (06_Baseline): TOA 1.531, RMSD 1.717, bias -0.939')
