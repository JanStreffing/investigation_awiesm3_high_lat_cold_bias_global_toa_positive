"""Paired radiation change maps: run A minus run B over the same years.

Panels: (a) net TOA, (b) SW cloud radiative effect, (c) LW cloud radiative effect, (d) net energy
flux into the surface (ssr+str+slhf+sshf minus snowfall latent heat, sf*1000*333550), all in
W/m2 with the cos-lat weighted global mean in the title.  Years come from outdata or a leg's work
directory (moved-aside legs skipped), so a control leg can be given through its own ROOT_B.
Default: PI200 with RCL_INPPMIN 50000 (S4) from 1500 against the cancelled leg that ran the
same years with 70000.

Usage:  YEARS=1500:1505 python3 scripts/figures/radiation_delta_maps.py
        ROOT_A=... ARM_A=PI200 ROOT_B=... ARM_B=PI200 LABEL="..." python3 ...
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCR = '/work/bb1469/a270092/runtime/PI200_noS4_view'   # symlink view of the cancelled no-S4 leg
RA = os.environ.get('ROOT_A', '/work/bb1469/a270092/runtime/awiesm3-v3.4'); AA = os.environ.get('ARM_A', 'PI200')
RB = os.environ.get('ROOT_B', SCR); AB = os.environ.get('ARM_B', 'PI200')
LABEL = os.environ.get('LABEL', 'PI200 with S4 minus without S4, same years')
Y0, Y1 = (int(x) for x in os.environ.get('YEARS', '1500:1505').split(':'))
LSMF = ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
        'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')


def path(root, arm, var, y):
    work = [q for q in sorted(glob.glob(f'{root}/{arm}/run_*/work/atm_remapped_1m_{var}_{y}-{y}.nc'))
            if os.path.basename(os.path.dirname(os.path.dirname(q))).count('.') == 0]
    for p in [f'{root}/{arm}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'] + work:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f'{root}/{arm} {var} {y}')


def fields(root, arm):
    acc = None
    for y in range(Y0, Y1 + 1):
        f = {}
        for v in ('tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str', 'slhf', 'sshf', 'sf'):
            with xr.open_dataset(path(root, arm, v, y), decode_times=False) as d:
                k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
                a = np.squeeze(d[k].values).astype('f8'); assert a.shape[0] == 12, (arm, v, y)
                f[v] = a.mean(0) / 3600.0; lat, lon = np.squeeze(d['lat'].values), np.squeeze(d['lon'].values)
        out = {'net TOA': f['tsr'] + f['ttr'], 'SW CRE': f['tsr'] - f['tsrc'], 'LW CRE': f['ttr'] - f['ttrc'],
               'net surface energy': f['ssr'] + f['str'] + f['slhf'] + f['sshf'] - f['sf'] * 1000 * 333550.0}
        acc = out if acc is None else {k: acc[k] + out[k] for k in acc}
    n = Y1 - Y0 + 1
    return {k: v / n for k, v in acc.items()}, lat, lon


FA, lat, lon = fields(RA, AA); FB, _, _ = fields(RB, AB)
with xr.open_dataset(LSMF, decode_times=False) as d:
    m = np.squeeze(d['lsm'].values); m = m[0] if m.ndim == 3 else m
o = np.argsort(lon % 360); LON = (lon % 360)[o]; land = (m > 0.5)[:, o]
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], land.shape)
fig, axes = plt.subplots(4, 1, figsize=(9, 13), dpi=150)
for ax, k, lim in zip(axes, ('net TOA', 'SW CRE', 'LW CRE', 'net surface energy'), (6, 6, 3, 6)):
    D = (FA[k] - FB[k])[:, o]; gm = float(np.average(D, weights=W))
    so = float(np.average(D[(lat >= -65) & (lat <= -45)], weights=W[(lat >= -65) & (lat <= -45)]))
    print(f'{k:<20} global {gm:+.3f}   45-65S {so:+.3f} W/m2')
    im = ax.pcolormesh(LON, lat, D, cmap='RdBu_r', vmin=-lim, vmax=lim, shading='auto')
    ax.contour(LON, lat, land.astype(float), levels=[0.5], colors='k', linewidths=0.4)
    ax.set_title(f'({"abcd"[["net TOA","SW CRE","LW CRE","net surface energy"].index(k)]}) {k}: global {gm:+.2f}, 45-65S {so:+.2f} W m$^{{-2}}$', loc='left', fontsize=10)
    ax.set_ylabel('lat'); fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label='W m$^{-2}$')
fig.suptitle(f'{LABEL}, {Y0}-{Y1} mean (positive = more energy into the system / surface)', x=0.01, ha='left', fontsize=11)
fig.tight_layout()
out = os.path.join(REPO, 'plots', f'radiation_delta_{AA}_S4_{Y0}-{Y1}.png'); fig.savefig(out); print('saved', out)
