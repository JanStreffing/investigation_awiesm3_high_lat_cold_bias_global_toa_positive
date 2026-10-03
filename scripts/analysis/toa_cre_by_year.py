"""Year-by-year TOA radiation of one run: net TOA, cloud radiative effect, Southern Ocean shortwave.

Per year (12 complete months of OpenIFS remapped output, from outdata or a running leg's work
directory): global net TOA, global SW and LW cloud radiative effect, SW CRE 45-65S, planetary
albedo, and at the surface: net radiation (ssr+str), net energy flux into the surface
(ssr+str+slhf+sshf minus the latent heat of snowfall, sf*3.3355e8), global and 45-65S.  The reference window REF gives the mean and the year-to-year standard deviation, so
single years after a change can be judged against the spread.  Accumulated J/m2 per hourly
step, divided by 3600; incoming solar asserted near 340 W/m2.

Usage:  ARM=PI200 REF=1490:1499 YEARS=1500:1509 python3 scripts/analysis/toa_cre_by_year.py
"""
import os, glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
R = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM = os.environ.get('ARM', 'PI200'); ACC = 3600.0
rng = lambda k, d: range(int(os.environ.get(k, d).split(':')[0]), int(os.environ.get(k, d).split(':')[1]) + 1)
REF, YEARS = rng('REF', '1490:1499'), rng('YEARS', '1500:1509')


def path(var, y):
    work = [q for q in sorted(glob.glob(f'{R}/{ARM}/run_*/work/atm_remapped_1m_{var}_{y}-{y}.nc'))
            if os.path.basename(os.path.dirname(os.path.dirname(q))).count('.') == 0]   # skip moved-aside legs (run_*.cancelled_*)
    for p in [f'{R}/{ARM}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'] + work:
        if os.path.exists(p):
            return p
    return None


def year(y):
    f = {}
    for v in ('tsr', 'ttr', 'tsrc', 'ttrc', 'tisr', 'ssr', 'str', 'slhf', 'sshf', 'sf'):
        p = path(v, y)
        if p is None:
            return None
        with xr.open_dataset(p, decode_times=False) as d:
            k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
            a = np.squeeze(d[k].values).astype('f8')
            if a.shape[0] != 12:
                return None
            f[v] = a.mean(0) / ACC; lat = np.squeeze(d['lat'].values)
    W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], f['tsr'].shape)
    g = lambda x, k=None: float(np.average(x if k is None else x[k], weights=W if k is None else W[k]))
    so = np.broadcast_to(((lat >= -65) & (lat <= -45))[:, None], W.shape)
    inc = g(f['tisr']); assert 320 < inc < 360, 'incoming solar off'
    return {'net TOA': g(f['tsr'] + f['ttr']), 'SW CRE glob': g(f['tsr'] - f['tsrc']), 'LW CRE glob': g(f['ttr'] - f['ttrc']),
            'SW CRE 45-65S': g(f['tsr'] - f['tsrc'], so), 'planet. albedo': 1 - g(f['tsr']) / inc,
            'sfc net rad': g(f['ssr'] + f['str']),
            'sfc net energy': g(f['ssr'] + f['str'] + f['slhf'] + f['sshf'] - f['sf'] * 1000 * 333550.0),
            'sfc energy 45-65S': g(f['ssr'] + f['str'] + f['slhf'] + f['sshf'] - f['sf'] * 1000 * 333550.0, so)}


ref = [r for r in (year(y) for y in REF) if r]
keys = list(ref[0])
mu = {k: np.mean([r[k] for r in ref]) for k in keys}; sd = {k: np.std([r[k] for r in ref], ddof=1) for k in keys}
print(f'{ARM}: reference {REF[0]}-{REF[-1]} ({len(ref)} yr), mean and year-to-year sd; then each year as a difference from the mean [W/m2]')
print(f'{"":>10}' + ''.join(f'{k:>16}' for k in keys))
print(f'{"ref mean":>10}' + ''.join(f'{mu[k]:16.3f}' for k in keys))
print(f'{"ref sd":>10}' + ''.join(f'{sd[k]:16.3f}' for k in keys))
got = []
for y in YEARS:
    r = year(y)
    if r is None:
        continue
    got.append(r)
    print(f'{y:>10}' + ''.join(f'{r[k] - mu[k]:+16.3f}' for k in keys))
if got:
    n = len(got)
    print(f'{"mean diff":>10}' + ''.join(f'{np.mean([r[k] for r in got]) - mu[k]:+16.3f}' for k in keys) + f'   ({n} yr; 1 sd of an n-yr mean = ref sd/sqrt(n))')
