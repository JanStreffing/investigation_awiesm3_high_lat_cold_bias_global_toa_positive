"""Does the S4 radiative change fall where PI200 was too warm?  Pattern correlation, three scales.

X: annual T2m bias of PI200 over BIAS years against the ERA5 monthly climatology (1990-2014).
Y: S4 minus no-S4 change over YEARS (paired, same years) in net TOA, SW CRE and net surface
energy flux (as scripts/figures/radiation_delta_maps.py).
Correlations are cos-lat weighted: at every gridpoint, on 10x10 degree ocean boxes, and on
5-degree zonal means.  A negative r means the forcing drops most where the model was warmest.
Also prints the zonal means side by side and saves a zonal-mean figure.

Usage:  BIAS=1480:1499 YEARS=1500:1505 python3 scripts/analysis/s4_forcing_vs_bias.py
"""
import os, sys
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'figures'))
os.environ.setdefault('YEARS', '1500:1505')
import radiation_delta_maps as rdm                       # paired S4 fields (runs on import)
R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
B0, B1 = (int(x) for x in os.environ.get('BIAS', '1480:1499').split(':'))

t2 = []
for y in range(B0, B1 + 1):
    with xr.open_dataset(f'{R}/PI200/outdata/oifs/atm_remapped_1m_2t_{y}-{y}.nc', decode_times=False) as d:
        k = [c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]
        t2.append(np.squeeze(d[k].values).astype('f8'))
T = np.stack(t2).mean(0); T = T - 273.15 if np.nanmean(T) > 100 else T
lat, lon = rdm.lat, rdm.lon
with xr.open_dataset('/work/ab0246/a270092/obs/era5/netcdf/T2M.nc') as d:
    c = d['t2m'].groupby('time.month').mean('time')
    la = [x for x in c.dims if 'lat' in x][0]; lo = [x for x in c.dims if 'lon' in x][0]
    O = np.asarray(c.interp({la: ('lat', lat), lo: ('lon', np.sort(lon % 360))}).values, float)
O = O - 273.15 if np.nanmean(O) > 100 else O
o = rdm.o
bias = T[:, :, o].mean(0) - O.mean(0)
ocean = ~rdm.land
W = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], bias.shape)
LAT = np.broadcast_to(lat[:, None], bias.shape); LON = np.broadcast_to(rdm.LON[None, :], bias.shape)


def wcorr(x, y, w):
    k = np.isfinite(x) & np.isfinite(y) & (w > 0); x, y, w = x[k], y[k], w[k]
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.sum(w * (x - mx) * (y - my)) / np.sqrt(np.sum(w * (x - mx) ** 2) * np.sum(w * (y - my) ** 2)))


def boxes(f, mask, dl=10):
    out, ww = [], []
    for a in np.arange(-90, 90, dl):
        for b in np.arange(0, 360, dl):
            k = (LAT >= a) & (LAT < a + dl) & (LON >= b) & (LON < b + dl) & mask
            if k.sum() > 20:
                out.append(np.average(f[k], weights=W[k])); ww.append(np.cos(np.deg2rad(a + dl / 2)))
    return np.array(out), np.array(ww)


ZB = np.arange(-90, 91, 5)
zon = lambda f, mask: np.array([np.average(f[(LAT >= a) & (LAT < a + 5) & mask], weights=W[(LAT >= a) & (LAT < a + 5) & mask])
                                if ((LAT >= a) & (LAT < a + 5) & mask).sum() else np.nan for a in ZB[:-1]])
zc = np.cos(np.deg2rad(ZB[:-1] + 2.5))
print(f'PI200 T2m bias {B0}-{B1} vs ERA5  against  S4 change {rdm.Y0}-{rdm.Y1} (paired):  weighted r')
print(f'{"quantity":<20}{"gridpoint all":>15}{"gridpoint ocean":>17}{"10deg ocean":>13}{"zonal all":>11}{"zonal ocean":>13}')
ZZ = {}
for q in ('net TOA', 'SW CRE', 'net surface energy'):
    D = (rdm.FA[q] - rdm.FB[q])[:, o]
    bx, bw = boxes(bias, ocean); dx, _ = boxes(D, ocean)
    za, zo = zon(D, np.ones_like(ocean)), zon(D, ocean); ZZ[q] = zo
    print(f'{q:<20}{wcorr(bias, D, W):15.2f}{wcorr(np.where(ocean, bias, np.nan), np.where(ocean, D, np.nan), W):17.2f}'
          f'{wcorr(bx, dx, bw):13.2f}{wcorr(zon(bias, np.ones_like(ocean)), za, zc):11.2f}{wcorr(zon(bias, ocean), zo, zc):13.2f}')
zb = zon(bias, ocean)
print('\nzonal means over ocean:   lat   T2m bias [K]   dNetTOA   dSWCRE   dSfcEnergy [W/m2]')
for i, a in enumerate(ZB[:-1]):
    if np.isfinite(zb[i]):
        print(f'{a + 2.5:6.1f} {zb[i]:+10.2f} {ZZ["net TOA"][i]:+10.2f} {ZZ["SW CRE"][i]:+9.2f} {ZZ["net surface energy"][i]:+11.2f}')
fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150); x = ZB[:-1] + 2.5
ax.plot(x, zb, color='k', lw=2, label='T2m bias vs ERA5 [K], ocean'); ax.axhline(0, color='grey', lw=0.6)
ax2 = ax.twinx()
for q, c_ in (('net TOA', '#2a78d6'), ('SW CRE', '#eb6834'), ('net surface energy', '#1baf7a')):
    ax2.plot(x, ZZ[q], color=c_, lw=1.4, label=f'S4 change: {q} [W m$^{{-2}}$]')
ax.set_xlabel('latitude'); ax.set_ylabel('K'); ax2.set_ylabel('W m$^{-2}$')
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, fontsize=7.5, frameon=False, loc='upper left')
ax.set_title(f'Ocean zonal means: PI200 bias {B0}-{B1} and S4 change {rdm.Y0}-{rdm.Y1}', loc='left', fontsize=10)
fig.tight_layout(); out = os.path.join(REPO, 'plots', f's4_forcing_vs_bias_zonal_{rdm.Y0}-{rdm.Y1}.png'); fig.savefig(out); print('saved', out)
