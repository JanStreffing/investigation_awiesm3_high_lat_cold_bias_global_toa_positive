"""Where would a GEOMETRIC eddy-energy closure build eddy energy in the current ocean state?

GEOMETRIC (Marshall et al. 2012; Mak et al. 2018) sets the GM coefficient from a depth-integrated
eddy energy E:
    K = alpha * E / int( M2 / N ) dz ,   M2 = |grad_h b| ,  N2 = db/dz
and E obeys  dE/dt = int( K M4 / N2 ) dz - lambda E  (+ advection, diffusion).
Substituting K, the local source is  E * sigma  with
    sigma = alpha * int( M4 / N2 ) dz / int( M2 / N ) dz      [1/s]
so eddy energy grows where sigma > lambda and decays where sigma < lambda. E itself, and with it K,
cannot be had offline: it equilibrates through its feedback on the density slopes. What can be
mapped from existing output is sigma, the growth rate the scheme would see at switch-on.

Here: annual-mean and September/March T, S of one year of PICAL_crunveg_tke_albsn082; b from the
in-situ density at the level's pressure (gradients along the level), N2 from the same density
referenced locally; M2 on elements from the mesh's gradient operator; slope |S| = M2/N2 capped at
0.01 and N2 floored at 1e-8 s-2 (as any implementation would have to); integrals from 100 m to the
bottom (below the mixed layer: GM is tapered there). alpha = 0.06, 1/lambda = 100 days.
Regional means are area-weighted over elements deeper than 1000 m.

Usage: geometric_growth_rate.py [year]   (reval environment, ~25 GB; writes data/clim/geometric_growth_rate.txt)
"""
import os, sys
for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import numpy as np
import xarray as xr
import gsw

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = '/albedo/work/projects/p_awiesm3_cmip7/jstreffi'
D = f'{P}/runtime/awiesm3-v3.4/PICAL_crunveg_tke_albsn082/outdata/fesom'
YEAR = int(sys.argv[1]) if len(sys.argv) > 1 else 2205
ALPHA, LAM = 0.06, 1.0 / (100 * 86400.0)
SMAX, N2MIN, ZTOP, G, RHO0 = 0.01, 1e-8, 100.0, 9.81, 1027.0

md = xr.open_dataset(f'{P}/reval_obs/mesh/core3/fesom.mesh.diag.nc')
tri = md['face_nodes'].values.T - 1                       # (elem, 3)
gx, gy = md['gradient_sca_x'].values.T, md['gradient_sca_y'].values.T
earea = md['elem_area'].values
zi = np.abs(md['nz'].values); zm = 0.5 * (zi[1:] + zi[:-1]); dz = np.diff(zi)
lon, lat = md['lon'].values, md['lat'].values
elon = np.rad2deg(np.arctan2(np.sin(np.deg2rad(lon[tri])).mean(1), np.cos(np.deg2rad(lon[tri])).mean(1)))
elat = lat[tri].mean(1)
edepth = np.abs(md['zbar_e_bottom'].values)
ecav = md['ulevels'].values > 1

REG = [('Labrador 56-62N 60-50W', (elat >= 56) & (elat < 62) & (elon >= -60) & (elon < -50)),
       ('SPNA 45-66N', (elat >= 45) & (elat < 66) & (elon >= -65) & (elon < 0)),
       ('Nordic Seas 66-80N', (elat >= 66) & (elat < 80) & (elon >= -20) & (elon < 20)),
       ('N Pacific 45-66N', (elat >= 45) & (elat < 66) & ((elon >= 140) | (elon < -120))),
       ('tropics 10S-10N', (np.abs(elat) < 10)),
       ('subtropics 20-40', (np.abs(elat) >= 20) & (np.abs(elat) < 40)),
       ('SO 40-55S', (elat >= -55) & (elat < -40)),
       ('SO 55-65S', (elat >= -65) & (elat < -55)),
       ('Antarctic <65S', (elat < -65)),
       ('Weddell 55-75S 60W-20E', (elat >= -75) & (elat < -55) & (elon >= -60) & (elon < 20))]


def sigma(T, S):
    """T, S (nod, nz) -> growth rate sigma [1/s] and int(M2/N) per element."""
    S = np.where(S > 1, S, np.nan); T = np.where(np.isfinite(S), T, np.nan)
    p = zm[None, :] * np.ones_like(T)
    rho = gsw.rho(S, T, p)                                   # in situ at level pressure (SA, CT approximated by S, theta)
    b = -G * rho / RHO0
    be = b[tri]                                              # (elem, 3, nz)
    bx = np.einsum('ej,ejk->ek', gx, be); by = np.einsum('ej,ejk->ek', gy, be)
    M2 = np.sqrt(bx ** 2 + by ** 2)                          # NaN where any node is dry
    # N2 at level interfaces from locally referenced density, then to mid-levels
    pm = 0.5 * (zm[1:] + zm[:-1])[None, :]
    r_up = gsw.rho(S[:, :-1], T[:, :-1], pm); r_dn = gsw.rho(S[:, 1:], T[:, 1:], pm)
    n2i = G / RHO0 * (r_dn - r_up) / np.diff(zm)[None, :]    # (nod, nz-1), positive when stable
    n2 = np.full_like(T, np.nan); n2[:, 1:-1] = 0.5 * (n2i[:, 1:] + n2i[:, :-1]); n2[:, 0] = n2i[:, 0]; n2[:, -1] = n2i[:, -1]
    N2 = np.nanmean(n2[tri], 1)
    N2 = np.where(np.isfinite(M2), np.maximum(N2, N2MIN), np.nan)
    slope = np.minimum(M2 / N2, SMAX)
    N = np.sqrt(N2)
    w = np.where((zm[None, :] >= ZTOP) & np.isfinite(slope), dz[None, :], 0.0)
    num = np.nansum(N2 * slope ** 2 * w, 1)                  # int M4/N2 dz with the capped slope
    den = np.nansum(N * slope * w, 1)                        # int M2/N dz
    with np.errstate(invalid='ignore', divide='ignore'):
        return ALPHA * num / den, den, np.nansum(slope * w, 1) / np.maximum(w.sum(1), 1e-9)


ds_t = xr.open_dataset(f'{D}/temp.fesom.{YEAR}.nc')['temp']; ds_s = xr.open_dataset(f'{D}/salt.fesom.{YEAR}.nc')['salt']
FIELDS = {'annual': (ds_t.mean('time').values, ds_s.mean('time').values),
          'March': (ds_t.isel(time=2).values, ds_s.isel(time=2).values),
          'September': (ds_t.isel(time=8).values, ds_s.isel(time=8).values)}
ok = (edepth > 1000) & ~ecav
out = [f'GEOMETRIC growth rate of eddy energy, sigma = alpha int(M4/N2)dz / int(M2/N)dz, from PICAL_crunveg_tke_albsn082 {YEAR}.',
       f'alpha = {ALPHA}, slope capped at {SMAX}, N2 floored at {N2MIN:g}, integrals below {ZTOP:.0f} m, elements deeper than 1000 m, no cavities.',
       'sigma in 1/(100 days): > 1 means eddy energy grows against a 100-day dissipation, < 1 means it decays.',
       'frac>1 = area fraction with sigma > lambda; slope = depth-mean isopycnal slope below 100 m (1e-4).', '',
       f"{'region':26s}" + ''.join(f"{k + ' sigma':>16s}{'frac>1':>8s}" for k in FIELDS) + f"{'slope ann':>11s}"]
res = {k: sigma(*v) for k, v in FIELDS.items()}
for name, m in REG:
    mm = m & ok
    line = f'{name:26s}'
    for k in FIELDS:
        s = res[k][0] / LAM
        g = mm & np.isfinite(s)
        line += f"{np.sum(s[g] * earea[g]) / np.sum(earea[g]):16.2f}{np.sum(earea[g & (s > 1)]) / np.sum(earea[g]):8.2f}"
    sl = res['annual'][2]; g = mm & np.isfinite(sl)
    out.append(line + f"{1e4 * np.sum(sl[g] * earea[g]) / np.sum(earea[g]):11.2f}")
txt = '\n'.join(out)
print(txt)
open(f'{REPO}/data/clim/geometric_growth_rate.txt', 'w').write(txt + '\n')
np.save(f'{REPO}/data/clim/geometric_growth_rate_annual_{YEAR}.npy', (res['annual'][0] / LAM).astype('f4'))
