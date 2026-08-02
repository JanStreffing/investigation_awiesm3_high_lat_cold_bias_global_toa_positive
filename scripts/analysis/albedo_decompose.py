"""Decompose the Siberian surface albedo into the three terms that can be wrong.

Round 13 killed the snow-COVER hypothesis: raising the critical snow depth moved
the albedo response into September-November and did nothing in June. So the June
albedo bias is set by something else. This script asks what, using output we
already have rather than another 48-year run.

HTESSEL builds the grid-box albedo out of three things:

    alpha_box  ~  f_snow * alpha_snow  +  (1 - f_snow) * alpha_snowfree

  f_snow          snow cover fraction, ZCVS = min(1, d_cm * RQSNCR), recomputable
                  here from sd (snow water equivalent) and rsn (snow density)
  alpha_snow      the prognostic snow albedo -- output directly as `asn`
  alpha_snowfree  the vegetation/soil blend -- not output

and the model's actual all-sky albedo is output directly as `fal`.

IMPORTANT -- the naive two-term split does NOT close, and the residual is the
finding, not a bug. Backing alpha_snowfree out of the identity gives NEGATIVE
values over the boreal box in May, June and October, because HTESSEL does not
apply `asn` to the whole snow-covered fraction: tile 7 is snow SHELTERED under
high vegetation and carries a canopy-shaded albedo far below `asn`, while tile 5
is exposed snow. So `fal` sits well below f_snow*asn wherever there is forest.
That the residual is large and negative in exactly the melt months is direct
evidence the canopy-masking term is already doing most of the work there.

This script therefore reports only measured quantities -- fal, asn, f_snow,
snow depth -- and compares them against ERA5, rather than inventing a
snow-free albedo the output cannot constrain.

Three candidate causes, and what would show them:
  1. too much snow            -> f_snow too high in June
  2. snow too bright          -> alpha_snow (asn) too high
  3. snow-free surface bright -> alpha_snowfree too high (this is the vegetation
                                 albedo table RVVEGALB, which the campaign has
                                 never touched)

Reference: ERA5 has all three fields (fal, asn, sd, rsn) on the DKRZ pool and
assimilates snow observations, so its snow EXTENT is observationally constrained
even though its albedo scheme is a relative of ours. Where ERA5 is unavailable
the script still prints the model's own decomposition, which localises the term
carrying the bias even without an absolute reference.

Period note: the model is 1870s and ERA5 is 1990-2014, so `amip_presentday` is
used for the ERA5 comparison whenever it is present (the campaign's standard
period-clean pairing), and the PI control is shown alongside for the offset.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, Y0, Y1

BOX = ((55, 75), (60, 180))
ACC = 3600.0
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
# Inverse critical snow depth, per run. As-released LESN09=T is 1/10
# (sussoil_mod.F90:157); the H-series rebuilt with 1/30. Using the wrong one
# silently mis-reconstructs f_snow, so it is keyed by run rather than assumed.
RQSNCR_DEFAULT = 1.0 / 10.0
RQSNCR_BY_RUN = {'amip_H1_snowcr30': 1.0 / 30.0, 'amip_H2_G1_snowcr30': 1.0 / 30.0}

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def box_mean(a2d, lat, lon, weights=None):
    ys = (lat >= BOX[0][0]) & (lat <= BOX[0][1])
    l180 = ((lon + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub, L = a2d[ii], lsm[ii]
    m = (L > 0.5) & np.isfinite(sub)
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape).copy()
    if weights is not None:
        w = w * weights[ii]
    return np.average(sub[m], weights=w[m]) if m.any() else np.nan


def load(run, var, years):
    """Monthly climatology [12, ny, nx] of `var`, plus lat/lon. None if absent."""
    acc, n, lat, lon = None, 0, None, None
    for y in years:
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        ds = xr.open_dataset(f)
        a = ds[var].values
        lat, lon = ds[var].lat.values, ds[var].lon.values
        ds.close()
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lon) if n else (None, None, None)


def decompose(run, years, tag):
    print(f'\n{"="*78}\n{tag}   ({run}, {years[0]}-{years[-1]})\n{"="*78}')
    fal, lat, lon = load(run, 'fal', years)
    asn, _, _ = load(run, 'asn', years)
    sd, _, _ = load(run, 'sd', years)
    rsn, _, _ = load(run, 'rsn', years)
    ssr, _, _ = load(run, 'ssr', years)
    ssrd, _, _ = load(run, 'ssrd', years)
    if fal is None or asn is None or sd is None or rsn is None:
        print('  missing fields, skipped'); return None

    rq = RQSNCR_BY_RUN.get(run, RQSNCR_DEFAULT)
    rows = {}
    for m in range(12):
        # snow cover fraction exactly as surfbc_ctl_mod.F90:317-322 computes it:
        # sd is snow water equivalent [m], rsn is density [kg/m3] -> depth_cm
        depth_cm = 100.0 * sd[m] * 1000.0 / np.maximum(rsn[m], 1.0)
        fsnow = np.clip(depth_cm * rq, 0.0, 1.0)

        # Residual of the naive two-term model, as a DIAGNOSTIC of how much
        # canopy masking (tile 7) is doing. Large negative = fal is far below
        # f_snow*asn = the sheltered-snow tile is heavily shading the snow.
        resid = box_mean(fal[m] - fsnow * asn[m] - (1.0 - fsnow) * 0.15, lat, lon)

        a_box = box_mean(fal[m], lat, lon)
        a_snow = box_mean(asn[m], lat, lon, weights=fsnow) if fsnow.max() > 0 else np.nan
        f_sn = box_mean(fsnow, lat, lon)

        # all-sky albedo actually seen by the radiation, as a closure check
        a_rad = (1.0 - box_mean(ssr[m] / ACC, lat, lon) /
                 max(box_mean(ssrd[m] / ACC, lat, lon), 1e-6)) if ssr is not None else np.nan
        rows[m] = (a_box, f_sn, a_snow, resid, a_rad, box_mean(depth_cm, lat, lon))

    print(f'  {"":5s} {"alb(fal)":>9s} {"f_snow":>8s} {"alb_snow":>9s} '
          f'{"canopy":>9s} {"alb(SW)":>8s} {"depth_cm":>9s}')
    for m in range(12):
        a_box, f_sn, a_snow, a_free, a_rad, dcm = rows[m]
        star = ' <<<' if m in (5, 6) else ''
        print(f'  {MON[m]:5s} {a_box:9.4f} {f_sn:8.3f} {a_snow:9.4f} '
              f'{a_free:9.4f} {a_rad:8.4f} {dcm:9.2f}{star}')
    print('  alb(fal) is the model surface albedo; alb(SW)=1-ssr/ssrd is what the')
    print('  radiation actually saw (they differ by the cloud/zenith weighting).')
    print('  canopy = fal - f_snow*asn - (1-f_snow)*0.15, the shortfall of the naive')
    print('  two-term model. Strongly negative = tile 7 canopy shading dominates.')
    return rows


years = list(range(Y0, Y1 + 1))
pi = decompose('amip_pi_base', years, 'PI control')

pd_years = list(range(1990, 2015))
if os.path.exists(f'{RT}/amip_presentday/outdata/oifs/atm_remapped_1m_fal_1m_1990-1990.nc'):
    decompose('amip_presentday', pd_years, 'Present day (period-clean vs ERA5)')

# H1 is the falsification control: it changed f_snow and nothing else, so it is a
# direct check that the f_snow reconstruction above tracks what the model did.
h1 = decompose('amip_H1_snowcr30', years, 'H1 snowcr30 (RQSNCR=1/30) -- reconstruction check')
if pi and h1:
    print('\n  Check: H1 changed only RQSNCR, so f_snow should fall and alb_snow should not.')
    for m in (5, 9):
        print(f'    {MON[m]}: f_snow {pi[m][1]:.3f} -> {h1[m][1]:.3f}, '
              f'alb_snow {pi[m][2]:.4f} -> {h1[m][2]:.4f}, '
              f'alb_box {pi[m][0]:.4f} -> {h1[m][0]:.4f}')


# ---------------------------------------------------------------------------
# ERA5 reference. Snow EXTENT in ERA5 is observationally constrained (IMS snow
# cover is assimilated), so f_snow and snow depth are a fair test. The albedo
# scheme is a relative of HTESSEL, so `asn` agreement is NOT independent
# evidence -- a matching snow albedo means "same scheme", not "correct".
# ---------------------------------------------------------------------------
E5 = '/work/ab0246/a270092/obs/era5/snow'
E5V = {'fal': ('era5_243_clim_1990-2014.nc', 'var243'),
       'asn': ('era5_032_clim_1990-2014.nc', 'var32'),
       'sd':  ('era5_141_clim_1990-2014.nc', 'var141'),
       'rsn': ('era5_033_clim_1990-2014.nc', 'var33')}


def era5_box():
    out = {}
    for k, (fn, vn) in E5V.items():
        p = os.path.join(E5, fn)
        if not os.path.exists(p):
            return None
        ds = xr.open_dataset(p)
        a = ds[vn].values
        la, lo = ds['lat'].values, ds['lon'].values
        ds.close()
        out[k] = (a, la, lo)
    return out


# ERA5's own land-sea mask, on ERA5's grid. Without it the box is diluted by the
# Arctic coast and the Sea of Okhotsk -- water has albedo ~0.07 and never any
# snow, which drags ERA5's fal and f_snow down and would manufacture exactly the
# model-too-bright / model-too-snowy signal this script is looking for.
_E5LSM = None
_p = os.path.join(E5, 'era5_lsm.nc')
if os.path.exists(_p):
    _d = xr.open_dataset(_p)
    _E5LSM = np.squeeze(_d['var172'].values)
    _d.close()


def e5_mean(a2d, la, lo, weights=None):
    """Area-weighted, LAND-MASKED Siberian-box mean on ERA5's own grid."""
    ys = (la >= BOX[0][0]) & (la <= BOX[0][1])
    l180 = ((lo + 180) % 360) - 180
    xs = (l180 >= BOX[1][0]) & (l180 <= BOX[1][1])
    ii = np.ix_(np.where(ys)[0], np.where(xs)[0])
    sub = a2d[ii]
    w = np.broadcast_to(np.cos(np.deg2rad(la[ys]))[:, None], sub.shape).copy()
    if weights is not None:
        w = w * weights[ii]
    m = np.isfinite(sub)
    if _E5LSM is not None:
        m = m & (_E5LSM[ii] > 0.5)
    return np.average(sub[m], weights=w[m]) if m.any() else np.nan


e5 = era5_box()
if e5 is None:
    print('\nERA5 files absent -- run albedo_decompose_prep.sh first.')
else:
    fal5, la, lo = e5['fal']
    asn5 = e5['asn'][0]; sd5 = e5['sd'][0]; rsn5 = e5['rsn'][0]
    print(f'\n{"="*78}\nERA5 1990-2014 vs model amip_presentday -- Siberian box'
          f'\n{"="*78}')
    print('  ERA5 is land-masked with its own lsm>0.5, matching the model box.')
    print('  f_snow for BOTH columns is min(1, depth_cm/10) -- the same HTESSEL')
    print('  formula applied to each dataset\'s snow mass, so a difference means')
    print('  different SNOW AMOUNT, not a different cover formula.')
    print('  Snow EXTENT in ERA5 is observation-constrained (IMS assimilated);')
    print('  its snow ALBEDO shares HTESSEL heritage and is NOT independent.')
    print(f'\n  {"":5s} {"fal_ERA5":>9s} {"fal_mod":>9s} {"d(fal)":>8s} '
          f'{"asn_ERA5":>9s} {"asn_mod":>8s} {"fsn_ERA5":>9s} {"fsn_mod":>8s}')
    md = decompose('amip_presentday', pd_years, 'model present day (for the ERA5 row)') \
        if os.path.exists(f'{RT}/amip_presentday/outdata/oifs/atm_remapped_1m_fal_1m_1990-1990.nc') \
        else None
    for m in range(12):
        d_cm5 = 100.0 * sd5[m] * 1000.0 / np.maximum(rsn5[m], 1.0)
        f5 = np.clip(d_cm5 * RQSNCR_DEFAULT, 0.0, 1.0)
        fa5 = e5_mean(fal5[m], la, lo)
        as5 = e5_mean(asn5[m], la, lo, weights=f5)
        fs5 = e5_mean(f5, la, lo)
        if md:
            fam, asm_, fsm = md[m][0], md[m][2], md[m][1]
            print(f'  {MON[m]:5s} {fa5:9.4f} {fam:9.4f} {fam-fa5:+8.4f} '
                  f'{as5:9.4f} {asm_:8.4f} {fs5:9.3f} {fsm:8.3f}')
        else:
            print(f'  {MON[m]:5s} {fa5:9.4f} {"--":>9s} {"--":>8s} '
                  f'{as5:9.4f} {"--":>8s} {fs5:9.3f} {"--":>8s}')
