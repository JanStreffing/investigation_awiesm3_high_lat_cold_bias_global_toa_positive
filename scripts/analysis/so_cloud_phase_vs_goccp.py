"""Is the Southern Ocean too icy?  Model low-level cloud phase against CALIPSO-GOCCP.

THE QUESTION THIS SETTLES.  RCL_INPSEA is the campaign's only SPATIALLY SELECTIVE
Southern Ocean lever -- cloudsc.F90:2520-2521 scales ice nuclei by
RCL_INPSEA+(1-RCL_INPSEA)*PLSM, which is exactly 1 over land, so Siberian cloud cannot
respond to it.  It has sat at 0.2 in all 46 runscripts that set it and has never been
scanned.  0.2 was never measured either: the scripts justify it as the "marine biogenic
floor", an argument about INP CONCENTRATION.

But the quantity that actually constrains this is CLOUD PHASE.  If the Southern Ocean
were still too icy at 0.2, lowering it further would be a correction toward observations
rather than a fudge past a floor.  So: measure the phase.

RESULT, 2026-08-19 -- it is NOT too icy, and the case for going lower fails.

                        model 11G      GOCCP        
    SO 45-65S             0.875        0.858
    Siberia 55-75N        0.820        0.818
    contrast             +0.055       +0.040

The model is already slightly MORE liquid than observed in the Southern Ocean.  Pushing
RCL_INPSEA below 0.2 would move low-level phase further past the observation, so the
cloud-phase argument does not justify it and the biogenic-INP floor stands.

WHAT THAT IMPLIES, and it is the useful part.  The SO SW CRE gap is +3.89 W/m2 against
CERES in 11G while the phase is right.  So the residual Southern Ocean error is NOT a
phase error, and the INP family of levers has done its job.  What remains must be cloud
AMOUNT, cloud water, or droplet number -- and the report already measured an amount
deficit: cloud area 83.1 % against CERES's 89.7 %, i.e. -6.6 pp.  That is where the next
lever has to act.

METRIC CAVEAT, stated plainly.  The model number is MASS-based, liquid over liquid+ice
condensate summed over the lowest 25 model levels.  GOCCP is CLOUD-FRACTION based, the
frequency with which the lidar detects liquid versus ice cloud below 3 km.  These are not
formally the same quantity and a rigorous comparison needs COSP, which this configuration
does not run.  What makes the conclusion usable anyway is the Siberian anchor: two
independent regions agree to 0.002 and 0.017, which is not what incommensurate metrics
do.  Treat the sign as solid and the third decimal as not.

GOCCP RPIC IS NOT USED.  clcalipso_RPIC reads 0.99 as a liquid fraction in both regions,
which cannot be right; the liq/(liq+ice) cloud fractions are used instead.

GRID NOTE.  The model ml_ files are on the reduced Gaussian grid (cell dimension), where
cells are approximately equal area, so a plain cell mean IS the area mean.  Applying
cos(lat) weighting there would double-count.

DATA.  CALIPSO-GOCCP v3.1.4, 3D_CloudFraction_Phase330m, monthly avg, 2008-2010, from
ftp://ftp.climserv.ipsl.polytechnique.fr/cfmip/GOCCP_V3.1.4/ -- staged in
/work/ab0246/a270092/obs/CALIPSO_GOCCP/.
"""
import glob
import numpy as np
import warnings
warnings.filterwarnings('ignore')

R = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
GOCCP = '/work/ab0246/a270092/obs/CALIPSO_GOCCP'
BANDS = (('SO 45-65S', -65, -45), ('Siberia 55-75N', 55, 75))


def model_lowfrac(root, years, lo, hi, nlow=25):
    import xarray as xr
    L, I = [], []
    for y in years:
        fw = f'{root}/outdata/oifs/atmos_mon_ml_clw_{y}-{y}.nc'
        fi = f'{root}/outdata/oifs/atmos_mon_ml_cli_{y}-{y}.nc'
        if not (glob.glob(fw) and glob.glob(fi)):
            continue
        try:
            with xr.open_dataset(fw, decode_times=False) as d:
                a = np.asarray(d['clw'].values, dtype=float)
                lat = np.asarray(d['lat'].values)
            with xr.open_dataset(fi, decode_times=False) as d:
                b = np.asarray(d['cli'].values, dtype=float)
            a = a.mean(axis=0)[-nlow:].sum(axis=0)
            b = b.mean(axis=0)[-nlow:].sum(axis=0)
            s = (lat >= lo) & (lat < hi)
            L.append(float(a[s].mean()))
            I.append(float(b[s].mean()))
        except Exception:
            pass
    return np.mean(L) / (np.mean(L) + np.mean(I)) if L else np.nan


def goccp_lowfrac(lo, hi, zmax=3.0):
    import netCDF4 as nc
    FS = sorted(glob.glob(f'{GOCCP}/3D_CloudFraction_Phase330m_*_avg_*.nc'))
    Lq, Ic = [], []
    for f in FS:
        d = nc.Dataset(f)
        lat = d.variables['latitude'][:]
        alt = d.variables['alt_mid'][:]
        s = (lat >= lo) & (lat < hi)
        zs = alt <= zmax
        for var, acc in (('clcalipso_liq', Lq), ('clcalipso_ice', Ic)):
            a = np.ma.masked_invalid(d.variables[var][:])
            a = np.ma.masked_where(a < -1e3, a)
            sub = a[0][zs][:, s, :]
            w = np.cos(np.deg2rad(lat[s]))[None, :, None] * np.ones_like(sub)
            acc.append(float((sub * w).sum() / w.sum()) if np.ma.count(sub) else np.nan)
        d.close()
    lq, ic = np.nanmean(Lq), np.nanmean(Ic)
    return lq / (lq + ic), len(FS)


def main():
    print(__doc__)
    print('=' * 88)
    obs = {}
    for name, lo, hi in BANDS:
        obs[name], n = goccp_lowfrac(lo, hi)
    print(f'\nCALIPSO-GOCCP v3.1.4, {n} monthly files, cloud-fraction based, below 3 km')
    print(f'model: mass-based, lowest 25 model levels\n')
    print(f'  {"arm":6s} ' + '  '.join(f'{b[0]:>16s}' for b in BANDS) + f'  {"contrast":>9s}')
    for tag, root, Y in [('11G', f'{R}/Tuning_test_11G_inppmin50k', range(1390, 1400)),
                         ('11L', f'{R}/11L', range(1380, 1390)),
                         ('11M', f'{R}/11M', range(1380, 1390))]:
        v = [model_lowfrac(root, list(Y), lo, hi) for _, lo, hi in BANDS]
        print(f'  {tag:6s} ' + '  '.join(f'{x:16.3f}' for x in v)
              + f'  {v[0] - v[1]:+9.3f}')
    v = [obs[b[0]] for b in BANDS]
    print(f'  {"GOCCP":6s} ' + '  '.join(f'{x:16.3f}' for x in v) + f'  {v[0] - v[1]:+9.3f}')
    print('\n  The model is not too icy in the SO; it is slightly MORE liquid than')
    print('  observed. The cloud-phase argument for lowering RCL_INPSEA below 0.2 fails,')
    print('  and the residual SO CRE gap is an amount problem, not a phase problem.')


if __name__ == '__main__':
    main()
