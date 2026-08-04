"""The land albedo excess, split into its snow half and its snow-free half.

WHY THIS EXISTS. `tropo_bias_section.py` traced the global tropospheric cold bias to a
CLEAR-SKY shortwave deficit (cloud is the wrong sign), and a first pass attributed part of
it to land surface albedo being +0.0154 too high. That pass classified land by dominant
VEGETATION type (tvh/cvh/tvl/cvl) -- and in HTESSEL snow is not a vegetation type, it is a
separate tile (5 = snow on low vegetation, 7 = snow under high vegetation) that OVERLIES
them. So every vegetation type in that table was part-time snow, and the ranking it
produced was mostly a map of who is snow-covered longest:

    Evergreen Needleleaf   53.9 % of the year snow-covered
    Tundra                 63.2 %
    bare soil              49.2 %   (Sahara and high-Arctic bare ground in one bucket)

An annual-mean T2m > 5 C filter does NOT fix this -- a cell averaging 6 C still has months
of snow. This script masks snow PER CELL PER MONTH (sd < 1 mm w.e.) and restricts to sunlit
months (SW down > 20 W/m2), where albedo is defined at all.

THE RESULT: the +0.0154 is almost exactly half snow.

    ALL LAND, sunlit months      +0.0154
      of which snow              +0.0074
      of which snow-free surface +0.0080

and the per-type ranking inverts once snow is removed (all-months -> snow-free):

    Evgr Needleleaf   +0.0250 -> +0.0032     entirely snow
    Bogs/Marshes      +0.0143 -> +0.0039     entirely snow
    Evgr Shrubs       +0.0150 -> +0.0058     entirely snow
    Evgr Broadleaf    +0.0017 -> -0.0003     was already right
    ...
    Crops             +0.0215 -> +0.0204     SURVIVES
    Semidesert        +0.0342 -> +0.0169     survives, halved
    bare soil         +0.0240 -> +0.0145     survives
    Irrig Crops       +0.0253 -> +0.0138     survives
    Tundra            +0.0064 -> +0.0105     gets WORSE -- snow was masking it

TWO HALVES, TWO DIFFERENT FIXES.

  * The SNOW half is what the round-15 snow-cover-fraction scheme acts on, and the
    model-internal comparison below shows I3 (mode 2, SDOR-scaled) removes ~70 % of it
    while I1 (mode 1) removes ~32 %. NOTE that I3 was calibrated to match I1 *in the
    Siberian box* so that any remaining difference would be spatial structure -- globally
    it is twice as strong. That is new information about the scheme.

  * The SNOW-FREE half sits on crops, irrigated crops, semidesert and bare soil, i.e.
    the sparse and cultivated surfaces where the soil background shows through, plus a
    genuine in-season tundra error. EVERY FOREST TYPE IS WITHIN +0.006 and tropical
    broadleaf is at -0.0003, so RVVEGALB's high-vegetation entries are sound -- that part
    of report sub:albreg holds up. Nothing in 41 runs has touched this half.

CAVEATS. CERES surface fluxes are a derived radiative-transfer product (Kato et al. 2018),
not direct observation, with several W/m2 of regional uncertainty; the TOA anchor is the
robust number and this is attribution. Albedos are assumption-free -- each dataset divided
by its OWN downward flux -- because an earlier pass that borrowed the CERES denominator for
the model inflated the land excess from +0.015 to +0.030. The G4-vs-I comparison is PI
runs against each other while the +0.0074 target is diagnosed from amip_presentday, so the
fraction-removed is indicative, not exact.
"""
import numpy as np, xarray as xr, os, glob, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, OBS

ACC = 3600.0
PD = list(range(1990, 2015))
SNOW_WE = 1e-3      # m water equivalent; below this a cell is effectively snow-free
SUN_MIN = 20.0      # W/m2 downward; below this albedo is meaningless

NAMES = {0: 'bare soil', 1: 'Crops', 2: 'Short Grass', 3: 'Evgr Needleleaf',
         4: 'Dec Needleleaf', 5: 'Dec Broadleaf', 6: 'Evgr Broadleaf', 7: 'Tall Grass',
         8: 'Desert', 9: 'Tundra', 10: 'Irrig Crops', 11: 'Semidesert', 12: 'Ice Caps',
         13: 'Bogs/Marshes', 16: 'Evgr Shrubs', 17: 'Dec Shrubs', 18: 'Mixed Forest',
         19: 'Interrupted Forest'}

lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values
land = lsm > 0.5


def clim(run, v, years=None):
    """Monthly climatology (12, ny, nx) for a run, or (None, None, None) if absent."""
    if years is None:
        fs = sorted(glob.glob(f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{v}_1m_*.nc'))
    else:
        fs = [f'{RT}/{run}/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc' for y in years]
        fs = [f for f in fs if os.path.exists(f)]
    acc, n, lat, lon = None, 0, None, None
    for f in fs:
        d = xr.open_dataset(f)
        a = d[v].values
        lat, lon = d[v].lat.values, d[v].lon.values
        d.close()
        acc = a if acc is None else acc + a
        n += 1
    return (acc / n, lat, lon) if n else (None, None, None)


def dominant_type():
    def g(v):
        d = xr.open_dataset(f'{RT}/amip_pi_base/outdata/oifs/atm_remapped_1d_{v}_1d_1900-1900.nc')
        a = d[v].values[0]
        d.close()
        return a
    tvh, cvh, tvl, cvl = g('tvh'), g('cvh'), g('tvl'), g('cvl')
    bare = np.clip(1.0 - cvh - cvl, 0.0, 1.0)
    return np.where(cvh >= np.maximum(cvl, bare), np.round(tvh),
                    np.where(cvl >= bare, np.round(tvl), 0)).astype(int)


# --- model present-day, and CERES on the model grid -------------------------
ssr, lat, lon = clim('amip_presentday', 'ssr', PD)
if ssr is None:
    raise SystemExit('amip_presentday ssr missing')
ssrd, _, _ = clim('amip_presentday', 'ssrd', PD)
sd, _, _ = clim('amip_presentday', 'sd', PD)

ds = xr.open_dataset(OBS)
ds = xr.concat([ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360.0), ds,
                ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360.0)], dim='lon')
tl = np.clip(lat, -89.5, 89.5)
tlo = np.where(lon < 0, lon + 360, lon)


def ip(v):
    return ds[v].interp(lat=xr.DataArray(tl, dims='y'),
                        lon=xr.DataArray(tlo, dims='x')).values


Cup, Cdn = ip('sfc_sw_up_all_clim'), ip('sfc_sw_down_all_clim')
Ms, Md = ssr / ACC, ssrd / ACC
dom = dominant_type()
w = np.cos(np.deg2rad(lat))[:, None] * np.ones_like(Ms[0])
W3 = np.broadcast_to(w, Ms.shape)
L3 = np.broadcast_to(land, Ms.shape)
sunlit = (Md > SUN_MIN) & L3
snowfree = sd < SNOW_WE


def alb(m3):
    """(model albedo, CERES albedo, mean downward flux) over a boolean 3-D mask."""
    dn = np.average(Md[m3], weights=W3[m3])
    am = 1.0 - np.average(Ms[m3], weights=W3[m3]) / dn
    ac = np.average(Cup[m3], weights=W3[m3]) / np.average(Cdn[m3], weights=W3[m3])
    return am, ac, dn


print('Model amip_presentday (1990-2014) vs CERES EBAF, land only, all-sky.')
print('Snow is a TILE, not a vegetation type -- it overlies every type below.\n')
print(f'  {"surface type":26s} {"snow-cov":>9s} {"model":>8s} {"CERES":>8s} '
      f'{"snow-free":>10s} {"W/m2":>7s} {"all-mo":>9s}')
print('  ' + '-' * 82)
rows = []
for t in sorted(NAMES):
    m = land & (dom == t)
    if not m.any() or 100 * w[m].sum() / w[land].sum() < 0.4:
        continue
    m3 = sunlit & np.broadcast_to(dom == t, Ms.shape)
    sf = m3 & snowfree
    if sf.sum() < 50:
        continue
    am, ac, dn = alb(sf)
    aa, ab, _ = alb(m3)
    cov = np.average((~snowfree[:, m]).mean(0), weights=w[m])
    rows.append((t, NAMES[t], cov, am, ac, am - ac, (am - ac) * dn, aa - ab))
for t, nm, cov, am, ac, d, lost, old in sorted(rows, key=lambda r: r[5]):
    print(f'  {t:2d} {NAMES[t]:23s} {100*cov:8.1f}% {am:8.4f} {ac:8.4f} '
          f'{d:+10.4f} {-lost:7.2f} {old:+9.4f}')
print('  ' + '-' * 82)
am, ac, _ = alb(sunlit & snowfree)
aa, ab, _ = alb(sunlit)
print(f'  {"ALL LAND snow-free":26s}           {am:8.4f} {ac:8.4f} {am-ac:+10.4f}')
print(f'  {"ALL LAND all sunlit months":26s}   {aa:8.4f} {ab:8.4f} {aa-ab:+10.4f}')
print(f'\n  => SNOW contributes {(aa-ab)-(am-ac):+.4f} of the {aa-ab:+.4f} land albedo excess,')
print(f'     the snow-free surface the other {am-ac:+.4f}.')

# --- does the round-15 snow scheme eat into the snow half? ------------------
print('\n\nLAND albedo over SNOW-COVERED sunlit months: does the new SCF scheme darken it?')
print('(model-internal; PI runs, so the absolute level is not comparable to CERES)\n')
print(f'  {"run":24s} {"nyr":>4s} {"snow-mo alb":>12s} {"vs G4":>9s} '
      f'{"all-land":>10s} {"vs G4":>9s}')
ref = None
for nm, r in (('G4 (adopted)', 'amip_G4_tundra'), ('I1 = G4+scf mode1', 'amip_I1_scf'),
              ('I2 = scf mode1 alone', 'amip_I2_scf_only'),
              ('I3 = G4+scf mode2', 'amip_I3_scf_sdor')):
    s, la, _ = clim(r, 'ssr')
    sdn, _, _ = clim(r, 'ssrd')
    sn, _, _ = clim(r, 'sd')
    if s is None or sn is None:
        print(f'  {nm:24s}   -- missing output')
        continue
    ny = len(glob.glob(f'{RT}/{r}/outdata/oifs/atm_remapped_1m_ssr_1m_*.nc'))
    A, D = s / ACC, sdn / ACC
    ww = np.cos(np.deg2rad(la))[:, None] * np.ones_like(A[0])
    w3 = np.broadcast_to(ww, A.shape)
    sun = (D > SUN_MIN) & np.broadcast_to(land, A.shape)
    snw = sun & (sn >= SNOW_WE)
    a_s = 1.0 - np.average(A[snw], weights=w3[snw]) / np.average(D[snw], weights=w3[snw])
    a_a = 1.0 - np.average(A[sun], weights=w3[sun]) / np.average(D[sun], weights=w3[sun])
    if ref is None:
        ref, d1, d2 = (a_s, a_a), 0.0, 0.0
    else:
        d1, d2 = a_s - ref[0], a_a - ref[1]
    print(f'  {nm:24s} {ny:4d} {a_s:12.4f} {d1:+9.4f} {a_a:10.4f} {d2:+9.4f}')

print("""
  Reading it: a negative "vs G4" is the scheme eating into the SNOW half (+0.0074).
  I3 removes ~70 % of it, I1 ~32 % -- yet I3 was calibrated to match I1 in the Siberian
  box, so globally it is twice as strong as intended. Evaluate I3 on its own merits
  rather than as an I1 equivalent.

  The SNOW-FREE half (+0.0080) -- crops, irrigated crops, semidesert, bare soil, and an
  in-season tundra error -- is untouched by every lever built so far.
""")
