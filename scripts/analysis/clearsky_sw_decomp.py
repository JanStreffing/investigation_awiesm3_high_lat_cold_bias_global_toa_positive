"""Where the clear-sky shortwave deficit comes from: three separable causes.

The chain that led here. The boreal surface budget is nearly spent (G4 took ~0.95 K of a
1.3-1.6 K boreal-specific envelope), leaving the global tropospheric cold bias of
report sub:vprof as the dominant remaining term for land 2 m temperature. That bias is
energetically consistent with a shortwave deficit, and the deficit is entirely CLEAR-SKY:
cloud radiative effect is the WRONG SIGN to explain it (model clouds reflect 0.90 W/m2 LESS
than observed). Every lever in 38 runs was a cloud or surface-vegetation lever; none touched
clear-sky shortwave.

METHOD NOTE, and it matters. CERES TOA fluxes are accurate to ~0.5 W/m2; the CERES SURFACE
product is a derived radiative-transfer calculation (Kato et al. 2018) with several W/m2
regional uncertainty. So every headline number here is anchored at TOA, and surface
quantities are used only for attribution, never for the size of the deficit. Albedos are
computed assumption-free -- each dataset divided by its OWN downward flux -- because an
earlier pass that borrowed the CERES denominator for the model inflated the land albedo
excess from +0.015 to +0.030.

THE RESULT: global clear-sky absorbed SW is -2.68 W/m2, and it separates into three
independent causes that need three different fixes.
"""
import numpy as np, xarray as xr, os, warnings
warnings.filterwarnings('ignore')

from runs import RT, LSMF, OBS

ACC = 3600.0
PD = list(range(1990, 2015))
lsm = xr.open_dataset(LSMF)['lsm'].isel(time_counter=0).values


def clim(v):
    acc, n, lat, lon = None, 0, None, None
    for y in PD:
        f = f'{RT}/amip_presentday/outdata/oifs/atm_remapped_1m_{v}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f)
        a = d[v].values
        lat, lon = d[v].lat.values, d[v].lon.values
        d.close()
        acc = a if acc is None else acc + a
        n += 1
    return acc / n, lat, lon


tsrc, lat, lon = clim('tsrc')
tisr, _, _ = clim('tisr')
ssr, _, _ = clim('ssr')
ssrd, _, _ = clim('ssrd')
ci, _, _ = clim('ci')

ds = xr.open_dataset(OBS)
ds = xr.concat([ds.isel(lon=[-1]).assign_coords(lon=ds.lon.values[-1:] - 360.0), ds,
                ds.isel(lon=[0]).assign_coords(lon=ds.lon.values[:1] + 360.0)], dim='lon')
tl = np.clip(lat, -89.5, 89.5)
tlo = np.where(lon < 0, lon + 360, lon)


def ip(v):
    return ds[v].interp(lat=xr.DataArray(tl, dims='y'),
                        lon=xr.DataArray(tlo, dims='x')).values.mean(0)


Ctoa_clr = ip('toa_sw_clr_t_clim')
Cs_dn_all, Cs_up_all = ip('sfc_sw_down_all_clim'), ip('sfc_sw_up_all_clim')
M_tsrc, M_isr = tsrc.mean(0) / ACC, tisr.mean(0) / ACC
M_ssr, M_ssrd = ssr.mean(0) / ACC, ssrd.mean(0) / ACC

w = np.cos(np.deg2rad(lat))[:, None] * np.ones_like(M_tsrc)
Wt = w.sum()
land = lsm > 0.5
ice = (np.nanmean(ci, 0) > 0.15) & (~land)
ocean = (~land) & (~ice)
l180 = ((lon + 180) % 360) - 180
diff = M_tsrc - (M_isr - Ctoa_clr)          # model minus CERES, clear-sky absorbed SW

print(f'GLOBAL clear-sky absorbed SW, model minus CERES: '
      f'{np.average(diff, weights=w):+.2f} W/m2   (TOA -- the robust anchor)\n')

print('CAUSE 1 -- ATMOSPHERIC over-scattering, seen over ocean where the surface is right')
print('  Ice-free ocean ALL-SKY surface albedo: model %.4f vs CERES %.4f (%+.4f) -- correct.'
      % (1 - np.average(M_ssr[ocean], weights=w[ocean]) / np.average(M_ssrd[ocean], weights=w[ocean]),
         np.average(Cs_up_all[ocean], weights=w[ocean]) / np.average(Cs_dn_all[ocean], weights=w[ocean]),
         (1 - np.average(M_ssr[ocean], weights=w[ocean]) / np.average(M_ssrd[ocean], weights=w[ocean]))
         - np.average(Cs_up_all[ocean], weights=w[ocean]) / np.average(Cs_dn_all[ocean], weights=w[ocean])))
print('  So the ocean TOA excess must be atmospheric. In the aerosol-clean SOUTHERN')
print('  hemisphere it grows steadily with latitude, over a nearly constant surface:\n')
print(f'    {"lat band":14s} {"d(TOA alb)":>11s}')
for a, b in ((-15, 0), (-30, -15), (-45, -30), (-60, -45)):
    m = ocean & np.broadcast_to((lat[:, None] >= a) & (lat[:, None] < b), M_tsrc.shape)
    if m.sum() < 15:
        continue
    isr = np.average(M_isr[m], weights=w[m])
    am = 1 - np.average(M_tsrc[m], weights=w[m]) / isr
    ac = np.average(Ctoa_clr[m], weights=w[m]) / isr
    print(f'    {a:4d} to {b:<6d} {am-ac:+11.4f}')
print('\n  A slant-path signature over a correct surface: clear-sky scattering (sea-salt or')
print('  background aerosol, or Rayleigh) is too strong. NOTE the ocean albedo knobs are at')
print('  their defaults in these runs -- RALBSEAD=0.06 diffuse (susrad_mod.F90:107) and the')
print('  Taylor direct-beam formula (surfrad_ctl_mod.F90:588-590); TMPRALBSEAD is NOT set in')
print('  the campaign runscripts, only in the older base one. So this is not a stale override.')

print('\nCAUSE 2 -- LAND surface albedo too high (all-sky, assumption-free):')
for nm, m in (('land', land), ('ocean ice-free', ocean), ('sea ice', ice)):
    am = 1 - np.average(M_ssr[m], weights=w[m]) / np.average(M_ssrd[m], weights=w[m])
    ac = np.average(Cs_up_all[m], weights=w[m]) / np.average(Cs_dn_all[m], weights=w[m])
    print(f'    {nm:16s} model {am:.4f}  CERES {ac:.4f}  {am-ac:+.4f}')

print('\nCAUSE 3 -- contributions to the global -2.68 W/m2:')
dust = np.broadcast_to(((lat[:, None] >= 0) & (lat[:, None] <= 35))
                       & ((l180[None, :] >= -45) & (l180[None, :] <= 80)), M_tsrc.shape)
tot = 0.0
for nm, m in (('dust belt 0-35N 45W-80E', dust), ('sea ice', ice),
              ('land outside dust belt', land & ~dust), ('ocean outside dust belt', ocean & ~dust)):
    c = (diff[m] * w[m]).sum() / Wt
    tot += c
    print(f'    {nm:26s} {100*w[m].sum()/Wt:5.1f}% area  {c:+7.3f}  '
          f'(local {np.average(diff[m], weights=w[m]):+6.2f})')
print(f'    {"sum":26s}                {tot:+7.3f}')

print("""
SUMMARY -- three separable causes, three different fixes:

  1. ATMOSPHERIC clear-sky over-scattering        ~-1.1 W/m2 global
     Ocean surface albedo is right, yet TOA clear-sky albedo is +0.0065 too high and
     grows with latitude in the clean Southern Hemisphere. Candidate: sea-salt or
     background aerosol optical depth, or Rayleigh. NEVER EXAMINED by this campaign.

  2. LAND surface albedo +0.0153 too high         ~-0.7 W/m2 global
     BUT SEE land_albedo_snow_split.py, which supersedes the attribution below. Snow is
     a TILE overlying every vegetation type, not a type itself, and masking it per cell
     per month splits the +0.0154 almost exactly in half: +0.0074 snow, +0.0080 snow-free
     surface. Once snow is removed EVERY FOREST TYPE is within +0.006 (tropical broadleaf
     -0.0003), so sub:albreg's claim about RVVEGALB's high-vegetation entries HOLDS. What
     survives is crops (+0.0204), semidesert (+0.0169), bare soil (+0.0145), irrigated
     crops (+0.0138) and an in-season tundra error (+0.0105) -- the sparse and cultivated
     surfaces where the soil background shows through. The regional numbers quoted here
     (Europe +0.0099, Great Plains +0.0153, Sahara +0.0207, Siberia +0.0165) are
     all-months values and are therefore snow-contaminated outside the subtropics.

  3. SEA-ICE albedo +0.0182 surface / +0.0327 TOA ~-0.45 W/m2 global
     The largest per-unit-area error of the three (-6.8 W/m2 locally). Directly relevant
     to the coupled sea-ice tuning already under way.

  The dust belt contributes only -0.40 of the -2.68, so dust is real but secondary.

CAVEATS. CERES surface fluxes are derived, not observed -- attribution between surface and
atmosphere carries several W/m2 of uncertainty, though the TOA total does not. The clean-SH
zenith-angle result is the strongest single piece of evidence because it needs only TOA
fluxes and a surface that is independently verified as correct.
""")
