"""Seasonal cycle of sea-ice VOLUME and EXTENT, both hemispheres, model against observations.

Two figures, because volume and extent are different measures and must never share an
axis: plots/ice_volume_seasonal_cycle.png and plots/ice_extent_seasonal_cycle.png.

Reading them together is the point.  Extent says how much ocean is covered, volume says
how much ice there is, and the Antarctic failure of this configuration is only visible in
the pair: roughly the right extent carrying roughly half the volume, i.e. a pack that is
far too thin rather than too small.

REFERENCES
  volume  GIOMAS 1989-2014 (global PIOMAS configuration; heff is effective thickness, so
          volume = sum(heff * dxt * dyt)).  A model with assimilation, not a measurement:
          ~10-20 % uncertainty in the Arctic and considerably more in the Antarctic, where
          thin snow-loaded ice is where altimetry is worst.
  extent  OSI-SAF monthly concentration on its 25 km grid, cells with ice_conc > 15 %.
  Both are PRESENT DAY (1989-2014) and the run is PRE-INDUSTRIAL, so the curves are not
  like-for-like - and the offset differs by hemisphere, which is the trap:
    NH  the record sits inside the satellite-era decline and PI is further above that
        again, so the model SHOULD exceed the observations by an amount nobody here has
        pinned down.  A model NH curve merely matching the blue one is probably too LOW.
    SH  Antarctic ice carried no significant trend over 1979-2015 (the decline begins
        ~2016, after this record ends), so PD ~ PI is a defensible first approximation
        and the SH panel can be read near face value.
  There is no observational PI sea-ice product; only proxy reconstructions, not held here.
  Bands are +-1 interannual sd.

Usage:  ARM=PICAL_ccnice NLAST=3 python3 scripts/figures/ice_seasonal_cycle.py
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

REPO = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive'
R    = os.environ.get('ROOT', '/work/bb1469/a270092/runtime/awiesm3-v3.4')
ARM  = os.environ.get('ARM', 'PICAL_ccnice')
NLAST= int(os.environ.get('NLAST', 3))
# Observation window.  The default is the whole record; OBS_Y0/OBS_Y1 narrow it.
# Pre-2000 is the better pre-industrial analogue in the NH: the Arctic decline steepens
# after ~2000, so 1989-1999 sits closer to PI while still being observed.  It makes
# little difference in the SH, which carries no significant trend over this record.
OY0, OY1 = int(os.environ.get('OBS_Y0', 1989)), int(os.environ.get('OBS_Y1', 2014))
GIO  = '/work/ab0246/a270092/obs/GIOMAS/GIOMAS_heff_miss.nc'
OSI  = {'NH': '/work/ab0246/a270092/obs/osisaf_nh.nc',
        'SH': '/work/ab0246/a270092/obs/osisaf_sh.nc'}
MON  = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()

OBS_C, MOD_C = '#2a78d6', '#eb6834'      # validated categorical slots 1 and 2
INK, MUTED   = '#0b0b0b', '#52514e'

# ---------------------------------------------------------------- model, last N years
def _load(path, var, scale):
    with xr.open_dataset(path, decode_times=False) as d:
        a = np.squeeze(np.asarray(d[var].values, float))/scale
    return a if a.size == 12 else None

def model_years():
    """Completed years from outdata, plus any already written by a running leg."""
    cand = {}
    for y in range(1850, 2101):
        s = {v: f'{R}/{ARM}/outdata/fesom/{v}.fesom.{y}.nc'
             for v in ('sivoln', 'sivols', 'siextentn', 'siextents')}
        if all(os.path.exists(p) for p in s.values()): cand[y] = s
    import glob
    for w in sorted(glob.glob(f'{R}/{ARM}/run_*/work')):
        for y in range(1850, 2101):
            s = {v: f'{w}/{v}.fesom_{y}-{y}.nc'
                 for v in ('sivoln', 'sivols', 'siextentn', 'siextents')}
            if all(os.path.exists(p) for p in s.values()): cand[y] = s
    out = {}
    for y, s in cand.items():
        a = {v: _load(p, v, 1e3 if v.startswith('sivol') else 1.0) for v, p in s.items()}
        if all(x is not None for x in a.values()): out[y] = a
    return out

MY = model_years()
yrs = sorted(MY)[-NLAST:]
print(f'{ARM}: using model years {yrs}')
mod = {k: np.array([MY[y][k] for y in yrs]) for k in ('sivoln', 'sivols', 'siextentn', 'siextents')}

# ---------------------------------------------------------------- GIOMAS volume
with xr.open_dataset(GIO, decode_times=False) as d:
    heff = np.nan_to_num(np.asarray(d['heff'].values, float))
    glat = np.asarray(d['lat_scaler'].values, float)
    garea= np.asarray(d['dxt'].values, float)*1e3*np.asarray(d['dyt'].values, float)*1e3
    gmon = np.asarray(d['month'].values, int)
gyr = 1989 + np.arange(heff.shape[0])//12          # verified contiguous Jan 1989 start
gsel = (gyr >= OY0) & (gyr <= OY1)
gv = {}
for h, m in (('NH', glat > 0), ('SH', glat < 0)):
    series = np.array([float((heff[t]*garea*m).sum())/1e12 for t in range(heff.shape[0])])
    gv[h] = np.array([[series[(gmon == k) & gsel].mean(), series[(gmon == k) & gsel].std()]
                      if ((gmon == k) & gsel).any() else [np.nan, np.nan] for k in range(1, 13)])

# ---------------------------------------------------------------- OSI-SAF extent
ov = {}
for h, path in OSI.items():
    with xr.open_dataset(path, decode_times=False) as d:
        conc = np.asarray(d['ice_conc'].values, np.float32)
        t    = np.asarray(d['time'].values, float)
        dx   = abs(float(d['xc'].values[1]-d['xc'].values[0]))   # km
        dy   = abs(float(d['yc'].values[1]-d['yc'].values[0]))
    cell = dx*dy                                                  # km2
    dates  = np.datetime64('1978-01-01') + t.astype('timedelta64[s]')
    months = (dates.astype('datetime64[M]').astype(int) % 12) + 1
    oyr    = dates.astype('datetime64[Y]').astype(int) + 1970
    ext = np.array([float(np.nansum(conc[i] > 15.0))*cell/1e6 for i in range(conc.shape[0])])
    ext = np.where((oyr >= OY0) & (oyr <= OY1), ext, np.nan)
    # this file carries Jan-Nov only (286 = 26 yr x 11); December has NO samples, so it
    # must come out NaN and leave a gap rather than a fabricated point.
    ov[h] = np.array([[np.nanmean(ext[months == k]), np.nanstd(ext[months == k])]
                      if np.isfinite(ext[months == k]).any() else [np.nan, np.nan]
                      for k in range(1, 13)])

# ---------------------------------------------------------------- draw
def figure(obs, mkey, ylab, title, fname, note):
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5), sharex=True)
    x = np.arange(12)
    for ax, h in zip(axes, ('NH', 'SH')):
        o = obs[h]
        good = np.isfinite(o[:, 0])
        ax.fill_between(x[good], (o[:, 0]-o[:, 1])[good], (o[:, 0]+o[:, 1])[good],
                        color=OBS_C, alpha=0.16, lw=0)
        ax.plot(x, np.where(good, o[:, 0], np.nan), color=OBS_C, lw=2.0, marker='o', ms=5,
                label='observed', zorder=3)
        m = mod[mkey[h]]
        ax.fill_between(x, m.min(0), m.max(0), color=MOD_C, alpha=0.18, lw=0)
        ax.plot(x, m.mean(0), color=MOD_C, lw=2.0, marker='s', ms=5,
                label=f'model {yrs[0]}-{yrs[-1]}', zorder=4)
        ax.set_title(f'{h}', fontsize=11, color=INK, loc='left', fontweight='bold')
        ax.set_xticks(x); ax.set_xticklabels(MON, fontsize=8, color=MUTED)
        ax.tick_params(axis='y', labelsize=8, colors=MUTED)
        ax.grid(alpha=0.18, lw=0.6)
        for s in ('top', 'right'): ax.spines[s].set_visible(False)
        for s in ('left', 'bottom'): ax.spines[s].set_color('#d8d7d2')
        ax.set_ylim(bottom=0)
        # Selective direct labels: each series' own max and min only.  The offset is
        # chosen from which curve is higher AT THAT MONTH, so the two never overlap.
        ob, mb = o[:, 0], m.mean(0)
        for arr, col in ((ob, OBS_C), (mb, MOD_C)):
            if not np.isfinite(arr).any(): continue
            for j in (int(np.nanargmax(arr)), int(np.nanargmin(arr))):
                other = mb[j] if col == OBS_C else ob[j]
                up = (not np.isfinite(other)) or (arr[j] >= other)
                # a label placed below a near-zero minimum lands on the month ticks, so
                # flip it up and nudge it sideways instead
                lo, hi = ax.get_ylim()
                if not up and (arr[j]-lo)/(hi-lo) < 0.10:
                    up, dxp = True, 16
                else:
                    dxp = 0
                ax.annotate(f'{arr[j]:.1f}', (x[j], arr[j]), textcoords='offset points',
                            xytext=(dxp, 10 if up else -16), ha='center',
                            fontsize=8, color=col, zorder=5)
    axes[0].set_ylabel(ylab, fontsize=9, color=MUTED)
    axes[0].legend(frameon=False, fontsize=9, loc='upper right', labelcolor=MUTED)
    fig.suptitle(title, fontsize=12.5, color=INK, x=0.008, ha='left', y=0.99)
    import textwrap
    fig.text(0.008, 0.005, '\n'.join(textwrap.wrap(note, 158)),
             fontsize=7.6, color=MUTED, ha='left', va='bottom')
    fig.tight_layout(rect=[0, 0.085, 1, 0.94])
    p = f'{REPO}/plots/{fname}'
    fig.savefig(p, dpi=160, facecolor='#fcfcfb'); plt.close(fig)
    print('wrote', p)

figure(gv, {'NH': 'sivoln', 'SH': 'sivols'}, 'volume  [10$^3$ km$^3$]',
       f'Sea-ice volume, seasonal cycle - {ARM} against GIOMAS',
       'ice_volume_seasonal_cycle.png',
       f'Observed = GIOMAS {OY0}-{OY1} (assimilation, not measurement; Antarctic weakly constrained), band +-1 sd. '
       'Model band = min-max over the years shown. OBSERVATIONS ARE PRESENT DAY, model is pre-industrial: the NH offset is large (PI well above the satellite era), the SH offset is small (no significant Antarctic trend over this record), so only the SH panel is near like-for-like.')

figure(ov, {'NH': 'siextentn', 'SH': 'siextents'}, 'extent  [10$^6$ km$^2$]',
       f'Sea-ice extent, seasonal cycle - {ARM} against OSI-SAF',
       'ice_extent_seasonal_cycle.png',
       f'Observed = OSI-SAF {OY0}-{OY1} monthly concentration, cells > 15 %, 25 km grid, band +-1 sd; December absent from the record. '
       'Model band = min-max over the years shown. OBSERVATIONS ARE PRESENT DAY, model is pre-industrial: large NH offset, small SH offset. Extent saturates geometrically - read it beside the volume figure.')
