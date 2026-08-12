"""Draw the ocp-tool land-sea pipeline: what each stage did to lsm, slt and cl.

Reads the npz written by lake_pipeline_extract.py.  Two figures:

  plots/lake_pipeline_stages.png   the Caspian, through all five stages, three
                                   fields each -- where the story is visible
  plots/lake_pipeline_masks.png    the three land-sea masks that leave the tool,
                                   globally, and where they disagree

WHY THE CASPIAN.  It is the cleanest case in the domain: 27 contiguous cells,
all of them flagged cl = 1 by ECMWF, all of them flipped to land by the FESOM
reconciliation, and large enough that a wrong surface type there is a regional
climate signal rather than a rounding error.  Everything the pipeline does to a
water body it does to the Caspian, in one place, at a readable size.

RENDERING.  The reduced Gaussian grid has a different number of longitudes in
every latitude row, so it is not an array and pcolormesh cannot draw it.  Cells
are built as native quads from the row structure and drawn with a PolyCollection,
which is the model's own geometry rather than an interpolation of it.

COLOUR, and why soil type is shown as three classes and not eight.  A choropleth
puts every pair of categories side by side, and eight categorical hues cannot be
told apart under that condition -- the palette's own validator fails them.  Three
can.  The three chosen are the distinction that actually matters here: no soil at
all (an ocean cell), mineral (types 1-5), and organic (types 6-7), because it is
organic soil on a cell whose state was spun as mineral that produces the
"Illegal Frac_air" abort.  The exact type is printed on the panel, so nothing is
hidden by the binning.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

ROOT = ('/work/ab0246/a270092/postprocessing/'
        'investigation_awiesm3_high_lat_cold_bias_global_toa_positive')
MESH = os.environ.get('OCP_MESH', 'CORE3')
NPZ = f'{ROOT}/scripts/analysis/lake_pipeline_{MESH}.npz'
# The LPJ-GUESS soil-type files. The ICMGG "_v2" is a no-op on CORE3 (byte
# identical), so the soil change the campaign calls "_v2" lives here instead:
# slt_TCO95.nc is the Dec-2025 file the runs still stage, slt_TCO95_CORE3*.nc
# is the corrected one.
SLT_DIR = f'/work/ab0246/a270092/software/ocp-tool/output/TCO95_{MESH}/lpj-guess'
SLT_OLD = f'{SLT_DIR}/slt_TCO95.nc'
SLT_NEW = f'{SLT_DIR}/slt_TCO95_{MESH}_v2.nc'

# ---------------------------------------------------------------- palette
# Validated with the dataviz skill's validator, light surface #fcfcfb,
# --pairs all (choropleth): all checks pass, worst CVD dE 9.2, worst
# normal-vision dE 24.0.  The aqua slot warns on contrast, so every class is
# also carried by a visible label -- the relief rule.
SURFACE = '#fcfcfb'
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#898781'
GRIDLINE, BASELINE = '#e1e0d9', '#c3c2b7'
SLOT1, SLOT2, SLOT3 = '#2a78d6', '#eb6834', '#1baf7a'   # blue, orange, aqua
CRITICAL = '#d03b3b'

# Sequential blue, steps 100..700 from the reference palette.
BLUES = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7',
         '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281',
         '#0d366b']
# Sequential orange -- the second sequential context takes the next slot's hue.
ORANGES = ['#fde8dd', '#fbcfb8', '#f8b593', '#f59a6e', '#f28150', '#eb6834',
           '#d4552a', '#b34521', '#8f3719', '#6b2912']
CMAP_LSM = LinearSegmentedColormap.from_list('lsm', BLUES)
CMAP_CL = LinearSegmentedColormap.from_list('cl', ORANGES)

# no soil / mineral / organic
CMAP_SLT = ListedColormap(['#e8e7e1', SLOT1, SLOT2])
NORM_SLT = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], 3)

CASPIAN = (35.0, 50.0, 45.0, 63.0)      # lat0, lat1, lon0, lon1


def build_quads(lat, lon):
    """Native reduced-Gaussian cell corners, as an (n, 4, 2) vertex array."""
    order = np.lexsort((lon, -lat))
    rows, start = [], 0
    slat = lat[order]
    for i in range(1, len(order) + 1):
        if i == len(order) or slat[i] != slat[start]:
            rows.append(order[start:i])
            start = i
    lats_of_row = [lat[r[0]] for r in rows]
    verts = np.zeros((len(lat), 4, 2))
    for k, r in enumerate(rows):
        y = lats_of_row[k]
        up = (y + lats_of_row[k - 1]) / 2 if k > 0 else 90.0
        dn = (y + lats_of_row[k + 1]) / 2 if k < len(rows) - 1 else -90.0
        half = 180.0 / len(r)
        x = lon[r]
        verts[r, 0] = np.column_stack([x - half, np.full(len(r), dn)])
        verts[r, 1] = np.column_stack([x + half, np.full(len(r), dn)])
        verts[r, 2] = np.column_stack([x + half, np.full(len(r), up)])
        verts[r, 3] = np.column_stack([x - half, np.full(len(r), up)])
    return verts


def draw(ax, verts, values, sel, cmap, norm=None, vmin=None, vmax=None,
         edge=None, lw=0.0):
    pc = PolyCollection(verts[sel], array=np.asarray(values)[sel], cmap=cmap,
                        norm=norm, edgecolors=edge or 'none', linewidths=lw)
    if norm is None:
        pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    return pc


def soil_class(slt):
    """0 = no soil, 1 = mineral (1-5), 2 = organic (6-7)."""
    out = np.zeros_like(slt)
    out[(slt >= 1) & (slt <= 5)] = 1
    out[slt >= 6] = 2
    return out


def note(ax, x, y, text, color=None, bold=False, ha='left'):
    """An annotation that stays readable on top of a dark fill."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=7.5, ha=ha,
            color=color or INK2, fontweight='bold' if bold else 'normal',
            bbox=dict(facecolor=SURFACE, alpha=0.88, edgecolor='none',
                      boxstyle='round,pad=0.22'))


def style(ax, extent, title=None):
    ax.set_xlim(extent[2], extent[3])
    ax.set_ylim(extent[0], extent[1])
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(BASELINE)
        s.set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=7, length=2, width=0.6)
    ax.grid(True, color=GRIDLINE, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(title, color=INK, fontsize=9, pad=4)


d = np.load(NPZ)
lat, lon = d['lat'], d['lon']
lon = np.where(lon > 180, lon - 360, lon)
verts = build_quads(lat, lon)
print(f'  {len(lat)} cells, quads built')

_ident = all(np.array_equal(d[f'v1_{k}'], d[f'v2_{k}'])
             for k in ('lsm', 'slt', 'cl', 'dl'))
print(f'  mesh {MESH}: ICMGG v1 and v2 identical = {_ident}')

STAGES = [
    ('ecmwf',  'ECMWF as shipped',
     'water, flagged lake'),
    ('legacy', 'legacy rule (2020-22, removed)',
     'promoted to land, soil 6'),
]
if _ident:
    # On CORE3 the two ICMGGs are byte-identical, so showing both would invent
    # a distinction that does not exist on this mesh.
    STAGES.append(('v2', 'current master  (v1 == v2)',
                   'flip + NN fill: lake destroyed'))
else:
    STAGES.append(('v1', 'v1  master before the NN fill',
                   'lake kept, soil 6 everywhere'))
    STAGES.append(('v2', 'v2  after the NN fill',
                   'soil fixed, lake destroyed'))
STAGES.append(('lake', 'proposed  restore_flipped_lakes',
               'soil fixed, lake kept'))

# ============================================================ figure 1: stages
fig, axes = plt.subplots(len(STAGES), 3, figsize=(10.5, 15.0))
fig.patch.set_facecolor(SURFACE)
box = ((lat >= CASPIAN[0] - 3) & (lat <= CASPIAN[1] + 3)
       & (lon >= CASPIAN[2] - 3) & (lon <= CASPIAN[3] + 3))
# The cells the panels are ABOUT: the Caspian proper, i.e. the ones ECMWF
# flagged as lake. Counting over the whole drawn box instead would mix in the
# surrounding land and make every tally meaningless.
casp = (box & (lat >= CASPIAN[0]) & (lat <= CASPIAN[1])
        & (lon >= CASPIAN[2]) & (lon <= CASPIAN[3])
        & (d['ecmwf_cl'] >= 0.5))
NCASP = int(casp.sum())
print(f'  Caspian cells (ecmwf cl >= 0.5 in the box): {NCASP}')

for r, (key, name, rownote) in enumerate(STAGES):
    lsm, slt, cl = d[f'{key}_lsm'], d[f'{key}_slt'], d[f'{key}_cl']

    draw(axes[r, 0], verts, lsm, box, CMAP_LSM, vmin=0, vmax=1,
         edge=BASELINE, lw=0.15)
    style(axes[r, 0], CASPIAN, 'land-sea mask  lsm' if r == 0 else None)
    axes[r, 0].text(-0.30, 0.5, name, transform=axes[r, 0].transAxes,
                    rotation=90, va='center', ha='center', fontsize=9,
                    color=INK, fontweight='bold')
    axes[r, 0].text(-0.20, 0.5, rownote, transform=axes[r, 0].transAxes,
                    rotation=90, va='center', ha='center', fontsize=7.5,
                    color=INK2)
    note(axes[r, 0], 0.03, 0.04,
         f'land: {int((lsm[casp] >= 0.5).sum())}/{NCASP}')

    draw(axes[r, 1], verts, soil_class(slt), box, CMAP_SLT, norm=NORM_SLT,
         edge=BASELINE, lw=0.15)
    style(axes[r, 1], CASPIAN, 'soil type  slt' if r == 0 else None)
    vals = np.unique(slt[casp].astype(int))
    n_org = int((slt[casp] >= 6).sum())
    note(axes[r, 1], 0.03, 0.04,
         'types: ' + ','.join(str(v) for v in vals)
         + f'   organic: {n_org}/{NCASP}',
         color=CRITICAL if n_org == NCASP else INK2, bold=n_org == NCASP)

    draw(axes[r, 2], verts, cl, box, CMAP_CL, vmin=0, vmax=1,
         edge=BASELINE, lw=0.15)
    style(axes[r, 2], CASPIAN, 'lake cover  cl' if r == 0 else None)
    n_lake = int((cl[casp] >= 0.5).sum())
    note(axes[r, 2], 0.03, 0.04, f'lake: {n_lake}/{NCASP}',
         color=CRITICAL if n_lake == 0 else INK2, bold=n_lake == 0)

for ax in axes[:-1, :].ravel():
    ax.set_xticklabels([])

cax = fig.add_axes([0.13, 0.055, 0.20, 0.008])
fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP_LSM), cax=cax,
             orientation='horizontal').set_label(
                 'lsm   0 = ocean, 1 = land', color=INK2, fontsize=7.5)
cax.tick_params(colors=MUTED, labelsize=7)

lax = fig.add_axes([0.40, 0.030, 0.22, 0.035])
lax.axis('off')
lax.legend(handles=[Patch(facecolor='#e8e7e1', edgecolor=BASELINE,
                          label='no soil (ocean)'),
                    Patch(facecolor=SLOT1, edgecolor=BASELINE,
                          label='mineral, types 1-5'),
                    Patch(facecolor=SLOT2, edgecolor=BASELINE,
                          label='organic, types 6-7')],
           loc='center', frameon=False, fontsize=7.5, labelcolor=INK2,
           handlelength=1.2, borderpad=0)

cax2 = fig.add_axes([0.70, 0.055, 0.20, 0.008])
fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP_CL), cax=cax2,
             orientation='horizontal').set_label(
                 'cl   lake fraction of the box', color=INK2, fontsize=7.5)
cax2.tick_params(colors=MUTED, labelsize=7)

fig.text(0.5, 0.988,
         f'The Caspian through the ocp-tool pipeline ({MESH})',
         ha='center', va='top', fontsize=13, color=INK, fontweight='bold')
fig.text(0.5, 0.971,
         'ECMWF ships the Caspian as water flagged cl = 1 with lsm < 0.5. FESOM has no wet nodes there, so the flip must make it land - that part is right.\n'
         'The nearest-neighbour fill then overwrites cl from a dry land '
         'neighbour, so all 27 stop being lakes -- and nothing puts it back: '
         'the shipped file has cl = 0.',
         ha='center', va='top', fontsize=8.5, color=INK2, linespacing=1.5)
fig.subplots_adjust(left=0.09, right=0.98, top=0.938, bottom=0.085,
                    hspace=0.10, wspace=0.13)
p1 = f'{ROOT}/plots/lake_pipeline_stages.png'
fig.savefig(p1, dpi=170, facecolor=SURFACE)
print(f'  wrote {p1}')
plt.close(fig)

# ============================================================= figure 2: masks
v2lsm = d['v2_lsm']
fesom = d['fesom']
oa, ol, orr = d['oasis_A'], d['oasis_L'], d['oasis_R']

oland = d['oasis_land']
# Dateline land fraction: the CORE2 defect metric, computed before it is used
# both by the panel callout and by the header.
_dl = np.abs(np.abs(lon) - 180) < 0.5
_dlf = 100 * (fesom[_dl] >= 0.5).mean()
PANELS = [
    ('OpenIFS lsm, pristine', d['ecmwf_lsm'],
     'land fraction as ECMWF ships it'),
    ('OpenIFS lsm, after the flip', v2lsm,
     'reconciled with the FESOM mesh'),
    ('FESOM mesh coverage', fesom,
     'what the flip is driven BY, not a product of it'),
    ('OASIS A096.msk', oa.astype(float),
     'measured: identical to the post-flip lsm'),
    ('OASIS L096.msk', ol.astype(float),
     'measured: identical to the PRE-flip lsm'),
    ('OASIS TCO95-land.msk', oland.astype(float),
     'measured: identical to 1 - A096 (1 = masked out)'),
]

fig, axes = plt.subplots(3, 3, figsize=(15.5, 8.6))
fig.patch.set_facecolor(SURFACE)
allsel = np.ones(len(lat), bool)
GLOBE = (-90, 90, -180, 180)

for k, (name, field, panelnote) in enumerate(PANELS):
    ax = axes.ravel()[k]
    draw(ax, verts, field, allsel, CMAP_LSM, vmin=0, vmax=1)
    style(ax, GLOBE, name)
    ax.text(0.01, -0.16, panelnote, transform=ax.transAxes, fontsize=7.5,
            color=INK2)
    note(ax, 0.985, 0.05, f'{int(np.round(field).sum())} set', ha='right')

# The disagreements worth naming. A096 minus L096 is NOT among them: it was
# measured to be exactly the flip, so it would have duplicated the first panel.
flipped_signed = np.where(
    d['flipped_to_land'], 1.0,
    np.where(np.round(d['ecmwf_lsm']) > np.round(v2lsm), -1.0, 0.0))
# The two populations the lake branch has to tell apart: cells the input already
# called lake (restoring them invents nothing) and cells it called open ocean
# (making those lakes would invent both cover and depth).
was_lake = d['ecmwf_cl'] >= 0.5
populations = np.where(d['flipped_to_land'] & was_lake, 1.0,
                       np.where(d['flipped_to_land'], -1.0, 0.0))
_nl = int(d['flipped_to_land'].sum())
_no = int(((np.round(d['ecmwf_lsm']) > np.round(v2lsm))).sum())
_wl = int((d['flipped_to_land'] & was_lake).sum())
DIFFS = [
    ('What the FESOM mesh flips', flipped_signed,
     f'red = made land ({_nl}), blue = made ocean ({_no}).  '
     f'A096 - L096 is exactly this'),
    (f'The two populations among the {_nl}', populations,
     f'red = input already called it lake ({_wl}), '
     f'blue = input called it ocean ({_nl - _wl})'),
    ('Lake cover destroyed by the NN fill', d['ecmwf_cl'] - d['v2_cl'],
     'orange = cl the nearest-neighbour fill overwrote from a dry neighbour'),
]

# The dateline defect, found while drawing this: the CORE2 polygon mask marks
# EVERY cell within 0.5 deg of 180 as land. Called out rather than left for the
# reader to notice, because it is 44 % of the flip.
DATELINE = np.abs(180 - np.abs(lon)) < 0.5
DIVERGE = LinearSegmentedColormap.from_list(
    'div', [SLOT1, '#f0efec', CRITICAL])

for k, (name, field, panelnote) in enumerate(DIFFS):
    ax = axes.ravel()[6 + k] if k < 3 else None
    cmap = CMAP_CL if k == 2 else DIVERGE
    if k == 2:
        draw(ax, verts, np.abs(field), allsel, cmap, vmin=0, vmax=1)
    else:
        draw(ax, verts, field, allsel, cmap, vmin=-1, vmax=1)
    style(ax, GLOBE, name)
    ax.text(0.01, -0.16, panelnote, transform=ax.transAxes, fontsize=7.5,
            color=INK2)
    nz = int((np.abs(field) > 0.5).sum())
    note(ax, 0.985, 0.05, f'{nz} cells', ha='right')
    if k == 0:
        if _dlf > 50:
            ax.annotate('CORE2 dateline defect:\n'
                        'all 192 cells at lon 180\n'
                        'marked land (latent)',
                        xy=(178, 10), xytext=(96, -52), fontsize=7.5,
                        color=CRITICAL, fontweight='bold', ha='center',
                        bbox=dict(facecolor=SURFACE, alpha=0.92,
                                  edgecolor=CRITICAL, linewidth=0.7,
                                  boxstyle='round,pad=0.3'),
                        arrowprops=dict(arrowstyle='-|>', color=CRITICAL,
                                        lw=1.0, shrinkB=1))

fig.text(0.5, 0.992,
         f'{MESH}: four named OASIS masks, two distinct land-sea partitions',
         ha='center', va='top', fontsize=13, color=INK, fontweight='bold')
_A, _L = int(np.round(oa).sum()), int(np.round(ol).sum())
fig.text(0.5, 0.966,
         f'Measured, not assumed: A096 is bit-identical to the post-flip '
         f'OpenIFS mask ({_A} = FESOM land), L096 to the PRE-flip one ({_L}), '
         f'R096 and TCO95-land are both 1 - A096.\n'
         f'So A096 - L096 is not a lake mask; it is exactly the flip. On this '
         f'mesh the dateline is clean: {_dlf:.0f} % of the 192 cells at lon '
         f'180 are land, against 100 % on CORE2, where the polygon pass aborts '
         f'with "Dateline triangulation failed".\n'
         f'That CORE2 defect never reached a run - every run checked stages '
         f'the {MESH} product.',
         ha='center', va='top', fontsize=8.5, color=INK2, linespacing=1.5)
fig.subplots_adjust(left=0.025, right=0.985, top=0.875, bottom=0.075,
                    hspace=0.40, wspace=0.08)
p2 = f'{ROOT}/plots/lake_pipeline_masks.png'
fig.savefig(p2, dpi=155, facecolor=SURFACE)
print(f'  wrote {p2}')
plt.close(fig)
