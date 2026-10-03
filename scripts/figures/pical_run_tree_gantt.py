"""The PICAL lineage as a Gantt chart: who branched from whom, when, and what changed.

Five runs, each a span of MODEL years.  Those years carry no forcing meaning -- CO2, Ndep
and land use are fixed at 1850 -- with one exception that cost this campaign a great deal:
OpenIFS read the calendar year for its MACv2-SP anthropogenic aerosol until the fix of
2026-09-20, so 1850-2019 of the momixoff line carries a transient aerosol ramp.

Solid bar  = years the run integrated itself.
Hatched    = years symlinked from the parent so the folder reads as one record.
Open bar   = queued or still to run.
Labels on the connectors say what changed AT the branch; labels in the gutters are changes
made in-line, to a run already going.

Palette validated with the dataviz skill's checker: all checks pass, the one CVD warning
(teal/magenta, deltaE 7.1 deutan) is covered by the per-row labels, which are the required
secondary encoding.

Out: plots/pical_run_tree_gantt.png
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

C = ['#1f6fb4', '#c07000', '#4a7d1f', '#00856b', '#a03a8f']
INK, MUTED, GRID, BG = '#1a1a1a', '#5a5a5a', '#dcdcd8', '#fcfcfb'

# name, own(start,end) or None, inherited or None, parent row, branch yr, status, queued(end)
RUNS = [
    ('PICAL',                (1850, 1929), None,         None, None, 'ended',    None),
    ('PICAL_momixoff',       (1920, 2049), (1850, 1919), 0,    1920, 'ended',    None),
    ('PICAL_snowalb',        None,         None,         1,    1940, 'running',  1989),
    ('PICAL_v35def',         (2020, 2049), None,         1,    2020, 'killed',   None),
    ('PICAL_v35def_darkice', (2050, 2099), None,         3,    2050, 'complete', None),
]
# what changed at each branch: child row -> text
BRANCH_TXT = {1: 'momix off', 2: 'albsn/albsnm\nback to 0.75/0.65',
              3: 'S4 + RSNOWLIN2\ndropped', 4: 'albi/albim\n-> 0.66/0.64'}
# in-line changes: row, year, text, gutter (-1 above the bar, +1 below), h-align
INLINE = [(0, 1880, 'LPJ-GUESS offline\nstate swapped in', -1, 'center', 0.40),
          (1, 1950, 'sea-ice albedos\n-> FESOM defaults',  -1, 'left',   0.36),
          (1, 2020, 'MACv2-SP aerosol\nfixed to 1850',     -1, 'center', 0.36)]
H, X0, X1, RGUT = 0.16, 1838, 2190, 2104

fig, ax = plt.subplots(figsize=(15.0, 7.2))
fig.patch.set_facecolor(BG); ax.set_facecolor(BG)

for i, (name, own, inh, par, byr, status, qend) in enumerate(RUNS):
    col = C[i]
    if par is not None:                                    # connector + its label
        ax.plot([byr, byr], [par + H, i - H], color=MUTED, lw=1.1, ls=(0, (2, 2)), zorder=2)
        ax.plot([byr], [i - H], marker='o', ms=4.5, mfc=BG, mec=MUTED, mew=1.2, zorder=5)
        ax.text(byr - 6, i - H - 0.30, BRANCH_TXT[i], ha='right', va='center',
                fontsize=8.4, color=INK, linespacing=1.3, zorder=6)
    if inh:
        ax.barh(i, inh[1] - inh[0] + 1, left=inh[0], height=H * 2, color=col, alpha=0.22,
                hatch='////', edgecolor=col, lw=0.7, zorder=3)
    n = 0
    if own:
        n = own[1] - own[0] + 1
        ax.barh(i, n, left=own[0], height=H * 2, color=col, zorder=3)
    if qend:                                               # queued remainder
        s = own[1] + 1 if own else byr
        ax.barh(i, qend - s + 1, left=s, height=H * 2, facecolor='none', edgecolor=col,
                lw=1.4, ls=(0, (3, 2)), zorder=3)
    end = qend or own[1]
    if status == 'killed':
        ax.plot([end + 5], [i], marker='x', ms=9, mew=2.3, color=col, zorder=5)
    elif status == 'running':
        ax.annotate('', xy=(end + 20, i), xytext=(end + 4, i), zorder=5,
                    arrowprops=dict(arrowstyle='-|>', color=col, lw=1.7))
    ax.text(RGUT, i, f'{n if n else qend-byr+1} yr · {status}', ha='left', va='center',
            fontsize=9.2, color=MUTED)

for r, yr, txt, gut, ha, off in INLINE:
    ax.plot([yr], [r + gut * H], marker='v' if gut > 0 else '^', ms=6, color=INK, zorder=6)
    dx = 6 if ha == 'left' else 0
    ax.text(yr + dx, r + gut * off, txt, ha=ha, va='center', fontsize=8.4,
            color=INK, linespacing=1.3, zorder=6)

ax.set_yticks(range(len(RUNS)))
ax.set_yticklabels([r[0] for r in RUNS], fontsize=10.5)
for lbl, col in zip(ax.get_yticklabels(), C):
    lbl.set_color(INK)
ax.set_ylim(len(RUNS) - 0.40, -0.72)
ax.set_xlim(X0, X1)
ax.set_xticks(range(1850, 2101, 50))
ax.set_xlabel('model year   (forcing fixed at 1850 throughout; year labels are arbitrary)',
              fontsize=10, color=MUTED, labelpad=8)
ax.tick_params(axis='x', colors=MUTED, labelsize=9.5)
ax.tick_params(axis='y', length=0, pad=8)
ax.grid(axis='x', color=GRID, lw=0.8, zorder=0)
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color(GRID)
from matplotlib.patches import Rectangle
ax.add_patch(Rectangle((X0, -0.72), 2019.5 - X0, 1.42 + 0.72, facecolor='#b03030',
                       alpha=0.055, edgecolor='none', zorder=0))
ax.plot([2019.5, 2019.5], [-0.72, 1.42], color='#b03030', alpha=0.35, lw=1.0, zorder=1)
ax.text(1928, 1.30, 'transient anthropogenic aerosol (the MACv2-SP bug) -- these years only',
        ha='center', va='center', fontsize=8.6, color='#8a3030', style='italic', zorder=6)

ax.legend(handles=[Patch(facecolor='#8a8a8a', label='integrated by this run'),
                   Patch(facecolor='#8a8a8a', alpha=0.22, hatch='////', edgecolor='#8a8a8a',
                         label='symlinked from parent'),
                   Patch(facecolor='none', edgecolor='#8a8a8a', ls=(0, (3, 2)), label='queued'),
                   Line2D([], [], color=MUTED, ls=(0, (2, 2)), marker='o', ms=4.5,
                          mfc=BG, mec=MUTED, label='branch point')],
          loc='upper center', bbox_to_anchor=(0.5, -0.11), ncol=4, frameon=False,
          fontsize=9, labelcolor=MUTED, handlelength=1.9, columnspacing=2.2)
ax.set_title('The PICAL lineage: branch points and configuration changes',
             fontsize=13.5, color=INK, loc='left', pad=18)
fig.tight_layout()
fig.savefig('plots/pical_run_tree_gantt.png', dpi=150, bbox_inches='tight', facecolor=BG)
print('wrote plots/pical_run_tree_gantt.png')
