#!/usr/bin/env python3
"""Level-range bars at the nodes that trigger FESOM/fesom2#960.

For each node: the node's own water column (draft to bottom), then one bar
per adjacent element spanning that element's draft to bottom. The GM solve
takes nzmin = max element draft and nzmax = min element bottom; where the
magenta line lies below the cyan one there is no interior range.

Usage: plot_gm_range_bars.py <mesh_dir> <out_png> <nodes 1-based...>
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mesh, out = sys.argv[1], sys.argv[2]
seeds = [int(s) - 1 for s in sys.argv[3:]]
zbar = np.loadtxt(os.path.join(mesh, "levels3d.txt"), skiprows=1)
ne = int(open(os.path.join(mesh, "elem2d.out")).readline().split()[0])
el = np.loadtxt(os.path.join(mesh, "elem2d.out"), skiprows=1, dtype=np.int64, max_rows=ne) - 1
cav = np.loadtxt(os.path.join(mesh, "cavity_nlvls.out"), dtype=np.int64)
nlv = np.loadtxt(os.path.join(mesh, "nlvls.out"), dtype=np.int64)
ecav = np.loadtxt(os.path.join(mesh, "cavity_elvls.out"), dtype=np.int64)
enlv = np.loadtxt(os.path.join(mesh, "elvls.out"), dtype=np.int64)

fig, axes = plt.subplots(1, len(seeds), figsize=(3.2 * len(seeds), 7), sharey=True)
axes = np.atleast_1d(axes)
for ax, n in zip(axes, seeds):
    adj = np.where((el == n).any(axis=1))[0]
    adj = adj[np.argsort(ecav[adj])]
    nzmin, nzmax = ecav[adj].max(), enlv[adj].min()
    # node column
    ax.bar(0, zbar[nlv[n] - 1] - zbar[cav[n] - 1], bottom=zbar[cav[n] - 1], width=0.7,
           color="red", alpha=0.6, label="node column")
    ax.text(0, zbar[cav[n] - 1] + 6, "node\n%d..%d" % (cav[n], nlv[n]), ha="center", fontsize=7)
    # element bars
    for k, e in enumerate(adj, start=1):
        ax.bar(k, zbar[enlv[e] - 1] - zbar[ecav[e] - 1], bottom=zbar[ecav[e] - 1], width=0.7,
               color=plt.get_cmap("viridis")((ecav[e] - 1) / 16.0), alpha=0.8, edgecolor="k")
        ax.text(k, zbar[ecav[e] - 1] + 6, "%d..%d" % (ecav[e], enlv[e]), ha="center", fontsize=7)
    ax.axhline(zbar[nzmin - 1], color="magenta", lw=2.5, label="nzmin = max elem draft (%d)" % nzmin)
    ax.axhline(zbar[nzmax - 1], color="cyan", lw=2.5, label="nzmax = min elem bottom (%d)" % nzmax)
    inv = nzmin >= nzmax
    ax.set_title("node %d\n%s" % (n + 1, "INVERTED: no interior range" if inv else "ok (%d levels)" % (nzmax - nzmin - 1)),
                 color="magenta" if inv else "green", fontsize=9)
    ax.set_xticks(range(len(adj) + 1)); ax.set_xticklabels(["node"] + ["e%d" % k for k in range(1, len(adj) + 1)], fontsize=7)
    ax.legend(fontsize=6, loc="lower right")
    ax.grid(axis="y", alpha=0.3)
axes[0].set_ylabel("depth (m); bars span draft level .. bottom level")
plt.suptitle("GM solve range at the crash nodes: nzmin = max over adjacent elements of the element draft level, nzmax = min of the element bottom level", fontsize=9)
plt.tight_layout()
plt.savefig(out, dpi=140); print("wrote", out)
