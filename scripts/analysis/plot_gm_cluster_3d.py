#!/usr/bin/env python3
"""3-D view of the vertical structure at the cavity-front cluster that triggers
the uninitialised GM streamfunction (FESOM/fesom2#960) on the corrected CORE3
mesh.

Nodes are drawn as vertical water columns from their ice draft
(cavity_nlvls) to their bottom (nlvls). Elements are translucent triangular
prisms from their own draft (cavity_elvls) to their bottom (elvls). At the
inverted node the GM solve uses nzmin = max over adjacent elements of the
element draft level and nzmax = min over adjacent elements of the element
bottom level; both are drawn as horizontal bars, and where nzmin sits below
nzmax the solve has no interior range.

Usage: plot_gm_cluster_3d.py <mesh_dir> <out_png> <seed nodes, 1-based...>
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

mesh, out = sys.argv[1], sys.argv[2]
seeds = np.array([int(s) for s in sys.argv[3:]]) - 1

zbar = np.loadtxt(os.path.join(mesh, "levels3d.txt"), skiprows=1)
nod = np.loadtxt(os.path.join(mesh, "nod2d.out"), skiprows=1)
ne = int(open(os.path.join(mesh, "elem2d.out")).readline().split()[0])
el = np.loadtxt(os.path.join(mesh, "elem2d.out"), skiprows=1, dtype=np.int64, max_rows=ne) - 1
cav = np.loadtxt(os.path.join(mesh, "cavity_nlvls.out"), dtype=np.int64)
nlv = np.loadtxt(os.path.join(mesh, "nlvls.out"), dtype=np.int64)
ecav = np.loadtxt(os.path.join(mesh, "cavity_elvls.out"), dtype=np.int64)
enlv = np.loadtxt(os.path.join(mesh, "elvls.out"), dtype=np.int64)

elems = np.where(np.isin(el, seeds).any(axis=1))[0]
nodes = np.unique(el[elems])
lon0, lat0 = nod[nodes, 1].mean(), nod[nodes, 2].mean()
kmx = lambda lo: (lo - lon0) * 111.32 * np.cos(np.radians(lat0))
kmy = lambda la: (la - lat0) * 111.32

fig = plt.figure(figsize=(13, 9))
ax = fig.add_subplot(111, projection="3d")

# element prisms: draft to bottom
cmap = plt.get_cmap("viridis")
for e in elems:
    p = el[e]
    x, y = kmx(nod[p, 1]), kmy(nod[p, 2])
    zt, zb = zbar[ecav[e] - 1], zbar[enlv[e] - 1]
    col = cmap((ecav[e] - 1) / 16.0)
    top = [list(zip(x, y, [zt] * 3))]
    bot = [list(zip(x, y, [zb] * 3))]
    sides = [[(x[i], y[i], zt), (x[j], y[j], zt), (x[j], y[j], zb), (x[i], y[i], zb)]
             for i, j in ((0, 1), (1, 2), (2, 0))]
    ax.add_collection3d(Poly3DCollection(top, facecolor=col, edgecolor="k", lw=0.6, alpha=0.55))
    ax.add_collection3d(Poly3DCollection(sides, facecolor=col, edgecolor="k", lw=0.3, alpha=0.18))
    ax.add_collection3d(Poly3DCollection(bot, facecolor="0.4", edgecolor="k", lw=0.4, alpha=0.35))

# node columns: draft to bottom
for n in nodes:
    x, y = kmx(nod[n, 1]), kmy(nod[n, 2])
    zt, zb = zbar[cav[n] - 1], zbar[nlv[n] - 1]
    is_seed = n in seeds
    ax.plot([x, x], [y, y], [zt, zb], color="red" if is_seed else "k", lw=3.5 if is_seed else 1.5)
    ax.scatter([x], [y], [zt], color="red" if is_seed else "k", s=30 if is_seed else 12)
    ax.text(x, y, zt + 4, str(n + 1), fontsize=7, color="red" if is_seed else "k")

# GM solve range at each seed node
for n in seeds:
    adj = np.where((el == n).any(axis=1))[0]
    nzmin, nzmax = ecav[adj].max(), enlv[adj].min()
    x, y = kmx(nod[n, 1]), kmy(nod[n, 2])
    for lvl, lab, c in ((nzmin, "nzmin", "magenta"), (nzmax, "nzmax", "cyan")):
        z = zbar[lvl - 1]
        ax.plot([x - 3, x + 3], [y, y], [z, z], color=c, lw=3)
        ax.text(x + 3.5, y, z, "%s=%d" % (lab, lvl), fontsize=7, color=c)
    inv = "INVERTED" if nzmin >= nzmax else "ok"
    ax.text(x, y, zbar[nlv[n] - 1] - 15, "%d: %s" % (n + 1, inv), fontsize=8,
            color="magenta" if inv == "INVERTED" else "green", weight="bold")

ax.set_xlabel("km east"); ax.set_ylabel("km north"); ax.set_zlabel("depth (m)")
ax.set_zlim(zbar[max(enlv[elems].max(), nlv[nodes].max()) - 1] - 20, 5)
ax.view_init(elev=22, azim=-50)
ax.set_title("Cavity-front cluster, %d nodes / %d elements: node columns (draft to bottom), element prisms coloured by element draft level;\n"
             "magenta/cyan bars are the GM solve range nzmin = max elem draft, nzmax = min elem bottom at each red node" % (nodes.size, elems.size), fontsize=9)
plt.tight_layout()
plt.savefig(out, dpi=140)
print("wrote", out)
