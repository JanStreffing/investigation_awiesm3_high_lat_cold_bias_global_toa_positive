#!/usr/bin/env python3
"""Check a FESOM2 dist_<n> partitioning for pathologies.

Written to diagnose why a locally generated dist_1792 for the corrected CORE3
mesh blew up at mstep=1 while the shipped dist_512 and dist_1024 ran clean.

Reports, per partition: owned-node count, halo count, and the number of
connected components of the owned-node subgraph. A partition split into
several disconnected pieces is the classic METIS failure mode that FESOM's
halo exchange does not tolerate.

Usage:
    python3 check_fesom_partition.py <mesh_dir> <dist_dir_name> [...]
"""
import sys
import os
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def read_elem2d(mesh_dir):
    path = os.path.join(mesh_dir, "elem2d.out")
    with open(path) as fh:
        nelem = int(fh.readline().split()[0])
        elems = np.loadtxt(fh, dtype=np.int64, max_rows=nelem)
    return elems - 1  # to 0-based


def node_adjacency(elems, nnodes):
    """Sparse node-node adjacency from triangle connectivity."""
    pairs = np.vstack([elems[:, [0, 1]], elems[:, [1, 2]], elems[:, [2, 0]]])
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    data = np.ones(rows.size, dtype=np.int8)
    return coo_matrix((data, (rows, cols)), shape=(nnodes, nnodes)).tocsr()


def read_my_list(path):
    """Return (owned, halo) global 0-based node indices for one rank."""
    vals = np.fromstring(open(path).read(), dtype=np.int64, sep=" ")
    my_dim, e_dim = int(vals[1]), int(vals[2])
    nodes = vals[3:3 + my_dim + e_dim] - 1
    return nodes[:my_dim], nodes[my_dim:]


def check(mesh_dir, dist_name):
    dist_dir = os.path.join(mesh_dir, dist_name)
    nnodes = int(open(os.path.join(mesh_dir, "nod2d.out")).readline().split()[0])
    adj = node_adjacency(read_elem2d(mesh_dir), nnodes)

    npes = int(open(os.path.join(dist_dir, "rpart.out")).readline().split()[0])
    broken, owned_sizes, halo_sizes = [], [], []
    for pe in range(npes):
        owned, halo = read_my_list(os.path.join(dist_dir, "my_list%05d.out" % pe))
        owned_sizes.append(owned.size)
        halo_sizes.append(halo.size)
        sub = adj[owned][:, owned]
        ncomp, _ = connected_components(sub, directed=False)
        if ncomp > 1:
            broken.append((pe, ncomp, owned.size))

    owned_sizes = np.array(owned_sizes)
    halo_sizes = np.array(halo_sizes)
    print("%s/%s" % (os.path.basename(mesh_dir.rstrip("/")), dist_name))
    print("  npes            %d" % npes)
    print("  owned per rank  min %d  max %d  mean %.1f  sum %d"
          % (owned_sizes.min(), owned_sizes.max(), owned_sizes.mean(), owned_sizes.sum()))
    print("  halo per rank   min %d  max %d  mean %.1f"
          % (halo_sizes.min(), halo_sizes.max(), halo_sizes.mean()))
    print("  disconnected    %d of %d ranks" % (len(broken), npes))
    for pe, ncomp, size in broken[:20]:
        print("      rank %5d: %d components over %d owned nodes" % (pe, ncomp, size))
    if len(broken) > 20:
        print("      ... %d more" % (len(broken) - 20))
    print()


if __name__ == "__main__":
    mesh = sys.argv[1]
    for name in sys.argv[2:]:
        check(mesh, name)
