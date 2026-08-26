#!/usr/bin/env python3
"""Validate the node communication tables of a FESOM2 dist_<n> partitioning.

For every receive edge B -> A the check is exact, elementwise and in order:
the global IDs A expects in its halo segment for B (positions rptr[B] of the
halo section of A's my_list, since com_global2local makes rlist sequential)
must equal the global IDs of the owned nodes B sends (slist local indices
mapped through B's owned list).

Also checks that every slist entry is a valid owned-node index of the sender:
com_global2local maps send-list globals through a table built only over owned
nodes, so a node the sender does not own becomes local index 0 and would make
the runtime exchange read arr(0).

Usage: check_fesom_cominfo.py <dist_dir> <npes>
"""
import sys, os
import numpy as np


def toks(path):
    return np.fromstring(open(path).read(), dtype=np.int64, sep=" ")


def parse_my_list(d, pe):
    v = toks(os.path.join(d, "my_list%05d.out" % pe))
    md, ed = int(v[1]), int(v[2])
    nodes = v[3:3 + md + ed]
    return nodes[:md], nodes[md:]          # owned, halo (both global, halo in receive order)


def parse_com_nod(d, pe):
    v = toks(os.path.join(d, "com_info%05d.out" % pe))
    i = 1
    rPEnum = int(v[i]); i += 1
    rPE = v[i:i + rPEnum]; i += rPEnum
    rptr = v[i:i + rPEnum + 1]; i += rPEnum + 1
    nr = int(rptr[-1] - 1)
    rlist = v[i:i + nr]; i += nr
    sPEnum = int(v[i]); i += 1
    sPE = v[i:i + sPEnum]; i += sPEnum
    sptr = v[i:i + sPEnum + 1]; i += sPEnum + 1
    ns = int(sptr[-1] - 1)
    slist = v[i:i + ns]; i += ns
    return dict(rPE=rPE, rptr=rptr, rlist=rlist, sPE=sPE, sptr=sptr, slist=slist)


def main(d, npes):
    owned = {}; halo = {}; com = {}
    for pe in range(npes):
        owned[pe], halo[pe] = parse_my_list(d, pe)
        com[pe] = parse_com_nod(d, pe)

    bad_slist = bad_pairs = pairs = 0
    for a in range(npes):
        ca = com[a]
        for k, b in enumerate(ca["rPE"]):
            b = int(b)
            seg = slice(int(ca["rptr"][k]) - 1, int(ca["rptr"][k + 1]) - 1)
            expect = halo[a][seg]                      # globals A expects, in order
            cb = com[b]
            ks = np.where(cb["sPE"] == a)[0]
            pairs += 1
            if ks.size != 1:
                print("PAIR MISSING: rank %d receives from %d but %d has no matching send" % (a, b, b))
                bad_pairs += 1
                continue
            ks = int(ks[0])
            sseg = cb["slist"][int(cb["sptr"][ks]) - 1:int(cb["sptr"][ks + 1]) - 1]
            if sseg.min() < 1 or sseg.max() > owned[b].size:
                print("BAD SLIST: rank %d -> %d, local index out of owned range [1,%d]: min %d max %d"
                      % (b, a, owned[b].size, sseg.min(), sseg.max()))
                bad_slist += 1
                continue
            sent = owned[b][sseg - 1]                  # globals B sends, in order
            if sent.size != expect.size or not np.array_equal(sent, expect):
                nbad = int((sent != expect).sum()) if sent.size == expect.size else -1
                print("MISMATCH: %d -> %d, len %d vs %d, differing entries %s"
                      % (b, a, sent.size, expect.size, nbad))
                bad_pairs += 1

    print("%s: %d receive edges checked, %d slist violations, %d pair mismatches"
          % (os.path.basename(d.rstrip('/')), pairs, bad_slist, bad_pairs))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
