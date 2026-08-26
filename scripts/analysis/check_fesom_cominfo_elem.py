#!/usr/bin/env python3
"""Validate the element communication tables (com_elem2D and com_elem2D_full)
of a FESOM2 dist_<n> partitioning, elementwise and in order.

Local element numbering per rank is the my_list order: owned (myDim), then the
eDim halo in com_elem2D receive order, then the eXDim complement in
com_elem2D_full receive order; com_global2local assigns eXDim slots in exactly
that iteration order, so the two agree by construction and rlist entries in the
file are already local slots.

slist entries are mapped through a table holding owned AND eDim-halo elements,
so a send may reference the sender's halo copy (value in myDim+1..myDim+eDim,
stale by one exchange) and an element the sender does not know at all maps to
0. Both are counted.

Usage: check_fesom_cominfo_elem.py <dist_dir> <npes>
"""
import sys, os
import numpy as np


def toks(path):
    return np.fromstring(open(path).read(), dtype=np.int64, sep=" ")


def parse_my_list(d, pe):
    v = toks(os.path.join(d, "my_list%05d.out" % pe))
    i = 0
    _pe = int(v[i]); i += 1
    mdn, edn = int(v[i]), int(v[i+1]); i += 2
    i += mdn + edn
    mde, ede, exde = int(v[i]), int(v[i+1]), int(v[i+2]); i += 3
    elems = v[i:i + mde + ede + exde]
    return mde, ede, exde, elems


def parse_com(d, pe):
    v = toks(os.path.join(d, "com_info%05d.out" % pe))
    i = 1
    out = []
    for blk in range(3):                       # nod2D, elem2D, elem2D_full
        rPEnum = int(v[i]); i += 1
        rPE = v[i:i+rPEnum]; i += rPEnum
        rptr = v[i:i+rPEnum+1]; i += rPEnum+1
        nr = int(rptr[-1]-1); rlist = v[i:i+nr]; i += nr
        sPEnum = int(v[i]); i += 1
        sPE = v[i:i+sPEnum]; i += sPEnum
        sptr = v[i:i+sPEnum+1]; i += sPEnum+1
        ns = int(sptr[-1]-1); slist = v[i:i+ns]; i += ns
        out.append(dict(rPE=rPE, rptr=rptr, rlist=rlist, sPE=sPE, sptr=sptr, slist=slist))
    return out[1], out[2]                      # elem2D, elem2D_full


def check(d, npes):
    dims = {}; loc = {}; ce = {}; cf = {}
    for pe in range(npes):
        mde, ede, exde, elems = parse_my_list(d, pe)
        dims[pe] = (mde, ede, exde); loc[pe] = elems
        ce[pe], cf[pe] = parse_com(d, pe)

    for tag, com in (("elem2D", ce), ("elem2D_full", cf)):
        pairs = mism = zero = halo_fwd = 0
        for a in range(npes):
            ca = com[a]; mde_a = dims[a][0]
            for k, b in enumerate(ca["rPE"]):
                b = int(b); pairs += 1
                seg = slice(int(ca["rptr"][k])-1, int(ca["rptr"][k+1])-1)
                rslots = ca["rlist"][seg]
                expect = loc[a][rslots - 1]
                cb = com[b]; mde_b = dims[b][0]; ede_b = dims[b][1]
                ks = np.where(cb["sPE"] == a)[0]
                if ks.size != 1:
                    print("%s PAIR MISSING %d->%d" % (tag, b, a)); mism += 1; continue
                ks = int(ks[0])
                sseg = cb["slist"][int(cb["sptr"][ks])-1:int(cb["sptr"][ks+1])-1]
                zero += int((sseg < 1).sum())
                halo_fwd += int((sseg > mde_b).sum())
                ok = sseg >= 1
                sent = np.where(ok, loc[b][np.clip(sseg, 1, None) - 1], -1)
                if sent.size != expect.size or not np.array_equal(sent, expect):
                    mism += 1
        print("%s %s: %d edges, %d mismatched, %d zero-index sends, %d halo-forwarded sends"
              % (os.path.basename(d.rstrip('/')), tag, pairs, mism, zero, halo_fwd))


if __name__ == "__main__":
    check(sys.argv[1], int(sys.argv[2]))
