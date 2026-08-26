"""Assert that a mesh and a model output file describe the same mesh.

WHY THIS EXISTS.  On 2026-08-26 the shipped CORE3 mesh was replaced in place at
/work/ab0246/a270092/input/fesom2/core3 and the previous one renamed core3_beta.
The two have DIFFERENT node counts, 220509 against 211567, and the cavity depths
differ.  Every runscript, every reval config and every analysis script in this
campaign referred to the bare path `core3`, so after the swap they all silently
pointed at a mesh none of the completed runs had used.

Node count is the cheap discriminator and it is checked here.  A mismatch raises
instead of broadcasting: numpy will happily area-weight a 211567-node field with
the first 211567 entries of a 220509-node area array and return a number that
looks entirely reasonable.

RULE.  Runs completed before 2026-08-26 are core3_beta.  Runs started after are
core3.  Never mix them.
"""
import os
import xarray as xr

BETA = '/work/ab0246/a270092/input/fesom2/core3_beta'
NEW = '/work/ab0246/a270092/input/fesom2/core3'
N_BETA, N_NEW = 211567, 220509


def nodes_of(meshpath):
    """Node count of a mesh directory.

    Prefers fesom.mesh.diag.nc, but the freshly shipped core3 does not carry one
    yet, so fall back to the first line of nod2d.out, which every FESOM mesh has.
    """
    diag = f'{meshpath}/fesom.mesh.diag.nc'
    if os.path.exists(diag):
        with xr.open_dataset(diag, decode_times=False) as d:
            return d.sizes.get('nod2')
    with open(f'{meshpath}/nod2d.out') as fh:
        return int(fh.readline().split()[0])


def check(meshpath, nod2_from_output, context=''):
    """Raise unless the mesh at meshpath has nod2_from_output nodes."""
    n = nodes_of(meshpath)
    if n != nod2_from_output:
        want = 'core3_beta' if nod2_from_output == N_BETA else \
               'core3' if nod2_from_output == N_NEW else 'an unknown mesh'
        raise SystemExit(
            f'MESH MISMATCH{" in " + context if context else ""}: '
            f'{meshpath} has {n} nodes but the model output has '
            f'{nod2_from_output}. That output belongs to {want}.')
    return n


def mesh_for(nod2_from_output):
    """The mesh directory matching an output file's node count."""
    if nod2_from_output == N_BETA:
        return BETA
    if nod2_from_output == N_NEW:
        return NEW
    raise SystemExit(f'unknown mesh: output has {nod2_from_output} nodes')
