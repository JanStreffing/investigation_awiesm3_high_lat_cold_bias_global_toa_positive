"""Arctic sea-ice concentration with tripyview, as ESM_auto_tripy.yaml's hslice_np a_ice.

    python tripy_aice_np.py EXPID Y0 Y1 MONTH
"""
import os
import sys

import dask
import tripyview as tpv

# serial reads: the XIOS output is Blosc-compressed (HDF5 filter 32001) and needs
# HDF5_PLUGIN_PATH; threaded netCDF4 opens crash in this environment
dask.config.set(scheduler="synchronous")

expid, y0, y1, mon = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
base = "/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4"
mesh_path = "/albedo/work/projects/p_awiesm3_cmip7/jstreffi/input/fesom2/core3/"
data_path = f"{base}/{expid}/outdata/fesom/"
save_dir = f"{base}/{expid}/viz/fesom"
os.makedirs(save_dir, exist_ok=True)

mesh = tpv.load_mesh_fesom2(mesh_path, do_rot="None", focus=0, do_info=False, do_pickle=True)
data = tpv.load_data_fesom2(mesh, data_path, vname="a_ice", year=[y0, y1], mon=[mon],
                            descript=expid, do_info=False,
                            # one chunk per file: automatic chunks differ between variables
                            chunks={"time": -1, "nod2": -1, "elem": -1, "edg_n": -1, "nz": -1, "nz1": -1, "ndens": -1})
# XIOS also writes bounds_lon/bounds_lat; tripyview plots the first data variable
data = data[["a_ice"]]
mname = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][mon - 1]
out = f"{save_dir}/a_ice_nps_{y0}-{y1}_{mname}.png"
tpv.plot_hslice(mesh, [data], cinfo={"cstr": "wbgyr", "crange": [0, 1.0, 0.5]},
                box=[-180, 180, 40, 90], proj="nps", do_save=out, save_dpi=400)
print("wrote", out)
