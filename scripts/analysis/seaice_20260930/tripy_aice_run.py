"""NH sea-ice concentration map with tripyview (ESM_auto_tripy hslice_np a_ice settings).
    python tripy_aice_run.py DATADIR MESHDIR Y0 Y1 MONTH TITLE OUT.png
"""
import sys
import dask
import tripyview as tpv
dask.config.set(scheduler="synchronous")
data_path, mesh_path, y0, y1, mon, title, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7]
mesh = tpv.load_mesh_fesom2(mesh_path, do_rot="None", focus=0, do_info=False, do_pickle=True)
data = tpv.load_data_fesom2(mesh, data_path, vname="a_ice", year=[y0, y1], mon=[mon], descript=title, do_info=False,
                            chunks={"time": -1, "nod2": -1, "elem": -1, "edg_n": -1, "nz": -1, "nz1": -1, "ndens": -1})
assert data["a_ice"].sizes["nod2"] == mesh.n2dn, (data["a_ice"].sizes, mesh.n2dn)
data = data[["a_ice"]]      # XIOS output also carries bounds_lon/bounds_lat
tpv.plot_hslice(mesh, [data], cinfo={"cstr": "wbgyr", "crange": [0, 1.0, 0.5]},
                box=[-180, 180, 40, 90], proj="nps", do_save=out, save_dpi=400)
print("wrote", out)
