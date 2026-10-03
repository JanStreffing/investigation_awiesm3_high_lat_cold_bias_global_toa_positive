"""PI200 March a_ice map: its years (<1678) break tripyview's time handling, so average March
from the raw files and hand the field to tripyview inside a data object from a same-mesh run."""
import dask, numpy as np, xarray as xr, tripyview as tpv
dask.config.set(scheduler="synchronous")
R = "/work/bb1469/a270092/runtime/awiesm3-v3.4"
Y0, Y1 = 1615, 1619
mesh = tpv.load_mesh_fesom2("/work/ab0246/a270092/input/fesom2/core3/", do_rot="None", focus=0, do_info=False, do_pickle=True)
tmpl = tpv.load_data_fesom2(mesh, f"{R}/PICAL_ccnice/outdata/fesom/", vname="a_ice", year=[2095, 2099], mon=[3],
                            descript="awiesm3-v3.4_PI200", do_info=False,
                            chunks={"time": -1, "nod2": -1, "elem": -1, "edg_n": -1, "nz": -1, "nz1": -1, "ndens": -1})[["a_ice"]]
fields = []
for y in range(Y0, Y1 + 1):
    d = xr.open_dataset(f"{R}/PI200/outdata/fesom/a_ice.fesom.{y}.nc", decode_times=False)
    if y == Y0:
        dl = np.abs(((d.lon.values - mesh.n_x + 180) % 360) - 180).max(); dp = np.abs(d.lat.values - mesh.n_y).max()
        print("node order vs mesh: max dlon", dl, "max dlat", dp); assert dl < 1e-2 and dp < 1e-2  # polar nodes differ by rotation roundoff (<0.006 deg)
    a = d["a_ice"]
    fields.append((a.isel(time=2) if a.sizes["time"] == 12 else a.isel(time=slice(59, 90)).mean("time")).values)
tmpl["a_ice"].values[:] = np.mean(fields, axis=0)
tmpl["a_ice"].attrs.update(year=f"[{Y0}, {Y1}]", str_ltim=f"y:{Y0}-{Y1}, m:Mar.", str_lsave=f"y{Y0}-{Y1}_mMar.",
                           datapath=f"{R}/PI200/outdata/fesom/")
out = f"ice_maps/awiesm3-v3.4_PI200_a_ice_nps_{Y0}-{Y1}_Mar.png"
tpv.plot_hslice(mesh, [tmpl], cinfo={"cstr": "wbgyr", "crange": [0, 1.0, 0.5]}, box=[-180, 180, 40, 90], proj="nps",
                do_save=out, save_dpi=400)
print("wrote", out)
