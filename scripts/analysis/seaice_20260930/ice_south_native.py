"""March ice edge in the Bering/North Pacific and Labrador boxes on the native FESOM mesh.
Mesh-independent: southernmost node with a_ice >= 0.15. Coordinates from the file (XIOS output)
or, for older FESOM I/O, from nod2d.out picked by node count.
    python ice_south_native.py RUNDIR [STEP]      one line per sampled year:
    run year bering_southmost labrador_southmost
"""
import glob, sys
import numpy as np, xarray as xr
MESH = {126858: "/work/ab0246/a270092/input/fesom2/core2/nod2d.out",
        204875: "/work/ab0246/a270092/input/fesom2/core3/nod2d.out"}
run, step = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1
files = sorted(glob.glob(f"{run}/outdata/fesom/a_ice.fesom.[0-9][0-9][0-9][0-9].nc"))[::step]
lon = lat = None
for f in files:
    y = int(f[-7:-3])
    try:
        d = xr.open_dataset(f, decode_times=False)
        a = d["a_ice"]
        n = a.sizes["nod2"]
        if lon is None or lon.size != n:
            if "lon" in d.coords or "lon" in d:
                lon, lat = d["lon"].values.astype(float), d["lat"].values.astype(float)
            else:
                xy = np.loadtxt(MESH[n], skiprows=1, usecols=(1, 2)); lon, lat = xy[:, 0], xy[:, 1]
            lon = ((lon + 180) % 360) - 180
            ber = ((lon >= 163) | (lon < -157)) & (lat >= 40) & (lat < 66)
            lab = (lon >= -70) & (lon < -44) & (lat >= 40) & (lat < 66)
        nt = a.sizes["time"]
        if nt == 12: m = a.isel(time=2).values
        elif nt in (365, 366): m = a.isel(time=slice(59, 90)).mean("time").values
        else: print(run.split("/")[-1], y, "nt", nt); continue
        ice = m >= 0.15
        sb = lat[ice & ber].min() if (ice & ber).any() else np.nan
        sl = lat[ice & lab].min() if (ice & lab).any() else np.nan
        print(run.split("/")[-1], y, f"{sb:.2f}", f"{sl:.2f}", flush=True)
    except Exception as e:
        print(run.split("/")[-1], y, "ERR", str(e)[:80], flush=True)
