"""March state of the North Atlantic deep-convection regions, native FESOM mesh.
    python convection_boxes.py LABEL DATADIR Y0 Y1 [NOD2D.OUT]
Per box: mean a_ice, mean and max mixed-layer depth (MLD2), % of nodes deeper than 1000 m, mean SSS."""
import sys
import numpy as np, xarray as xr
lab, dd, y0, y1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
nod = sys.argv[5] if len(sys.argv) > 5 else None
BOX = {"Labrador": (56, 62, -60, -50), "Irminger": (58, 63, -40, -30), "Greenland Sea": (72, 77, -10, 5)}
def march(v, y):
    d = xr.open_dataset(f"{dd}/{v}.fesom.{y}.nc", decode_times=False)
    a = d[v]; nt = a.sizes["time"]
    m = a.isel(time=2).values if nt == 12 else a.isel(time=slice(59, 90)).mean("time").values
    if "lon" in d: return m, d["lon"].values.astype(float), d["lat"].values.astype(float)
    return m, None, None
acc = {v: [] for v in ("a_ice", "MLD2", "sss")}
lon = lat = None
for y in range(y0, y1 + 1):
    for v in acc:
        m, lo, la = march(v, y); acc[v].append(m)
        if lon is None and lo is not None: lon, lat = lo, la
if lon is None:
    xy = np.loadtxt(nod, skiprows=1, usecols=(1, 2)); lon, lat = xy[:, 0], xy[:, 1]
lon = ((lon + 180) % 360) - 180
f = {v: np.mean(acc[v], axis=0) for v in acc}
mld = np.abs(f["MLD2"])
for b, (a0, a1, o0, o1) in BOX.items():
    k = (lat >= a0) & (lat < a1) & (lon >= o0) & (lon < o1)
    print(f"{lab:26s} {b:14s} n={k.sum():5d}  ice {np.nanmean(f['a_ice'][k]):4.2f}  MLD mean {np.nanmean(mld[k]):6.0f} m  max {np.nanmax(mld[k]):6.0f} m"
          f"  >1000m {100*np.mean(mld[k] > 1000):4.0f} %  SSS {np.nanmean(f['sss'][k]):6.3f}")
