"""March sea-ice extent in the marginal seas, per year, for every levante AWI-ESM3 run.
Extent = area of 0.5-degree cells with a_ice >= 0.15 (a_ice.fesom.gr.YYYY.nc).
Prints one line per run and year: bering, labrador, okhotsk, south-of-50N, NH total (1e6 km2),
and the southernmost ice latitude in the Bering box and in the Pacific 150E-130W.
    python ice_south_levante.py RUN [RUN ...]    (years taken from what exists)
"""
import glob, os, sys
import numpy as np, xarray as xr

R = "/work/bb1469/a270092/runtime/awiesm3-v3.4"
RE = 6.371e6
BOX = {"bering": (50, 66, 163, -157), "labrador": (40, 66, -70, -44), "okhotsk": (42, 62, 135, 163),
       "s50": (0, 50, -180, 180), "pac": (30, 66, 150, -130)}
geo = None
def setup(d):
    lat = d["lat"].values; lon = ((d["lon"].values + 180) % 360) - 180
    A = (RE**2 * np.deg2rad(0.5)**2 * np.cos(np.deg2rad(lat)))[:, None] * np.ones(lon.size)[None, :]
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    m = {}
    for k, (a0, a1, o0, o1) in BOX.items():
        lo = (LO >= o0) & (LO < o1) if o0 <= o1 else (LO >= o0) | (LO < o1)
        m[k] = (LA >= a0) & (LA < a1) & lo
    return A, LA, m
for run in sys.argv[1:]:
    for f in sorted(glob.glob(f"{R}/{run}/outdata/fesom/a_ice.fesom.gr.[0-9][0-9][0-9][0-9].nc")):
        y = int(f[-7:-3])
        try:
            d = xr.open_dataset(f)
            v = [k for k in d.data_vars if k.startswith("a_ice")][0]
            a = d[v].sel(time=d.time.dt.month == 3).mean("time").values
        except Exception as e:
            print(run, y, "ERR", str(e)[:60], flush=True); continue
        if geo is None: geo = setup(d)
        A, LA, m = geo
        ice = (a >= 0.15) & np.isfinite(a)
        ext = {k: (A * (ice & m[k])).sum() / 1e12 for k in m}
        nh = (A * (ice & (LA > 0))).sum() / 1e12
        sb = LA[ice & m["bering"]].min() if (ice & m["bering"]).any() else np.nan
        sp = LA[ice & m["pac"]].min() if (ice & m["pac"]).any() else np.nan
        print(f"{run} {y} {ext['bering']:.3f} {ext['labrador']:.3f} {ext['okhotsk']:.3f} {ext['s50']:.3f} {nh:.3f} {sb:.2f} {sp:.2f}", flush=True)
