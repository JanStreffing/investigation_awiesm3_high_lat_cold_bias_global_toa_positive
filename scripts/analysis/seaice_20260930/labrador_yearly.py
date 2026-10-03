import sys, numpy as np, xarray as xr
N = "/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4"
lon = lat = None
for exp, y0, y1 in (("PICAL_crunveg", 2100, 2119), ("PICAL_crunveg_ihf0", 2120, 2139)):
    for y in range(y0, y1 + 1):
        out = {}
        for v in ("MLD2", "sss", "a_ice"):
            d = xr.open_dataset(f"{N}/{exp}/outdata/fesom/{v}.fesom.{y}.nc", decode_times=False)
            if lon is None:
                lon = ((d["lon"].values + 180) % 360) - 180; lat = d["lat"].values
                k = (lat >= 56) & (lat < 62) & (lon >= -60) & (lon < -50)
            a = d[v]; nt = a.sizes["time"]
            m = a.isel(time=2).values if nt == 12 else a.isel(time=slice(59, 90)).mean("time").values
            out[v] = m[k]
        mld = np.abs(out["MLD2"])
        print(f"{exp:20s} {y}  Labrador March: MLD {np.nanmean(mld):6.0f} m  >1000m {100*np.mean(mld > 1000):3.0f} %  SSS {np.nanmean(out['sss']):6.3f}  ice {np.nanmean(out['a_ice']):4.2f}", flush=True)
