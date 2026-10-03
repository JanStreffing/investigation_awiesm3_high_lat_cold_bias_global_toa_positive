import glob, numpy as np, xarray as xr
N = "/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4"
lon = lat = k = None
for exp, lab in (("TEST_perfE", "1200 s, old build"), ("TEST_perfG", "1200 s, old build"),
                 ("TEST_perfF", "1800 s, old build"), ("TEST_perfH", "1800 s, old build"),
                 ("TEST_perfI", "1800 s, new build"), ("TEST_nolockB1", "1800 s, new build + stack"),
                 ("TEST_nolockB4", "1800 s, new build + stack")):
    out = {}
    try:
        for v in ("MLD2", "sss", "a_ice"):
            d = xr.open_dataset(f"{N}/{exp}/outdata/fesom/{v}.fesom.2100.nc", decode_times=False)
            if k is None:
                lon = ((d["lon"].values + 180) % 360) - 180; lat = d["lat"].values
                k = (lat >= 56) & (lat < 62) & (lon >= -60) & (lon < -50)
            a = d[v]; nt = a.sizes["time"]
            out[v] = (a.isel(time=slice(nt - 10, nt)).mean("time").values if nt > 12 else a.isel(time=nt - 1).values)[k]
        mld = np.abs(out["MLD2"])
        print(f"{exp:14s} {lab:26s} Labrador, last 10 days of Feb: MLD {np.nanmean(mld):5.0f} m  >1000m {100*np.mean(mld > 1000):3.0f} %  SSS {np.nanmean(out['sss']):6.3f}  ice {np.nanmean(out['a_ice']):4.2f}")
    except Exception as e:
        print(exp, "ERR", str(e)[:80])
