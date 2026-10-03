"""Where is the extra winter sea ice? Regional NH ice extent from the 0.5-degree FESOM output.

Part 1: March extent by basin, decade means, crunveg 2110s / ihf0 2120s / ihf0 2130s.
Part 2: late-February extent in the 2-month tests that separate the FESOM time step from
        the FESOM build (all from the same 2100-01-01 restart).
Extent = area of 0.5-degree cells with concentration >= 0.15.
"""
import glob

import numpy as np
import xarray as xr

RT = "/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4"
R_E = 6.371e6

REG = {  # name: (lat0, lat1, lon0, lon1), lon in -180..180; lon0 > lon1 wraps the dateline
    "Okhotsk": (42, 62, 135, 163),
    "Bering": (50, 66, 163, -157),
    "Labrador+Newfoundland": (40, 66, -70, -44),
    "GIN seas": (60, 80, -44, 20),
    "Barents+Kara": (66, 82, 20, 90),
    "Hudson Bay": (50, 66, -96, -70),
    "Baltic": (53, 66, 9, 31),
    "south of 50N (all)": (0, 50, -180, 180),
}


def cell_area(lat, lon):
    dlat = np.deg2rad(abs(float(lat[1] - lat[0])))
    dlon = np.deg2rad(abs(float(lon[1] - lon[0])))
    return (R_E**2 * dlat * dlon * np.cos(np.deg2rad(lat)))[:, None] * np.ones(lon.size)[None, :]


def mask(name, lat, lon):
    la0, la1, lo0, lo1 = REG[name]
    LA, LO = np.meshgrid(lat, ((np.asarray(lon) + 180.0) % 360.0) - 180.0, indexing="ij")  # file lon is 0..360
    inlat = (LA >= la0) & (LA < la1)
    inlon = (LO >= lo0) & (LO < lo1) if lo0 <= lo1 else (LO >= lo0) | (LO < lo1)
    return inlat & inlon


def month_mean(exp, years, month):
    fields = []
    for y in years:
        f = f"{RT}/{exp}/outdata/fesom/a_ice.fesom.gr.{y}.nc"
        d = xr.open_dataset(f)
        v = [k for k in d.data_vars if k.startswith("a_ice")][0]
        fields.append(d[v].sel(time=d.time.dt.month == month).mean("time").values)
        lat = d[[c for c in d.coords if "lat" in c][0]].values
        lon = d[[c for c in d.coords if "lon" in c][0]].values
    return np.mean(fields, axis=0), lat, lon


def regions(a, lat, lon):
    A = cell_area(lat, lon)
    ice = (a >= 0.15) & np.isfinite(a)
    nh = lat[:, None] > 0
    out = {"NH total": float((A * (ice & nh)).sum()) / 1e12}
    for r in REG:
        out[r] = float((A * (ice & mask(r, lat, lon))).sum()) / 1e12
    # southernmost ice-covered latitude per basin
    for r in ("Okhotsk", "Bering", "Labrador+Newfoundland"):
        m = ice & mask(r, lat, lon)
        out[r + " southmost lat"] = float(np.broadcast_to(lat[:, None], m.shape)[m].min()) if m.any() else float("nan")
    return out


print("PART 1: March, decade means, extent in 1e6 km2")
rows = [("crunveg 2111-19 (leadclose, 1200 s, albsn .80)", "PICAL_crunveg", range(2111, 2120)),
        ("ihf0 2121-29 (main+#1060, 1800 s, albsn .80)", "PICAL_crunveg_ihf0", range(2121, 2130)),
        ("ihf0 2131-39 (main+#1060, 1800 s, albsn .82)", "PICAL_crunveg_ihf0", range(2131, 2140))]
res = {}
for lab, exp, yrs in rows:
    a, lat, lon = month_mean(exp, yrs, 3)
    res[lab] = regions(a, lat, lon)
keys = list(next(iter(res.values())).keys())
print(f"{'':30s}" + "".join(f"{lab.split(' (')[0]:>18s}" for lab, _, _ in rows))
for k in keys:
    print(f"{k:30s}" + "".join(f"{res[lab][k]:18.2f}" for lab, _, _ in rows))

print("\nPART 2: 2-month tests from the same 2100-01-01 restart, extent over Feb 19-28, 1e6 km2")
for exp, lab in (("TEST_perfE", "1200 s, leadclose"), ("TEST_perfG", "1200 s, leadclose"),
                 ("TEST_perfF", "1800 s, leadclose"), ("TEST_perfH", "1800 s, leadclose"),
                 ("TEST_perfI", "1800 s, main+#1060"), ("TEST_nolockB1", "1800 s, main+#1060+stack")):
    fs = sorted(glob.glob(f"{RT}/{exp}/outdata/fesom/a_ice.fesom.gr.2100*.nc"))
    if not fs:
        print(f"{exp:14s} {lab:28s} no regridded a_ice"); continue
    d = xr.open_dataset(fs[0])
    v = [k for k in d.data_vars if k.startswith("a_ice")][0]
    t = d.time.values
    sel = (d.time.dt.month == 2) & (d.time.dt.day >= 19)
    if int(sel.sum()) == 0:
        print(f"{exp:14s} {lab:28s} no late-Feb records ({len(t)} records)"); continue
    a = d[v].sel(time=sel).mean("time").values
    lat = d[[c for c in d.coords if "lat" in c][0]].values
    lon = d[[c for c in d.coords if "lon" in c][0]].values
    r = regions(a, lat, lon)
    print(f"{exp:14s} {lab:28s} NH {r['NH total']:6.2f}  Okhotsk {r['Okhotsk']:5.2f}  Bering {r['Bering']:5.2f}  "
          f"Labrador {r['Labrador+Newfoundland']:5.2f}  <50N {r['south of 50N (all)']:5.2f}")
