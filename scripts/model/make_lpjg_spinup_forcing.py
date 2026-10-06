#!/usr/bin/env python3
"""
Build the daily forcing file for an offline LPJ-GUESS spin-up from coupled AWI-ESM3 output.

LPJ-GUESS (framework.cpp, OIFS48r1spinup / getdailyIFSforcing_fast) reads eight daily
fields on the native grid, (time_counter, cell): tas, tasmin, tasmax, pr, rsns, rlns,
hurs, sfcWind. It reads calc_total_days(spin-up start year, nyearifs) records per cell
and fills year by year with its own calendar, so the file must carry the leap days of the
SPIN-UP calendar, not of the source years. Where the two differ, 29 February is dropped
or 28 February is repeated, so every forcing year starts on 1 January.

  tas      daily mean of hourly-mean 2t
  tasmin   atmos_day_minmax_tasmin
  tasmax   atmos_day_minmax_tasmax
  pr       daily mean of hourly-mean pr (kg m-2 s-1)
  rsns     atmos_day_rad_rss (net surface shortwave, W m-2)
  rlns     atmos_day_rad_rls (net surface longwave, W m-2, negative)
  hurs     daily mean of hourly 100 * es(2d) / es(2t), saturation over water
  sfcWind  daily mean of hourly instantaneous sqrt(10u^2 + 10v^2)

Storage is contiguous and uncompressed, as in the CRUNCEP file this replaces: the reader
pulls one cell's whole time series per call.
"""
import argparse
import calendar

import numpy as np
import netCDF4

RUNTIME = "/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4"

ATTRS = [
    ("rsns", "Surface Net Solar Radiation", "surface_net_downward_shortwave_flux", "W m-2"),
    ("rlns", "Surface Net Thermal Radiation", "surface_net_downward_longwave_flux", "W m-2"),
    ("pr", "Precipitation", "precipitation_flux", "kg m-2 s-1"),
    ("tas", "Near-Surface Air Temperature", "air_temperature", "K"),
    ("tasmin", "Daily Minimum Near-Surface Air Temperature", "air_temperature", "K"),
    ("tasmax", "Daily Maximum Near-Surface Air Temperature", "air_temperature", "K"),
    ("sfcWind", "Near-Surface Wind Speed", "wind_speed", "m s-1"),
    ("hurs", "Near-Surface Relative Humidity", "relative_humidity", "%"),
]


def esat(t):
    """Saturation vapour pressure over water (Pa), Magnus form used by IFS."""
    return 611.21 * np.exp(17.502 * (t - 273.16) / (t - 32.19))


def read(outdir, stream, var, year, nrec):
    with netCDF4.Dataset(f"{outdir}/{stream}_{year}-{year}.nc") as d:
        v = d[var]
        if v.shape != (nrec, 40320):
            raise SystemExit(f"{stream} {year}: shape {v.shape}, expected {(nrec, 40320)}")
        a = np.asarray(v[:], dtype=np.float32)
    if not np.isfinite(a).all():
        raise SystemExit(f"{stream} {year}: non-finite values")
    return a


def daily(a, nd):
    return a.reshape(nd, 24, -1).mean(axis=1, dtype=np.float64).astype(np.float32)


def year_fields(outdir, year):
    nd = 366 if calendar.isleap(year) else 365
    f = {}
    f["tasmin"] = read(outdir, "atmos_day_minmax_tasmin", "tasmin", year, nd)
    f["tasmax"] = read(outdir, "atmos_day_minmax_tasmax", "tasmax", year, nd)
    f["rsns"] = read(outdir, "atmos_day_rad_rss", "rss", year, nd)
    f["rlns"] = read(outdir, "atmos_day_rad_rls", "rls", year, nd)
    f["pr"] = daily(read(outdir, "atmos_1h_pr_pr", "pr", year, nd * 24), nd)
    t = read(outdir, "atmos_1h_sfc_2t", "2t", year, nd * 24)
    f["tas"] = daily(t, nd)
    td = read(outdir, "atmos_1h_sfc_2d", "2d", year, nd * 24)
    f["hurs"] = daily(np.clip(100.0 * esat(td) / esat(t), 0.0, 100.0), nd)
    del t, td
    u = read(outdir, "atmos_1h_pt_10u", "10u", year, nd * 24)
    v = read(outdir, "atmos_1h_pt_10v", "10v", year, nd * 24)
    f["sfcWind"] = daily(np.hypot(u, v), nd)
    return f


def to_calendar(a, src_leap, dst_leap):
    if src_leap == dst_leap:
        return a
    if src_leap:                       # drop 29 February
        return np.delete(a, 59, axis=0)
    return np.insert(a, 59, a[58], axis=0)   # repeat 28 February


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exp", default="PICAL_crunveg_tke_albsn082")
    p.add_argument("--first", type=int, default=2190)
    p.add_argument("--nyears", type=int, default=20)
    p.add_argument("--spinup-start", type=int, default=1900)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    outdir = f"{RUNTIME}/{a.exp}/outdata/oifs"
    ntime = sum(366 if calendar.isleap(a.spinup_start + i) else 365 for i in range(a.nyears))

    with netCDF4.Dataset(f"{outdir}/atmos_day_rad_rss_{a.first}-{a.first}.nc") as d:
        lat, lon = d["lat"][:], d["lon"][:]

    root = netCDF4.Dataset(a.out, "w", format="NETCDF4")
    root.createDimension("time_counter", ntime)
    root.createDimension("cell", 40320)
    tv = root.createVariable("time_counter", "f8", ("time_counter",))
    tv.units = f"days since {a.spinup_start}-01-01 12:00:00"
    tv.calendar = "proleptic_gregorian"
    tv.axis = "T"
    tv.standard_name = "time"
    tv[:] = np.arange(ntime)
    la = root.createVariable("lat", "f4", ("cell",))
    la.units = "degrees_north"
    la.standard_name = "latitude"
    la[:] = lat
    lo = root.createVariable("lon", "f4", ("cell",))
    lo.units = "degrees_east"
    lo.standard_name = "longitude"
    lo[:] = lon
    for name, long_name, std, units in ATTRS:
        v = root.createVariable(name, "f4", ("time_counter", "cell"), contiguous=True)
        v.long_name, v.standard_name, v.units = long_name, std, units
    root.title = (f"LPJ-GUESS spin-up forcing from {a.exp} {a.first}-{a.first + a.nyears - 1}, "
                  f"on the calendar of {a.spinup_start}-{a.spinup_start + a.nyears - 1}")
    root.source = f"{outdir}"
    root.history = ("make_lpjg_spinup_forcing.py; leap days follow the spin-up calendar "
                    "(29 Feb dropped or 28 Feb repeated where the source year differs)")

    t0 = 0
    for i in range(a.nyears):
        ys, yd = a.first + i, a.spinup_start + i
        f = year_fields(outdir, ys)
        n = 366 if calendar.isleap(yd) else 365
        for name in f:
            root[name][t0:t0 + n] = to_calendar(f[name], calendar.isleap(ys), calendar.isleap(yd))
        print(f"{ys} -> {yd}: {n} days, tas {f['tas'].mean():.3f} K, pr {f['pr'].mean() * 86400:.3f} mm/d, "
              f"hurs {f['hurs'].mean():.2f} %, wind {f['sfcWind'].mean():.3f} m/s", flush=True)
        t0 += n
    assert t0 == ntime
    root.close()


if __name__ == "__main__":
    main()
