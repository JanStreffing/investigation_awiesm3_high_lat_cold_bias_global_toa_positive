"""Labrador interior convection in the build x time-step test, scored on persistence.

Single collapsed winters happen in the old configuration too (PICAL_ccnice, old build
1200 s: 2108, 2112, 2117, 2121, 2126, 2129, each recovering within 1-2 winters), so a weak
winter is not the signal.  What separates the new-build 1800 s runs is that they never
recover (ihf0 from 2128, ihf1 from 2127).  Per arm, March of 2121-2129 (2120 is the
OpenIFS cold start): MLD2 mean, % nodes > 1000 m, SSS, ice; then the count of collapsed
winters (MLD < 200 m), the longest consecutive run of them, and whether the arm ends
convecting (last two winters > 500 m).

Usage: python labrador_2x2.py [y0 y1]
"""
import sys, numpy as np, xarray as xr
N = "/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4"
Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (2121, 2129)
ARMS = [("PICAL_crunveg_ob1200", "old build 1200 s"), ("PICAL_crunveg_ts1200", "new build 1200 s"),
        ("PICAL_crunveg_ob1800", "old build 1800 s"), ("PICAL_crunveg_ihf0", "new build 1800 s"),
        ("PICAL_crunveg_ihf1", "new 1800 s + McPhee")]
lon = lat = k = None
summary = []
for exp, lab in ARMS:
    mlds = []
    for y in range(Y0, Y1 + 1):
        out = {}
        try:
            for v in ("MLD2", "sss", "a_ice"):
                d = xr.open_dataset(f"{N}/{exp}/outdata/fesom/{v}.fesom.{y}.nc", decode_times=False)
                if lon is None:
                    lon = ((d["lon"].values + 180) % 360) - 180; lat = d["lat"].values
                    k = (lat >= 56) & (lat < 62) & (lon >= -60) & (lon < -50)
                a = d[v]; nt = a.sizes["time"]
                m = a.isel(time=2).values if nt == 12 else a.isel(time=slice(59, 90)).mean("time").values
                out[v] = m[k]
        except (OSError, KeyError):
            break
        mld = np.abs(out["MLD2"]); mlds.append(np.nanmean(mld))
        print(f"{exp:22s} {y}  MLD {mlds[-1]:6.0f} m  >1000m {100*np.mean(mld > 1000):3.0f} %  "
              f"SSS {np.nanmean(out['sss']):6.3f}  ice {np.nanmean(out['a_ice']):4.2f}", flush=True)
    if not mlds:
        summary.append(f"{lab:22s} no output yet"); continue
    c = [m < 200 for m in mlds]
    run = best = 0
    for x in c:
        run = run + 1 if x else 0; best = max(best, run)
    ends = len(mlds) >= 2 and all(m > 500 for m in mlds[-2:])
    summary.append(f"{lab:22s} years {len(mlds)}  mean MLD {np.mean(mlds):5.0f} m  collapsed winters {sum(c)}  "
                   f"longest run {best}  ends convecting: {'yes' if ends else 'no'}")
print("\n" + "\n".join(summary))
