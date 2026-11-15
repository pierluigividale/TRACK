"""Convert one month of EERIE (dz876/ec100) raw pp to the monthly netCDF
files TRACK's pipeline needs. Usage: pp_to_nc_iris_WCSI_EERIE.py RUN YEAR MON

Raw files are read from their original location or, once pp files have been
removed from disk, from a copy staged out of NLDS beforehand by
`tc_workflow.py stage` (run on a sci node - NLDS rejects calls from LOTUS
compute nodes, see eerie_nlds.py). This script never talks to NLDS itself.
"""
import iris
import os
import sys
import time

import eerie_nlds as E

name = sys.argv[1]
year = int(sys.argv[2])
mon_num = int(sys.argv[3])
if mon_num not in E.MONTH_ABBR:
    raise ValueError(f"month must be 1-12, got {mon_num}")

os.makedirs(E.dirout(name), exist_ok=True)
out = E.outputs(name, year, mon_num)

# Skip-if-already-done (ported 2026-11-12 from pp_to_nc_iris_WCSI.py's
# 2026-11-07 NH/SH-overlap fix, made per-variable/per-stream). NH and SH
# windows overlap by 6 months, so the other hemisphere's job (or a retry of
# this one) may already have converted this month. SLURMscript_n1280_para_WCSI
# holds a per-month lock around this script + ll_to_gg for EERIE runs, so a
# concurrent job for the same month waits and then lands here.
need = E.needed_vars(name, year, mon_num)
streams = E.needed_streams(name, year, mon_num)
if not streams:
    print(f"Converted output for {name} {year}{mon_num:02d} already exists - skipping pp conversion.")
    sys.exit(0)
print("still to convert: " + ", ".join(k for k, v in need.items() if v))

src = {}
for stream, original in streams.items():
    p = E.locate(name, year, mon_num, original)
    if p is None:
        raise FileNotFoundError(
            f"{original} is neither on disk nor staged+verified under "
            f"{E.staging_month_dir(name, year, mon_num)} - run on a sci node:\n"
            f"  python tc_workflow.py stage --experiments {name} --years <year> --yes --wait")
    print(f"{stream}: reading {p}")
    src[stream] = p


def save_atomic(cube, final_path):
    # write under a temporary name and rename only once complete, so a job
    # killed mid-save (walltime, quota) never leaves a truncated file under
    # the final name that the skip check above would then trust
    tmp = final_path + ".part.nc"
    iris.save(cube, tmp)
    os.replace(tmp, final_path)


start = time.time()
print("start ", start)

levelsw = iris.Constraint(pressure=[250, 500, 600, 700, 850, 925])
if need["u"]:
    save_atomic(iris.load_cube(src["p8"], iris.Constraint(name='x_wind') & levelsw), out["u"])
if need["v"]:
    save_atomic(iris.load_cube(src["p8"], iris.Constraint(name='y_wind') & levelsw), out["v"])
if need["geo"]:
    # EERIE Hart levels are 925/600/250 hPa (not 900/600/300) - see
    # track_n1280_para_WCSI.py's EERIE_EXPERIMENTS dispatch
    levelsg = iris.Constraint(pressure=[250, 600, 925])
    zgcon = iris.AttributeConstraint(STASH='m01s30i207')
    save_atomic(iris.load_cube(src["p8"], zgcon & levelsg), out["geo"])
if need["mslp"]:
    save_atomic(iris.load_cube(src["p8"], iris.Constraint(name='air_pressure_at_sea_level')), out["mslp"])
if need["wsp"]:
    save_atomic(iris.load_cube(src["pc"], iris.Constraint(name='wind_speed')), out["wsp"])

print('time taken ', time.time() - start)

# Only reached on success: free the staged copy (~450GB/month on a GWS at
# ~96%). On failure it is deliberately kept, so a retry doesn't need another
# tape recall - `tc_workflow.py stage` reports leftover staged months.
E.remove_staged_month(name, year, mon_num)
