import iris
import os
import re
import shutil
import subprocess
import sys
import time

MONTH_ABBR = {
    1: "jan", 2: "feb", 3: "mar", 4: "apr", 5: "may", 6: "jun",
    7: "jul", 8: "aug", 9: "sep", 10: "oct", 11: "nov", 12: "dec",
}

# Holding id of the NLDS archive covering each EERIE run's early (now
# off-disk) months - see tc_workflow.py's EERIE comment block for the real
# usable date ranges. Add an entry here once a run's early months are
# archived and confirmed retrievable; a run missing from this dict is
# assumed to still be fully on disk (checked below regardless).
NLDS_HOLDING_ID = {
    "dz876": 5455,
}
NLDS_BIN = os.path.expanduser("~/.local/bin/nlds")
STAGING_ROOT = "/gws/ssde/j25a/hrcm/plvidale/nlds_staging/"
NLDS_TERMINAL_STATES = {"COMPLETE", "FAILED", "COMPLETE_WITH_ERRORS", "COMPLETE_WITH_WARNINGS"}
NLDS_POLL_SECONDS = 20
NLDS_TIMEOUT_SECONDS = 7200

name = sys.argv[1]
year = sys.argv[2]
mon_num = int(sys.argv[3])

if mon_num not in MONTH_ABBR:
    raise ValueError(f"month must be 1-12, got {mon_num}")
monname = MONTH_ABBR[mon_num]
mon = f"{mon_num:02d}"

#PLV's runs - each month's data is archived under a directory dated one month
#ahead of the data itself, e.g. 1950jan data sits under u-dz876/19500201T0000Z/
dir_year, dir_mon_num = int(year), mon_num + 1
if dir_mon_num == 13:
    dir_mon_num = 1
    dir_year += 1
cycle_date = f"{dir_year}{dir_mon_num:02d}01T0000Z"

if name == "dz876" and year == "1959" and mon_num == 12:
    # this one cycle bundled Nov+Dec 1959 output into the Nov cycle dir
    # instead of giving Dec its own - the normal +1-month dir (19600101)
    # has no atmosphere pp files at all; dz876a.p81959dec/pc1959dec etc.
    # actually sit in 19591201T0000Z alongside the 1959nov files.
    cycle_date = "19591201T0000Z"

dirin='/gws/ssde/j25a/hiresgw/aosprey/u-'+name+'/'+cycle_date+'/'
dirin_lev='/gws/ssde/j25a/hiresgw/aosprey/u-'+name+'/'+cycle_date+'/'

dirout = '/work/scratch-pw4/plvidale/u-'+name+'/'
os.makedirs(dirout, exist_ok=True)

namein=name+'a.pc'+year+monname
namein_lev=name+'a.p8'+year+monname

filenamein       = dirin + namein
filenamein_lev   = dirin_lev +namein_lev

print(filenamein)

#TO FIND OUT WHAT IS IN A PP FILE:
#cubes=iris.load(filenamein)
#print(cubes[6])
#print(cubes[7])
#print(cubes[8])
#print(cubes.standard_name)
#print(cubes.long_name)

filenameout  = dirout+name+'_'+year+mon+'_uwinds.nc'
filenameout1 = dirout+name+'_'+year+mon+'_vwinds.nc'
filenameout2 = dirout+name+'_'+year+mon+'_geo.nc'
filenameout3 = dirout+name+'_'+year+mon+'_mslp.nc'
filenameout4 = dirout+name+'_'+year+mon+'_10m_wspeed.nc'

_STATE_RE = re.compile(r"state\s*:\s*(\S+)")
_staged_dirs = []  # staging dirs created this run - always cleaned up on exit


def _nlds_get(original_path, staging_dir, job_label):
    subprocess.run(
        [NLDS_BIN, "get", original_path, "-A", "-i", str(NLDS_HOLDING_ID[name]),
         "-r", staging_dir, "-b", job_label],
        check=True, capture_output=True, text=True,
    )
    deadline = time.time() + NLDS_TIMEOUT_SECONDS
    first = True
    while time.time() < deadline:
        if not first:
            time.sleep(NLDS_POLL_SECONDS)
        first = False
        out = subprocess.run([NLDS_BIN, "stat", "-b", job_label],
                              check=True, capture_output=True, text=True).stdout
        m = _STATE_RE.search(out)
        state = m.group(1) if m else None
        if state in NLDS_TERMINAL_STATES:
            if state != "COMPLETE":
                raise RuntimeError(f"NLDS get of {original_path} ended in state "
                                    f"{state!r} (job_label={job_label}):\n{out}")
            return
    raise TimeoutError(f"NLDS get of {original_path} did not finish within "
                        f"{NLDS_TIMEOUT_SECONDS}s (job_label={job_label})")


def ensure_local(original_path, stream_tag):
    """Return a locally-readable path for original_path, retrieving it from
    NLDS to a scratch-of-permanent-disk staging dir if it's no longer on
    disk at its normal location. Any staging dir this creates is deleted
    again before the process exits (see cleanup at the bottom)."""
    if os.path.exists(original_path):
        return original_path
    if name not in NLDS_HOLDING_ID:
        raise FileNotFoundError(
            f"{original_path} is missing on disk and no NLDS holding is "
            f"configured for run {name!r} in NLDS_HOLDING_ID - add one once "
            f"its early months are archived.")
    staging_dir = STAGING_ROOT + f"u-{name}/{cycle_date}/"
    os.makedirs(staging_dir, exist_ok=True)
    if staging_dir not in _staged_dirs:
        _staged_dirs.append(staging_dir)
    job_label = f"{name}_{year}{mon}_{stream_tag}_auto"
    print(f"{original_path} not on disk - retrieving from NLDS "
          f"(holding {NLDS_HOLDING_ID[name]}, job_label={job_label})")
    _nlds_get(original_path, staging_dir, job_label)
    staged_path = staging_dir + os.path.basename(original_path)
    if not os.path.exists(staged_path):
        raise RuntimeError(f"NLDS reported success but {staged_path} is not there")
    return staged_path


try:
    filenamein_lev = ensure_local(filenamein_lev, "p8")
    filenamein = ensure_local(filenamein, "pc")

    start = time.time()
    print("start ", start)

    levelsw = iris.Constraint(pressure=[250,500,600,700,850,925])
    windu = iris.Constraint(name='x_wind')
    cubeu  = iris.load_cube(filenamein_lev, windu&levelsw)
    iris.save(cubeu, filenameout)
    windv = iris.Constraint(name='y_wind')
    cubev  = iris.load_cube(filenamein_lev, windv&levelsw)
    iris.save(cubev, filenameout1)

    levelsg = iris.Constraint(pressure=[250,600,925])
    zgcon = iris.AttributeConstraint(STASH='m01s30i207')
    cubeg   = iris.load_cube(filenamein_lev, zgcon&levelsg)
    iris.save(cubeg, filenameout2)

    levelsp = iris.Constraint(name='air_pressure_at_sea_level')
    cubep   = iris.load_cube(filenamein_lev, levelsp)
    iris.save(cubep, filenameout3)

    levelwsp  = iris.Constraint(name='wind_speed')
    cubewsp   = iris.load_cube(filenamein, levelwsp)
    iris.save(cubewsp, filenameout4)

    end = time.time()
    print('time taken ', end - start)
finally:
    # Staged pp copies live on permanent gws disk, not scratch - always
    # delete them again (success or failure) so a month's worth (~500GB)
    # never lingers; gws_hrcm is already at 94% of its 1PB quota.
    for d in _staged_dirs:
        print(f"removing staged NLDS copy {d}")
        shutil.rmtree(d, ignore_errors=True)
