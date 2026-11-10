import iris
import subprocess
import time
import sys
import os

name = sys.argv[1]

#dirin ='/gws/nopw/j04/hrcm/n1280run/u-'+name+'/6hrly'
dirin ='/gws/ssde/j25a/hrcm/plvidale/u-'+name+'/6hrly'
#dirin ='/gws/nopw/j04/hrcm/n1280run/u-'+name+'/6hrly/'
#dirout = '/work/scratch-pw3/plvidale/u-'+name+'/'
dirout = '/work/scratch-pw4/plvidale/u-'+name+'/'

year=sys.argv[2]
mon =sys.argv[3]

if len(mon)<2:
    mon = '0'+str(mon)
else:
    mon = str(mon)

namein=name+'a.pc'+year+mon+'*.pp'

filenamein   = dirin +namein

filenameout  = dirout+name+'_'+year+mon+'_uwinds.nc'
filenameout1 = dirout+name+'_'+year+mon+'_vwinds.nc'
filenameout2 = dirout+name+'_'+year+mon+'_geo.nc'
filenameout3 = dirout+name+'_'+year+mon+'_mslp.nc'
filenameout4 = dirout+name+'_'+year+mon+'_10m_wspeed.nc'

tempout  = dirout+name+'_'+year+mon+'_temp_uwinds.nc'
tempout1 = dirout+name+'_'+year+mon+'_temp_vwinds.nc'
tempout2 = dirout+name+'_'+year+mon+'_temp_geo.nc'
tempout3 = dirout+name+'_'+year+mon+'_temp_mslp.nc'
tempout4 = dirout+name+'_'+year+mon+'_temp_10m_wspeed.nc'

filenameout_gg  = dirout+name+'_'+year+mon+'_uwinds_gg.nc'
filenameout1_gg = dirout+name+'_'+year+mon+'_vwinds_gg.nc'

def file_ntime(filename):
    try:
        out = subprocess.check_output(["cdo", "-s", "ntime", filename], stderr=subprocess.DEVNULL)
        return int(out.decode().split()[0])
    except Exception:
        return -1

def file_days(filename):
    try:
        out = subprocess.check_output(["cdo", "-s", "showdate", filename], stderr=subprocess.DEVNULL)
        return len(set(out.decode().split()))
    except Exception:
        return -1

def already_6hourly(filename):
    # True only if this file is already down to a 6-hourly-or-coarser
    # cadence (<=4 steps/day, +2 slack for the stray edge-of-month
    # timestamp cdo cat sometimes pulls in) - regardless of WHICH hours,
    # so this also correctly recognises wspeed's legitimate 3/9/15/21
    # offset cadence (see WSPEED_SELHOUR below) as already-done.
    #
    # Found 2026-11-11: dd828's Jan-Mar 2013 mslp sat genuinely 3-hourly (8
    # steps/day) under its FINAL filename because pp_to_nc_iris_WCSI.py had
    # re-extracted it fresh from raw .pp (its own skip-if-exists check
    # isn't per-variable - one missing wind intermediate left over from an
    # interrupted earlier attempt was enough to trigger redoing all 5
    # outputs) while this script's OLD top-of-file gate only checked that
    # all 5 final files existed, not that they'd actually been resampled -
    # so it exited early and the fresh 3-hourly mslp was never downsampled.
    nt = file_ntime(filename)
    nd = file_days(filename)
    return nt >= 0 and nd > 0 and nt <= nd * 4 + 2

start = time.time()
print("start ", start)

sixhour_resample = True
winds = True

# 10m_wspeed is natively 6-hourly but on a 03/09/15/21 UTC cycle for dd828/
# dc009/ec100 (confirmed directly for each: every delta is exactly 6h, never
# a 0/6/12/18 timestamp) - 3h offset from everything else (geo/mslp/winds,
# all 0/6/12/18 for every experiment including these). selhour,0,6,12,18
# then finds zero matching timesteps and aborts, leaving the original
# unresampled file in place under its final filename (harmless but noisy -
# "Abort: No timesteps selected!" plus a failed mv). Kevin confirmed
# (2026-10-22, re dd828/dc009) the 3h offset is fine to live with for now -
# addwind doesn't need it fixed immediately - so just select the hours this
# field actually has instead of erroring. ec100's EERIE conversion hit the
# identical pattern (found 2026-10-22 from its first real test job,
# 61572175) - same fix applies. dz876 not yet checked; verify its wspeed
# hours the same way before assuming 0,6,12,18 works for it. Revisit
# alongside the ap9-hourly-derived-windspeed plan noted in memory once
# that's scheduled.
WSPEED_SELHOUR = {"dd828": "3,9,15,21", "dc009": "3,9,15,21", "ec100": "3,9,15,21"}.get(name, "0,6,12,18")

# Fast-path skip for months that are genuinely fully done: both regrids
# already produced, and geo/mslp/wspeed already at a 6-hourly-or-coarser
# cadence (not just present - see already_6hourly's docstring above for why
# existence alone isn't enough). winds/filenameout(1) aren't checked here
# since a successful regrid already deletes them, so gg existing is the
# right signal for them specifically.
if (os.path.exists(filenameout_gg) and os.path.exists(filenameout1_gg)
        and already_6hourly(filenameout2) and already_6hourly(filenameout3)
        and already_6hourly(filenameout4)):
    print("Regridded output for "+name+" "+year+mon+" already exists and is correctly resampled - skipping.")
    sys.exit(0)

# Per-file, not all-or-nothing (replaces the old single os.path.exists-only
# gate that used to sit at the top of this script - see already_6hourly's
# docstring above for why that was wrong): each of winds/geo/mslp/wspeed is
# only resampled if it isn't ALREADY at a 6-hourly-or-coarser cadence, and
# only if its input file is actually present (it may already have been
# deleted by a prior successful regrid, for winds specifically).
RESAMPLE_TARGETS = [
    (filenameout,  tempout,  '0,6,12,18',    winds),
    (filenameout1, tempout1, '0,6,12,18',    winds),
    (filenameout2, tempout2, '0,6,12,18',    True),
    (filenameout3, tempout3, '0,6,12,18',    True),
    (filenameout4, tempout4, WSPEED_SELHOUR, True),
]

if sixhour_resample:
    for fname, tname, hours, enabled in RESAMPLE_TARGETS:
        if not enabled or not os.path.exists(fname) or already_6hourly(fname):
            continue
        command = 'cdo -P 16 selhour,' + hours + ' ' + fname + ' ' + tname
        os.system(command)
        if os.path.exists(tname):
            command = 'mv ' + tname + ' ' + fname
            os.system(command)

# Only delete the pre-regrid intermediate once remapcon has actually
# produced the final "_gg" file - found 2026-11-10, diagnosing a gap where
# dd828 had uwinds_gg for only 198 of ~455 months despite vwinds_gg being
# complete: the old unconditional rm here deleted filenameout/filenameout1
# even when cdo remapcon failed (e.g. a scratch quota hit mid-regrid),
# leaving neither the intermediate nor the final file - a dead end that
# could only be recovered by re-running pp_to_nc_iris_WCSI.py from the raw
# pp file for that month. Same verify-before-delete pattern already used in
# track_n1280_para_WCSI.py's link_variables().
if os.path.exists(filenameout):
    command = 'cdo -P 16 remapcon,n640 '+filenameout+' '+filenameout_gg
    os.system(command)
    if os.path.exists(filenameout_gg):
        command = 'rm '+filenameout
        os.system(command)
    else:
        print('remapcon failed to produce '+filenameout_gg+' - keeping '+filenameout+' for retry')

if os.path.exists(filenameout1):
    command = 'cdo -P 16 remapcon,n640 '+filenameout1+' '+filenameout1_gg
    os.system(command)
    if os.path.exists(filenameout1_gg):
        command = 'rm '+filenameout1
        os.system(command)
    else:
        print('remapcon failed to produce '+filenameout1_gg+' - keeping '+filenameout1+' for retry')

end = time.time()
print('time taken ', end - start)

