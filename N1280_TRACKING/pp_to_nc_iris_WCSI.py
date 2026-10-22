import iris
import time
import sys
import os

name = sys.argv[1]
year = sys.argv[2]
mon  = sys.argv[3]

# Raw .pp source group per experiment - keep in sync with tc_workflow.py's
# SOURCE_DIRS/EXPERIMENTS tables, which already encode this same split.
SOURCE_GROUP = {
    "ch330": "n1280run", "ck777": "n1280run", "ck778": "n1280run",
    "dd828": "plvidale",   # PLV's more recent runs
    "dc009": "torau",      # Torsten's runs
}
group = SOURCE_GROUP[name]

# 10m winds, geopotential and mslp are always in stream apc. Level winds
# (u/v) are also in apc for the n1280run group (ch330/ck777/ck778), but in a
# separate ap7 stream for plvidale/torau runs (dd828/dc009) - see
# filenamein_lev below.
if group == "n1280run":
    dirin = ('/gws/ssde/j25a/hrcm/n1280run/u-%s/%s/'
             % (name, "6hrly" if int(year) < 2000 else "apc"))
    dirin_lev = dirin
elif group == "plvidale":
    dirin     = '/gws/ssde/j25a/hrcm/plvidale/u-%s/apc/' % name
    dirin_lev = '/gws/ssde/j25a/hrcm/plvidale/u-%s/ap7/' % name
elif group == "torau":
    dirin     = '/gws/ssde/j25a/hrcm/torau/u-%s/apc/' % name
    dirin_lev = '/gws/ssde/j25a/hrcm/torau/u-%s/ap7/' % name
else:
    raise ValueError("no raw-data source configured for experiment %r" % name)

dirout = '/work/scratch-pw4/plvidale/u-' + name + '/'
if not os.path.exists(dirout):
    os.system('mkdir ' + dirout)

if len(mon) < 2:
    mon = '0' + str(mon)
else:
    mon = str(mon)

namein     = name + 'a.pc' + year + mon + '*.pp'
namein_lev = name + 'a.p7' + year + mon + '*.pp'

filenamein = dirin + namein
# n1280run keeps level winds in the same apc-stream file as everything else
# (filenamein); plvidale/torau runs keep them in a separate ap7-stream file.
filenamein_lev = filenamein if group == "n1280run" else dirin_lev + namein_lev

filenameout  = dirout+name+'_'+year+mon+'_uwinds.nc'
filenameout1 = dirout+name+'_'+year+mon+'_vwinds.nc'
filenameout2 = dirout+name+'_'+year+mon+'_geo.nc'
filenameout3 = dirout+name+'_'+year+mon+'_mslp.nc'
filenameout4 = dirout+name+'_'+year+mon+'_10m_wspeed.nc'

start = time.time()
print("start ", start)

levelsw = iris.Constraint(pressure=[250,500,600,700,850,925])

windu = iris.Constraint(name='x_wind')
cubeu  = iris.load_cube(filenamein_lev, windu&levelsw)
iris.save(cubeu, filenameout)

windv = iris.Constraint(name='y_wind')
cubev  = iris.load_cube(filenamein_lev, windv&levelsw)
iris.save(cubev, filenameout1)

# mslp is NOT in the apc-stream file for plvidale/torau runs (dd828/dc009) -
# confirmed by direct inspection: apc only has surface_air_pressure (a
# different, non-reduced-to-sea-level field), while air_pressure_at_sea_level
# lives in the ap7 stream alongside the level winds. For n1280run
# (ch330/ck777/ck778), filenamein_lev == filenamein, so this is a no-op there
# - the 170-job dd828/dc009 batch that hit this (2026-10-22) got a Python
# ConstraintMismatchError crash on this exact line, after u/v winds had
# already succeeded, which is why only uwinds/vwinds ever got produced.
levelsp = iris.Constraint(name='air_pressure_at_sea_level')
cubep   = iris.load_cube(filenamein_lev, levelsp)
iris.save(cubep, filenameout3)

levelwsp  = iris.Constraint(name='wind_speed')
cubewsp   = iris.load_cube(filenamein, levelwsp)
iris.save(cubewsp, filenameout4)

levelsg = iris.Constraint(pressure=[300,600,900])
zgcon = iris.AttributeConstraint(STASH='m01s30i207')
cubeg   = iris.load_cube(filenamein, zgcon&levelsg)
iris.save(cubeg, filenameout2)

end = time.time()
print('time taken ', end - start)
