import iris
import time
import sys
import os

MONTH_ABBR = {
    1: "jan", 2: "feb", 3: "mar", 4: "apr", 5: "may", 6: "jun",
    7: "jul", 8: "aug", 9: "sep", 10: "oct", 11: "nov", 12: "dec",
}

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
