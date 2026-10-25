import iris
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

if sixhour_resample:
    if winds:
        command = 'cdo -P 16 selhour,0,6,12,18 ' + filenameout + ' ' + tempout
        os.system(command)
        command = 'mv ' + tempout + ' ' +filenameout
        os.system(command)

        command = 'cdo -P 16 selhour,0,6,12,18 ' + filenameout1 + ' ' + tempout1
        os.system(command)
        command = 'mv ' + tempout1 + ' ' +filenameout1
        os.system(command)

    command = 'cdo -P 16 selhour,0,6,12,18 ' + filenameout2 + ' ' + tempout2
    os.system(command)
    command = 'mv ' + tempout2 + ' ' +filenameout2
    os.system(command)

    command = 'cdo -P 16 selhour,0,6,12,18 ' + filenameout3 + ' ' + tempout3
    os.system(command)
    command = 'mv ' + tempout3 + ' ' +filenameout3
    os.system(command)

    command = 'cdo -P 16 selhour,' + WSPEED_SELHOUR + ' ' + filenameout4 + ' ' + tempout4
    os.system(command)
    command = 'mv ' + tempout4 + ' ' +filenameout4
    os.system(command)

command = 'cdo -P 16 remapcon,n640 '+filenameout+' '+filenameout_gg
os.system(command)
command = 'rm '+filenameout
os.system(command)

command = 'cdo -P 16 remapcon,n640 '+filenameout1+' '+filenameout1_gg
os.system(command)
command = 'rm '+filenameout1
os.system(command)

end = time.time()
print('time taken ', end - start)

