import os
import sys

run     = sys.argv[1]
year    = sys.argv[2]
tr      = sys.argv[3] # CAREFUL: there is currently a dependence in vorticity generation that is done at T63 and is needed at T42. Can re-code later.
home    = sys.argv[4]
scrdir  = sys.argv[5]
link_var_log = sys.argv[6]
add_names    = [sys.argv[7]]

levs=['850','700','600','500','250']

def cmd_do(command,printflag='False'):
    os.system(command)
    if (printflag):
        print(command)

def link_variables(varnames, run, scrdir, year, lev):
    for var in varnames:
        filename=scrdir+run+"_"+year+var
        linkname="indat/"+run+"_"+year+var
        if not os.path.exists(filename):
            cmd_do("cdo cat "+scrdir+run+"_"+year+"{01..12}"+var+" "+filename)
            cmd_do("ln -sf " + filename + " " + linkname)

def link_lev_full_vort(vortname_nc):
    # select and LINK full resolution level vorticity
    if not os.path.exists("indat/"+vortname_nc):
        #cmd_do("converters/bin2nc indat/" + vortname_dat + " " +scrdir+vortname)
        #cmd_do("cdo selmon,5,6,7,8,9,10,11 "+ scrdir + "/annual/"+vortname_nc + " "+scrdir+vortname_nc,printflag=True)
        cmd_do("cdo selmon,1,2,3,4,5,6,7,8,9,10,11,12 "+ scrdir + "/annual/"+vortname_nc + " "+scrdir+vortname_nc,printflag=True)
        cmd_do("ln -sf " + scrdir+vortname_nc + " indat/"+vortname_nc,printflag=True)

def names(run,year,tr):
    specname              = run+"_"+year+"_"+tr
    specname_1            = "specfil."+specname+"_band001"
    specname_0            = "specfil."+specname+"_band000"
    return specname,specname_0,specname_1
           
def filter_lev_vort(lev,run,year,tr):
    # Filter each LEVEL to specified TRUNCATION
    vort_tr_name    = "vor_"+year+"_"+lev+"_"+tr+"_filt.dat"
    vort_tr_name_nc = "vor_"+year+"_"+lev+"_"+tr+"_filt.nc"
    specname,specname_0,specname_1=names(run,year,tr)
    #print("inside filter_levl_vort: ",lev,run,year,tr)
    if not os.path.exists("indat/"+vort_tr_name):
        if tr == 'T63':
            cmd_do("bin/track.pl -i " + vortname_nc + " -f "+specname+" <  specfilt_era5_"+tr+".in", printflag=True)
            cmd_do("mv outdat/"+specname_0+" "+scrdir+vort_tr_name,  printflag=True)
        if tr == 'T42':
            cmd_do("bin/track.pl -i " + vortname_nc + " -f "+specname+" <  specfilt_era5_track_"+tr+".in", printflag=True)
            cmd_do("mv outdat/"+specname_1+" "+scrdir+vort_tr_name,  printflag=True)
        cmd_do("ln -sf " + scrdir+vort_tr_name + " indat/"+vort_tr_name, printflag=True)

def avg_filter_lev_vort(run,year,tr):
    # Filter the full resolution level vorticities to T63 average vertical vorticity
    #vort_T63_avg_name =  "vor_"+year+"_avg"+"_T63_filt.dat"
    vort_tr_avg_filt_name =  "vor_"+year+"_avg"+"_"+tr+"_filt.nc"
    vortname_avg_nc       =  "vor_"+year+"_avg.nc"
    specname,specname_0,specname_1=names(run,year,tr)
    if ( (not os.path.exists("indat/"+vort_tr_avg_filt_name)) and (os.path.exists("indat/vor_"+year+"_850.nc") and os.path.exists("indat/vor_"+year+"_700.nc") and os.path.exists("indat/vor_"+year+"_600.nc")) ):
        cmd_do("cdo -O ensavg indat/vor_"+year+"_850.nc indat/vor_"+year+"_700.nc indat/vor_"+year+"_600.nc " + scrdir+vortname_avg_nc)
        cmd_do("ln -sf " + scrdir+vortname_avg_nc + " indat/"+vortname_avg_nc)
        cmd_do("bin/track.pl -i " + vortname_avg_nc + " -f "+specname+" <  specfilt_era5_track_"+tr+".in")
        cmd_do("mv outdat/"+specname_1+" "+scrdir+vort_tr_avg_filt_name)
        cmd_do("ln -sf " + scrdir+vort_tr_avg_filt_name + " indat/"+vort_tr_avg_filt_name)
        cmd_do("rm outdat/"+specname_0)
    return vort_tr_avg_filt_name      
           
def master_tracking(run,year,tr,vort_tracking_name):
    ## TRACK each year at truncation tr
    if tr == 'T63':
        trackdir=run+'/VOR_VERTAVG_'+tr+'filt_'+year
    else:
        trackdir=run+'/VOR_850_'+tr+'filt_'+year
    trackname=trackdir+'/tr_trs_pos'
    trackfile=trackname+'.gz'
    specname,specname_0,specname_1=names(run,year,tr)
    if not os.path.exists(trackfile):
        print('TRACKING AT '+tr+' using vertical vorticity file called: '+vort_tracking_name)
        #MAY-NOV setting
        #cmd_do('master -c='+trackdir+' -e=track.pl -f=y' +year+' -i='+vort_T63_avg_name+' -j=RUN_AT.in -n=1,64,13 -o='+home+' -r=RUN_AT_ -s=RUNDATIN.vor_global_hilat_NH_T63',printflag=False)
        #JAN-DEC setting
        cmd_do('master -c='+trackdir+' -e=track.pl -f=' +specname+' -i='+vort_tracking_name+' -j=RUN_AT.in -n=1,62,24 -o='+home+' -r=RUN_AT_ -s=RUNDATIN.vor_global_hilat_NH_'+tr,printflag=False)
    return trackdir,trackname

def add_fields(add_names,run,year,tr,trackname,two_day_file):
    trackfile_added=two_day_file+'_addvor'
    if not os.path.exists(trackfile_added):
        cmd_do('gunzip '+trackname+'.gz')
        nml_ext='_'+run
        cmd_do('cp addvor_NH_'+tr+'_era5.in addvor'+nml_ext, printflag=True)
        outfile=two_day_file
        for addv in add_names:
            addname=addv+nml_ext
            waddname='work_'+year+'_'+addname
            outfile=outfile+'_'+addv
            #print("NOW DOING NAMELIST:   ",addname,'    ',outfile)
            cmd_do('cp '+addname+' '+waddname)
            if tr == 'T42':
                cmd_do("sed -i 's/YYYY/"+year+"/g; s/RUN/"+run+"/g; s/VORLEV/VOR_850/g; s/TTT/"+tr+"/g' "+waddname)
            else:
                cmd_do("sed -i 's/YYYY/"+year+"/g; s/RUN/"+run+"/g; s/VORLEV/VOR_VERTAVG/g; s/TTT/"+tr+"/g' "+waddname)
            cmd_do('bin/track.pl -f y'+year+' < '+waddname)
            cmd_do('mv outdat/ff_trs.y'+year+'_addfld '+outfile)   
        cmd_do('gzip '+trackname)
    return trackfile_added
    
def tcident(trackfile_tcident,home,trackdir):
    if not os.path.exists(trackdir+'/'+trackfile_tcident):
        cmd_do('ln -fs '+home+'/utils/TC/lmask_linux_wrap_hires.dat '+home+trackdir+'/lmask.dat')
        os.chdir(trackdir)   # cd into each track directory
        cmd_do(home+'utils/bin/tcident < '+home+'tcident_era5.in')
        cmd_do('mv track.dat '         + home + trackfile_tcident+'.wc' , printflag=True)
        cmd_do(home+'utils/bin/tr2nc ' + home + trackfile_tcident       + ' s '+home+'utils/TR2NC/tr2nc.meta', printflag=True)
        cmd_do(home+'utils/bin/tr2nc ' + home + trackfile_tcident+'.wc' + ' s '+home+'utils/TR2NC/tr2nc.meta', printflag=True)
        os.chdir(home)

########    MAIN PART, THIS IS WHERE WE DO THE FULL SEQUENCE OF TRACKING AT EACH TRUNCATION
os.chdir(home)
for il,lev in enumerate(levs):
    # Create full resolution vorticities
    vortname     = "vor_"+year+"_"+lev
    vortname_dat = vortname+".dat"
    vortname_nc  = vortname+".nc"
    if tr == 'T63':
        link_lev_full_vort(vortname_nc)
        filter_lev_vort(lev,run,year,tr)
    if tr == 'T42' and lev == '850':
        link_lev_full_vort(vortname_nc)
        filter_lev_vort(lev,run,year,tr)
            
if tr == 'T63':
    vort_tracking_name=avg_filter_lev_vort(run,year,tr)
else:
    vort_tracking_name="vor_"+year+"_850"+"_"+tr+"_filt.dat"

# MASTER TRACKING
trackdir, trackname = master_tracking(run,year,tr,vort_tracking_name)

# DONE TRACKING, now adding fields
two_day_file = trackname+'.2day'
trackfile_added=add_fields(add_names,run,year,tr,trackname,two_day_file)

#   now do the tcident
trackfile_tcident=trackfile_added+'.tcident'
tcident(trackfile_tcident,home,trackdir)

# DONE TRACKING
#cmd_do('rm /home/users/plvidale/TRACK-1.5.2/outdat/*')
quit()
