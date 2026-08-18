import os
import sys

run     = sys.argv[1]
year    = sys.argv[2]
tr      = sys.argv[3] # CAREFUL: there is currently a dependence in vorticity generation that is done at T63 and is needed at T42. Can re-code later.
home    = '/home/users/plvidale/TRACK-1.5.2/'
scrdir  = sys.argv[4]
permdir  = sys.argv[5]
test = sys.argv[6]
link_var_log = [False, True][test.lower()[0] == 't']
test = sys.argv[8]
only_tcident = [False, True][test.lower()[0] == 't']
hemisphere = sys.argv[9]
if hemisphere == '_SH':
    hemi = '_SH'
else:
    hemi = ''
test = sys.argv[10] if len(sys.argv) > 10 else 'False'
skip_to_hart = [False, True][test.lower()[0] == 't']
print('start of script: HEMI is: ',hemi,' tc_ident is: ',only_tcident,' skip_to_hart is: ',skip_to_hart)

add_names = ['addvor','addmslp','addwind']
levs      = ['850','700','600','500','250']

def cmd_do(command,printflag=False):
    os.system(command)
    if (printflag):
        print(command)

def link_variables(varnames, run, scrdir, year, hemi):
    for var in varnames:
        filename=  scrdir+run+"_"+year+var+hemi+".nc"
        linkname="indat/"+run+"_"+year+var+hemi+".nc"
        if not os.path.exists(filename):
            if hemi == '':
                if year == '1979':
                    cmd_do("cdo cat "+scrdir+run+"_"+year+"{03..12}"+var+".nc "+filename,printflag=True)
                else:
                    cmd_do("cdo cat "+scrdir+run+"_"+year+"{01..12}"+var+".nc "+filename,printflag=True)
                cmd_do("ln -sf " + filename + " " + linkname)
                cmd_do("rm "+scrdir+run+"_"+year+"??"+var+'.nc',printflag=True)
            elif hemi == '_SH':
                yearp1 = str(int(year)+1)
                cmd_do("cdo cat "+scrdir+run+"_"+year+"{07..12}"+var+".nc "+scrdir+run+"_"+yearp1+"{01..06}"+var+".nc " + filename,printflag=True)
                cmd_do("ln -sf " + filename + " " + linkname)
                cmd_do("rm "+scrdir+run+"_"+year+"{07..12}"+var+".nc "+scrdir+run+"_"+yearp1+"{01..06}"+var+".nc",printflag=True)

def link_lev_full_vort(vortname_nc):
    # select and LINK full resolution level vorticity
    if not os.path.exists("indat/"+vortname_nc):
        cmd_do("ln -sf " + scrdir+vortname_nc + " indat/"+vortname_nc,printflag=True)

def names(run,year,tr):
    specname              = run+"_"+year+"_"+tr
    specname_1            = "specfil."+specname+"_band001"
    specname_0            = "specfil."+specname+"_band000"
    return specname,specname_0,specname_1
           
def compute_vort(lev,run,year,scrdir,hemi):
    # Create full resolution vorticities
    if hemi == '':
        vortname     = "vor_"+run+"_"+year+"_"+lev
    elif hemi == '_SH':
        vortname     = "vor_"+run+"_"+year+"_"+lev+"_SH"
    vortname_dat = vortname+".dat"
    vortname_nc  = vortname+".nc"
    work_vorcalc = "work_vorcalc_"+run+"_"+year+".in"
    if not os.path.exists("indat/"+vortname_dat):
        cmd_do('cp vorcalc.in '+work_vorcalc)
        if hemi == '':
            cmd_do("sed -i 's/YYYY/"+year+"/g; s/RUN/"+run+"/g; s/LEV/"+lev+"/g; s|indat|"+scrdir+"|g' "+work_vorcalc,printflag=True)
        elif hemi == '_SH':
            cmd_do("sed -i 's/YYYY/"+year+"/g; s/RUN/"+run+"/g; s/LEV.dat/LEV_SH.dat/g; s/LEV/"+lev+"/g; s/_gg.nc/_gg_SH.nc/g; s|indat|"+scrdir+"|g' "+work_vorcalc,printflag=True)
        cmd_do('bin/track.pl <  '+work_vorcalc,printflag=True)
        cmd_do("ln -sf " + scrdir + vortname_dat + " indat/" + vortname_dat,printflag=True)
    # Convert full resolution level vorticity to netcdf
    if not os.path.exists("indat/"+vortname_nc):
        cmd_do("converters/bin2nc indat/" + vortname_dat + " " + scrdir + vortname + " 0")
        cmd_do("ln -sf " + scrdir+vortname_nc + " indat/"+vortname_nc)
    return vortname_nc

def filter_lev_vort(lev,run,year,tr):
    # Filter each LEVEL to specified TRUNCATION
    vort_tr_name    = "vor_"+run+"_"+year+"_"+lev+"_"+tr+"_filt"+hemi+".dat"
    vort_tr_name_nc = "vor_"+run+"_"+year+"_"+lev+"_"+tr+"_filt"+hemi+".nc"
    specname,specname_0,specname_1=names(run,year,tr)
    if not os.path.exists("indat/"+vort_tr_name):
        cmd_do("bin/track.pl -i " + vortname_nc + " -f "+specname+" <  specfilt_n1280_"+tr+".in", printflag=True) # as we do in ERA5: for the netcdf file
        cmd_do("mv outdat/"+specname_0+" "+scrdir+vort_tr_name,  printflag=True)
        cmd_do("ln -sf " + scrdir+vort_tr_name + " indat/"+vort_tr_name, printflag=True)

def avg_filter_lev_vort(run,year,tr):
    # Filter the full resolution level vorticities to T63 average vertical vorticity
    vort_tr_avg_filt_name =  "vor_"+run+"_"+year+"_avg"+"_"+tr+"_filt"+hemi+".dat"
    vortname_avg_nc       =  "vor_"+run+"_"+year+"_avg"+hemi+".nc"
    specname,specname_0,specname_1=names(run,year,tr)
    if ( (not os.path.exists("indat/"+vort_tr_avg_filt_name)) and (os.path.exists("indat/vor_"+run+"_"+year+"_850"+hemi+".nc") and os.path.exists("indat/vor_"+run+"_"+year+"_700"+hemi+".nc") and os.path.exists("indat/vor_"+run+"_"+year+"_600"+hemi+".nc")) ):
        cmd_do("cdo -O ensavg indat/vor_"+run+"_"+year+"_850"+hemi+".nc indat/vor_"+run+"_"+year+"_700"+hemi+".nc indat/vor_"+run+"_"+year+"_600"+hemi+".nc " + scrdir+vortname_avg_nc)
        cmd_do("ln -sf " + scrdir+vortname_avg_nc + " indat/"+vortname_avg_nc)
        cmd_do("bin/track.pl -i " + vortname_avg_nc + " -f "+specname+" <  specfilt_n1280_track_"+tr+".in")
        cmd_do("mv outdat/"+specname_1+" "+scrdir+vort_tr_avg_filt_name)
        cmd_do("ln -sf " + scrdir+vort_tr_avg_filt_name + " indat/"+vort_tr_avg_filt_name)
        cmd_do("rm outdat/"+specname_0)
    return vort_tr_avg_filt_name      
           
def master_tracking(run,year,tr,hemi,vort_tracking_name,trackdir,trackdir_full,trackname):
    trackfile=trackdir_full+trackname
    print('Inside MASTER TRACKING, the trackfile is: ',trackfile)
    specname,specname_0,specname_1=names(run,year,tr)
    if not os.path.exists(trackfile+'.gz'):
        print('TRACKING AT '+tr+' using vertical vorticity file called: '+vort_tracking_name)
        if hemi == '':
            #JAN-DEC setting for NH
            cmd_do('master -c='+trackdir+' -e=track.pl -f=' +specname+' -i='+vort_tracking_name+' -j=RUN_AT.in -n=1,62,24 -o='+home+' -r=RUN_AT_ -s=RUNDATIN.vor_global_hilat_NH_'+tr,printflag=True)
        elif hemi == '_SH':
            #JUL-JUN setting for SH
            cmd_do('master -c='+trackdir+' -e=track.pl -f=' +specname+' -i='+vort_tracking_name+' -j=RUN_AT.in -n=1,62,24 -o='+home+' -r=RUN_AT_ -s=RUNDATIN.vor_global_hilat_SH_'+tr,printflag=True)  
    return trackfile

def add_fields(add_names,run,year,tr,trackfile,two_day_file,hemi):
    trackfile_added=two_day_file+'_addvor_addmslp_addwind' # for HadGEM3 with 3 added fields
    if not os.path.exists(trackfile_added):
        cmd_do('gunzip '+trackfile+'.gz')
        outfile=two_day_file
        for addv in add_names:
            outfile=outfile+'_'+addv
            waddname='work_'+run+'_'+year+'_'+addv+hemi
            if hemi == '':
                cmd_do('cp '+addv+'_NH_'+tr+'.in '+waddname, printflag=True)
            elif hemi == '_SH':
                cmd_do('cp '+addv+'_SH_'+tr+'.in '+waddname, printflag=True)
            if tr == 'T42':
                cmd_do("sed -i 's/YYYY/"+year+"/g; s/RUN/"+run+"/g; s/VORLEV/VOR_850/g; s/TTT/"+tr+"/g' "+waddname)
            else:
                cmd_do("sed -i 's/YYYY/"+year+"/g; s/RUN/"+run+"/g; s/VORLEV/VOR_VERTAVG/g; s/TTT/"+tr+"/g' "+waddname)
            cmd_do('bin/track.pl -f y'+year+' < '+waddname, printflag=True)
            cmd_do('mv outdat/ff_trs.y'+year+'_addfld '+outfile, printflag=True)   
        cmd_do('gzip '+trackfile)
        cmd_do('rm '+two_day_file+'_addvor' )
        cmd_do('rm '+two_day_file+'_addvor_addmslp')
    return trackfile_added
    
def tcident(trackfile_tcident,home,trackdir,year):
    # THE FOLLOWING LINE IS WRONG, but I temporarily leave it as is, so that all tcident is done no matter
    if not os.path.exists(trackdir+'/'+trackfile_tcident):
        os.chdir(trackdir)   # cd into each track directory
        cmd_do('ln -fs '+home+'utils/TC/lmask_linux_wrap_hires.dat ' + 'lmask.dat')
        if hemi == '':
            startdate=year+'010100'
            cmd_do(home+'utils/bin/tcident < '+home+'tcident_n1280.in', printflag=True)
        elif hemi == '_SH':
            startdate=year+'070100'
            cmd_do(home+'utils/bin/tcident < '+home+'tcident_n1280_SH.in', printflag=True)
        cmd_do('mv track.dat '         + trackfile_tcident+'.wc' , printflag=True)
        cmd_do('/home/users/mjrobert/track_148b/TRACK/programs/count_30days '  + trackfile_tcident+' 0 0 5 4 0 ' + startdate + ' 6', printflag=True)
        cmd_do(home+'utils/bin/tr2nc ' + trackfile_tcident+'.new' + ' s '+home+'utils/TR2NC/tr2nc.meta.addmslp.addwind', printflag=True)       
        cmd_do(home+'utils/bin/tr2nc ' + trackfile_tcident        + ' s '+home+'utils/TR2NC/tr2nc.meta.addmslp.addwind', printflag=True)
        cmd_do(home+'utils/bin/tr2nc ' + trackfile_tcident+'.wc'  + ' s '+home+'utils/TR2NC/tr2nc.meta.addmslp.addwind', printflag=True)
        cmd_do('/home/users/mjrobert/track_148b/TRACK/programs/count_30days '  + trackfile_tcident+'.wc 0 0 5 4 0 ' + startdate + ' 6', printflag=True)        
        cmd_do(home+'utils/bin/tr2nc ' + trackfile_tcident+'.wc.new' + ' s '+home+'utils/TR2NC/tr2nc.meta.addmslp.addwind', printflag=True)
        os.chdir(home)
        
def Hart(run, hemi, year, hart_output):
    if os.path.exists(hart_output):
        print('Hart output already exists, skipping Hart step: '+hart_output)
        return
    if hemi == '':
        regin=home+'regZ_NH.in'
        hartin=home+'hart_NH.in'
    elif hemi == '_SH':
        regin=home+'regZ_SH.in'
        hartin=home+'hart_SH.in'
    work_regZ=home+'work_'+run+'_'+year+'_regZ'+hemi+'.in'
    work_Hart=home+'work_'+run+'_'+year+'_hart'+hemi+'.in'
    cmd_do('cp '+regin+' '+work_regZ)
    cmd_do('cp '+hartin+' '+work_Hart)
    cmd_do("sed -i 's/YYYY/"+year+"/g; s/EXPT/"+run+"/g' "+work_regZ, printflag=True)
    cmd_do("sed -i 's/YYYY/"+year+"/g; s/EXPT/"+run+"/g' "+work_Hart, printflag=True)
    regZ_input = home+'outdat/ff_trs.'+run+'_'+year+hemi+'_addfld_reg'
    if not os.path.exists(regZ_input):
        cmd_do(home+'bin/track.pl -f ' +run+'_'+year+hemi+' ' + ' < ' + work_regZ,  printflag=True)
    cmd_do('mv '+regZ_input+' '+home+'outdat/regZ_'+run+'_'+year+hemi+'.dat', printflag=True)
    cmd_do(home+'utils/bin/hart < ' + work_Hart, printflag=True)
    cmd_do('rm '+home+'outdat/regZ_'+run+'_'+year+hemi+'.dat', printflag=True)

########    MAIN PART, THIS IS WHERE WE DO THE FULL SEQUENCE OF TRACKING AT EACH TRUNCATION
os.chdir(home)
tdir_name='u-'+run+'/'
if not os.path.exists(tdir_name):
    os.mkdir(tdir_name)
trackdir=tdir_name+'/VOR_VERTAVG_'+tr+'filt_'+year+hemi
trackdir_full=home+trackdir+'/'

print(run, ' ', year, ' ', tr, ' ', home, ' ', scrdir, ' ', permdir, ' ', link_var_log, ' ', add_names, ' ', only_tcident, ' ', hemi)
print("AT start, TRACKDIR_FULL is: ",trackdir_full)

if hemi == '':
    trackname         = 'tr_trs_pos'
elif hemi == '_SH':
    trackname         = 'tr_trs_neg'

# Concatenate the monthly wind/mslp/geo/10m_wspeed files into annual files and
# link them into indat/ - the Hart regZ step needs the annual geo file even
# when skip_to_hart is True, so this must run unconditionally rather than only
# as a side effect of the (now possibly skipped) per-level vorticity loop below
varnames = ["_uwinds_gg", "_vwinds_gg", "_mslp", "_geo", "_10m_wspeed"]
if link_var_log:
    link_variables(varnames, run, scrdir, year, hemi)

## if only_tcident is true, only carry out the tcident portion of the code
## if only_tcident is false, compute vorticities and track
## if skip_to_hart is true, skip vorticity computation and tracking entirely -
## for experiments already tracked previously, where trackfile_added (the
## .2day_addvor_addmslp_addwind file) already exists and only the Hart step
## (and the final count/tr2nc conversion) still needs to be run

if not skip_to_hart:
    # PREPARATORY STEP: check if vorticities already exist etc.
    if not only_tcident:
        if not os.path.exists(trackdir_full):
            os.mkdir(trackdir_full)
        if hemi == '':
            computed_vortname=permdir+"vor_"+run+"_"+year+"_avg_"+tr+"_filt"+hemi+".dat"
        elif hemi == '_SH':
            computed_vortname=permdir+"vor_"+run+"_"+year+"_avg_"+tr+"_filt"+hemi+".dat"
        if os.path.exists(computed_vortname):
            cmd_do("ln -sf " + permdir + "vor_"+run+"_"+year+"_*filt"+hemi+".dat "     + scrdir, printflag=True)
            cmd_do("ln -sf " + permdir + run+"_"+year+"_10m_wspeed"+hemi+".nc " + scrdir, printflag=True)
            cmd_do("ln -sf " + permdir + run+"_"+year+"_mslp"+hemi+".nc "       + scrdir, printflag=True)
            cmd_do("ln -sf " + permdir + run+"_"+year+"_geo"+hemi+".nc "       + scrdir, printflag=True)
        else:
            for il,lev in enumerate(levs):
                # Create and link full resolution vorticities
                print('Creating vorticities')
                vortname_nc=compute_vort(lev,run,year,scrdir,hemi)
                link_lev_full_vort(vortname_nc)
                filter_lev_vort(lev,run,year,tr)

    # VORTICITY COMPUTATION STEP
    vort_tracking_name=avg_filter_lev_vort(run,year,tr)

    # TRACKING COMPUTATION STEP
    # MASTER TRACKING
    trackfile = master_tracking(run,year,tr,hemi,vort_tracking_name,trackdir,trackdir_full,trackname)
    # DONE TRACKING
else:
    print('skip_to_hart is True: skipping vorticity computation and tracking phase')

#CPL need this repeated, in case we only do tc_ident (or skip straight to Hart)
trackfile=trackdir_full+trackname
print("just after MASTER, the trackfile is: ",trackfile)

# ADD FIELDS STEP
#two_day_file = trackdir_full+'/'+trackname+'.2day'
two_day_file = trackfile+'.2day'
trackfile_added=add_fields(add_names,run,year,tr,trackfile,two_day_file,hemi)

# List here the name of the final file
final_file=trackfile_added+'.hart'

# HART STEP
# The output of this step is:
# EXPT/VOR_VERTAVG_T63filt_YYYY/tr_trs_pos.2day_addvor_addmslp_addwind.hart
hart_phase=Hart(run, hemi, year, final_file)

# with WCSI, the TCIDENT pahse is no longer needed
#   now do the tcident
#trackfile_tcident = trackfile_added+'.tcident'
#tcident(trackfile_tcident,home,trackdir_full,year)

if hemi == '':
    startdate=year+'010100'
elif hemi == '_SH':
    startdate=year+'070100'

# Now inserts the proper dates and convert to NetCDF
final_nc = final_file+'.new.nc'
if not os.path.exists(final_nc):
    cmd_do(home+'utils/bin/count ' + final_file + ' 0 0 5 4 0 ' + startdate + ' 6', printflag=True)
    cmd_do(home+'utils/bin/tr2nc ' + final_file+'.new' + ' s '+home+'utils/TR2NC/tr2nc.meta.addmslp.addwind', printflag=True)
else:
    print('Final NetCDF already exists, skipping count/tr2nc step: '+final_nc)

if not only_tcident and os.path.exists(final_nc):
    cmd_do("rm " + home + "outdat/*"+run+"_"+year+hemi+"*", printflag=True)
    cmd_do("mv " + scrdir + "vor_"+run+"_"+year+"_*filt*.dat " + permdir, printflag=True)
    cmd_do("mv " + scrdir + run+"_"+year+"_10m_wspeed*.nc " + permdir, printflag=True)
    cmd_do("mv " + scrdir + run+"_"+year+"_mslp*.nc "       + permdir, printflag=True)
    cmd_do("mv " + scrdir + run+"_"+year+"_geo*.nc "        + permdir, printflag=True)
    cmd_do("rm " + scrdir + run + "_" + year + "??_*.nc", printflag=True)
    cmd_do("rm " + scrdir + "vor_"+ run + "_" + year+"_???.*", printflag=True)
    cmd_do("rm " + scrdir + run + "_" + year + "??_uwinds*.nc", printflag=True)
    cmd_do("rm " + scrdir + run + "_" + year + "??_vwinds*.nc", printflag=True)
    cmd_do("rm " + home + "work*" + run + "_" + year + "*.in", printflag=True)
    print("Done tracking")
quit()