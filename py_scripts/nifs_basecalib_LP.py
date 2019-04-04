"""

This routine reduces the baseline calibrations using the Gemini IRAF
reduction tasks (or modified Gemini IRAF reduction tasks).

"""

def nifs_basecalib_LP(workdir, date, flatlist, flatdarklist, arclist,
                      arcdarklist, ronchilist, ronchidarklist, darklist,
                      refdir, flinter_arc, flinter_sdist):

    ###########################################################################
    #  STEP 1: Prepare IRAF  		                                          #
    ###########################################################################

    import sys
    import getopt
    import os
    import time
    import shutil
    from pyraf import iraf
    iraf.gemini(_doprint=0)
    iraf.nifs(_doprint=0)
    iraf.gnirs(_doprint=0)
    iraf.gemtools(_doprint=0)
    from pyraf import iraffunctions
    import astropy.io.fits as pyfits

    #unlearn the used tasks
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs,
                 'nfsdist_jonelle')
    iraf.set(stdimage='imt2048')
    
    #create a log file and back up the previous one if it already
    #exists
    log = 'basecalib_'+date+'.log'
    if os.path.exists(log):
        t = time.localtime()
        app = '_'+str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)+'_'+ \
        str(t[3]).zfill(2)+':'+str(t[4]).zfill(2)+':'+str(t[5]).zfill(2)
        shutil.move(log,log+app)

    #change to the workdir within pyraf
    iraffunctions.chdir(workdir)

    #prepare the package for NIFS
    iraf.nsheaders('nifs',logfile=log)

    #set clobber to 'yes' for the script. this still does not make the
    #gemini tasks overwrite files, so you will likely have to remove
    #files if you re-run the script.
    user_clobber=iraf.envget('clobber')
    iraf.reset(clobber='yes')

    #set the file names for the main calibration outputs. just use the
    #first name in the list of relevant files, which we get from the
    #file lists
    calflat=str(open(flatlist, 'r').readlines()[0]).strip()
    flatdark=str(open(flatdarklist, 'r').readlines()[0]).strip()
    arc = []
    for i in range(len(arclist)):
        arc.append(str(open(arclist[i], 'r').readlines()[0]).strip())
    arcdark=str(open(arcdarklist, 'r').readlines()[0]).strip()
    ronchiflat=str(open(ronchilist, 'r').readlines()[0]).strip()
    ronchiflatdark=str(open(ronchidarklist, 'r').readlines()[0]).strip()
    if len(darklist) > 0:
        darkobj = []
        for i in range(len(darklist)):
            darkobj.append(str(open(darklist[i], 'r').readlines()[0]).strip())

        
    ###########################################################################
    #  STEP 2: Determine the shift to the MDF file  		                  #
    ###########################################################################
    
    iraf.nfprepare(calflat, rawpath='', outpref='s', shiftx='INDEF',
                   shifty='INDEF', fl_vardq='no', fl_corr='no', fl_nonl='no',
                   logfile=log)
    
    
    ###########################################################################
    #  STEP 3: Make the Flat Field and BPM  		                          #
    ###########################################################################
        
    iraf.nfprepare('@'+flatlist, rawpath='', shiftim='s'+calflat,
                    fl_vardq='yes', fl_inter='yes', fl_corr='no',
                    fl_nonl='no', logfile=log)
                   
    iraf.nfprepare('@'+flatdarklist, rawpath='', shiftim='s'+calflat,
                    fl_vardq='yes', fl_inter='yes', fl_corr='no',
                    fl_nonl='no', logfile=log)
    
    iraf.gemcombine('n//@'+flatlist, output='gn'+calflat, fl_dqpr='yes',
                    fl_vardq='yes', masktype='none', logfile=log)
    iraf.gemcombine('n//@'+flatdarklist, output='gn'+flatdark, fl_dqpr='yes',
                    fl_vardq='yes', masktype='none', logfile=log)
    
    iraf.nsreduce('gn'+calflat, fl_cut='yes', fl_nsappw='yes', fl_vardq='yes',
                   fl_sky='no', fl_dark='no', fl_flat='no', logfile=log)
    iraf.nsreduce('gn'+flatdark, fl_cut='yes', fl_nsappw='yes', fl_vardq='yes',
                   fl_sky='no', fl_dark='no', fl_flat='no', logfile=log)
    
    #creating flat image, final name = rgnN....._sflat.fits
    iraf.nsflat('rgn'+calflat, darks='rgn'+flatdark,
                 flatfile='rgn'+calflat+'_sflat', darkfile='rgn'+flatdark+\
                 '_dark', fl_save_dark='yes', process='fit', thr_flo=0.15,
                 thr_fup=1.55, fl_vardq='yes', logfile=log)
    
    #rectify the flat for slit function differences - make the final flat
    iraf.nsslitfunction('rgn'+calflat, 'rgn'+calflat+'_flat',
                         flat='rgn'+calflat+'_sflat',
                         dark='rgn'+flatdark+'_dark', combine='median', order=3,
                         fl_vary='no', logfile=log)
    
    
    ###########################################################################
    # STEP 4: Reduce the Object Darks If They Exist                           #
    ###########################################################################
    
    if len(darklist) > 0:
        for i in range(len(darklist)):
            iraf.nfprepare('@'+darklist[i], rawpath='', shiftimage='s'+calflat,
                            bpm='rgn'+calflat+'_sflat_bpm.pl', fl_vardq='yes',
                            fl_corr='no', fl_nonl='no', logfile=log)
            if len(darklist[i]) > 1:
                iraf.gemcombine('n//@'+darklist[i], output='gn'+darkobj[i],
                                fl_dqpr='yes', fl_vardq='yes', masktype='none',
                                logfile=log)
            else:
                iraf.copy('n'+darkobj[i]+'.fits','gn'+darkobj[i]+'.fits')
    
    
    ###########################################################################
    # STEP 5: Reduce the Arcs and determine the wavelength solution           #
    ###########################################################################
    
    iraf.nfprepare('@'+arcdarklist, rawpath='', shiftimage='s'+calflat,
                    bpm='rgn'+calflat+'_sflat_bpm.pl', fl_vardq='yes',
                    fl_corr='no', fl_nonl='no', logfile=log)
    
    nfiles = len(open(arcdarklist).readlines())
    if nfiles > 1:
        iraf.gemcombine('n//@'+arcdarklist, output='gn'+arcdark,
                         fl_dqpr='yes', fl_vardq='yes', masktype='none',
                         logfile=log)
    else:
        iraf.copy('n'+arcdark+'.fits','gn'+arcdark+'.fits')
            
    #account for multiple arclists
    for i in range(len(arclist)):
        iraf.nfprepare('@'+arclist[i], rawpath='', shiftimage='s'+calflat,
                        bpm='rgn'+calflat+'_sflat_bpm.pl', fl_vardq='yes',
                        fl_corr='no', fl_nonl='no', logfile=log)
    
        nfiles = len(open(arclist[i]).readlines())
        if nfiles > 1:
            iraf.gemcombine('n//@'+arclist[i], output='gn'+arc[i],
                            fl_dqpr='yes', fl_vardq='yes', masktype='none',
                            logfile=log)
        else:
            iraf.copy('n'+arc[i]+'.fits','gn'+arc[i]+'.fits')
    
        iraf.nsreduce('gn'+arc[i], outpr='r', darki='gn'+arcdark,
                      flati='rgn'+calflat+'_flat', fl_vardq='no', fl_cut='yes',
                      fl_nsappw='yes', fl_sky='no', fl_dark='yes',fl_flat='yes',
                      logfile=log)
    
        #determine the wavelength of the observation and set the arc
        #coordinate file. if the user wishes to change the coordinate
        #file to a different one, they need only to change the "clist"
        #variable to their line list in the coordli= parameter in the
        #nswavelength call.
    
        hdulistobj = pyfits.open('rgn'+arc[i]+'.fits')
        bandobj = hdulistobj[0].header['GRATING'][0:1]
    
        #only tested for K
        if bandobj == "Z":
            clistobj="nifs$data/ArXe_Z.dat"
            my_threshobj=100.0
        elif bandobj == "K":
            clistobj='anil_ArXe_K.dat'
            my_threshobj=50.0
        else:
            clistobj="gnirs$data/argon.dat"
            my_threshobj=100.0    
    
        iraf.nswavelength('rgn'+arc[i], coordli=clistobj, nsum=10,
                          thresho=my_threshobj, fwidth=2.0, match=-6,
                          cradius=8.0, fl_inter=flinter_arc, nfound=10,
                          nlost=10, logfile=log, order=4)

    
    ############################################################################
    # STEP 6: Trace the spatial curvature/spectral distortion in the Ronchi    #
    # flat                                                                     #
    ############################################################################

    iraf.nfprepare('@'+ronchilist, rawpath='', shiftimage='s'+calflat,
                    bpm='rgn'+calflat+'_sflat_bpm.pl', fl_vardq='yes',
                    fl_corr='no', fl_nonl='no', logfile=log)
    iraf.nfprepare('@'+ronchidarklist, rawpath='', shiftimage='s'+calflat,
                    bpm='rgn'+calflat+'_sflat_bpm.pl', fl_vardq='yes',
                    fl_corr='no', fl_nonl='no', logfile=log)

    #determine the number of input Ronchi calibration mask files and
    #Ronchi dark files so that the routine runs automatically for
    #single or multiple files.
    nfiles = len(open(ronchilist).readlines())
    if nfiles > 1:
        iraf.gemcombine('n//@'+ronchilist, output='gn'+ronchiflat,
                        fl_dqpr='yes', masktype='none', fl_vardq='yes',
                        logfile=log)
    else:
        iraf.copy('n'+ronchiflat+'.fits','gn'+ronchiflat+'.fits')

    nfiles = len(open(ronchidarklist).readlines())
    if nfiles > 1:
        iraf.gemcombine('n//@'+ronchidarklist, output='gn'+ronchiflatdark,
                        fl_dqpr='yes', masktype='none', fl_vardq='yes',
                        logfile=log)
    else:
        iraf.copy('n'+ronchiflatdark+'.fits','gn'+ronchiflatdark+'.fits')

    iraf.nsreduce('gn'+ronchiflat, outpref='r', dark='gn'+ronchiflatdark,
                  flatimage='rgn'+calflat+'_flat', fl_cut='yes',
                  fl_nsappw='yes', fl_flat='yes', fl_sky='no', fl_dark='yes',
                  fl_vardq='no', logfile=log)

    iraf.nfsdist_jonelle('rgn'+ronchiflat, fwidth=6.0, cradius=8.0, glshift=2.8,
                         minsep=6.5, thresh=2000.0, nlost=3,
                         fl_inter=flinter_sdist, logfile=log, order=4)

    #apply result from nfsdist above to the (average) ronchi mask as a
    #check

    #derive the 2D to 3D spatial/spectral transformation. for some reason,
    #setting lyorder=5 below shows up as lyorder=4 when running an interactive
    #fit (and setting lyorder=4 below comes up as lyorder=3 during the
    #interactive fit). want lyorder=4, so keeping lyorder=5 set below.
    iraf.nsfitcoords('rgn'+ronchiflat, outpref='f', fl_int='no',
                     lamptr='wrgn'+arc[0], sdisttr='rgn'+ronchiflat,
                     logfile=log, lxorder=4, lyorder=5, sxorder=4,
                     syorder=4)

    #apply the transformation determined in the nffitcoords step
    iraf.nstransform('frgn'+ronchiflat, outpref='t', logfile=log)

    #create a cube
    iraf.nifcube('tfrgn'+ronchiflat, outcubes='cube_tfrgn'+ronchiflat,
                 sscale=0.05, verbose='no')

                  
    ###########################################################################
    # Reset to user defaults                                                  #
    ###########################################################################
    if user_clobber == "no":
        iraf.set(clobber='no')


###########################################################################
# End of the Baseline Calibration reduction                               #
###########################################################################
#	                                                                      #
#  The final output files created from this script for later science      #
#  reduction have prefixes and file names of:                             #
#     1. Shift reference file:  "s"+calflat                               #
#     2. Flat field:  "rn"+calflat+"_flat"                                #
#     3. Flat BPM (for DQ plane generation):  "rn"+calflat+"_flat_bpm.pl" #
#     4. Combined, prepared object dark file: "gn"+darkobj                # 
#     5. Wavelength referenced Arc:  "wrgn"+arc                           #
#     6. Spatially referenced Ronchi Flat:  "rgn"+ronchiflat              #
#     7. Database directory with wavelength solution and spatial          #
#        rectification information for each slice.                        #
#	                                                                      #
###########################################################################
