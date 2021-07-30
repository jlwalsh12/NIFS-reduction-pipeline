"""

This routine reduces the telluric star observations using the Gemini
IRAF reduction tasks (or modified Gemini IRAF reduction tasks).

"""

def nifs_telluric_LP(workdir, caldir, date, flatlist, arclist, ronchilist,
                     telluriclist, skylist, skylistshort, flinter_nsfitcoords,
                     flinter_extract, obs_setup):

    ###########################################################################
    #  STEP 1: Prepare IRAF  		                                      #
    ###########################################################################

    import sys
    import getopt
    import os
    import time
    import shutil
    import astropy.io.fits as pyfits
    import numpy as np
    import subprocess
    from pyraf import iraf
    iraf.gemini(_doprint=0)
    iraf.nifs(_doprint=0)
    iraf.gnirs(_doprint=0)
    iraf.gemtools(_doprint=0)
    iraf.onedspec(_doprint=0)
    from pyraf import iraffunctions
    
    #unlearn the used tasks
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs,
                 'nffixbad_anil')
    iraf.set(stdimage='imt2048')

    #create a log file and back up the previous one if it already exists
    log = 'telluric_'+date+'.log'
    if os.path.exists(log):
        t = time.localtime()
        app = '_'+str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)+'_'+ \
        str(t[3]).zfill(2)+':'+str(t[4]).zfill(2)+':'+str(t[5]).zfill(2)
        shutil.move(log,log+app)

    #change to the workdir within pyraf
    iraffunctions.chdir(workdir)
        
    #prepare the package for NIFS
    iraf.nsheaders('nifs',logfile=log)

    #set clobber to 'yes' for the script. this still does not make the gemini
    #tasks overwrite files, so you will likely have to remove files if you
    #re-run the script.
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')

    #get the file names of the reduced flat, ronci mask, and arcs
    calflat=str(open(caldir+flatlist, 'r').readlines()[0]).strip()
    ronchiflat=str(open(caldir+ronchilist, 'r').readlines()[0]).strip()
    arc = []
    for i in range(len(arclist)):
        arc.append(str(open(caldir+arclist[i], 'r').readlines()[0]).strip())
    #use the first telluric frame as the base name for the combined telluric
    #spectrum
    telluric = []
    for i in range(len(telluriclist)):
        telluric.append(str(open(telluriclist[i], 'r').readlines()[0]).strip())

    
    ############################################################################
    # STEP 2:  Get the Calibrations for the Reduction                          #
    ############################################################################

    #copy required files and transformation database into the current
    #working directory
    iraf.copy(caldir+'rgn'+ronchiflat+'.fits',output='./')
    iraf.copy(caldir+'rgn'+calflat+'_sflat_bpm.pl',output='./')
    iraf.copy(caldir+'rgn'+calflat+'_flat.fits',output='./')
    for i in range(len(arc)):
        iraf.copy(caldir+'wrgn'+arc[i]+'.fits',output='./')
    if not os.path.isdir('./database'):
        os.mkdir('./database/')
    iraf.copy(caldir+'database/*',output='./database/')


    ###########################################################################
    # STEP 3:  Reduce the Telluric Standard                                   #
    ###########################################################################

    for i in range(len(telluriclist)):
    
        #prepare the data
        iraf.nfprepare('@'+telluriclist[i], rawpath='',
                       shiftim=caldir+'s'+calflat,
                       bpm='rgn'+calflat+'_sflat_bpm.pl', fl_vardq='yes',
                       fl_int='yes', fl_corr='no', fl_nonl='no', logfile=log)

        iraf.nfprepare('@'+skylistshort[i], rawpath='',
                       shiftim=caldir+'s'+calflat,
                       bpm='rgn'+calflat+'_sflat_bpm.pl', fl_vardq='yes',
                       fl_int='yes', fl_corr='no', fl_nonl='no', logfile=log)

        #do the sky subtraction on all the individual frames. read the
        #list and get rid of '\n' character returns first.
        telluricexps=open(telluriclist[i], 'r').readlines()
        telluricexps=[word.strip() for word in telluricexps]
        skyexps=open(skylist[i], 'r').readlines()
        skyexps=[word.strip() for word in skyexps]
        for j in range(len(telluricexps)):
            iraf.gemarith('n'+telluricexps[j], '-', 'n'+skyexps[j],
                          'sn'+telluricexps[j], fl_vardq='yes', logfile=log)

        #cut the slices and flat field the telluric data
        
        #use different calls to nsreduce depending on the
        #observational setup. the pipeline will have already quit if
        #one of these two setups were not used.
        if obs_setup == 'hk_2.20':
            iraf.nsreduce('sn@'+telluriclist[i], outpref='r',
                          flatim='rgn'+calflat+'_flat', fl_cut='yes',
                          fl_nsappw='no', fl_vardq='yes', fl_sky='no',
                          fl_dark='no', fl_flat='yes', logfile=log)
        if obs_setup == 'hk_2.30':
            iraf.nsreduce('sn@'+telluriclist[i], outpref='r',
                          flatim='rgn'+calflat+'_flat', fl_cut='yes',
                          fl_nsappw='no', fl_vardq='yes', fl_sky='no',
                          fl_dark='no', fl_flat='yes', crval=23000.000000,
                          cdelt=-2.115000, logfile=log)

        #fix bad pixels from the DQ plane
        iraf.nffixbad_anil('rsn@'+telluriclist[i], outpref='b', logfile=log)

        #read in the arc2uselist file to get the arc to use for this
        #telluric set (arc2uselist will only have 1 entry)
        arc2use=str(open('arc2uselist_'+str(i+1), 'r').readlines()[0]).strip()

        #derive the 2D to 3D spatial/spectral transformation. for some
        #reason, setting lyorder=5 below shows up as lyorder=4 when
        #running an interactive fit (and setting lyorder=4 below comes
        #up as lyorder=3 during the interactive fit). want lyorder=4,
        #so keeping lyorder=5 set below.
        iraf.nsfitcoords('brsn@'+telluriclist[i], outpref='f',
                         fl_int=flinter_nsfitcoords,
                         lamptr=arc2use, sdisttr='rgn'+ronchiflat,
                         logfile=log, lxorder=4, lyorder=5, sxorder=4,
                         syorder=4)

        #apply the transformation determined in the nffitcoords step
        iraf.nstransform('fbrsn@'+telluriclist[i], outpref='t', logfile=log)

        tmp = open(telluriclist[i], 'r').readlines()
        for j in range(len(tmp)):
            telluric_now = \
              str(open(telluriclist[i], 'r').readlines()[j]).strip()
        
            #in order to extract a 1D spectrum from the 2D data
            #non-interactively, create a temporary cube
            iraf.nifcube('tfbrsn'+telluric_now, outcubes='tmp_tfbrsn'+\
                         telluric_now, sscale=0.043, verbose='no')
        
            tmp_fitsfile = pyfits.open('tmp_tfbrsn'+telluric_now+'.fits')
            tmp_tel_cube = tmp_fitsfile[1].data
            #tmp_tel_cube is in the form: wavelength, y, x. sum along
            #wavelength to get collapsed image of telluric star
            tmp_tel_image = np.sum(tmp_tel_cube,axis=0)
            #determine the y and x values of the maximum
            max_tmp_tel_image = tmp_tel_image.argmax()
            k,m = np.unravel_index(tmp_tel_image.argmax(), tmp_tel_image.shape)
            tmp_fitsfile.close()

            #the temporary cube doesn't have the same size as
            #nfextract (x has a length f 29 and y a length of
            #69). also iraf starts indexing a 1 not 0. calculate the x
            #and y center values to feed into nfextract here.
            xsize_tel_img = tmp_tel_image.shape[1]
            ysize_tel_img = tmp_tel_image.shape[0]
            x_nfextract = round( (m / (xsize_tel_img/29.)) + 1., 2)
            y_nfextract = round( (k / (ysize_tel_img/69.)) + 1., 2)
        
            #delete the temporary cube and the nifs.log file that was
            #automatically created
            tmp=subprocess.call(['rm',workdir+'tmp_tfbrsn'+\
                                 telluric_now+'.fits'],
                                stderr=open(os.devnull,'w'))
            tmp=subprocess.call(['rm',workdir+'nifs.log'],
                                stderr=open(os.devnull,'w'))

            if flinter_extract == 'no':
                #extract 1D spectra from the 2D data
                iraf.nfextract('tfbrsn'+telluric_now, outpref='x',
                               diameter=0.5, fl_int='no', xc=x_nfextract,
                               yc=y_nfextract, logfile=log)
                
            if flinter_extract == 'yes':
                #print the (x,y) pixel with the maximum value to the
                #screen
                print('The pixel with the maximum value '+\
                      'is: x=%1.1f y=%1.1f.' % (x_nfextract, y_nfextract))
                #extract 1D spectra from the 2D data. need to have ds9
                #open.
                iraf.nfextract('tfbrsn'+telluric_now, outpref='x',
                               diameter=0.5, fl_int='yes', logfile=log)
                
            #for flux calibration later on will want a cube of one
            #telluric star
            if i == 0 and j == 0:
                iraf.nifcube('tfbrsn'+telluric_now, outcubes='ctfbrsn'+\
                             telluric_now, sscale=0.05, verbose='no')

        #combine all the 1D spectra to one final output file
        iraf.gemcombine('xtfbrsn//@'+telluriclist[i], output='gxtfbrsn'+\
                        telluric[i], statsec='[*]', combine='median',
                        logfile=log, masktype='none', fl_vardq='yes')

    
    ###########################################################################
    # Reset to user defaults                                                  #
    ###########################################################################
    if user_clobber == "no":
        iraf.set(clobber='no')

###########################################################################
#          End of the Telluric Calibration Data Reduction                 #
#                                                                         #
# The output of this reduction script is a 1D spectrum used for           #
# telluric calibration of NIFS science data. Also created is a datacube   #
# of one of the telluric stars, which is used later during the flux       #
# calibration step.                                                       #
#                                                                         #
# The output files are: gxtfbrsn + xxx --> combined 1D telluric           #
#                       ctfbrsn + xxx --> cube of a telluric star         #
#                                                                         #
# The file prefixes are described below.                                  #
#                                                                         #
# g = gemcombined/gemarithed   n=nfprepared  s=skysubtracted              #
# r=nsreduced  b = bad pixel corrected  f= run through nffitcoords        # 
# t = nftransformed   x = extracted to a 1D spectrum                      #
#                                                                         #
###########################################################################
