"""

This routine reduces the PSF star observations using the Gemini IRAF
reduction tasks (or modified Gemini IRAF reduction tasks).

"""

def nifs_psf_LP(workdir, caldir, teldir, date, flatlist, arclist, ronchilist,
                telcorrlist, galname, psflist, skylist, skylistshort,
                sigclip_in, sigfrac_in, objlim_in, nsigneg_in,
                rootdir, reducedir, obs_setup, flinter_nsfitcoords,
                flinter_telluric):

                
    ###########################################################################
    #  STEP 1: Prepare IRAF  		                                          #
    ###########################################################################

    #import some useful python utilities
    import os
    import sys
    import time
    import shutil
    import astropy.io.fits as pyfits
    import numpy as np
    import subprocess
    from nifs_lacosmic_LP import nifs_lacosmic_LP
    from nifs_xzap_LP import nifs_xzap_LP
    #import the pyraf module and relevant packages
    from pyraf import iraf
    iraf.gemini(_doprint=0)
    iraf.nifs(_doprint=0)
    iraf.gnirs(_doprint=0)
    iraf.gemtools(_doprint=0)
    from pyraf import iraffunctions

    #unlearn the used tasks
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs,\
                 'nffixbad_anil','nftelluric_anil','nftelluric_aniljonelle')
    iraf.set(stdimage='imt2048')

    #create a log file and back up the previous one if it already exists
    log = 'psf_'+galname+'_'+date+'.log'
    if os.path.exists(log):
        t = time.localtime()
        app = "_"+str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)+'_'+ \
        str(t[3]).zfill(2)+':'+str(t[4]).zfill(2)+':'+str(t[5]).zfill(2)
        shutil.move(log,log+app)

    #change to the workdir within pyraf
    iraffunctions.chdir(workdir)
        
    #prepare the package for NIFS
    iraf.nsheaders('nifs',logfile=log)

    #set clobber to 'yes' for the script. this still does not make the
    #gemini tasks overwrite files, so you will likely have to remove
    #files if you re-run the script.
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')

    #get the file names of the reduced flat, ronci mask, and arc
    calflat=str(open(caldir+flatlist, 'r').readlines()[0]).strip()
    ronchiflat=str(open(caldir+ronchilist, 'r').readlines()[0]).strip()
    arc = []
    for i in range(len(arclist)):
        arc.append(str(open(caldir+arclist[i], 'r').readlines()[0]).strip())

    
    ###########################################################################
    # STEP 2:  Get the Calibrations for the Reduction                         #
    ###########################################################################

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
    # STEP 3:  Reduce the Science Data                                        #
    ###########################################################################
    
    iraf.nfprepare('@'+psflist, rawpath='', shiftimage=caldir+'s'+calflat,
                    fl_vardq='yes', bpm='rgn'+calflat+'_sflat_bpm.pl',
                    logfile=log)
    
    iraf.nfprepare('@'+skylistshort, rawpath='', shiftimage=caldir+'s'+calflat,
                    fl_vardq='yes', bpm='rgn'+calflat+'_sflat_bpm.pl',
                    logfile=log)

    #read in the frame lists (removing '\n' line breaks from the strings)
    psfexps=open(psflist, 'r').readlines()
    psfexps=[word.strip() for word in psfexps]
    skyexps=open(skylist, 'r').readlines()
    skyexps=[word.strip() for word in skyexps]
    #subtract the sky
    for i in range(len(psfexps)):
        print('n%s, n%s'% (psfexps[i],skyexps[i]))
        iraf.gemarith('n'+psfexps[i], '-', 'n'+skyexps[i], 'sn'+psfexps[i],
                       fl_vardq='yes', logfile=log)

    #flat field and cut the psf
    iraf.nsreduce('sn@'+psflist, fl_cut='yes', fl_nsappw='yes', fl_dark='no',
                   fl_sky='no', fl_flat='yes', flatimage='rgn'+calflat+'_flat',
                   fl_vardq='yes', logfile=log)
    
    #interpolate over bad pixels flagged in the DQ plane for the psf
    iraf.nffixbad_anil('rsn@'+psflist,logfile=log)

    #do further cleaning
    
    #start with LA COSMIC to identify and fixpix to correct the cosmic rays

    #need gain and readnoise from header to pass into LA COSMIC
    tmp_file = pyfits.open('n'+psfexps[i]+'.fits')
    gain_in = tmp_file[0].header['GAIN']
    readn_in = tmp_file[0].header['RDNOISE']
    tmp_file.close()
    xorder_in = 4         #LA COSMIC xorder
    yorder_in = 3         #LA COSMIC yorder
    nifs_lacosmic_LP(workdir, psfexps, gain_in, readn_in, xorder_in,
                     yorder_in, sigclip_in, sigfrac_in, objlim_in)
    
    #LA COSMIC isn't great at removing negative bad pixels, so using
    #xzap to identify negative pixels and fixpix to correct
    nifs_xzap_LP(workdir, psfexps, nsigneg_in)

    #read in the arc2uselist file to get the arc to use for this psf
    #set (arc2uselist will only have 1 entry)
    arc2use=str(open('arc2uselist', 'r').readlines()[0]).strip()


    #derive the 2D to 3D spatial/spectral transformation. for some
    #reason, setting lyorder=5 below shows up as lyorder=4 when
    #running an interactive fit (and setting lyorder=4 below comes up
    #as lyorder=3 during the interactive fit). want lyorder=4, so
    #keeping lyorder=5 set below.
    iraf.nsfitcoords('brsn@'+psflist, outpref='f',
                     fl_int=flinter_nsfitcoords, lamptr=arc2use,
                     sdisttr='rgn'+ronchiflat, logfile=log, lxorder=4,
                     lyorder=5, sxorder=4, syorder=4)
    
    #apply the transformation determined in the nffitcoords step for
    #the psf
    iraf.nstransform('fbrsn@'+psflist, outpref='t', logfile=log)
    
    #read in the telluric star list
    telcorr=open(telcorrlist, 'r').readlines()
    telcorr=[word.strip() for word in telcorr]
    
    for i in range(len(psfexps)):
    
        #copy the telluric star to this directory (iraf will truncate
        #long path names and the file may not be found).
        tmp=subprocess.call(['cp',teldir+'/'+telcorr[i],
                                 workdir+(telcorr[i].split('/'))[1]],
                                 stderr=open(os.devnull,'w'))
        
        #in order to extract a 1D spectrum from the 2D data
        #non-interactively, create a temporary cube
        iraf.nifcube('tfbrsn'+psfexps[i], outcubes='tmp_tfbrsn'+psfexps[i],
                     sscale=0.043, verbose='no', logfile=log)
        
        tmp_fitsfile = pyfits.open('tmp_tfbrsn'+psfexps[i]+'.fits')
        tmp_psf_cube = tmp_fitsfile[1].data
        #tmp_psf_cube is in the form: wavelength, y, x. sum along
        #wavelength to get collapsed image of telluric star
        tmp_psf_image = np.sum(tmp_psf_cube,axis=0)
        #determine the y and x values of the maximum
        max_tmp_psf_image = tmp_psf_image.argmax()
        j,k = np.unravel_index(tmp_psf_image.argmax(), tmp_psf_image.shape)
        tmp_fitsfile.close()
        
        #the temporary cube doesn't have the same size as nfextract (x
        #has a length f 29 and y a length of 69). also iraf starts
        #indexing a 1 not 0. calculate the x and y center values to
        #feed into nfextract here.
        xsize_psf_img = tmp_psf_image.shape[1]
        ysize_psf_img = tmp_psf_image.shape[0]
        x_nfextract = round( (k / (xsize_psf_img/29.)) + 1., 2)
        y_nfextract = round( (j / (ysize_psf_img/69.)) + 1., 2)
        
        #delete the temporary cube and the nifs.log file that was
        #automatically created
        tmp=subprocess.call(['rm',workdir+'tmp_tfbrsn'+psfexps[i]+'.fits'],
                            stderr=open(os.devnull,'w'))
        tmp=subprocess.call(['rm',workdir+'nifs.log'],
                            stderr=open(os.devnull,'w'))

        #correct the data for telluric absorption features
        
        #the telluric correction. below the extraction of a 1D psf
        #spectrum and the determination of the shift/scale parameters
        #(applied to the correction spectrum) are not done
        #interactively. for the 1D psf spectrum extraction, use an
        #aperture centered on the spaxel with the most
        #counts. fl_flux='no' means no multiplication by a blackbody
        #is done (as the telluric already has the blackbody shape
        #removed).
        if flinter_telluric == 'no':
            iraf.nftelluric_aniljonelle('tfbrsn'+psfexps[i],
                                        (telcorr[i].split('/'))[1],
                                        fl_inter='no', fl_flux='no',
                                        sample='23000:24100',
                                        xc=x_nfextract, yc=y_nfextract,
                                        logfile=log)

        #the telluric correction. below the extraction of a 1D psf
        #spectrum and the determination of the shift/scale parameters
        #(applied to the correction spectrum) are done
        #interactively. for the 1D psf spectrum extraction, the (x,y)
        #location of the spaxel with the most counts is printed to the
        #screen to guide the user. fl_flux='no' means no
        #multiplication by a blackbody is done (as the telluric
        #already has the blackbody shape removed). note, need to have
        #a DS9 window open.
        if flinter_telluric == 'yes':
            
            #print the (x,y) pixel with the maximum value to the screen
            print('The pixel with the maximum value '+\
                  'is: x=%1.1f y=%1.1f' % (x_nfextract,y_nfextract))
              
            iraf.nftelluric_anil('tfbrsn'+psfexps[i],(telcorr[i].split('/'))[1],
                                 fl_inter='yes', fl_flux='no',
                                 sample='23000:24100', logfile=log)
        
        #remove the copy of the telluric star now that the psf
        #exposure has been reduced
        tmp=subprocess.call(['rm',workdir+(telcorr[i].split('/'))[1]],
                            stderr=open(os.devnull,'w'))

        #create the 3D cube of the psf exposures only for the purposes
        #of collapsing and figuring out spatial offsets. therefore,
        #just use gemini's nifcube for now even though the variance/dq
        #extensions are not properly dealt with.
        print('')
        print('Creating temporary cube of psf exposure and '+\
                  'collapsing to form an image. Please be patient.')
        print('')
        iraf.nifcube(inimages='atfbrsn'+psfexps[i], sscale=0.05,
                     verbose='no', logfile=log)
                         
        #collapse the psf cube to get an image
        tmp_cube = pyfits.open('catfbrsn'+psfexps[i]+'.fits')
        tmp_psf_cube = tmp_cube[1].data
        #tmp_psf_cube is in the form: wavelength, y, x. sum along
        #wavelength to get collapsed image of galaxy
        tmp_psf_image = np.sum(tmp_psf_cube,axis=0)
        #will use the primary header (with wcs info) from the psf
        #cube and make some adjustments
        hdu = pyfits.PrimaryHDU(tmp_psf_image)
        hdu.header = tmp_cube[0].header
        hdu.header['BITPIX'] = -32
        hdu.header['NAXIS'] = 2
        hdu.header['NAXIS1'] = (np.shape(tmp_psf_image))[1]
        hdu.header['NAXIS2'] = (np.shape(tmp_psf_image))[0]
        hdu.header['PIXSCALE'] = 0.05
        #write to the fits file
        hdu.writeto('image_catfbrsn'+psfexps[i]+'.fits',\
                    output_verify='silentfix')
        tmp_cube.close()
        #remove the temporary cube
        tmp=subprocess.call(['rm',workdir+'catfbrsn'+psfexps[i]+'.fits'],
                            stderr=open(os.devnull,'w'))

        #copy the collapsed psf cubes to the merged directory. don't
        #simply move the files (like was done for the galaxy
        #reduction) because we will also create a merged psf cube each
        #night as well.
        tmp=subprocess.call(['cp','image_catfbrsn'+psfexps[i]+'.fits',
                             rootdir+reducedir+'psfs/merged/'+obs_setup+\
                             '/psf_'+galname],stderr=open(os.devnull,'w'))
        #copy the telluric corrected psf exposures to the merged
        #directory and make a copy for header modification both in the
        #current and merged directories.
        tmp=subprocess.call(['cp',workdir+'atfbrsn'+psfexps[i]+'.fits',
                             rootdir+reducedir+'psfs/merged/'+obs_setup+\
                             '/psf_'+galname],stderr=open(os.devnull,'w'))
        tmp=subprocess.call(['cp',workdir+'atfbrsn'+psfexps[i]+'.fits',
                             rootdir+reducedir+'psfs/merged/'+obs_setup+\
                             '/psf_'+galname+'/hatfbrsn'+psfexps[i]+'.fits'],
                             stderr=open(os.devnull,'w'))
        tmp=subprocess.call(['cp',workdir+'atfbrsn'+psfexps[i]+'.fits',
                             workdir+'hatfbrsn'+psfexps[i]+'.fits'],
                             stderr=open(os.devnull,'w'))

        
    ###########################################################################
    # Reset to user defaults                                                  #
    ###########################################################################
    if user_clobber == "no":
        iraf.set(clobber='no')

        
###########################################################################
#                    End of the PSF Data Reduction                        #
#                                                                         #
# The output of this reduction script are prepared, sky-subtracted,       #
# reduced, bad pixel corrected, geometrically rectified, telluric         #
# corected PSF slices, ready to be merged and turned into a final data    #
# cube. Also, a temporary cube of each PSF exposure is created and summed #
# along the wavelength axis to produce an image. The images will later be #
# used to determine the spatial offsets between individual psf exposures. #
# As for the merging of psf exposures, we will do it for each individual  #
# night and then all nights together (in a "../merged" directory. In the  #
# case of this reduction, the final output files are called: atfbrsn+psf. #
#                                                                         #
# The meaning of the output prefixes are described below:                 #
#                                                                         #
# n=nfprepared   s=sky subtracted   r=nsreduced                           #
# b = bad pixel corrected /LA COSMIC cleaned/ XZAP cleaned                #
# f= run through nffitcoords   t = nftransformed                          #
# a = corrected for telluric absorption features                          #
# some baseline calibs are copied to the galaxy working directory as well #
#                                                                         #
###########################################################################
