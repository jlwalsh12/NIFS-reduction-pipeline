"""

This routine reduces the galaxy observations using the Gemini IRAF
reduction tasks (or modified Gemini IRAF reduction tasks).

"""

def nifs_galaxy_LP(workdir, caldir, teldir, date, flatlist, arclist,
                   ronchilist, darklist, telcorrlist, galname, gallist,
                   skylist, skylistshort, sigclip_in, sigfrac_in,
                   objlim_in, nsigneg_in, rootdir, reducedir, obs_setup,
                   flinter_nsfitcoords, flinter_telluric):

    ###########################################################################
    #  STEP 1: Prepare IRAF  		                                      #
    ###########################################################################

    import os
    import sys
    import time
    import shutil
    import astropy.io.fits as pyfits
    import numpy as np
    import subprocess
    from nifs_lacosmic_LP import nifs_lacosmic_LP
    from nifs_xzap_LP import nifs_xzap_LP
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
    log = galname+'_'+date+'.log'
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

    #get the file names of the reduced flat, ronci mask, arc, and dark
    calflat=str(open(caldir+flatlist, 'r').readlines()[0]).strip()
    ronchiflat=str(open(caldir+ronchilist, 'r').readlines()[0]).strip()
    arc = []
    for i in range(len(arclist)):
        arc.append(str(open(caldir+arclist[i], 'r').readlines()[0]).strip())
    if len(darklist) > 0:
        dark = []
        for i in range(len(darklist)):
            dark.append(str(open(caldir+darklist[i],
                                 'r').readlines()[0]).strip())

    
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
    
    iraf.nfprepare('@'+gallist, rawpath='', shiftimage=caldir+'s'+calflat,
                   fl_vardq='yes', bpm='rgn'+calflat+'_sflat_bpm.pl',
                   logfile=log)
    
    iraf.nfprepare('@'+skylistshort, rawpath='', shiftimage=caldir+'s'+calflat,
                   fl_vardq='yes', bpm='rgn'+calflat+'_sflat_bpm.pl',
                   logfile=log)

    #read in the frame lists (removing '\n' line breaks from the strings)
    galexps=open(gallist, 'r').readlines()
    galexps=[word.strip() for word in galexps]
    skyexps=open(skylist, 'r').readlines()
    skyexps=[word.strip() for word in skyexps]
    for i in range(len(galexps)):
        print('n%s, n%s'% (galexps[i],skyexps[i]))
        iraf.gemarith('n'+galexps[i], '-', 'n'+skyexps[i], 'sn'+galexps[i],
                       fl_vardq='yes', logfile=log)

    #want to also reduce the sky exposures, so subtract off the
    #dark. need to find darks with the same exposure time as the sky
    #exposures.
    if len(darklist) > 0:
        
        skyshortexps=open(skylistshort, 'r').readlines()
        skyshortexps=[word.strip() for word in skyshortexps]
        sky_exptime = \
          ((pyfits.open(skyshortexps[0]+'.fits'))[0].header)['EXPTIME']
        
        for i in range(len(dark)):
            dark_exptime = \
              ((pyfits.open(caldir+'/'+dark[i]+'.fits'))[0].header)['EXPTIME']
            if sky_exptime == dark_exptime:
                for j in range(len(skyshortexps)):
                    iraf.gemarith('n'+skyshortexps[j], '-',
                                  caldir+'/'+'gn'+dark[i],
                                  'sn'+skyshortexps[j], fl_vardq='yes',
                                  logfile=log)
    
    #flat field and cut the object data

    #use different nsreduce calls depending on the observational
    #setup. the pipeline will have already quit if one of these two
    #setups wasn't used.
    if obs_setup == 'hk_2.20':
        iraf.nsreduce('sn@'+gallist, fl_cut='yes', fl_nsappw='yes',
                      fl_dark='no', fl_sky='no', fl_flat='yes',
                      flatimage='rgn'+calflat+'_flat',
                      fl_vardq='yes', logfile=log)
        #flat field and cut the sky data
        if len(darklist) > 0:
            iraf.nsreduce('sn@'+skylistshort, fl_cut='yes', fl_nsappw='yes',
                          fl_dark='no', fl_sky='no', fl_flat='yes',
                          flatimage='rgn'+calflat+'_flat', fl_vardq='yes',
                          logfile=log)
    if obs_setup =='hk_2.30':
        iraf.nsreduce('sn@'+gallist, fl_cut='yes', fl_nsappw='yes',
                      fl_dark='no', fl_sky='no', fl_flat='yes',
                      flatimage='rgn'+calflat+'_flat',
                      fl_vardq='yes', crval=23000.000000, cdelt=-2.115000,
                      logfile=log)
        #flat field and cut the sky data
        if len(darklist) > 0:
            iraf.nsreduce('sn@'+skylistshort, fl_cut='yes', fl_nsappw='yes',
                          fl_dark='no', fl_sky='no', fl_flat='yes',
                          flatimage='rgn'+calflat+'_flat', fl_vardq='yes',
                          crval=23000.000000, cdelt=-2.115000, logfile=log)
    
    #interpolate over bad pixels flagged in the DQ plane for both the
    #data and the sky
    iraf.nffixbad_anil('rsn@'+gallist,logfile=log)
    if len(darklist) > 0:
        iraf.nffixbad_anil('rsn@'+skylistshort,logfile=log)

    #do further cleaning
    
    #start with LA COSMIC to identify and fixpix to correct the cosmic
    #rays
    
    #need gain and readnoise from header to pass into LA COSMIC
    tmp_file = pyfits.open('n'+galexps[i]+'.fits')
    gain_in = tmp_file[0].header['GAIN']
    readn_in = tmp_file[0].header['RDNOISE']
    tmp_file.close()
    xorder_in = 4         #LA COSMIC xorder
    yorder_in = 3         #LA COSMIC yorder
    nifs_lacosmic_LP(workdir, galexps, gain_in, readn_in, xorder_in,
                     yorder_in, sigclip_in, sigfrac_in, objlim_in)
    if len(darklist) > 0:
        nifs_lacosmic_LP(workdir, skyshortexps, gain_in, readn_in,
                         xorder_in, yorder_in, sigclip_in, sigfrac_in,
                         objlim_in)

    #LA COSMIC isn't great at removing negative bad pixels, so using
    #xzap to identify negative pixels and fixpix to correct
    nifs_xzap_LP(workdir, galexps, nsigneg_in)
    if len(darklist) > 0:
        nifs_xzap_LP(workdir, skyshortexps, nsigneg_in)
    
    #read in the arc2uselist file to get the arc to use for this
    #galaxy/sky set (arc2uselist will only have 1 entry)
    arc2use=str(open('arc2uselist', 'r').readlines()[0]).strip()

    #derive the 2D to 3D spatial/spectral transformation. for some
    #reason, setting lyorder=5 below shows up as lyorder=4 when
    #running an interactive fit (and setting lyorder=4 below comes up
    #as lyorder=3 during the interactive fit). want lyorder=4, so
    #keeping lyorder=5 set below.
    iraf.nsfitcoords('brsn@'+gallist, outpref='f',
                     fl_int=flinter_nsfitcoords, lamptr=arc2use,
                     sdisttr='rgn'+ronchiflat, logfile=log, lxorder=4,
                     lyorder=5, sxorder=4, syorder=4)
    #repeat the 2D to 3D spatial/spectral transformation for the sky.
    if len(darklist) > 0:
        iraf.nsfitcoords('brsn@'+skylistshort, outpref='f',
                         fl_int=flinter_nsfitcoords, lamptr=arc2use,
                         sdisttr='rgn'+ronchiflat, logfile=log, lxorder=4,
                         lyorder=5, sxorder=4, syorder=4)
    
    #apply the transformation determined in the nffitcoords step for
    #both the galaxy and sky
    iraf.nstransform('fbrsn@'+gallist, outpref='t', logfile=log)
    if len(darklist) > 0:
        iraf.nstransform('fbrsn@'+skylistshort, outpref='t', logfile=log)
    
    #read in the telluric star list
    telcorr=open(telcorrlist, 'r').readlines()
    telcorr=[word.strip() for word in telcorr]
    
    for i in range(len(galexps)):
    
        #copy the telluric star to this directory (iraf will truncate
        #long path names and the file may not be found).
        tmp=subprocess.call(['cp',teldir+'/'+telcorr[i],
                            workdir+(telcorr[i].split('/'))[1]],
                            stderr=open(os.devnull,'w'))
        
        #in order to extract a 1D spectrum from the 2D data
        #non-interactively, create a temporary cube.
        iraf.nifcube('tfbrsn'+galexps[i], outcubes='tmp_tfbrsn'+galexps[i],
                     sscale=0.043, verbose='no', logfile=log)
        
        tmp_fitsfile = pyfits.open('tmp_tfbrsn'+galexps[i]+'.fits')
        tmp_gal_cube = tmp_fitsfile[1].data
        #tmp_gal_cube is in the form: wavelength, y, x. sum along
        #wavelength to get collapsed image of telluric star
        tmp_gal_image = np.sum(tmp_gal_cube,axis=0)
        #determine the y and x values of the maximum
        max_tmp_gal_image = tmp_gal_image.argmax()
        j,k = np.unravel_index(tmp_gal_image.argmax(), tmp_gal_image.shape)
        tmp_fitsfile.close()
        
        #the temporary cube doesn't have the same size as nfextract (x
        #has a length f 29 and y a length of 69). also iraf starts
        #indexing a 1 not 0. calculate the x and y center values to
        #feed into nfextract here.
        xsize_gal_img = tmp_gal_image.shape[1]
        ysize_gal_img = tmp_gal_image.shape[0]
        x_nfextract = round( (k / (xsize_gal_img/29.)) + 1., 2)
        y_nfextract = round( (j / (ysize_gal_img/69.)) + 1., 2)
        
        #delete the temporary cube and the nifs.log file that was
        #automatically created
        tmp=subprocess.call(['rm',workdir+'tmp_tfbrsn'+galexps[i]+'.fits'],
                            stderr=open(os.devnull,'w'))
        tmp=subprocess.call(['rm',workdir+'nifs.log'],
                            stderr=open(os.devnull,'w'))
    
        #correct the data for telluric absorption features
        
        #the telluric correction. below the extraction of a 1D galaxy
        #spectrum and the determination of the shift/scale parameters
        #(applied to the correction spectrum) are not done
        #interactively. for the 1D galaxy spectrum extraction, use an
        #aperture centered on the spaxel with the most
        #counts. fl_flux='no' means no multiplication by a blackbody
        #is done (as the telluric already has the blackbody shape
        #removed).
        if flinter_telluric == 'no':
            iraf.nftelluric_aniljonelle('tfbrsn'+galexps[i],
                                        (telcorr[i].split('/'))[1],
                                        fl_inter='no', fl_flux='no',
                                        xc=x_nfextract, yc=y_nfextract,
                                        logfile=log)
    
        #the telluric correction. below the extraction of a 1D galaxy
        #spectrum and the determination of the shift/scale parameters
        #(applied to the correction spectrum) are done
        #interactively. for the 1D galaxy spectrum extraction, the
        #(x,y) location of the spaxel with the most counts is printed
        #to the screen to guide the user. fl_flux='no' means no
        #multiplication by a blackbody is done (as the telluric
        #already has the blackbody shape removed). note, need to have
        #a DS9 window open.
        if flinter_telluric == 'yes':
            
            #print the (x,y) pixel with the maximum value to the
            #screen
            print('The pixel with the maximum value '+\
                  'is: x=%1.1f y=%1.1f' % (x_nfextract,y_nfextract))
              
            iraf.nftelluric_anil('tfbrsn'+galexps[i],(telcorr[i].split('/'))[1],
                                 fl_inter='yes', fl_flux='no',
                                 logfile=log)
        
        if len(darklist) > 0:
            
            #correct the sky for telluric absorption features. since
            #we will eventually need one sky exposure for each galaxy
            #exposure, create duplicates here. adopt the scale and
            #shift from the galaxy cube above.
            tmp = pyfits.open('atfbrsn'+galexps[i]+'.fits')
            tel_shift = tmp[0].header['NFTELSHF']
            tel_scale = tmp[0].header['NFTELSCL']
            tmp.close()
            iraf.nftelluric_aniljonelle('tfbrsn'+skyexps[i],
                                        (telcorr[i].split('/'))[1],
                                        outimages='atfbrsn'+skyexps[i]+'_'+\
                                        (galexps[i])[10:]+'.fits',
                                        fl_inter='no', fl_xcorr='no',
                                        fl_twea='no', fl_flux='no',
                                        shift=tel_shift, scale=tel_scale,
                                        dshift=0., dscale=0., xc=x_nfextract,
                                        yc=y_nfextract, logfile=log)
        
        #remove the copy of the telluric star now that the galaxy and
        #sky exposure (if applicable) has been reduced
        tmp=subprocess.call(['rm',workdir+(telcorr[i].split('/'))[1]],
                            stderr=open(os.devnull,'w'))
    
        #create the 3D cube of the galaxy exposures only for the
        #purposes of collapsing and figuring out spatial
        #offsets. therefore, just use gemini's nifcube for now even
        #though the variance/dq extensions are not properly dealt
        #with.
        print('')
        print('Creating temporary cube of galaxy exposure and '+\
              'collapsing to form an image. Please be patient.')
        print('')
        iraf.nifcube(inimages='atfbrsn'+galexps[i], sscale=0.05,
                     verbose='no', logfile=log)
                         
        #collapse the galaxy cube to get an image
        tmp_cube = pyfits.open('catfbrsn'+galexps[i]+'.fits')
        tmp_gal_cube = tmp_cube[1].data
        #tmp_gal_cube is in the form: wavelength, y, x. sum along
        #wavelength to get collapsed image of galaxy
        tmp_gal_image = np.sum(tmp_gal_cube,axis=0)
        #will use the primary header (with wcs info) from the galaxy
        #cube and make some adjustments
        hdu = pyfits.PrimaryHDU(tmp_gal_image)
        hdu.header = tmp_cube[0].header
        hdu.header['BITPIX'] = -32
        hdu.header['NAXIS'] = 2
        hdu.header['NAXIS1'] = (np.shape(tmp_gal_image))[1]
        hdu.header['NAXIS2'] = (np.shape(tmp_gal_image))[0]
        hdu.header['PIXSCALE'] = 0.05
        #write to the fits file
        hdu.writeto('image_catfbrsn'+galexps[i]+'.fits',
                    output_verify='silentfix')
        tmp_cube.close()
        #remove the temporary cube
        tmp=subprocess.call(['rm',workdir+'catfbrsn'+galexps[i]+'.fits'],
                            stderr=open(os.devnull,'w'))
    
        #move the collapsed galaxy cubes to the merged directory.
        tmp=subprocess.call(['mv','image_catfbrsn'+galexps[i]+'.fits',
                            rootdir+reducedir+galname+'/merged/'+\
                            obs_setup+'/'],stderr=open(os.devnull,'w'))
        #copy the telluric corrected galaxy and sky exposures and make
        #a copy for header modification.
        tmp=subprocess.call(['cp',workdir+'atfbrsn'+galexps[i]+'.fits',
                            rootdir+reducedir+galname+'/merged/'+\
                            obs_setup+'/'],stderr=open(os.devnull,'w'))
        tmp=subprocess.call(['cp',workdir+'atfbrsn'+galexps[i]+'.fits',
                            rootdir+reducedir+galname+'/merged/'+\
                            obs_setup+'/hatfbrsn'+galexps[i]+'.fits'],
                            stderr=open(os.devnull,'w'))
        if len(darklist) > 0:
            tmp=subprocess.call(['cp','atfbrsn'+skyexps[i]+'_'+\
                                (galexps[i])[10:]+'.fits',rootdir+reducedir+\
                                galname+'/merged/'+obs_setup+'/'],
                                stderr=open(os.devnull,'w'))
            tmp=subprocess.call(['cp','atfbrsn'+skyexps[i]+'_'+\
                                (galexps[i])[10:]+'.fits',rootdir+reducedir+\
                                galname+'/merged/'+obs_setup+'/hatfbrsn'+\
                                skyexps[i]+'_'+(galexps[i])[10:]+'.fits'],
                                stderr=open(os.devnull,'w'))

        
    ###########################################################################
    # Reset to user defaults                                                  #
    ###########################################################################
    if user_clobber == "no":
        iraf.set(clobber='no')
        
###########################################################################
#                  End of the Science Data Reduction                      #
#                                                                         #
# The output of this reduction script are prepared, sky-subtracted,       #
# reduced, bad pixel corrected, geometrically rectified, telluric         #
# corected galaxy slices, ready to be merged and turned into a final data #
# cube. Copies of those files are placed in a "../merged" directory.      #
# Also, a temporary cube of each galaxy exposure is created and summed    #
# along the wavelength axis to produce an image. The images are also      #
# saved in the "../merged" directory, and will later be used to determine #
# the spatial offsets between individual galaxy exposures. Finally,       #
# prepared, bias-subtracted, reduced, bad pixel corrected, geometrically  #
# rectified, telluric corrected sky slices are created, ready to be       #
# merged (in the same way as the galaxy exposures) and turned into a      #
# final data cube. In the case of this reduction, the final output files  #
# are called: atfbrsn+science and atfbrsn+sky+_+science (this is the sky  #
# with the corresponding science exposure identified).                    #
#                                                                         #
# The meaning of the output prefixes are described below:                 #
#                                                                         #
# n=nfprepared   s=sky subtracted/bias subtracted   r=nsreduced           #
# b = bad pixel corrected /LA COSMIC cleaned/ XZAP cleaned                #
# f= run through nffitcoords   t = nftransformed                          #
# a = corrected for telluric absorption features                          #
# some baseline calibs are copied to the galaxy working directory as well #
#                                                                         #
###########################################################################
