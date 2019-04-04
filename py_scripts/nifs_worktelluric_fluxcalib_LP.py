"""

This routine corrects an A0 V telluric star cube for telluric
features, leaving just A0 V spectra. If such a cube already exists,
then this routine will be skipped.

"""

def nifs_worktelluric_fluxcalib_LP(telstar_calib_dir, telstar_calib_file,
                                   galname, date, flinter_telluric):

    #import some useful python utilities
    import os
    import time
    import glob
    import shutil
    import astropy.io.fits as pyfits
    import numpy as np
    import subprocess

    #check to see if there is already a telluric corrected telluric
    #cube. will skip this entire function if there is already a cube.
    check_telcube = glob.glob(telstar_calib_dir+'catfbrsn'+\
                              telstar_calib_file+'.fits')
                              
    if len(check_telcube) == 0:
        
        print('')
        print('Creating a cube of the telluric that will be used for '+\
                  'flux calibration.')
        print('')
    
        ###################################################################
        #  STEP 1: Prepare IRAF  		                                  #
        ###################################################################
    
        #import the pyraf module and relevant packages
        from pyraf import iraf
        iraf.gemini(_doprint=0)
        iraf.nifs(_doprint=0)
        iraf.gnirs(_doprint=0)
        iraf.gemtools(_doprint=0)
        from pyraf import iraffunctions

        #unlearn the used tasks
        iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs,\
                     'nftelluric_anil','nftelluric_aniljonelle')
        iraf.set(stdimage='imt2048')

        #create a log file and back up the previous one if it already exists
        log = 'telluric_'+date+'.log'
        if os.path.exists(log):
            t = time.localtime()
            app = "_"+str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)+'_'+ \
            str(t[3]).zfill(2)+':'+str(t[4]).zfill(2)+':'+str(t[5]).zfill(2)
            shutil.move(log,log+app)

        #change to the telluric star directory within pyraf
        iraffunctions.chdir(telstar_calib_dir)
        
        #prepare the package for NIFS
        iraf.nsheaders('nifs',logfile=log)

        #set clobber to 'yes' for the script. this still does not make the
        #gemini tasks overwrite files, so you will likely have to remove
        #files if you re-run the script.
        user_clobber=iraf.envget("clobber")
        iraf.reset(clobber='yes')

    
        ####################################################################
        #  STEP 2: Telluric Correct the Telluric Star  		               #
        ####################################################################

        #in order to extract a 1D spectrum from the 2D data
        #non-interactively, create a temporary cube.
        iraf.nifcube('tfbrsn'+telstar_calib_file, outcubes='tmp_tfbrsn'+\
                     telstar_calib_file, sscale=0.043, verbose='no', \
                     logfile=log)
        
        tmp_fitsfile = pyfits.open('tmp_tfbrsn'+telstar_calib_file+'.fits')
        tmp_tel_cube = tmp_fitsfile[1].data
        #tmp_tel_cube is in the form: wavelength, y, x. sum along
        #wavelength to get collapsed image of telluric star
        tmp_tel_image = np.sum(tmp_tel_cube,axis=0)
        #determine the y and x values of the maximum
        max_tmp_tel_image = tmp_tel_image.argmax()
        j,k = np.unravel_index(tmp_tel_image.argmax(), tmp_tel_image.shape)
        tmp_fitsfile.close()

        #the temporary cube doesn't have the same size as nfextract (x
        #has a length f 29 and y a length of 69). also iraf starts
        #indexing a 1 not 0. calculate the x and y center values to
        #feed into nfextract here.
        xsize_tel_img = tmp_tel_image.shape[1]
        ysize_tel_img = tmp_tel_image.shape[0]
        x_nfextract = round( (k / (xsize_tel_img/29.)) + 1., 2)
        y_nfextract = round( (j / (ysize_tel_img/69.)) + 1., 2)
        
        #delete the temporary cube and the nifs.log file that was
        #automatically created
        tmp=subprocess.call(['rm',telstar_calib_dir+'tmp_tfbrsn'+\
                             telstar_calib_file+'.fits'],
                             stderr=open(os.devnull,'w'))
        tmp=subprocess.call(['rm',telstar_calib_dir+'nifs.log'],
                             stderr=open(os.devnull,'w'))

        #correct the telluric star for telluric absorption features
        telcorr = glob.glob('*final.fits')
        
        #the telluric correction. below the extraction of a 1D
        #telluric spectrum and the determination of the shift/scale
        #parameters (applied to the correction spectrum) are not done
        #interactively. for the 1D telluric spectrum extraction, use
        #an aperture centered on the spaxel with the most
        #counts. fl_flux='no' means no multiplication by a blackbody
        #is done (as the telluric already has the blackbody shape
        #removed).
        if flinter_telluric == 'no':
            iraf.nftelluric_aniljonelle('tfbrsn'+telstar_calib_file,
                                        telcorr[0], fl_inter='no',
                                        fl_flux='no',
                                        sample='23000:24100',
                                        xc=x_nfextract, yc=y_nfextract,
                                        logfile=log)

        #the telluric correction. below the extraction of a 1D
        #telluric spectrum and the determination of the shift/scale
        #parameters (applied to the correction spectrum) are done
        #interactively. for the 1D telluric spectrum extraction, the
        #(x,y) location of the spaxel with the most counts is printed
        #to the screen to guide the user. fl_flux='no' means no
        #multiplication by a blackbody is done (as the telluric
        #already has the blackbody shape removed). note, need to have
        #a DS9 window open.
        if flinter_telluric == 'yes':
            
            #print the (x,y) pixel with the maximum value to the screen
            print('The pixel with the maximum value is: '+\
                      'x=%1.1f y=%1.1f' % (x_nfextract,y_nfextract))
              
            iraf.nftelluric_anil('tfbrsn'+telstar_calib_file,telcorr[0],
                                 fl_inter='yes', fl_flux='no',
                                 sample='23000:24100', logfile=log)

    
        #######################################################################
        #  STEP 3: Create Cube of Telluric Corrected Telluric Star  		  #
        #######################################################################

        iraf.nifcube(inimages='atfbrsn'+telstar_calib_file, sscale=0.05,
                     verbose='no', logfile=log)

