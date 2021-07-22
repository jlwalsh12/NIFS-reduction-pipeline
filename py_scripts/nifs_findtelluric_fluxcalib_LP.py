"""

This routine looks in the header of the merged galaxy or sky cube for
the location and file name of one of the telluric stars. This telluric
star will be used for a rough flux calibration.

"""

def nifs_findtelluric_fluxcalib_LP(begindir, workdir, galname, obs_setup):

    import os
    import sys
    import glob
    import astropy.io.fits as pyfits

    os.chdir(workdir)

    #read in the header of the merged galaxy cube to find the file
    #name of one of the telluric stars. this star will be used to do
    #the flux calibration.
    checkfile = sorted(glob.glob(galname+'_combined.fits'),key=os.path.basename)
    if len(checkfile) == 1:
        tmp = pyfits.open(galname+'_combined.fits')
        telfile_calib = tmp[0].header['NFTELCAL']
        tmp.close()
    else:
        sys.exit('Cannot find merged galaxy cube to flux calibrate.')

    #extract the date and frame number of the telluric star file
    index_tmp = telfile_calib.find('N')
    date_tel = telfile_calib[index_tmp+1:index_tmp+9]
    index_tmp = telfile_calib.find('_')
    frame_tel = telfile_calib[index_tmp-5:index_tmp]

    #track down which telluric star the file name belongs to
    checkteldir = begindir+'tellurics/'+date_tel+'/'+obs_setup+'/'
    os.chdir(checkteldir)
    checktelstar = sorted(glob.glob('*'),key=os.path.basename)

    telstar_calib = ''
    for i in range(len(checktelstar)):
        
        os.chdir(begindir+'tellurics/'+date_tel+'/'+obs_setup+\
                 '/'+checktelstar[i])
        checktel = sorted(glob.glob('N'+date_tel+frame_tel+'.fits'),
                          key=os.path.basename)
        if len(checktel) == 1:
            telstar_calib = begindir+'tellurics/'+date_tel+'/'+\
              obs_setup+'/'+checktelstar[i]+'/'

    if len(telstar_calib) == 0:
        sys.exit('Cannot find the telluric star to use for flux calibration.')

    return [telstar_calib, date_tel, frame_tel]
