"""

This routine automatically creates the telluric correction list. The
telluric star exposure taken closest in time to the galaxy or psf star
exposure is used here. (After calling this routine, the pipeline asks
the user to check the telluric correction list and make manual
modifications, e.g., if a star closer in airmass should be used
instead.)

"""

def nifs_create_telcorrlist_LP(workdir, objlist, rootdir, reducedir, date,
                               obs_setup):

    import os
    import glob
    import astropy.io.fits as pyfits
    import numpy as np
    
    os.chdir(workdir)

    objfile=open(objlist, 'r').readlines()
    objfile=[word.strip() for word in objfile]

    #get the middle UT time for each of the galaxy/psf exposures
    objfile_ut = []
    for ii in range(len(objfile)):
        objfile_header = (pyfits.open(objfile[ii]+'.fits'))[0].header
        decimal_ut = ( float((objfile_header['UT'])[0:2]) + \
                       float((objfile_header['UT'])[3:5])/60. + \
                       float((objfile_header['UT'])[6:])/(60.*60.) ) + \
                       ( (float(objfile_header['EXPTIME'])/2.)/(60.*60.) )
        objfile_ut.append(decimal_ut)

    #get the names of the telluric stars observed on the same day and
    #in the same setup as the galaxy/psf
    teldir = rootdir+reducedir+'tellurics/'+date+'/'+obs_setup+'/'
    os.chdir(teldir)
    telluric_star = glob.glob('*')

    telfile = []
    telfile_ut = []
    for ii in range(len(telluric_star)):

        os.chdir(teldir+telluric_star[ii])
        telluriclist = glob.glob('telluriclist*')

        for jj in range(len(telluriclist)):

            telfile_all=open(telluriclist[jj], 'r').readlines()
            telfile_all=[word.strip() for word in telfile_all]
            telfile.append(telluric_star[ii] + '/'+ \
                           str(open(telluriclist[jj],
                               'r').readlines()[0]).strip())

            #get the UT time for each of the telluric exposures in
            #this group
            telfile_all_ut = []
            for kk in range(len(telfile_all)):
                telfile_all_header = \
                  (pyfits.open(telfile_all[kk]+'.fits'))[0].header
                telfile_all_ut.append(float((telfile_all_header['UT'])[0:2])+\
                                float((telfile_all_header['UT'])[3:5])/60.+\
                                float((telfile_all_header['UT'])[6:])/(60.*60.))

            #use the average UT time for this grouping of telluric
            #star observations
            telfile_ut.append(np.mean(telfile_all_ut))

    #for each galaxy/psf exposure find the closest telluric set in
    #time. use that as a guess for what the user wants for the
    #telluric correction
    telcorr_files = []
    for ii in range(len(objfile)):
        idx = telfile_ut.index( min(telfile_ut,
                                    key=lambda x:abs(x-objfile_ut[ii])) )
        tmp = (telfile[idx].split('/'))[0] + '/cgxtfbrsn' + \
          (telfile[idx].split('/'))[1] + '_final.fits'
        telcorr_files.append(tmp)

    #write the telluric correction file
    os.chdir(workdir)
    output = open('telcorrlist', 'w')
    for item in telcorr_files:
        output.write("%s\n" % item)
    output.close()
