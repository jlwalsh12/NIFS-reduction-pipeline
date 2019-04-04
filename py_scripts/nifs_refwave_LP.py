"""

This routine creates the wavelength array that the merged cube will
have. If requested, a barycentric or heliocentric correction is
applied.

"""

def nifs_refwave_LP(workdir, objexp, idlpath, corr_type):

    import astropy.io.fits as pyfits
    import numpy as np
    import pidly
    idl = pidly.IDL(idlpath)

    #if needed, determine the barycentric or heliocentric correction
    #for a reference galaxy exposure (the first file) in the merged
    #directory.
    if corr_type.lower() == 'bary' or corr_type.lower() == 'helio':
        
        objfile = pyfits.open(workdir+objexp)
        ra = (objfile[0].header['RA'])/15.
        dec = objfile[0].header['DEC']
        dateobs = objfile[0].header['DATE-OBS']
        dateobs_yr = dateobs[0:4]
        dateobs_month = dateobs[5:7]
        dateobs_day = dateobs[8:10]
        ut = objfile[0].header['UT']
        ut_hrs = ut[0:2]
        ut_mins = ut[3:5]
        ut_sec = ut[6:10]
        objfile.close()

        v_corr = idl.func('vcorr_anil', ra, dec, dateobs_yr,
                          dateobs_month, dateobs_day, ut_hrs,
                          ut_mins, ut_sec, corr_type)

    #set the reference wavelength in units of angstroms - need to
    #minimum wavelength, the maximum wavelength, and the wavelength
    #scale for nifcube_jonelle.
    objfile = pyfits.open(workdir+objexp)
    lambda0 = objfile[1].header['CRVAL1']
    dlambda = objfile[1].header['CDELT1']
    #iraf/the header values start at 1 whereas our array will start at
    #0
    crpix = objfile[1].header['CRPIX1'] - 1.
    wavelen = (np.shape(objfile[1].data))[1]
    objfile.close()
    refwave = \
      (((np.linspace(0.,wavelen-1,num=wavelen,dtype=int))-crpix)*dlambda)+\
      lambda0
    #make the barycentric or heliocentric correction if requested
    if corr_type.lower() == 'bary' or corr_type.lower() == 'helio':
        refwave = refwave*(1. + (v_corr/2.99792458e5))
    #get the minimum and maximum wavelength
    min_refwave_vcorr = refwave[0]
    max_refwave_vcorr = refwave[wavelen-1]
    dlambda_vcorr = refwave[100] - refwave[99]

    return [min_refwave_vcorr, max_refwave_vcorr, dlambda_vcorr]
