def nifs_xzap_LP(workdir, exps, nsigneg_in):

    #import some useful python utilities
    import os
    import sys
    import glob
    import subprocess
    import astropy.io.fits as pyfits
    import numpy as np
    #import the pyraf module and relevant packages
    from pyraf import iraf
    iraf.xdimsum(_doprint=0)
    iraf.images(_doprint=0)
    iraf.imutil(_doprint=0)
    iraf.proto(_doprint=0)
    from pyraf import iraffunctions

    #unlearn the used tasks
    iraf.unlearn('xzap_jonelle', 'imcopy', 'fixpix')

    #create a subdirectory to hold the bad pixel masks, if it doesn't
    #already exist.
    os.chdir(workdir)
    if len(glob.glob('clean')) == 0:
        tmp=subprocess.call(['mkdir','clean'],stderr=open(os.devnull,'w'))
    os.chdir(workdir+'clean')
    if len(glob.glob('masks_negpix')) == 0:
        tmp=subprocess.call(['mkdir','masks_negpix'],\
                            stderr=open(os.devnull,'w'))
    os.chdir(workdir)

    #set the science, variance, and dq extensions
    sci_exten = [2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,\
                 59,62,65,68,71,74,77,80,83,86]
    var_exten = [x+1 for x in sci_exten]
    dq_exten = [x+2 for x in sci_exten]

    for ii in range(len(exps)):
        
        exp_now = exps[ii]
        
        for jj in range(len(sci_exten)):

            #set high value for nsigzap - have identified cosmic rays
            #with la cosmic already. this is for identifying bad
            #negative pixels.
            iraf.xzap_jonelle(inlist='brsn'+exp_now+'['+str(sci_exten[jj])+']',
                              omasks='', outlist='brsn'+exp_now+'_exten'+\
                              str(sci_exten[jj]), crmasks='brsn'+exp_now+\
                              '_exten'+str(sci_exten[jj])+'_mask',
                              nsigzap=10000.0,
                              nsigneg=nsigneg_in, checklimits='no')

            #copy the bad pixel mask (.pl) to a fits file
            iraf.imcopy(input='brsn'+exp_now+'_exten'+str(sci_exten[jj])+\
                        '_mask.pl', output='brsn'+exp_now+'_exten'+\
                        str(sci_exten[jj])+'_mask.fits', verbose='no')

            #the bad pixel mask generated above has good values as 1
            #and bad values as 0. need to swap this for fixpix. also
            #update the dq extension.
            maskfile = pyfits.open('brsn'+exp_now+'_exten'+str(sci_exten[jj])+\
                                   '_mask.fits', mode='update')
            mask_tmp = maskfile[0].data
            expfile = pyfits.open('brsn'+exp_now+'.fits', mode='update')
            dq_tmp = expfile[dq_exten[jj]].data
            #find the pixels in the mask with values of 1 (good) and 0
            #(bad)
            index_good = np.where(mask_tmp == 1.)
            index_bad = np.where(mask_tmp == 0.)
            #check that there are good pixels in the mask, otherwise stop
            #and tell the user.
            if (np.shape(index_good))[1] == 0:
                print 'Problem finding good pixels in XZAP mask.'
                sys.exit('Quiting the pipeline.')
            #if there are good pixels in the mask, label them as 0
            #instead of 1. if there are bad pixels in the mask label
            #them as 1. also update the dq.
            if (np.shape(index_good))[1] > 0:
                mask_tmp[index_good] = 0.
                dq_tmp[index_good] = 0. + dq_tmp[index_good]
            if (np.shape(index_bad))[1] > 0:
                mask_tmp[index_bad] = 1.
                dq_tmp[index_bad] = 1.
            #modify some headers.
            #del maskfile[0].header['FIXPIX']
            #del maskfile[0].header['FIXPIX02']
            #del maskfile[0].header['SKYMED,S']
            #del expfile[0].header['CRMASK']
            #update the fits file
            maskfile.close(output_verify='silentfix')
            expfile.close(output_verify='silentfix')

            #clean up files
            tmp=subprocess.call(['mv','brsn'+exp_now+'_exten'+\
                                 str(sci_exten[jj])+'_mask.fits',\
                                 'clean/masks_negpix/'],\
                                 stderr=open(os.devnull,'w'))
            tmp=subprocess.call(['rm','brsn'+exp_now+'_exten'+\
                                 str(sci_exten[jj])+'_mask.pl'],\
                                 stderr=open(os.devnull,'w'))
            tmp=subprocess.call(['rm','brsn'+exp_now+'_exten'+\
                                 str(sci_exten[jj])+'.fits'],\
                                 stderr=open(os.devnull,'w'))

            #use the bad pixel mask and iraf's fixpix to correct the
            #bad pixels and do the same to the corresponding variance
            #extension
            iraf.fixpix(images='brsn'+exp_now+'['+str(sci_exten[jj])+']',
                        masks='clean/masks_negpix/brsn'+exp_now+'_exten'+\
                        str(sci_exten[jj])+'_mask.fits', linterp='INDEF',
                        cinterp='INDEF')
            iraf.fixpix(images='brsn'+exp_now+'['+str(var_exten[jj])+']',
                        masks='clean/masks_negpix/brsn'+exp_now+'_exten'+\
                        str(sci_exten[jj])+'_mask.fits', linterp='INDEF',
                        cinterp='INDEF')

            #remove some headers
            expfile = pyfits.open('brsn'+exp_now+'.fits', mode='update')
            del expfile[sci_exten[jj]].header['FIXPIX']
            #del expfile[sci_exten[jj]].header['FIXPIX02']
            del expfile[sci_exten[jj]].header['CRMASK']
            del expfile[var_exten[jj]].header['FIXPIX']
            #del expfile[var_exten[jj]].header['FIXPIX02']
            expfile.close(output_verify='silentfix')
            
