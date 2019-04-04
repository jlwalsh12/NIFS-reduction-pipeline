def nifs_lacosmic_LP(workdir, exps, gain_in, readn_in, xorder_in, \
                     yorder_in, sigclip_in, sigfrac_in, objlim_in):

    #import some useful python utilities
    import os
    import sys
    import glob
    import subprocess
    import astropy.io.fits as pyfits
    import numpy as np
    #import the pyraf module and relevant packages
    from pyraf import iraf
    iraf.proto(_doprint=0)
    iraf.stsdas(_doprint=0)
    from pyraf import iraffunctions

    #unlearn the used tasks
    iraf.unlearn('lacos_spec', 'fixpix')

    #create a subdirectory to hold the bad pixel masks, if it doesn't
    #already exist.
    os.chdir(workdir)
    if len(glob.glob('clean')) == 0:
        tmp=subprocess.call(['mkdir','clean'],stderr=open(os.devnull,'w'))
    os.chdir(workdir+'clean')
    if len(glob.glob('masks_lacosmic')) == 0:
        tmp=subprocess.call(['mkdir','masks_lacosmic'],\
                             stderr=open(os.devnull,'w'))
    os.chdir(workdir)

    #set the science, variance, and dq extensions
    sci_exten = [2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,\
                 59,62,65,68,71,74,77,80,83,86]
    var_exten = [x+1 for x in sci_exten]
    dq_exten = [x+2 for x in sci_exten]

    for ii in range(len(exps)):
        
        exp_now = exps[ii]
        print ''
        print 'Running LA COSMIC on '+exp_now
        
        for jj in range(len(sci_exten)):

            #call LA COSMIC to identify cosmic rays in each science
            #extension of each galaxy/sky exposure
            iraf.lacos_spec(input='brsn'+exp_now+'['+str(sci_exten[jj])+\
                            '].fits', output='brsn'+exp_now+'_exten'+\
                            str(sci_exten[jj])+'_clean.fits', outmask='brsn'+\
                            exp_now+'_exten'+str(sci_exten[jj])+\
                            '_mask.fits', gain=gain_in, readn=readn_in,\
                            xorder=xorder_in, yorder=yorder_in,\
                            sigclip=sigclip_in, sigfrac=sigfrac_in,\
                            objlim=objlim_in, niter=4, verbose='no')

            #update the dq extension
            maskfile = pyfits.open('brsn'+exp_now+'_exten'+str(sci_exten[jj])+\
                                   '_mask.fits')
            mask_tmp = maskfile[0].data
            expfile = pyfits.open('brsn'+exp_now+'.fits', mode='update')
            dq_tmp = expfile[dq_exten[jj]].data
            #find the pixels in the mask with values of 0 (good) and 1
            #(bad)
            index_good = np.where(mask_tmp == 0.)
            index_bad = np.where(mask_tmp == 1.)
            #check that there are good pixels in the mask, otherwise stop
            #and tell the user.
            if (np.shape(index_good))[1] == 0:
                print 'Problem finding good pixels in LA COSMIC mask.'
                sys.exit('Quiting the pipeline.')
            #update the dq
            if (np.shape(index_good))[1] > 0:
                dq_tmp[index_good] = 0. + dq_tmp[index_good]
            if (np.shape(index_bad))[1] > 0:
                dq_tmp[index_bad] = 1.
            #update the fits file
            maskfile.close(output_verify='silentfix')
            expfile.close(output_verify='silentfix')
            
            #move the bad pixel mask to a subdirectory. clean up
            #files.
            tmp=subprocess.call(['mv','brsn'+exp_now+'_exten'+\
                                 str(sci_exten[jj])+'_mask.fits',\
                                 'clean/masks_lacosmic/'],\
                                 stderr=open(os.devnull,'w'))
            tmp=subprocess.call(['rm','brsn'+exp_now+'_exten'+\
                                 str(sci_exten[jj])+'_clean.fits'],\
                                 stderr=open(os.devnull,'w'))

            #use the bad pixel mask and iraf's fixpix to correct the
            #bad pixels and do the same to the corresponding variance
            #extension
            iraf.fixpix(images='brsn'+exp_now+'['+str(sci_exten[jj])+']',\
                        masks='clean/masks_lacosmic/brsn'+exp_now+'_exten'+\
                        str(sci_exten[jj])+'_mask.fits', linterp='INDEF', \
                        cinterp='INDEF')
            iraf.fixpix(images='brsn'+exp_now+'['+str(var_exten[jj])+']',\
                        masks='clean/masks_lacosmic/brsn'+exp_now+'_exten'+\
                        str(sci_exten[jj])+'_mask.fits', linterp='INDEF', \
                        cinterp='INDEF')

            #remove some headers
            maskfile = pyfits.open('clean/masks_lacosmic/brsn'+exp_now+\
                                   '_exten'+str(sci_exten[jj])+\
                                   '_mask.fits', mode='update')
            del maskfile[0].header['FIXPIX']
            maskfile.close(output_verify='silentfix')

            expfile = pyfits.open('brsn'+exp_now+'.fits', mode='update')
            del expfile[sci_exten[jj]].header['FIXPIX']
            del expfile[sci_exten[jj]].header['FIXPIX02']
            del expfile[var_exten[jj]].header['FIXPIX']
            del expfile[var_exten[jj]].header['FIXPIX02']
            expfile.close(output_verify='silentfix')
            
