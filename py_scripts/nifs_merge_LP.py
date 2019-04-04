"""

Create a single data cube with science, variance and data quality
extensions from individual exposures. Each individual exposure is
composed of slices that have been geometrically rectified and
wavelength calibrated. Thus, the individual exposures are like "data
cubes", but have not yet been assembled into a data cube. The headers
of the individual exposures will be modified in order to trick PyRAF's
nifcube (actually a modified version of nifcube called
nifcube_jonelle) into accounting for the spatial offsets between the
individual exposures. The spatial offsets have been previously
determined and written to a file (which is provided as input into this
routine; see nifs_getoffmethod_LP.py). A wavelength vector is also
given as input into this routine (see nifs_refwave_LP.py), and will be
the third axis of the final data cube. This method uses less
interpolation compared to generating data cubes for each individual
exposure then merging to form a final data cube; the method can deal
with non-integer spatial offsets between individual exposures; and the
method incorporates the use of an input mask for the slices so that
the bad edges can be excluded. The output cube has units of counts/s.

"""

def nifs_merge_LP(workdir, obj_exps, galname, offsetfile, headerinfo, refwave,
                  intype):

    ###########################################################################
    #  STEP 1: Prepare IRAF  		                                          #
    ###########################################################################

    #import some useful python utilities
    import os
    import sys
    import time
    import shutil
    import subprocess
    import glob
    import astropy.io.fits as pyfits
    import numpy as np
    #import the pyraf module and relevant packages
    from pyraf import iraf
    iraf.gemini(_doprint=0)
    iraf.nifs(_doprint=0)
    iraf.gnirs(_doprint=0)
    iraf.gemtools(_doprint=0)
    from pyraf import iraffunctions

    #unlearn the used tasks
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs,\
                 'nifcube_jonelle')
    iraf.set(stdimage='imt2048')

    #create a log file and back up the previous one if it already
    #exists
    if intype == 'gal' or intype == 'sky':
        log = galname+'_merged'+'.log'
    if intype == 'psf':
        log = galname+'_psf_merged'+'.log'
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

    
    ###########################################################################
    # STEP 2:  Update Header Keywords and Create a Mask                       #
    ###########################################################################

    #will trick nifcube_jonelle into accounting for spatial offsets by
    #modifying some headers. will also create a mask for merging with
    #nifcube_jonelle simultaneously.

    #read in the offsets in pixels relative to the reference file.
    offsets=np.loadtxt(offsetfile, dtype='float', comments='#')
    comment_offset=str(open(offsetfile, 'r').readlines()[0]).strip()
    tmp_comment_offset = (comment_offset.split('='))[1]
    pixscale_offset = float( (tmp_comment_offset.split('"'))[0] )
    #set the science, variance, and dq extension numbers
    nslice=29
    sciexts = np.linspace(0.,nslice-1,num=nslice,dtype=int)*3+1
    varexts = np.linspace(0.,nslice-1,num=nslice,dtype=int)*3+2
    dqexts = np.linspace(0.,nslice-1,num=nslice,dtype=int)*3+3

    #make a directory to hold the masks for merging cubes
    if len(glob.glob('mask')) == 0:
        tmp=subprocess.call(['mkdir','mask'],
                            stderr=open(os.devnull,'w'))
                            
    print('')
    for i in range(len(obj_exps)):

        if len(obj_exps) > 1:

            print('Preparing file %i of %i for merging.' % (i+1,len(obj_exps)))

            #convert offset for this exposure from pixels to arcsec
            yoff_arc = (offsets[i,1])*pixscale_offset
            xoff_arc = (offsets[i,0])*pixscale_offset
                
        else:
            print('Creating a cube for the single object exposure.')

            #convert offset for this exposure from pixels to arcsec
            yoff_arc = (offsets[1])*pixscale_offset
            xoff_arc = (offsets[0])*pixscale_offset

        #open the hatfbrsn*.fits file to modify
        objfile = pyfits.open(workdir+'h'+obj_exps[i],mode='update')

        header_check = 'pass'
        for j in range(nslice):

            #for the science extensions
            objfile[sciexts[j]].header['CRVAL2'] = \
              objfile[sciexts[j]].header['CRVAL2'] + yoff_arc
            objfile[sciexts[j]].header['CRVAL3'] = \
              objfile[sciexts[j]].header['CRVAL3'] + xoff_arc
            objfile[sciexts[j]].header['GEMMASK'] = 'mask/'+ \
              (obj_exps[i].split('.'))[0]+'_sci'+str(sciexts[j])+ \
              '_var'+str(varexts[j])+'_dq'+str(dqexts[j])+'.pl'
            #for flux calibration want counts/s
            objfile[sciexts[j]].data = \
              objfile[sciexts[j]].data / float(objfile[0].header['EXPTIME'])

            if objfile[sciexts[j]].header['CTYPE2'] != 'LINEAR':
                header_check = 'fail'
            if objfile[sciexts[j]].header['CTYPE3'] != 'LINEAR':
                header_check = 'fail'
            if objfile[sciexts[j]].header['WAT2_001'] != \
              'wtype=linear axtype=eta':
                header_check = 'fail'
            if objfile[sciexts[j]].header['WAT3_001'] != \
              'wtype=linear axtype=xi':
                header_check = 'fail'
            if objfile[sciexts[j]].header['WCSDIM'] != 3:
                header_check = 'fail'
            if objfile[sciexts[j]].header['NAXIS'] != 2:
                header_check = 'fail'

            #for the variance extensions
            objfile[varexts[j]].header['CRVAL2'] = \
                objfile[varexts[j]].header['CRVAL2'] + yoff_arc
            objfile[varexts[j]].header['CRVAL3'] = \
                objfile[varexts[j]].header['CRVAL3'] + xoff_arc
            objfile[varexts[j]].header['GEMMASK'] = 'mask/'+ \
              (obj_exps[i].split('.'))[0]+'_sci'+str(sciexts[j])+ \
              '_var'+str(varexts[j])+'_dq'+str(dqexts[j])+'.pl'
            #for flux calibration want counts/s
            objfile[varexts[j]].data = objfile[varexts[j]].data / \
              (float(objfile[0].header['EXPTIME']))**2
              
            if objfile[varexts[j]].header['CTYPE2'] != 'LINEAR':
                header_check = 'fail'
            if objfile[varexts[j]].header['CTYPE3'] != 'LINEAR':
                header_check = 'fail'
            if objfile[varexts[j]].header['WAT2_001'] != \
              'wtype=linear axtype=eta':
                header_check = 'fail'
            if objfile[varexts[j]].header['WAT3_001'] != \
              'wtype=linear axtype=xi':
                header_check = 'fail'
            if objfile[varexts[j]].header['WCSDIM'] != 3:
                header_check = 'fail'
            if objfile[varexts[j]].header['NAXIS'] != 2:
                header_check = 'fail'
            if (objfile[varexts[j]].header['CRVAL2'] != \
                    objfile[sciexts[j]].header['CRVAL2']) and \
                    (objfile[varexts[j]].header['CRPIX2'] == \
                         objfile[sciexts[j]].header['CRPIX2']):
                header_check = 'fail'
            if objfile[varexts[j]].header['CRVAL3'] != \
              objfile[sciexts[j]].header['CRVAL3']:
                header_check = 'fail'

            #for the dq extensions. the original header isn't quite
            #the same as the science and variance extensions, so
            #modify/add other keywords too.
            objfile[dqexts[j]].header['WCSDIM'] = 3
            objfile[dqexts[j]].header['CTYPE2'] = 'LINEAR'
            objfile[dqexts[j]].header['CTYPE3'] = 'LINEAR'
            objfile[dqexts[j]].header['CRPIX2'] = \
              objfile[varexts[j]].header['CRPIX2']
            objfile[dqexts[j]].header['CRVAL2'] = \
              objfile[varexts[j]].header['CRVAL2']
            objfile[dqexts[j]].header['CD2_2'] = \
              objfile[varexts[j]].header['CD2_2']
            objfile[dqexts[j]].header['CRVAL3'] = \
              objfile[varexts[j]].header['CRVAL3']
            objfile[dqexts[j]].header['CD3_3'] = \
              objfile[varexts[j]].header['CD3_3']
            objfile[dqexts[j]].header['WAT1_001'] = 'wtype=linear axtype=wave'
            objfile[dqexts[j]].header['WAT2_001'] = 'wtype=linear axtype=eta'
            objfile[dqexts[j]].header['WAT3_001'] = 'wtype=linear axtype=xi'
            objfile[dqexts[j]].header['WAXMAP01'] = \
              objfile[varexts[j]].header['WAXMAP01']
            objfile[dqexts[j]].header['GEMMASK'] = 'mask/'+ \
              (obj_exps[i].split('.'))[0]+'_sci'+str(sciexts[j])+'_var'+\
              str(varexts[j])+'_dq'+str(dqexts[j])+'.pl'
            del objfile[dqexts[j]].header['LTM1_1']
            del objfile[dqexts[j]].header['LTM2_2']
            del objfile[dqexts[j]].header['WAT0_001']

            #create a mask based on the DQ extension. just want to
            #know the pixels on the edges that are bad.
            dq_img = objfile[dqexts[j]].data
            mask = np.zeros((dq_img.shape[0],dq_img.shape[1]))
            #examine columns first
            for k in range(dq_img.shape[1]):
                index = np.where(dq_img[:,k] == 1)
                #if more than 50% of the pixels in this column are
                #flagged as bad, then copy to gemcube mask
                if (np.shape(index))[1] > (0.50*dq_img.shape[0]):
                    mask[index,k] = dq_img[index,k]
            #examine rows now
            for k in range(dq_img.shape[0]):
                index = np.where(dq_img[k,:] == 1)
                #if more than 50% of the pixels in this row are
                #flagged as bad, then copy to gemcube mask
                if (np.shape(index))[1] > (0.50*dq_img.shape[1]):
                    mask[k,index] = dq_img[k,index]

            #create the .pl file.
            pyfits.writeto('mask_tmp.fits',mask)
            iraf.imcopy('mask_tmp.fits', 'mask/'+ \
                        (obj_exps[i].split('.'))[0]+'_sci'+str(sciexts[j])+ \
                        '_var'+str(varexts[j])+'_dq'+str(dqexts[j])+'.pl',
                        verbose='no')
            tmp=subprocess.call(['rm','mask_tmp.fits'],
                                stderr=open(os.devnull,'w'))
        
        #update the hatfbrsn*.fits file
        if header_check == 'fail':
            sys.exit('Problem with header for '+'h'+obj_exps[i])
        objfile.flush(output_verify='silentfix')
        objfile.close(output_verify='silentfix')

    
    ###########################################################################
    # STEP 3:  Create the Merged Cube                                         #
    ###########################################################################

    #create a list of object hatfbrsn*.fits files
    if intype == 'gal' or intype == 'psf':
        hlistname = 'hlist'
    if intype == 'sky':
        hlistname = 'hlist_sky'

    hlist = open(workdir+hlistname,'w')
    for i in range(len(obj_exps)):
        outstr = '%s\n' % ('h'+(obj_exps[i].split('.'))[0])
        hlist.write(outstr)
    hlist.close()
        
    #nifcube_jonelle takes only 1 set of extensions at time (e.g., the
    #science, variance, or dq)

    if intype == 'gal' or intype == 'psf':
        filenameend1 = '_sci.fits'
        filenameend2 = '_var.fits'
        filenameend3 = '_dq.fits'
        filenameend_weight = '_sci_weight'
        filenameend2_weight = '_var_weight'
        filenameend3_weight = '_dq_weight'
    if intype =='sky':
        filenameend1 = '_sky_sci.fits'
        filenameend2 = '_sky_var.fits'
        filenameend3 = '_sky_dq.fits'
        filenameend_weight = '_sky_sci_weight'
        filenameend2_weight = '_sky_var_weight'
        filenameend3_weight = '_sky_dq_weight'
    if intype == 'gal' or intype == 'sky':
        outfilebegin = galname
    if intype == 'psf':
        outfilebegin = galname + '_psf'
        
    
    #science
    iraf.nifcube_jonelle(inimages='@'+hlistname, inext='SCI',
                         outphu='h'+(obj_exps[0].split('.'))[0],
                         outcubes=outfilebegin+filenameend1,
                         outweight=outfilebegin+filenameend_weight+'.fits',
                         waxis=3, saxis=2, taxis=1, wmin=refwave[0],
                         wmax=refwave[1], wscale=refwave[2], sscale=0.05,
                         blank=0.0, logfile=log, verbose='yes')
    
    #variance
    iraf.nifcube_jonelle(inimages='@'+hlistname, inext='VAR',
                         outphu='h'+(obj_exps[0].split('.'))[0],
                         outcubes=outfilebegin+filenameend2,
                         outweight=outfilebegin+filenameend2_weight+'.fits',
                         waxis=3, saxis=2, taxis=1, wmin=refwave[0],
                         wmax=refwave[1], wscale=refwave[2], sscale=0.05,
                         blank=0.0, logfile=log, verbose='yes')
    #variance cube should be scaled as 1/weight^2, whereas gemcube
    #(assuming science) has scaled as 1/weight so divide by weight one
    #more time to get correct variance cube.
    weightfile = pyfits.open(outfilebegin+filenameend2_weight+'.fits')
    cube_weight = weightfile[0].data
    outvarfile = pyfits.open(outfilebegin+filenameend2, mode='update')
    cube_var = outvarfile[1].data
    #loop over each element of cube
    for i in range(cube_weight.shape[2]):
        for j in range(cube_weight.shape[1]):
            for k in range(cube_weight.shape[0]):
                if cube_weight[k,j,i] != 0:
                    cube_var[k,j,i] = cube_var[k,j,i] / cube_weight[k,j,i]
                else:
                    cube_var[k,j,i] = 0.
    weightfile.close()
    outvarfile.flush(output_verify='silentfix')
    outvarfile.close()
    
    #dq
    iraf.nifcube_jonelle(inimages='@'+hlistname, inext='DQ',
                         outphu='h'+(obj_exps[0].split('.'))[0],
                         outcubes=outfilebegin+filenameend3,
                         outweight=outfilebegin+filenameend3_weight+'.fits',
                         waxis=3, saxis=2, taxis=1, wmin=refwave[0],
                         wmax=refwave[1], wscale=refwave[2], sscale=0.05,
                         blank=0.0, logfile=log, verbose='yes')

    
    ###########################################################################
    # STEP 4:  Join Science, Variance, DQ Merged Cubes and Fix Header         #
    ###########################################################################

    #create a single fits file with a science, variance, and dq
    #extension.

    if intype == 'gal' or intype == 'sky':
        filenamebegin = galname
    if intype == 'psf':
        filenamebegin = galname+'_psf'
    if intype == 'gal' or intype =='psf':
        filenameend = '_combined.fits'
    if intype == 'sky':
        filenameend = '_sky_combined.fits'

        
    tmp=subprocess.call(['cp',filenamebegin+filenameend1,
                        filenamebegin+filenameend], stderr=open(os.devnull,'w'))
    
    tmp = pyfits.open(filenamebegin+filenameend2)
    pyfits.append(filenamebegin+filenameend,tmp[1].data,tmp[1].header)
    tmp.close()
    tmp = pyfits.open(filenamebegin+filenameend3)
    pyfits.append(filenamebegin+filenameend,tmp[1].data,tmp[1].header)
    tmp.close()

    #get the x and y pix reference and RA and Dec value
    ref_headinfo = np.loadtxt(headerinfo,dtype='float',comments='#')

    #correct headers
    objfile = pyfits.open(filenamebegin+filenameend, mode='update')
    #primary header
    objfile[0].header['CRPIX1'] = ref_headinfo[0]
    objfile[0].header['CRVAL1'] = ref_headinfo[2]
    objfile[0].header['CRPIX2'] = ref_headinfo[1]
    objfile[0].header['CRVAL2'] = ref_headinfo[3]
    objfile[0].header['PIXSCALE'] = 0.05

    #science header
    objfile[1].header['PIXSCALE'] = 0.05
    del objfile[1].header['DATE']
    del objfile[1].header['IRAF-TLM']
    del objfile[1].header['NFPAD']
    del objfile[1].header['GEM-TLM']
    del objfile[1].header['NSCUTSEC']
    del objfile[1].header['NSCUTSPC']
    #del objfile[1].header['FIXPIX']
    #del objfile[1].header['FIXPIX02']
    #del objfile[1].header['CRMASK']
    #del objfile[1].header['FIXPIX03']
    del objfile[1].header['GEMMASK']

    #variance header
    objfile[2].header['PIXSCALE'] = 0.05
    del objfile[2].header['DATE']
    del objfile[2].header['IRAF-TLM']
    del objfile[2].header['NFPAD']
    del objfile[2].header['GEM-TLM']
    del objfile[2].header['NSCUTSEC']
    del objfile[2].header['NSCUTSPC']
    #del objfile[2].header['FIXPIX']
    #del objfile[2].header['FIXPIX02']
    #del objfile[2].header['CRMASK']
    #del objfile[2].header['FIXPIX03']
    del objfile[2].header['GEMMASK']

    #dq header
    objfile[3].header['PIXSCALE'] = 0.05
    del objfile[3].header['DATE']
    del objfile[3].header['IRAF-TLM']
    del objfile[3].header['NFPAD']
    del objfile[3].header['GEM-TLM']
    del objfile[3].header['NSCUTSEC']
    del objfile[3].header['NSCUTSPC']
    #del objfile[3].header['FIXPIX']
    #del objfile[3].header['FIXPIX02']
    #del objfile[3].header['CRMASK']
    #del objfile[3].header['FIXPIX03']
    del objfile[3].header['GEMMASK']
    
    objfile.flush(output_verify='silentfix')
    objfile.close(output_verify='silentfix')

    #remove the individual sci, var, dq fits files
    tmp=subprocess.call(['rm',filenamebegin+filenameend1],
                            stderr=open(os.devnull,'w'))
    tmp=subprocess.call(['rm',filenamebegin+filenameend2],
                            stderr=open(os.devnull,'w'))
    tmp=subprocess.call(['rm',filenamebegin+filenameend3],
                            stderr=open(os.devnull,'w'))
    tmp=subprocess.call(['rm',outfilebegin+filenameend_weight+'.fits'],
                            stderr=open(os.devnull,'w'))
    tmp=subprocess.call(['rm',outfilebegin+filenameend2_weight+'.fits'],
                            stderr=open(os.devnull,'w'))
    tmp=subprocess.call(['rm',outfilebegin+filenameend3_weight+'.fits'],
                            stderr=open(os.devnull,'w'))
