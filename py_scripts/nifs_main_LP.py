"""

Created on June 11 2018
Author: Jonelle Walsh, Anil Seth, Richard McDermid, etc.

See setup_pipeline.pdf for installation instructions, and xx.pdf for
xx useage instructions.

To run, go to the rootdir (see below), open pyraf, then type:
---> pyexecute("[pyscriptpath (see below)]/nifs_main_LP.py")

Pieces of the code cannot yet deal with more than 1 observational
setup at a time. If needing to reduce more than 1 observational setup,
create a datadir (see below) holding raw files for 1 observational
setup, and a separate datadir (see below) holding raw files for a
different observational setup. Run the code separately on each
datadir.

The code calls PyRAF/IRAF and tasks within the Gemini data reduction
package. Thus this code must be run with Python 2.7. The code has been
tested with Gemini IRAF package version 1.14 and IRAF 2.16 on NIFS
data with the K grating and HK filter centered on 2.2 microns.

"""

################################################################################
################################################################################
# EDIT BELOW

#give the location of where IDL lives on your computer
idlpath = '/Applications/exelis/idl85/bin/idl'
#enter the location of the data reduction python scripts
#pyscriptspath = '/Users/jlwalsh/Data/LP_2016/nifs_reduction_info/'+\
#  'full_reduce/nifs_pipeline_v4/py_scripts/'
pyscriptspath = '/Users/jlwalsh/Dropbox/LP_reduction_pipeline/py_scripts/'
#directory where anil_ArXe_K.dat and kurucz star files are located.
#refdir = '/Users/jlwalsh/Data/LP_2016/nifs_reduction_info/'+\
#  'full_reduce/nifs_pipeline_v4/ref_files/'
refdir = '/Users/jlwalsh/Dropbox/LP_reduction_pipeline/ref_files/'

#location of the main directory. make sure to include the / at the end
#of the string. NEEDS TO EXIST before running pipeline.
rootdir = '/Users/jlwalsh/Data/LP_2016/nifs_data/2016b/'
#directory where the raw data is stored. make sure to include the / at
#the end of the string. inside of this directory, must have
#subdirectories holding the raw data broken up by date
#(datadir+20160920/). datadir and the subdirectories NEED TO EXIST
#before running pipeline.
datadir = 'raw_data/'

#directory that will hold the reduced data. make sure to include the /
#at the end of the string. will be created by the pipeline if doesn't
#already exist.
reducedir = 'nifs_pipeline_v4/reduced_data/'

#provide the dates for which you want to reduce data. enter as a
#string and use a yyyymmdd format. a single date or multiple dates can
#be entered. (e.g., dates = ['20121227'] or dates =
#['20121227','20121230'])
dates = ['20160921']
#galaxies observed during ANY of the dates listed above. use full name
#and lower case (e.g., ngc, ugc, pgc, ic, mrk). this will be used to
#name the directoires.
#galaxies = ['ngc7242','ugc11537','pgc12557'] #for 20160921
galaxies = ['ugc11537']
#galaxies = ['ugc11537','ngc772'] #for 20160920
#galaxies = ['ngc7242','ngc1589'] #for 20160919
#galaxies = ['ngc1589'] #for 20160918
#galaxies = ['ngc7242'] #for 20160928
#galaxies = ['ngc7242','ngc1022'] #for 20160925
#tellurics observed during ANY of the dates listed above. use full
#name and lower case (e.g., hip, hd, hr). will be used to name the
#directories.
#tellurics = ['hip105437','hip95793','hd203769','hip10559','hip18769'] #for 20160921
#tellurics = ['hip95793','hd203769','hip1603','hd19600'] #for 20160920
#tellurics = ['hip105437','hip114714','hip13917','hip24607'] #for 20160919
#tellurics = ['hip13917','hip24607']  #for 20160918
#tellurics = ['hip105437','hip114714'] #for 20160928
tellurics = ['hd123233','*8her']

#you can chose to do some parts of the reduction and not others. to do
#the step enter 'yes', to skip the step enter 'no'.
reduce_sortfiles = 'no'      #create directory structure and sort files
reduce_daycals = 'no'        #reduce daycals
reduce_tellurics = 'no'      #reduce telluric stars
reduce_galaxies = 'no'       #basic reduction of galaxy exposures
reduce_combine_gal = 'no'    #merge galaxy cubes
reduce_psf = 'no'            #reduce psf star observations
reduce_combine_psf = 'yes'    #merge psf cubes
reduce_fluxcal = 'no'        #flux calibrate the merged galaxy and sky cubes
reduce_measurelsf = 'no'     #measure the line-spread function from the
                              #merged sky cube

#do you want to do the wavelength calibration interactively, or accept
#the default fits? Set to 'yes' to do interactively.
flinter_arc = 'no'
#do you want to do the s-distortion calibration interactively, or
#accept the default fits? Set to 'yes' to do interactively?
flinter_sdist = 'no'
#do you want to determine the 2D to 3D spatial/spectral transformation
#interactively, or accept the default fitting function? Set to 'yes'
#to do interactively.
flinter_nsfitcoords = 'no'
#do you want to extract 1D spectra interactively (for tellurics), or
#accept the (x,y) pixel with the maximum value in the collapsed data
#cube? Set to 'yes' to do interactively.
flinter_extract = 'no'
#do you want to apply the telluric correction interactively
#(extracting a 1D spectrum, searching for values of the shift/scale
#applied to the telluric)? Set to 'yes' to do interactively.
flinter_telluric = 'no'

#clean cosmic rays and bad negataive pixels. probably can leave the
#numbers below unchanged, but be sure to check the masks in
#rootdir+reducedir+[galaxies]+[dates]+[setup]+'clean/masks_lacosmic/'
#and +'clean/masks_negpix/'. cosmic rays/bad pixels should not have
#been identified on top of any real signal. for parameters below,
#larger numbers correspond to being more conservative
sigclip_in = 10.0  #LA COSMIC sigclip
sigfrac_in = 3.0   #LA COSMIC sigfrac
objlim_in = 5.0    #LA COSMIC objlim
nsigneg_in = 5.0   #XZAP nsigneg (sigma for finding negative pixels)

#do you want to generate barycentric or heliocentric wavelengths, or
#no correction? for corr_type, options are 'none', 'bary', or
#'helio'. note, if no correction is applied, spectra are still sampled
#to the same starting/ending wavelength with the same delta lambda
#when a single merged cube is created.
corr_type = 'none'


# STOP EDITTING
################################################################################
################################################################################

# Import some useful python utilities.

import sys
sys.path.append(pyscriptspath)
import os
import glob
import astropy.io.fits as pyfits
import numpy as np
import time
import subprocess
import pidly
idl = pidly.IDL(idlpath)
from nifs_checkdaycal_LP import nifs_checkdaycal_LP
from nifs_checktel_LP import nifs_checktel_LP
from nifs_checkgal_LP import nifs_checkgal_LP
from nifs_checkpsf_LP import nifs_checkpsf_LP
from nifs_checkdata_LP import nifs_checkdata_LP
from nifs_basecalib_LP import nifs_basecalib_LP
from nifs_telluric_LP import nifs_telluric_LP
from nifs_galaxy_LP import nifs_galaxy_LP
from nifs_checkdata_merged_LP import nifs_checkdata_merged_LP
from nifs_refwave_LP import nifs_refwave_LP
from nifs_getoffmethod_LP import nifs_getoffmethod_LP
from nifs_merge_LP import nifs_merge_LP
from nifs_psf_LP import nifs_psf_LP
from nifs_findtelluric_fluxcalib_LP import nifs_findtelluric_fluxcalib_LP
from nifs_worktelluric_fluxcalib_LP import nifs_worktelluric_fluxcalib_LP

################################################################################

# Check that the user has selected at least one data reduction task to
# complete.

if reduce_sortfiles.lower() == 'no' and reduce_daycals.lower() == 'no' and \
    reduce_tellurics.lower() == 'no' and reduce_galaxies.lower() == 'no' and \
    reduce_combine_gal.lower() == 'no' and reduce_psf.lower() == 'no' and \
    reduce_combine_psf.lower() == 'no' and reduce_fluxcal.lower() == 'no' and \
    reduce_measurelsf.lower() == 'no':
    sys.exit('You need to set one of the reduction steps to "yes".')

################################################################################

# Run the IDL tasks to write important info to text files, create
# directory structure, and sort files.

if reduce_sortfiles.lower() == 'yes':

    print('')
    print('Gathering info...')
    print('')
    idl.pro('setup_writeinfo', rootdir, rootdir+datadir, dates)
    print('')
    print('Setting up directories...')
    print('')
    idl.pro('setup_makedirectories', rootdir, rootdir+datadir,
            rootdir+reducedir, dates, galaxies, tellurics)

################################################################################

# Check the daycal/telluric/galaxy/sky/psf lists that were created
# above and allow user to modify the lists if needed. At this time,
# also create a telluric correction list for each of the galaxy and
# psf exposures. Create a file that lists each telluric star and their
# flux density for flux calibration.

if reduce_daycals.lower() == 'yes':
    nifs_checkdaycal_LP(rootdir, reducedir, datadir, dates)
    
if reduce_tellurics.lower() == 'yes':
    nifs_checktel_LP(rootdir, reducedir, datadir, dates)

if reduce_galaxies.lower() == 'yes':
    nifs_checkgal_LP(rootdir, reducedir, datadir, galaxies, dates)

if reduce_psf.lower() == 'yes':
    nifs_checkpsf_LP(rootdir, reducedir, datadir, galaxies, dates)

################################################################################

# Start the real work.
    
#loop over each night to reduce the daycals and the tellurics
for a in range(len(dates)):

    if reduce_daycals.lower() == 'yes' or reduce_tellurics.lower() == 'yes': 

        os.chdir(rootdir+reducedir+'daycals/'+dates[a])
        obs_setups = glob.glob('*')

        for b in range(len(obs_setups)):

            #work on the baseline calibrations
            workdir = rootdir+reducedir+'daycals/'+dates[a]+'/'+\
              obs_setups[b]+'/'
            os.chdir(workdir)
        
            #the names of the daycal lists
            flatlist = 'flatlist'
            flatdarklist = 'flatdarklist'
            arcdarklist = 'arcdarklist'
            ronchilist = 'ronchilist'
            ronchidarklist = 'ronchidarklist'
            arclist = glob.glob('arclist*')
            darklist = glob.glob('darklist*')

            ####################################################################
        
            if reduce_daycals.lower() == 'yes':

                print('')
                print('Starting the baseline calibrations '+\
                    'for %s and %s' % (dates[a], obs_setups[b]))
                print('')

                #check to see if outputs of data reduction already
                #exist. if so, ask the user whether to proceed. if the
                #user does want to proceed, delete the previous
                #outputs and start over.
                nifs_checkdata_LP(workdir,'basecalib*log*')

                #if refdir is too long, pyraf won't be able to find
                #the reference file, so temporarily copy the arc line
                #list to the working directory
                tmp=subprocess.call(['cp',refdir+'/anil_ArXe_K.dat',workdir],
                                    stderr=open(os.devnull,'w'))
        
                #reduce the baseline calibrations
                nifs_basecalib_LP(workdir, dates[a], flatlist, flatdarklist,
                                arclist, arcdarklist, ronchilist,
                                ronchidarklist, darklist, refdir, flinter_arc,
                                flinter_sdist)

                #remove the copy of the arc line list now that the
                #daycals have been reduced
                tmp=subprocess.call(['rm',workdir+'/anil_ArXe_K.dat'],
                                    stderr=open(os.devnull,'w'))

            ####################################################################

            if reduce_tellurics.lower() == 'yes':
            
                os.chdir(rootdir+reducedir+'tellurics/'+dates[a]+'/'+\
                        obs_setups[b]+'/')
                telluric_stars = glob.glob('*')
        
                for c in range(len(telluric_stars)):
        
                    workdir = rootdir+reducedir+'tellurics/'+dates[a]+'/'+\
                    obs_setups[b]+'/'+telluric_stars[c]+'/'
                    caldir = rootdir+reducedir+'daycals/'+dates[a]+'/'+\
                      obs_setups[b]+'/'

                    print('')
                    print('Starting the '+\
                        'telluric %s for %s and %s' % (telluric_stars[c], \
                                                           dates[a], \
                                                           obs_setups[b]))
                    print('')

                    os.chdir(workdir)

                    #check to see if outputs of data reduction already
                    #exist. if so, ask the user whether to proceed. if
                    #the user does want to proceed, delete the
                    #previous outputs and start over.
                    nifs_checkdata_LP(workdir,'telluric*log*')

                    #the names of the telluric and sky lists
                    telluriclist = glob.glob('telluriclist*')
                    skylist = glob.glob('skylist_*')
                    skylistshort = glob.glob('skylistshort_*')

                    #remind the user to open a ds9 window if needed
                    if flinter_extract == 'yes':
                        response = \
                        raw_input('Need to have a DS9 window open. Hit '+\
                                    'any key to continue. ')
                
                    #reduce the telluric stars
                    nifs_telluric_LP(workdir, caldir, dates[a], flatlist,
                                     arclist, ronchilist, telluriclist, skylist,
                                     skylistshort, flinter_nsfitcoords,
                                     flinter_extract)

                    #remove absorption lines from the telluric star
                    #and correct for the blackbody shape. the output
                    #1D spectrum is ready to be fed into nftelluric to
                    #make the telluric correction.
                    telluric_file = glob.glob('gxtfbrsn*.fits')
                    idl.pro('fit_telluric_richard',workdir,telluric_file,refdir)

################################################################################

#loop over each galaxy to reduce the galaxy exposures and merge cubes
for a in range(len(galaxies)):
    
    for b in range(len(dates)):

        #dates contains all the days data are to be reduced, but this
        #galaxy may not have been observed on this date. check here.
        if os.path.isdir(rootdir+reducedir+galaxies[a]+'/'+dates[b]):

            os.chdir(rootdir+reducedir+galaxies[a]+'/'+dates[b])
            obs_setups = glob.glob('*')

            for c in range(len(obs_setups)):

                workdir = rootdir+reducedir+galaxies[a]+'/'+dates[b]+'/'+\
                  obs_setups[c]+'/'
                caldir = rootdir+reducedir+'daycals/'+dates[b]+'/'+\
                  obs_setups[c]+'/'
                teldir = rootdir+reducedir+'tellurics/'+dates[b]+'/'+\
                  obs_setups[c]+'/'

                #the names of the daycal lists. re-search for arclist
                #and darklist in case the user made changes.
                flatlist = 'flatlist'
                ronchilist = 'ronchilist'
                os.chdir(caldir)
                arclist = glob.glob('arclist*')
                darklist = glob.glob('darklist*')
                #the names of the galaxy, sky, and telluric correction
                #lists
                os.chdir(workdir)
                gallist = 'gallist'
                skylist = 'skylist'
                skylistshort = 'skylistshort'
                telcorrlist = 'telcorrlist'

                ################################################################
                
                if reduce_galaxies.lower() == 'yes':

                    print('')
                    print('Starting the '+\
                          'galaxy %s for %s and %s' % (galaxies[a], \
                                                       dates[b], obs_setups[c]))
                    print('')

                    #check to see if outputs of data reduction already
                    #exist. if so, ask the user whether to proceed. if
                    #the user does want to proceed, delete the
                    #previous outputs and start over.
                    nifs_checkdata_LP(workdir,galaxies[a]+'*log*')

                    #remind the user to open a ds9 window if needed
                    if flinter_telluric == 'yes':
                        response = \
                        raw_input('Need to have a DS9 window open. Hit '+\
                                  'any key to continue. ')
                
                    #reduce the galaxy exposures
                    nifs_galaxy_LP(workdir, caldir, teldir, dates[b], flatlist,
                                   arclist, ronchilist, darklist, telcorrlist,
                                   galaxies[a], gallist, skylist, skylistshort,
                                   sigclip_in, sigfrac_in, objlim_in,
                                   nsigneg_in, rootdir, reducedir,
                                   obs_setups[c], flinter_nsfitcoords,
                                   flinter_telluric)
                    
                ################################################################

    if reduce_combine_gal.lower() == 'yes':

        os.chdir(rootdir+reducedir+galaxies[a]+'/merged/')
        obs_setups = glob.glob('*')

        for b in range(len(obs_setups)):

            #set the working directory
            workdir = rootdir+reducedir+galaxies[a]+'/merged/'+obs_setups[b]+'/'
        
            #may have a combined cube from another night, but now have
            #an additional night and want to combine
            #everything. remove any files that are not
            #'image_catfbrsn*.fits' and atfbrsn*.fits and start the
            #process over, but check with the user first.
            nifs_checkdata_merged_LP(workdir,galaxies[a]+'*log*')

            #check that at least one atbrsn*.fits file exists
            os.chdir(workdir)
            gal_images = glob.glob('image_catfbrsn*.fits')
            if len(gal_images) == 0:
                sys.exit('Cannot find any image_catfbrsn*.fits '+\
                         'files in '+workdir+'.')
            #get the galaxy exposure file names (atfbrsn*.fits) from
            #the image file names
            gal_exps = []
            for c in range(len(gal_images)):
                gal_exps.append( (gal_images[c])[7:] )

            #get the sky exposure names, if they exist
            sky_exps = glob.glob('atfbrsn*_*.fits')

            #get the reference wavelength, which includes making a
            #barycentric or heliocentric correction (if the user
            #wants). even if no correction is done, the spectra are
            #all sampled to the same starting/ending wavelength with
            #the same delta lambda when creating a single merged
            #cube. after a correction is applied (if any) determine
            #the minimum, maximum, and wavelength scale for the
            #spectra.
            refwave_gal = nifs_refwave_LP(workdir, gal_exps[0], idlpath,
                                          corr_type)

            #determine the spatial offsets from the collapsed galaxy
            #cubes using the first file as the reference
            idl.pro('find_offsets_lp', workdir, gal_images)

            #set the method for determining offsets between cubes
            #(spaxel with maximum value, fitting 2D Gaussian,
            #cross-correlating the collapsed cubes)
            offsetresult = nifs_getoffmethod_LP(workdir, len(gal_images))
    
            #merge galaxy exposures
            nifs_merge_LP(workdir, gal_exps, galaxies[a], offsetresult[0],
                          offsetresult[1], refwave_gal, 'gal')
            #merge sky exposures using the galaxy spaital offsets and
            #the galaxy reference wavelength array
            if len(sky_exps) > 0:
                nifs_merge_LP(workdir, sky_exps, galaxies[a], offsetresult[0],
                              offsetresult[1], refwave_gal, 'sky')

################################################################################

#reduce psf star 
if reduce_psf.lower() == 'yes':

    for a in range(len(dates)):

        #dates contains all the days data are to be reduced, but this
        #psf may not have been observed on this date. check that here.
        if os.path.isdir(rootdir+reducedir+'psfs/'+dates[a]):

            os.chdir(rootdir+reducedir+'psfs/'+dates[a])
            obs_setups = glob.glob('*')
    
            for b in range(len(obs_setups)):

                    for c in range(len(galaxies)):

                        #a psf may or may not exist for each galaxy, so check.
                        if os.path.isdir(rootdir+reducedir+'psfs/'+dates[a]+\
                                             '/'+obs_setups[b]+'/psf_'+\
                                             galaxies[c]):

                            #set the psf work directory
                            workdir = rootdir+reducedir+'psfs/'+dates[a]+'/'+\
                              obs_setups[b]+'/psf_'+galaxies[c]+'/'
                            caldir = rootdir+reducedir+'daycals/'+dates[a]+'/'+\
                              obs_setups[b]+'/'
                            teldir = rootdir+reducedir+'tellurics/'+dates[a]+\
                              '/'+obs_setups[b]+'/'

                            #the names of the daycal lists; ; re-search for
                            #arclist and darklist in case the user made changes
                            flatlist = 'flatlist'
                            ronchilist = 'ronchilist'
                            os.chdir(caldir)
                            arclist = glob.glob('arclist*')
                            #the names of the psf, sky, and telluric correction
                            #lists
                            os.chdir(workdir)
                            psflist = 'psflist'
                            skylist = 'skylist'
                            skylistshort = 'skylistshort'
                            telcorrlist = 'telcorrlist'

                            gpr = galaxies[c]
                            opr = obs_setups[b]
                            dpr = dates[a]
                            print('Starting the psf for %s for %s and %s' \
                                      % (gpr,dpr,opr))

                            #check to see if outputs of data reduction
                            #already exist. if so, ask the user
                            #whether to proceed. if the user does want
                            #to proceed, delete the previous outputs
                            #and start over.
                            nifs_checkdata_LP(workdir,'psf_'+galaxies[c]+\
                                                  '*log*')
                
                            #reduce the psf exposures
                            nifs_psf_LP(workdir, caldir, teldir, dates[a],
                                        flatlist, arclist, ronchilist,
                                        telcorrlist, galaxies[c], psflist,
                                        skylist, skylistshort, sigclip_in,
                                        sigfrac_in, objlim_in, nsigneg_in,
                                        rootdir, reducedir, obs_setups[b],
                                        flinter_nsfitcoords, flinter_telluric)

                            #for this single night/observational setup
                            #combine into a single psf cube

                            #check that at least one atbrsn*.fits file exists
                            os.chdir(workdir)
                            psf_images = glob.glob('image_catfbrsn*.fits')
                            if len(psf_images) == 0:
                                sys.exit('Cannot find any '+\
                                         'image_catfbrsn*.fits '+\
                                         'files in '+workdir+'.')

                            #get the psf exposure file names (atfbrsn*.fits)
                            #from the image file names
                            psf_exps = []
                            for d in range(len(psf_images)):
                                psf_exps.append( (psf_images[d])[7:] )

                            #get the reference wavelength, which
                            #includes making a barycentric or
                            #heliocentric correction (if the user
                            #wants). even if no correction is done,
                            #the spectra are all sampled to the same
                            #starting/ending wavelength with the same
                            #delta lambda when creating a single
                            #merged cube. after a correction is
                            #applied (if any) determine the minimum,
                            #maximum, and wavelength scale for the
                            #spectra.
                            refwave_psf = nifs_refwave_LP(workdir, psf_exps[0],
                                                          idlpath, corr_type)

                            #determine the spatial offsets from the collapsed
                            #psf cubes using the first file as the reference
                            idl.pro('find_offsets_lp', workdir, psf_images)
                
                            #set the method for determining offsets between
                            #cubes (spaxel with maximum value, fitting 2D
                            #Gaussian, cross-correlating the collapsed cubes)
                            offsetresult = nifs_getoffmethod_LP(workdir,
                                                                len(psf_images))
    
                            #merge psf exposures
                            nifs_merge_LP(workdir, psf_exps, galaxies[c],
                                          offsetresult[0], offsetresult[1],
                                          refwave_psf, 'psf')

################################################################################

#merge psf star observations
if reduce_combine_psf.lower() == 'yes':

    os.chdir(rootdir+reducedir+'psfs/merged')
    obs_setups = glob.glob('*')

    for a in range(len(obs_setups)):
            
        for b in range(len(galaxies)):

            #a psf may or may not exist for each galaxy, so check.
            if os.path.isdir(rootdir+reducedir+'psfs/merged/'+\
                             obs_setups[a]+'/psf_'+galaxies[b]):

                #set the working directory
                workdir = rootdir+reducedir+'psfs/merged/'+obs_setups[a]+\
                  '/psf_'+galaxies[b]+'/'
        
                #may have a combined cube from another night, but now
                #have an additional night and want to combine
                #everything. remove any files that are not
                #'image_catfbrsn*.fits' and atfbrsn*.fits and start
                #the process over, but check with the user first.
                nifs_checkdata_merged_LP(workdir,galaxies[b]+'_psf*log*')

                #check that at least one atbrsn*.fits file exists
                os.chdir(workdir)
                psf_images = glob.glob('image_catfbrsn*.fits')
                if len(psf_images) == 0:
                    sys.exit('Cannot find any image_catfbrsn*.fits '+\
                                 'files in '+workdir+'.')
                      
                #get the psf exposure file names (atfbrsn*.fits) from
                #the image file names
                psf_exps = []
                for c in range(len(psf_images)):
                    psf_exps.append( (psf_images[c])[7:] )

                #get the reference wavelength, don't make a
                #barycentric or heliocentric correction. determine the
                #minimum, maximum, and wavelength scale for the
                #spectra.
                refwave_psf = nifs_refwave_LP(workdir, psf_exps[0], idlpath,
                                              corr_type)

                #determine the spatial offsets from the collapsed psf
                #cubes using the first file as the reference
                idl.pro('find_offsets_lp', workdir, psf_images)
                
                #set the method for determining offsets between cubes
                #(spaxel with maximum value, fitting 2D Gaussian,
                #cross-correlating the collapsed cubes)
                offsetresult = nifs_getoffmethod_LP(workdir, len(psf_images))
    
                #merge psf exposures
                nifs_merge_LP(workdir, psf_exps, galaxies[b], offsetresult[0],
                                offsetresult[1], refwave_psf, 'psf')


################################################################################

#flux calibrate the merged galaxy cube and the sky cube
if reduce_fluxcal.lower() == 'yes':

    for a in range(len(galaxies)):

        os.chdir(rootdir+reducedir+galaxies[a]+'/merged')
        obs_setups = glob.glob('*')

        for b in range(len(obs_setups)):
            
            workdir = rootdir+reducedir+galaxies[a]+'/merged/'+obs_setups[b]+'/'

            #need to track down one of the telluric stars used for the
            #galaxy.
            telstar_result = nifs_findtelluric_fluxcalib_LP(rootdir+reducedir,
                                                            workdir,
                                                            galaxies[a],
                                                            obs_setups[b])
            telstar_calib_dir = telstar_result[0]
            date_telstar = telstar_result[1]
            telstar_calib_file = 'N'+date_telstar+telstar_result[2]
            telstar_name = (telstar_calib_dir.split('/'))[-2]
                                                            
            #work with the telluric star. need to telluric correct the
            #telluric star and then generate a cube.
            nifs_worktelluric_fluxcalib_LP(telstar_calib_dir,
                                           telstar_calib_file, galaxies[a],
                                           date_telstar, flinter_telluric)

            #determine the multiplicative conversion factor to convert
            #from counts/s to ergs/s/cm^2/A. then apply to the merged
            #galaxy and sky cubes.
            idl.pro('measurephot_calibcube', telstar_calib_dir, 'catfbrsn'+\
                    telstar_calib_file+'.fits', rootdir, reducedir,
                    galaxies[a], obs_setups[b], telstar_name)

                    
################################################################################

#measure the LSF from the merged sky cube for each
#galaxy/observational setup
if reduce_measurelsf.lower() == 'yes':

    for a in range(len(galaxies)):

        os.chdir(rootdir+reducedir+galaxies[a]+'/merged')
        obs_setups = glob.glob('*')

        for b in range(len(obs_setups)):
            
            workdir = rootdir+reducedir+galaxies[a]+'/merged/'+obs_setups[b]+'/'
            os.chdir(workdir)

            sky_flux = glob.glob(galaxies[a]+'_sky_combined_flux.fits')
            sky_noflux = glob.glob(galaxies[a]+'_sky_combined.fits')

            if len(sky_flux) == 1 and len(sky_noflux) == 1:
                skycube = galaxies[a]+'_sky_combined_flux.fits'
            if len(sky_flux) == 0 and len(sky_noflux) == 1:
                skycube = galaxies[a]+'_sky_combined.fits'
            if len(sky_flux) == 0 and len(sky_noflux) == 0:
                sys.exit('There is no merged sky cube in '+workdir+\
                         ' from which to measure the LSF.')

            print('')
            print('Measuring the LSF as a function of spatial location from'+\
                      ' the merged sky cube...')

            #fit Gaussians to the sky lines
            idl.pro('measure_lsf_lp',workdir,skycube,refdir,galaxies[a])

            print('The LSF map is called '+\
                  '%s_fullskycube_fwhm_med.fits in %s ' % (galaxies[a],workdir))

################################################################################
                
sys.exit('Done.')
