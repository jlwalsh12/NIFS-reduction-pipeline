"""

This routine has the user check the automatically generated PSF star
file lists and telluric correction list (specifying which telluric
star exposure to use to correct each PSF star exposure). The user has
the ability to modify the lists manually at this point (if
needed). The user is also asked which arc exposure should be used to
wavelength calibrate each PSF star observation.

"""

def nifs_checkpsf_LP(rootdir, reducedir, datadir, galaxies, dates):

    import os
    import glob
    from nifs_checklists_LP import nifs_checklists_LP
    from nifs_create_telcorrlist_LP import nifs_create_telcorrlist_LP
    from nifs_writearc2use_LP import nifs_writearc2use_LP

    for i in range(len(dates)):
            
        #dates contains all the days data are to be reduced, but this
        #psf may not have been observed on this date. check that here.
        if os.path.isdir(rootdir+reducedir+'psfs/'+dates[i]):

            os.chdir(rootdir+reducedir+'psfs/'+dates[i])
            obs_setups = sorted(glob.glob('*'), key=os.path.basename)

            for j in range(len(obs_setups)):
                for k in range(len(galaxies)):

                    #a psf may or may not exist for each galaxy, so
                    #check.
                    if os.path.isdir(rootdir+reducedir+'psfs/'+dates[i]+'/'+\
                                     obs_setups[j]+'/psf_'+galaxies[k]):

                        #change to the psf work directory
                        workdir = rootdir+reducedir+'psfs/'+dates[i]+'/'+\
                          obs_setups[j]+'/psf_'+galaxies[k]+'/'
                        os.chdir(workdir)

                        #the names of the psf exposures and sky lists
                        psflist = 'psflist'
                        skylist = 'skylist'
                        skylistshort = 'skylistshort'

                        #have the user check/modify the psf and sky
                        #lists
                        listnames = [psflist,skylist,skylistshort]
                        nifs_checklists_LP(workdir,rootdir+datadir,
                                           listnames,'psf',dates[i])

                        #create a telluric correction list the same
                        #size as the psf list
                        nifs_create_telcorrlist_LP(workdir, psflist, rootdir,
                                                   reducedir, dates[i],
                                                   obs_setups[j])

                        #have the user check the telluric correction
                        #list that was made
                        telcorrlist = 'telcorrlist'
                        listnames = [telcorrlist]
                        nifs_checklists_LP(workdir,rootdir+datadir,listnames,
                                           'telluric_correction',dates[i])

                        #create a file that specifies which arc to use
                        #for the psflist
                        nifs_writearc2use_LP(workdir,rootdir,rootdir+datadir,
                                             reducedir, [psflist], dates[i],
                                             obs_setups[j],'psf')
