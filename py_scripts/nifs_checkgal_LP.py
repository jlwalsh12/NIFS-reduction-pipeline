"""

This routine has the user check the automatically generated galaxy
file lists and telluric correction list (specifying which telluric
star exposure to use to correct each galaxy exposure). The user has
the ability to modify the lists manually at this point (if
needed). The user is also asked which arc exposure should be used to
wavelength calibrate each galaxy observation.


"""

def nifs_checkgal_LP(rootdir, reducedir, datadir, galaxies, dates):

    import os
    import glob
    from nifs_checklists_LP import nifs_checklists_LP
    from nifs_create_telcorrlist_LP import nifs_create_telcorrlist_LP
    from nifs_writearc2use_LP import nifs_writearc2use_LP

    for i in range(len(galaxies)):
        for j in range(len(dates)):
            
            #dates contains all the days data are to be reduced, but
            #this galaxy may not have been observed on this
            #date. check that here.
            if os.path.isdir(rootdir+reducedir+galaxies[i]+'/'+dates[j]):

                os.chdir(rootdir+reducedir+galaxies[i]+'/'+dates[j])
                obs_setups = sorted(glob.glob('*'), key=os.path.basename)

                for k in range(len(obs_setups)):

                    #change to the galaxy work directory
                    workdir = rootdir+reducedir+galaxies[i]+'/'+dates[j]+\
                      '/'+obs_setups[k]+'/'
                    os.chdir(workdir)

                    #the names of the galaxy and sky lists
                    gallist = 'gallist'
                    skylist = 'skylist'
                    skylistshort = 'skylistshort'

                    #have the user check/modify the galaxy and sky lists
                    listnames = [gallist,skylist,skylistshort]
                    nifs_checklists_LP(workdir,rootdir+datadir,listnames,
                                       'galaxy',dates[j])

                    #create a telluric correction list the same size as the
                    #galaxy list names of the galaxy files
                    nifs_create_telcorrlist_LP(workdir, gallist, rootdir,
                                               reducedir, dates[j],
                                               obs_setups[k])

                    #have the user check the telluric correction list that was
                    #made
                    telcorrlist = 'telcorrlist'
                    listnames = [telcorrlist]
                    nifs_checklists_LP(workdir,rootdir+datadir,listnames,
                                       'telluric_correction',dates[j])

                    #create a file that specifies which arc to use for the
                    #gallist
                    nifs_writearc2use_LP(workdir,rootdir,rootdir+datadir,
                                         reducedir, [gallist], dates[j],
                                         obs_setups[k], 'gal')
                        
