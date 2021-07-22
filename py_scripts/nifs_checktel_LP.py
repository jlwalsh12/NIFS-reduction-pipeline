"""

This routine has the user check the telluric file lists that were
automatically generated. The user has the ability to modify the lists
manually at this point (if needed). The user is also asked which arc
exposure should be used to wavelength calibrate the telluric star
observation.

"""

def nifs_checktel_LP(rootdir, reducedir, datadir, dates):

    import os
    import glob
    from nifs_checklists_LP import nifs_checklists_LP
    from nifs_writearc2use_LP import nifs_writearc2use_LP

    for i in range(len(dates)):

        os.chdir(rootdir+reducedir+'tellurics/'+dates[i])
        obs_setups = sorted(glob.glob('*'), key=os.path.basename)

        for j in range(len(obs_setups)):

            os.chdir(rootdir+reducedir+'tellurics/'+dates[i]+'/'+\
                     obs_setups[j]+'/')
            telluric_stars = sorted(glob.glob('*'), key=os.path.basename)

            for k in range(len(telluric_stars)):

                #change to the telluric work directory
                workdir = rootdir+reducedir+'tellurics/'+dates[i]+'/'+\
                  obs_setups[j]+'/'+telluric_stars[k]+'/'
                os.chdir(workdir)

                #the names of the telluric and sky lists
                telluriclist = sorted(glob.glob('telluriclist*'),
                                      key=os.path.basename)
                skylist = sorted(glob.glob('skylist_*'),
                                      key=os.path.basename)
                skylistshort = sorted(glob.glob('skylistshort_*'),
                                      key=os.path.basename)

                #have the user check/modify the telluric and sky lists
                listnames = [telluriclist,skylist,skylistshort]
                nifs_checklists_LP(workdir,rootdir+datadir,listnames,
                                   'telluric',dates[i])

                #create a file that specifies which arc to use for
                #each telluriclist
                nifs_writearc2use_LP(workdir,rootdir,rootdir+datadir,reducedir,
                                     telluriclist, dates[i],obs_setups[j],
                                     'tel')
