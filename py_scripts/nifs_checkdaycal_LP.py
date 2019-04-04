"""

This routine has the user check the day calibration file lists that
were automatically generated. The user has the ability to modify the
lists manually at this point (if needed).

"""

def nifs_checkdaycal_LP(rootdir, reducedir, datadir, dates):

    import os
    import glob
    from nifs_checklists_LP import nifs_checklists_LP

    for i in range(len(dates)):
    
        os.chdir(rootdir+reducedir+'daycals/'+dates[i])
        obs_setups = glob.glob('*')

        for j in range(len(obs_setups)):
        
            #change the to daycals work directory
            workdir = rootdir+reducedir+'daycals/'+dates[i]+'/'+\
              obs_setups[j]+'/'
            os.chdir(workdir)
        
            #the names of the daycal lists
            flatlist = 'flatlist'
            flatdarklist = 'flatdarklist'
            arcdarklist = 'arcdarklist'
            ronchilist = 'ronchilist'
            ronchidarklist = 'ronchidarklist'
            arclist = glob.glob('arclist*')
            darklist = glob.glob('darklist*')

            #have the user check/modify the calibration lists
            if len(darklist) > 0:
                listnames = [flatlist,flatdarklist,arcdarklist,ronchilist,
                             ronchidarklist,arclist,darklist]
            else:
                listnames = [flatlist,flatdarklist,arcdarklist,ronchilist,
                             ronchidarklist,arclist]
                    
            nifs_checklists_LP(workdir,rootdir+datadir,listnames,'daycal',
                               dates[i])
