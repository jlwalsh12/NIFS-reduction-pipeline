"""

This routine creates a file that lists the arc exposure to use to
wavelength calibrate each of the telluric star observations, galaxy
observations, or PSF star observations.

"""

def nifs_writearc2use_LP(workdir,rootdir,rawdir,reducedir,objlist,date,
                         obs_setup,objtype):
    
    import os
    import glob

    #get the name of the arc files
    caldir = rootdir+reducedir+'daycals/'+date+'/'+obs_setup+'/'
    os.chdir(caldir)
    arclist = glob.glob('arclist*')
    arc = []
    for ii in range(len(arclist)):
        arc.append(str(open(caldir+arclist[ii], 'r').readlines()[0]).strip())
            
    os.chdir(workdir)

    print('')
    print('***********************************************************')
    print('Provide the following information for exposures in %s: ' % (workdir))

    #ask the user which arc file to use for each objlist.
    for ii in range(len(objlist)):

        print('')
        for jj in range(len(arc)):
            print('Type %i to use %s' % (jj+1, arc[jj]))

        objfile=str(open(objlist[ii],'r').readlines()[0]).strip()

        response = '5000'
        while not response.lower() in str( (range(len(arc)+1))[1:] ):
            response = raw_input('Which arc do you want to use to '+\
                                 'reduce '+objlist[ii]+' (starting '+\
                                 'with file '+objfile+')? It may be '+\
                                 'useful to look at datalist.txt in '+\
                                 rawdir+date+'/. ')
            if response.lower() == '':
                response = '5000'
            if not response.lower() in str( (range(len(arc)+1))[1:] ):
                print('That is not one of the choices. Enter again.')

        if objtype == 'gal' or objtype == 'psf':
            with open('arc2uselist', 'w') as arc2usefile:
                arc2usefile.write('%s\n' % ('wrgn'+arc[int(response)-1]))
        else:
            with open('arc2uselist_'+str(ii+1), 'w') as arc2usefile:
                arc2usefile.write('%s\n' % ('wrgn'+arc[int(response)-1]))

    print('***********************************************************')


