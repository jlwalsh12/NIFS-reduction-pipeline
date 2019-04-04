"""

Previously, the spatial offsets between galaxy data cubes (relative to
the first cube) have been determined based on the location of the
spaxel with the largest flux value, the location of the center of a
fitted 2D gaussian, and by cross-correlating the collapsed data
cubes. This routine adopts cross-correlation as the default method,
but checks that the offsets determined by the other two methods are
similar. If the offsets are not similar, then the user is asked to
look at the offsets determined by the other methods, and select which
method to adopt instead. Alternatively, the user can provide another
file with offsets to use instead.

"""

def nifs_getoffmethod_LP(workdir, nfiles):

    import numpy as np
    import glob
    import os
    
    if nfiles > 1:

        #read in each offset method file.
        offsets_max = np.loadtxt('offset_max.list', dtype='float',
                                 comments='#')
        offsets_gauss = np.loadtxt('offset_gauss.list', dtype='float',
                                   comments='#')
        offsets_crosscorr = np.loadtxt('offset_crosscorr.list',
                                       dtype='float', comments='#')

        #default will be offsets determined via
        #cross-correlation. check that the offsets are at least within
        #2.5 pixels away from offsets determined using maximum or
        #gaussian methods.
        oktouse = 'yes'
        for i in range(len(offsets_crosscorr)):

            xoff_crosscorr = offsets_crosscorr[i,0]
            yoff_crosscorr = offsets_crosscorr[i,1]

            xoff_max = offsets_max[i,0]
            yoff_max = offsets_max[i,1]

            xoff_gauss = offsets_gauss[i,0]
            yoff_gauss = offsets_gauss[i,1]

            if abs(xoff_crosscorr - xoff_max) > 2.5 and \
              abs(xoff_crosscorr - xoff_gauss) > 2.5:
                oktouse = 'no'

            if abs(yoff_crosscorr - yoff_max) > 2.5 and \
              abs(yoff_crosscorr - yoff_gauss) > 2.5:
                oktouse = 'no'

        if oktouse == 'yes':
            offsetfile = 'offset_crosscorr.list'
            headerinfo = 'headervals_offset_crosscorr.dat'
            print('Using cross-correlation method to determine offsets'+\
                  ' between cubes.')

        #if the cross-correlation numbers are strange, ask the user to
        #look at the offsets from the other methods and select a
        #different method to adopt, or to provide the name of an
        #offset file they would like to use instead.
        if oktouse == 'no':

            print('')
            print('Offsets from cross-correlation method seem off. '+\
                  'Can try using maximum or gaussian methods to determine'+\
                  ' offsets.')
            print('Look at offset_*.list in %s.' % (workdir))
            print('Type 1 to use offset_max.list.')
            print('Type 2 to use offset_gauss.list')
            print('Type 3 to provide the name of the offset and '+\
                  'header information files that you have created yourself.')

            response = ''
            while not (response == '1' or response == '2' or response == '3'):
            
                response = raw_input('Which offset file do you want to use? ')
            
                if response[0:1] == '"' or response[0:1] == "'":
                    response = response[1:-1]
                
                if response.lower() == '1':
                    offsetfile = 'offset_max.list'
                    headerinfo = 'headervals_offset_max.dat'
                    
                if response.lower() == '2':
                    offsetfile = 'offset_gauss.list'
                    headerinfo = 'headervals_offset_gauss.dat'
                    
                if response.lower() == '3':
                    userfiles_good = ''
                    while not userfiles_good == 'yes':
                        offsetfile = raw_input('What is the name of the '+\
                                               'offset file you want to use '+\
                                               '(file must be located in '+\
                                               workdir+')? ')
                        headerinfo = raw_input('What is the name of the '+\
                                               'header information file you'+\
                                               ' want to use (file must be'+\
                                               ' located in '+ workdir +')? ')
                        if offsetfile[0:1] == '"' or offsetfile[0:1] == "'":
                            offsetfile = offsetfile[1:-1]
                        if headerinfo[0:1] == '"' or headerinfo[0:1] == "'":
                            headerinfo = headerinfo[1:-1]
                        #check that the above files exist
                        os.chdir(workdir)
                        if len(glob.glob(offsetfile)) == 1 and \
                          len(glob.glob(headerinfo)) == 1:
                            userfiles_good = 'yes'
                        else:
                            print('File(s) do not exist. Enter again.')
                        
                if (response.lower() != '1' and response.lower() != '2' and \
                        response.lower() != '3'):
                    print('That is not one of the choices. Enter again.')

    if nfiles == 1:

        #read in each offset method file.
        offsets_max = np.loadtxt('offset_max.list', dtype='float',
                                 comments='#')
        offsets_gauss = np.loadtxt('offset_gauss.list', dtype='float',
                                   comments='#')

        #default will be center determined via gaussian. check that
        #the center is at least within 2.5 pixels away from center
        #determined using maximum method.
        oktouse = 'yes'

        xoff_max = offsets_max[0]
        yoff_max = offsets_max[1]

        xoff_gauss = offsets_gauss[0]
        yoff_gauss = offsets_gauss[1]

        if abs(xoff_gauss - xoff_max) > 2.5:
            oktouse = 'no'

        if abs(yoff_gauss - yoff_max) > 2.5:
            oktouse = 'no'

        if oktouse == 'yes':
            offsetfile = 'offset_gauss.list'
            headerinfo = 'headervals_offset_gauss.dat'
            print('Using Gaussian method to determine center of cube.')

        #if the gaussian center is strange, ask the user if they would
        #rather use the maximum value, or enter the name of an offset
        #file they would like to use instead.
        if oktouse == 'no':

            print('')
            print('Center from gaussian method seem off.')
            print('Type 1 to use offset_max.list in %s instead.' % (workdir))
            print('Type 2 to provide the name of the offset and '+\
                  'header information files that you have created yourself.')

            response = ''
            while not (response == '1' or response == '2'):
            
                response = raw_input('What do you want to do? ')
            
                if response[0:1] == '"' or response[0:1] == "'":
                    response = response[1:-1]
                
                if response.lower() == '1':
                    offsetfile = 'offset_max.list'
                    headerinfo = 'headervals_offset_max.dat'
                    
                if response.lower() == '2':
                    userfiles_good = ''
                    while not userfiles_good == 'yes':
                        offsetfile = raw_input('What is the name of the '+\
                                               'offset file you want to use '+\
                                               '(file must be located in '+\
                                               workdir+')? ')
                        headerinfo = raw_input('What is the name of the '+\
                                               'header information file '+\
                                               'you want to use (file must'+\
                                               ' be located in '+\
                                               workdir+')? ')
                        if offsetfile[0:1] == '"' or offsetfile[0:1] == "'":
                            offsetfile = offsetfile[1:-1]
                        if headerinfo[0:1] == '"' or headerinfo[0:1] == "'":
                            headerinfo = headerinfo[1:-1]
                        #check that the above files exist
                        os.chdir(workdir)
                        if len(glob.glob(offsetfile)) == 1 and \
                          len(glob.glob(headerinfo)) == 1:
                            userfiles_good = 'yes'
                        else:
                            print('File(s) do not exist. Enter again.')

                if (response.lower() != '1' and response.lower() != '2'):
                    print('That is not one of the choices. Enter again.')

    return [offsetfile, headerinfo]
