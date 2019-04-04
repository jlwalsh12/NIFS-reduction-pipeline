"""

This routine checks whether there is already a merged cube in the
output directory. If there is, the user is asked whether the product
should be deleted and all individual cubes should be merged together
from scratch. If the user does not want this, the pipeline will exit.

"""

def nifs_checkdata_merged_LP(workdir,filecheck):
    
    import sys
    import os
    import glob
    import subprocess

    if len(glob.glob(workdir+filecheck)) > 0:

        print('')
        print('There are outputs from a previous reduction already.')

        response = ''
        while not response == 'yes' or response == 'no':
            
            response = raw_input('Do you want to delete these files and '+\
                                 'combine the cubes from scratch (yes/no)? ')

            if response[0:1] == '"' or response[0:1] == "'":
                response = response[1:-1]

            if response == 'yes':

                tmp = subprocess.call(['mkdir',workdir+'holdtmp/'],
                                      stderr=open(os.devnull,'w'))
                tmp = subprocess.call(['mv']+\
                                      glob.glob(os.path.join(workdir,
                                                             'atfbrsn*.fits'))+\
                                      [workdir+'holdtmp/'],
                                      stderr=open(os.devnull,'w'))
                tmp = subprocess.call(['mv']+\
                                      glob.glob(os.path.join(workdir,
                                                    'image_catfbrsn*.fits'))+\
                                      [workdir+'holdtmp/'],
                                      stderr=open(os.devnull,'w'))
                tmp = subprocess.call(['rm']+\
                                      glob.glob(os.path.join(workdir,'*')),
                                      stderr=open(os.devnull,'w'))
                
                if os.path.isdir(workdir+'mask/'):
                    tmp = subprocess.call(['rm','-r',workdir+'mask/'],
                                          stderr=open(os.devnull,'w'))
                    
                tmp = subprocess.call(['mv']+\
                                      glob.glob(os.path.join(workdir+
                                                             'holdtmp/','*'))+\
                                      [workdir], stderr=open(os.devnull,'w'))
                tmp = subprocess.call(['rm','-r',workdir+'holdtmp/'],
                                      stderr=open(os.devnull,'w'))

                os.chdir(workdir)
                filenames = glob.glob('atfbrsn*.fits')
                for i in range(len(filenames)):
                    tmp = subprocess.call(['cp',filenames[i],'h'+filenames[i]],
                                          stderr=open(os.devnull,'w'))

            if response == 'no':
                sys.exit('Okay, NOT deleting any files. Modify inputs to'+\
                         ' the pipeline and run again.')

            if response != 'yes' and response != 'no':
                print('That is not one of the choices. Enter again.')
