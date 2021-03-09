
from datetime import datetime,timedelta
from os import path,makedirs,system
from pathlib import Path

def split_file(dir_from,files,unzipflag):

    for dir_sub in ['ITRF_ID','ITRF_TS','ITRF_PV']:
        dir_to_sub = dir_from + dir_sub
        if not path.exists(dir_to_sub): makedirs(dir_to_sub)

    igs_trf,psd_gnss = files
    dir_from_igs_trf = dir_from + igs_trf
    dir_from_psd_gnss = dir_from + psd_gnss  

    dir_to_ITRF_ID = dir_from + 'ITRF_ID/' + igs_trf[:-4] + '.ID'
    dir_to_ITRF_TS = dir_from + 'ITRF_TS/' + igs_trf[:-4] + '.TS'
    dir_to_ITRF_PV = dir_from + 'ITRF_PV/' + igs_trf[:-4] + '.PV'
    dir_to_psd_gnss = dir_from_psd_gnss[:-4] + '.site'

    cmd_ITRF_ID = "awk '/\+SITE\/ID/,/\-SITE\/ID/{print}'" + " {:s} > {:s}".format(dir_from_igs_trf,dir_to_ITRF_ID)
    cmd_ITRF_TS = "awk '/\+SOLUTION\/EPOCHS/,/\-SOLUTION\/EPOCHS/{print}'" + " {:s} > {:s}".format(dir_from_igs_trf,dir_to_ITRF_TS)
    cmd_ITRF_PV = "awk '/\+SOLUTION\/ESTIMATE/,/\-SOLUTION\/ESTIMATE/{print}'" + " {:s} > {:s}".format(dir_from_igs_trf,dir_to_ITRF_PV)
    cmd_psd_gnss = "awk '/\+SITE\/ID/,/\-SITE\/ID/{print}'" + " {:s} > {:s}".format(dir_from_psd_gnss,dir_from_psd_gnss[:-4])

    if unzipflag:

        print('Extracting block SITE/ID from {:s} ... '.format(igs_trf),end='')
        system(cmd_ITRF_ID)
        print('Finished.')

        print('Extracting block SOLUTION/EPOCHS from {:s} ... '.format(igs_trf),end='')
        system(cmd_ITRF_TS)
        print('Finished.')

        print('Extracting block SOLUTION/ESTIMATE from {:s} ... '.format(igs_trf),end='')
        system(cmd_ITRF_PV)
        print('Finished.')

        print('Extracting block SITE/ID from {:s} ... '.format(psd_gnss),end='')
        system(cmd_psd_gnss) 
        print('Finished.')

    block_files = [dir_to_ITRF_ID,dir_to_ITRF_TS,dir_to_ITRF_PV,dir_to_psd_gnss]    

    return block_files  