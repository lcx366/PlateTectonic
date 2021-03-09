
from datetime import datetime,timedelta
from os import path,makedirs,remove
from pathlib import Path
from ftplib import FTP
from gzip import GzipFile
from zipfile import ZipFile

from ..utils.try_download import tqdm_ftp,tqdm_request

def download_itrf_snx(dir_to=None,out_days=90):
    '''
    Download or update the space weather file from www.celestrak.com
    Usage: 
    swfile = download_sw([direc])
    Inputs: 
    direc -> [str, optionanl, default = $HOME+'/src/sw-data/'] Directory for storing sw file
    
    Outputs: 
    swfile -> [str] Path of sw file
    '''
    
    if dir_to is None:
        home = str(Path.home())
        dir_to = home + '/src/platemotion-data/'
    
    files = ['ITRF2014-IGS-TRF.SNX.gz','ITRF2014-psd-gnss.snx'] 

    server = 'itrf-ftp.ign.fr'
    ftp = FTP(server,timeout=200) 
    ftp.login()  
    ftp.cwd('~/pub/itrf/itrf2014/') 

    if not path.exists(dir_to): makedirs(dir_to)

    unzipflag = False

    for file in files:
        dir_file = dir_to + file
        dir_file_unzip = dir_file.replace('.gz','')

        if not path.exists(dir_file_unzip):
            if not path.exists(dir_file):
                desc = 'Downloading the ITRF2014 SINEX file {:s} from ITRF'.format(file)
                tqdm_ftp(ftp,dir_to,file,desc)      
            else:
                modified_time = datetime.fromtimestamp(path.getmtime(dir_file))
                if datetime.now() > modified_time + timedelta(days=out_days):
                    remove(dir_file)
                    desc = 'Updating the ITRF2014 SINEX file {:s} from ITRF'.format(file)
                    tqdm_ftp(ftp,dir_to,file,desc)              
        else:
            modified_time = datetime.fromtimestamp(path.getmtime(dir_file_unzip))
            if datetime.now() > modified_time + timedelta(days=out_days):
                remove(dir_file_unzip)
                desc = 'Updating the ITRF2014 SINEX file {:s} from ITRF'.format(file)
                tqdm_ftp(ftp,dir_to,file,desc)
 
        if '.gz' in file and path.exists(dir_file):

            print('Unzip {:s} ... '.format(dir_file),end='')   
            g_file = GzipFile(dir_file)
            open(dir_file[:-3], "wb").write(g_file.read())
            g_file.close()
            print('Finished.')
            remove(dir_file) 
            unzipflag = True

    files = [file.replace('.gz','') for file in files]
    return dir_to,files,unzipflag    

def download_gsrm(dir_to=None):
    '''
    Download or update the space weather file from www.celestrak.com
    Usage: 
    swfile = download_sw([direc])
    Inputs: 
    direc -> [str, optionanl, default = $HOME+'/src/sw-data/'] Directory for storing sw file
    
    Outputs: 
    swfile -> [str] Path of sw file
    '''
    
    if dir_to is None:
        home = str(Path.home())
        dir_to = home + '/src/platemotion-data/'
    
    file = 'GSRM_gridded_strain_v2.1.zip'
    file_unzip = file.replace('.zip','.txt')

    if not path.exists(dir_to): makedirs(dir_to)

    dir_file = dir_to + file
    dir_file_unzip = dir_to + file_unzip


    url = 'http://ftp.globalquakemodel.org/strain-rate/' + file

    if not path.exists(dir_file_unzip):
        if not path.exists(dir_file):
            desc = 'Downloading the Global Strain Rate {:s} from GEM'.format(file)    
            tqdm_request(url,dir_to,file,desc) 
 
    if path.exists(dir_file):
        print('Unzip {:s} ... '.format(dir_file),end='') 
        zipdata = ZipFile(dir_file)
        zipinfo = zipdata.infolist()[0]
        zipinfo.filename = file_unzip
        zipdata.extract(zipinfo,dir_to)
        print('Finished.')
        remove(dir_file) 

    return dir_file_unzip

def download_gia(dir_to=None):
    '''
    Download or update the space weather file from www.celestrak.com
    Usage: 
    swfile = download_sw([direc])
    Inputs: 
    direc -> [str, optionanl, default = $HOME+'/src/sw-data/'] Directory for storing sw file
    
    Outputs: 
    swfile -> [str] Path of sw file
    '''
    
    if dir_to is None:
        home = str(Path.home())
        dir_to = home + '/src/platemotion-data/'
    
    file = 'drad.12mgrid_512.nc'

    if not path.exists(dir_to): makedirs(dir_to)

    dir_file = dir_to + file

    url = 'https://www.atmosp.physics.utoronto.ca/~peltier/datasets/Ice6G_D_VM5a_O512/' + file

    if not path.exists(dir_file):
        desc = 'Downloading the GIA rate of radial displacement {:s}'.format(file)    
        tqdm_request(url,dir_to,file,desc) 

    return dir_file   
