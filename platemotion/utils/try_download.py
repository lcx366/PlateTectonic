from os import remove
import requests
from tqdm import tqdm
from colorama import Fore
from time import sleep

def tqdm_ftp(ftp,dir_to,file,desc):
    '''
    Try to download files from remote server by ftp with a colored progress bar.

    uasge:
    tqdm_ftp(ftp,dir_to,file,desc)

    Inputs:
    ftp
    dir_to
    file
    desc

    Outputs:
    '''
    dir_file = dir_to + file
    total_size = int(ftp.size(file))

    bar_format = "{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET)

    for idownload in range(5):
        try:
            local_file =  open(dir_file, 'wb')
            pos = local_file.tell()
            pbar = tqdm(desc = desc,total=total_size,unit='B',unit_scale=True,position=0,initial=pos,bar_format=bar_format)
            def progressbar(chunk):
                local_file.write(chunk)
                pbar.update(len(chunk))
            res = ftp.retrbinary('RETR ' + file, progressbar) # RETR is an FTP command
            break       
        except:
            sleep(2)
            if idownload == 4:
                remove(dir_file)
                print('No response, skip {:s}'.format(file))
        finally:
            pbar.close()
            local_file.close()   


def tqdm_request(url,dir_to,file,desc):
    '''
    Try to download files from remote server by request with a colored progress bar.
    '''
    dir_file = dir_to + file
    block_size = 1024*10
    bar_format = "{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET)

    for idownload in range(5):
        try:
            local_file = open(dir_file, 'ab')
            pos = local_file.tell()
            resume_header = {'Range': f'bytes={pos}-'}
            res = requests.get(url,stream=True,timeout=200,headers=resume_header)
            total_size = int(res.headers.get('Content-Length',0))
            pbar = tqdm(desc = desc,total=total_size,unit='B',unit_scale=True,bar_format = bar_format,position=0,initial=pos)
            for chunk in res.iter_content(block_size):
                local_file.write(chunk)  
                pbar.update(len(chunk))   
            pbar.close()   
            res.close()  
            break
        except: 
            sleep(2)
            if idownload == 4:
                remove(dir_file)
                print('No response, skip {:s}'.format(file)) 

        finally:    
            local_file.close()      