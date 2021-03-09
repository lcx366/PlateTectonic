
from .data_download import download_itrf_snx,download_gsrm,download_gia
from .seperate_file import split_file

def data_prepare(dir_to=None,out_days=90):
    dir_to_itrf,files,unzipflag = download_itrf_snx(dir_to,out_days)
    block_files = split_file(dir_to_itrf,files,unzipflag)

    dir_file_gsrm = download_gsrm(dir_to)
    dir_file_gia = download_gia(dir_to)

    return block_files,dir_file_gsrm,dir_file_gia