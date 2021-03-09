import numpy as np
from scipy.interpolate import RectSphereBivariateSpline
from astropy import units as u
from astropy.coordinates import Angle
import xarray as xr
from os import path,walk,makedirs
from pathlib import Path
from pkg_resources import resource_filename
from sphericalpolygon import Sphericalpolygon
import pandas as pd

from ..utils.utils import year_day_second,geode2geocen
from ..utils.data_prepare import data_prepare
from ..classes.plate import Plate,PlateMotion

def long_than_3yr(ITRF_TS,ITRF_PV):

    # read the start epoch and end epoch of sites
    data_ts = np.genfromtxt(ITRF_TS,dtype=np.str,skip_header=2,skip_footer=1)

    # reject sites with the time-span of observations shorter than 3 years
    sites_of_reject = []
    n = len(data_ts)
    for i in range(n):
        time_span = year_day_second(data_ts[i,5]) - year_day_second(data_ts[i,4])
        if time_span < 3:
            sites_of_reject.append(data_ts[i,0]+data_ts[i,1]+data_ts[i,2])

    # read the positions and velocities
    data_pv = np.genfromtxt(ITRF_PV,dtype=np.str,skip_header=2,skip_footer=1)

    # build a pandas dataframe
    n = len(data_pv)
    CODE,PT,SOLN,REF_EPOCH = [[] for j in range(4)]
    STAX,STAY,STAZ,VELX,VELY,VELZ = [[] for j in range(6)]
    STAX_STD,STAY_STD,STAZ_STD,VELX_STD,VELY_STD,VELZ_STD = [[] for j in range(6)]
    for i in range(0,n,6):
        CODE.append(data_pv[i][2])
        PT.append(data_pv[i][3])
        SOLN.append(data_pv[i][4])
        REF_EPOCH.append(data_pv[i][5])
        STAX.append(data_pv[i][8])
        STAY.append(data_pv[i+1][8])
        STAZ.append(data_pv[i+2][8])
        VELX.append(data_pv[i+3][8])
        VELY.append(data_pv[i+4][8])
        VELZ.append(data_pv[i+5][8])
        STAX_STD.append(data_pv[i][9])
        STAY_STD.append(data_pv[i+1][9])
        STAZ_STD.append(data_pv[i+2][9])
        VELX_STD.append(data_pv[i+3][9])
        VELY_STD.append(data_pv[i+4][9])
        VELZ_STD.append(data_pv[i+5][9])
    df_array = np.stack((CODE,PT,SOLN,REF_EPOCH,STAX,STAY,STAZ,VELX,VELY,VELZ,STAX_STD,STAY_STD,STAZ_STD,VELX_STD,VELY_STD,VELZ_STD)).T
    df_column = ['CODE','PT','SOLN','REF_EPOCH','STAX','STAY','STAZ','VELX','VELY','VELZ','STAX_STD','STAY_STD','STAZ_STD','VELX_STD','VELY_STD','VELZ_STD']
    df = pd.DataFrame(df_array,columns=df_column)      

    CODE_PT_SOLN = df['CODE'] + df['PT'] + df['SOLN']
    _,drop_index,__ = np.intersect1d(CODE_PT_SOLN, sites_of_reject, return_indices=True)  
    sites_of_retain0_pv = df.drop(drop_index).reset_index(drop=True)

    return sites_of_retain0_pv   

def less_than_1mm(sites_of_retain0_pv,vel_std_max=1):

    # uncertainty less than 1mm/yr 

    vx_vy_vz_std = sites_of_retain0_pv.loc[:,['VELX_STD','VELY_STD','VELZ_STD']]
    flags = np.max(vx_vy_vz_std,axis=1)*1e3 < vel_std_max    
    sites_of_retain1_pv = sites_of_retain0_pv[flags]

    return sites_of_retain1_pv

def remove_psd(sites_of_retain1_pv,ITRF_psd_gnss):

    # read the positions and velocities
    width = (5,3,10,2,23,4,3,5,4,3,5,8)
    data_psd = np.genfromtxt(ITRF_psd_gnss,delimiter = width,dtype=np.str,skip_header=2,skip_footer=1,autostrip=True)
    sites_psd = np.core.defchararray.add(data_psd[:,0], data_psd[:,1])

    df = sites_of_retain1_pv
    CODE_PT = df['CODE'] + df['PT']
    n = len(CODE_PT)
    drop_index = []
    for i in range(n):
        if CODE_PT[i] in sites_psd: drop_index.append(i)   
    sites_of_retain2_pv = df.drop(drop_index).reset_index(drop=True)
    # remove duplicates
    sites_of_retain3_pv = sites_of_retain2_pv.drop_duplicates(subset='CODE',keep='last',ignore_index=True)

    return sites_of_retain3_pv 

def apply_gsrm(sites_of_retain3_pv,ITRF_ID,dir_file_gsrm):

    # read the lons, lats, and heights 
    width = (5,3,10,25,4,3,5,4,3,5,8)
    data_lonlat = np.genfromtxt(ITRF_ID,delimiter = width,dtype=np.str,skip_header=2,skip_footer=1,autostrip=True)
    n = len(data_lonlat)
    code_pt = np.core.defchararray.add(data_lonlat[:,0],data_lonlat[:,1])

    CODE_PT = sites_of_retain3_pv['CODE'] + sites_of_retain3_pv['PT']
    _,__,_index = np.intersect1d(CODE_PT, code_pt, return_indices=True)  
    sites_of_retain3_lonlat = data_lonlat[_index]

    # convert lons and lats of sites in [d,m,s] to lons and colats in radian 
    sites_of_retain3_lons = np.array(sites_of_retain3_lonlat[:,4:7],dtype = np.float)
    sites_of_retain3_geod_lats = np.array(sites_of_retain3_lonlat[:,7:10],dtype = np.float)

    sites_of_retain3_lons = [tuple(x) for x in sites_of_retain3_lons]
    sites_of_retain3_geod_lats = [tuple(x) for x in sites_of_retain3_geod_lats]

    sites_of_retain3_lons_deg = Angle(sites_of_retain3_lons, unit=u.deg)
    sites_of_retain3_lats_geod_deg = Angle(sites_of_retain3_geod_lats, unit=u.deg)
    sites_of_retain3_lats_deg = geode2geocen(sites_of_retain3_lats_geod_deg)

    sites_of_retain3_lons = sites_of_retain3_lons_deg.rad
    sites_of_retain3_lats = sites_of_retain3_lats_deg.to(u.rad).value

    df = pd.DataFrame(sites_of_retain3_lonlat[:,:4], columns = ['CODE','PT','DOMES','STATION DESCRIPTION'])
    df['Geod. Lat. °N'] = sites_of_retain3_lats_geod_deg
    df['Lat. °N'] = sites_of_retain3_lats_deg
    df['Lon. °E'] = sites_of_retain3_lons_deg
    df['H'] = sites_of_retain3_lonlat[:,-1]
    sites_of_retain3_lonlat = df

    sites_of_retain3_colats = np.pi/2 - sites_of_retain3_lats

    # read the GSRMv2.1 
    data_gsrmv21 = np.loadtxt(dir_file_gsrm)

    # calculate second invariant strain rate
    sisr = np.sqrt(data_gsrmv21[:,8]**2+data_gsrmv21[:,9]**2).reshape(1750,3600)

    # flip the lons and lats
    sisr = np.flip(sisr,0)
    sisr[:,:1800],sisr[:,1800:] = sisr[:,1800:].copy(), sisr[:,:1800].copy()

    # lons and lats after flipping
    lons = np.linspace(0.05,359.95,3600)
    lats = np.linspace(87.45,-87.45,1750)

    # colats and lons in radians
    colats_rad = np.deg2rad(90 - lats)
    lons_rad = np.deg2rad(lons)

    lut = RectSphereBivariateSpline(colats_rad, lons_rad, sisr)

    sisr_sites = lut.ev(sites_of_retain3_colats,sites_of_retain3_lons) 

    sites_of_retain4_lonlat = sites_of_retain3_lonlat[sisr_sites < 0.1].reset_index(drop=True)

    CODE = sites_of_retain4_lonlat['CODE']
    CODE_pv = sites_of_retain3_pv['CODE']

    _,_index,__ = np.intersect1d(CODE_pv, CODE, return_indices=True)  
    sites_of_retain4_pv = sites_of_retain3_pv.loc[_index].reset_index(drop=True)
    sites_of_retain4 = pd.merge(sites_of_retain4_lonlat,sites_of_retain4_pv, on=['CODE','PT'],validate="one_to_one")

    return sites_of_retain4

def apply_gia(sites_of_retain4,dir_file_gia):
    
    # read the rate of radial displacement (UP) on a 0.2x0.2 grid from ICE-6G_D(VM5a)
    data_ICE_6G_D = dir_file_gia

    # Open the dataset and print out metadeta
    ds = xr.open_dataset(data_ICE_6G_D)
    
    # convert lons and lats of sites in [d,m,s] to lons and colats in radian 

    sites_of_retain4_lons = sites_of_retain4['Lon. °E']*u.deg.to(u.rad)
    sites_of_retain4_lats = sites_of_retain4['Lat. °N']*u.deg.to(u.rad)
    sites_of_retain4_colats = np.pi/2 - sites_of_retain4_lats

    # colats and lons in radians
    lats,lons,up_rate = ds['Lat'],ds['Lon'],ds['Drad_250']

    colats_rad = np.deg2rad(90 - lats)
    lons_rad = np.deg2rad(lons)

    lut = RectSphereBivariateSpline(colats_rad, lons_rad, up_rate)

    up_rate_sites = lut.ev(sites_of_retain4_colats,sites_of_retain4_lons) 

    sites_of_retain5 = sites_of_retain4[up_rate_sites < 0.75]

    return sites_of_retain5

def assign_sites(sites_of_retain5,modelname):

    '''
    avaliable modelname: NNR-MORVEL56 and GSRMv2.1
    '''

    data_path = resource_filename('platemotion', 'data/')

    if modelname == 'NNR-MORVEL56':      
        bnds_dir = 'NnrMRVL_PltBndsLatLon/'
    elif modelname == 'GSRMv2.1':   
        bnds_dir = 'NnrGSRMv2.1_PltBndsLatLon/'
    else:
        raise Exception("Currently, only 'NNR-MORVEL56' and 'GSRMv2.1' are avaliable.")    

    boundary_files = []
    # Go through all boundary files
    for (dirname, dirs, files) in walk(data_path + bnds_dir):
        for filename in files:
            if len(filename) == 2: boundary_files.append(path.join(dirname,filename))        

    sites_of_retain5_lons = sites_of_retain5['Lon. °E']
    sites_of_retain5_lats = sites_of_retain5['Lat. °N']

    # lats and lons of sites
    points = np.vstack((sites_of_retain5_lats, sites_of_retain5_lons)).T

    plates = []
    for boundary_file in boundary_files:
        polygon = Sphericalpolygon.from_file(boundary_file,skiprows=1) 
        # identify whether the given points are in the region
        flag = polygon.contains_points(points)
        n = np.sum(flag)
        if n >= 3: # At least three sites are in a given area
            plate = Plate.from_file(boundary_file,skiprows=1) 
            plate.set_name(boundary_file[-2:])
            plate.set_sites(sites_of_retain5[flag].reset_index(drop=True))
            plates.append(plate)
    platemodel = PlateMotion.initialize(plates)        
    return platemodel        

def sites_plate(modelname,dir_to=None,out_days=90):
    block_files,dir_file_gsrm,dir_file_gia = data_prepare(dir_to,out_days)
    ITRF_ID,ITRF_TS,ITRF_PV,ITRF_psd_gnss = block_files

    print('Screen out sites with an observation time of more than 3 years ... ',end='')
    sites_of_retain0_pv = long_than_3yr(ITRF_TS,ITRF_PV)
    print('Finished')

    print('Eliminate sites with velocity uncertainty greater than 1mm/yr ... ',end='')
    sites_of_retain1_pv = less_than_1mm(sites_of_retain0_pv)
    print('Finished')

    print('Eliminate sites affected by the post-earthquake deformation ... ',end='')
    sites_of_retain3_pv = remove_psd(sites_of_retain1_pv,ITRF_psd_gnss)
    print('Finished')

    print('Exclude sites with the second invariant strain rate greater than 1e-10 based on gsrmv2.1 ... ',end='')
    sites_of_retain4 = apply_gsrm(sites_of_retain3_pv,ITRF_ID,dir_file_gsrm)
    print('Finished')

    print('Exclude sites with the vertical displacement rate greater than 0.75mm/yr based on ICE-6G_D(VM5a) ... ',end='')
    sites_of_retain5 = apply_gia(sites_of_retain4,dir_file_gia)
    print('Finished')

    print('Assign sites to plates ... ',end='')
    platemodel = assign_sites(sites_of_retain5,modelname)
    print('Finished')

    return platemodel   