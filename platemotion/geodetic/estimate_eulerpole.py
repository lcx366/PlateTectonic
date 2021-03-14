import numpy as np
from numpy.linalg import inv,norm
from astropy import units as u
from astropy.coordinates import cartesian_to_spherical,EarthLocation

from ..utils import Const

# Rotation matrix from cartesian earth-fixed ref. frame to east and morth topocentric ref. frame
def R_M(lon,lat):
    R = np.empty((2,3))
    R[0] = -np.sin(lon),np.cos(lon),0
    R[1] = -np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)
    #R[2] = np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)
    return R

# Design matrix of plate motion equations
def A_M(lon,lat):
    Re = Const.r0*1e3 # volumetric radius of the Earth in meters
    A = np.empty((2,3))
    A[0] = -np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)                                                                                                    
    A[1] = np.sin(lon),-np.cos(lon),0
    return A*Re

# Jacobian matrix of euler vector
def J_M(w,w_lon,w_lat):
    J = np.empty((3,3))
    J[0] = np.cos(w_lon)*np.cos(w_lat),-w*np.sin(w_lon)*np.cos(w_lat),-w*np.cos(w_lon)*np.sin(w_lat)
    J[1] = np.sin(w_lon)*np.cos(w_lat),w*np.cos(w_lon)*np.cos(w_lat),-w*np.sin(w_lon)*np.sin(w_lat)
    J[2] = np.sin(w_lat),0,w*np.cos(w_lat)
    return J

# convert velocity in cartesian to velocity in east and morth
def v_xyz2v_en(v_xyz,lon,lat): 
    v_en = R_M(lon,lat)@v_xyz
    return v_en

# convert velocity std in cartesian to velocity std in east and morth
def v_xyz_std2v_en_std(v_xyz_std,lon,lat): 
    cov_v_xyz = np.diag(v_xyz_std**2)
    RM = R_M(lon,lat)
    RM_T = RM.T
    cov_v_en = RM@cov_v_xyz@RM_T
    v_en_std = np.sqrt(np.diag(cov_v_en))
    return v_en_std

# convert euler vector and its std in cartesian to euler vector and its std in norm, lon and lat
def w_xyz_std2w_std(w_xyz,w_xyz_std):
    w,w_lat,w_lon = cartesian_to_spherical(w_xyz[0],w_xyz[1],w_xyz[2])
    inv_J = inv(J_M(w,w_lon,w_lat))
    cov_w_xyz = np.diag(w_xyz_std**2)
    cov_w_lonlat = inv_J@cov_w_xyz@inv_J.T
    w_std,w_lon_std,w_lat_std = np.sqrt(np.diag(cov_w_lonlat))
    return w,w_lon,w_lat,w_std,w_lon_std,w_lat_std

# Plate Motion Calculator
def PMC(sites_info):
    n = len(sites_info)

    v_xyz = np.array(sites_info.loc[:,['VELX','VELY','VELZ']],dtype=float)
    v_xyz_std = np.array(sites_info.loc[:,['VELX_STD','VELY_STD','VELZ_STD']],dtype=float)
    v_en,weight_M = np.empty((n,2)),np.empty((n,2,2))
    
    lats = np.array(sites_info['Lat. °N'])*u.deg
    lons = np.array(sites_info['Lon. °E'])*u.deg

    ATA,ATB = np.zeros((3,3)),np.zeros(3)

    # calculate the the east-north components of velocity
    for i in range(n):
        v_en[i] = v_xyz2v_en(v_xyz[i],lons[i],lats[i])
        v_en_std_i = v_xyz_std2v_en_std(v_xyz_std[i],lons[i],lats[i])
        B = v_en[i]
        cov_v_en_i = np.diag(v_en_std_i**2)
        A = A_M(lons[i],lats[i])
        A_T = A.T
        weight_M[i] = inv(cov_v_en_i)
        ATA += A_T@weight_M[i]@A
        ATB += A_T@weight_M[i]@B

    w_xyz = inv(ATA)@ATB # rad/yr
    w_xyz_std = np.sqrt(np.diag(inv(ATA)))
    w,w_lon,w_lat,w_std,w_lon_std,w_lat_std = w_xyz_std2w_std(w_xyz,w_xyz_std) 
    
    w_lon,w_lat = w_lon.to(u.deg),w_lat.to(u.deg)

    return w_lat,w_lon,w,w_lat_std,w_lon_std,w_std,w_xyz,w_xyz_std,v_en,weight_M

# Calculate residual
def res(v_en,A,w_xyz):
    return v_en-A@w_xyz

def sigma_1(v_en_res,weight_M):
    n = len(v_en_res) 
    v_en_res2,weight_M_tr = 0,0
    for j in range(n):
        v_en_res2 += v_en_res[j]@weight_M[j]@v_en_res[j]
        weight_M_tr += np.trace(weight_M[j])    
    sigma = np.sqrt(v_en_res2/weight_M_tr*n/(2*n-3))*2 
    return sigma  
  
def sigma_2(v_en_res):
    n = len(v_en_res) 
    sigma = np.sqrt(2*np.sum(v_en_res**2)/(2*n-3))   
    return sigma  

def w1_w2(w1_xyz,w1_xyz_std,w2_xyz,w2_xyz_std):
    w_xyz = w2_xyz - w1_xyz
    w_xyz_std = np.sqrt(w1_xyz_std**2+w2_xyz_std**2)
    w,w_lon,w_lat,w_std,w_lon_std,w_lat_std = w_xyz_std2w_std(w_xyz,w_xyz_std)
    
    w_xyz *= u.deg/u.Ma
    w_xyz_std *= u.deg/u.Ma
    w *= u.deg/u.Ma
    w_std *= u.deg/u.Ma
    w_lon,w_lat = w_lon.to(u.deg),w_lat.to(u.deg)
    w_lon_std,w_lat_std = (w_lon_std*u.rad).to(u.deg),(w_lat_std*u.rad).to(u.deg)
    
    omega_cartesian = w_xyz
    omega_cartesian_std = w_xyz_std
    omega_spherical = w_lat,w_lon,w
    omega_spherical_std = w_lat_std,w_lon_std,w_std
    
    return omega_cartesian,omega_cartesian_std,omega_spherical,omega_spherical_std      

def PMC_iterate(sites_info):
    n = len(sites_info)
    D_M = np.empty((n,2,3))

    lats = np.array(sites_info['Lat. °N'])*u.deg
    lons = np.array(sites_info['Lon. °E'])*u.deg

    # calculate the Design Matrix
    for i in range(n): D_M[i] = A_M(lons[i],lats[i])
    
    # First estimation
    w_lat,w_lon,w,w_lat_std,w_lon_std,w_std,w_xyz,w_xyz_std,v_en,weight_M = PMC(sites_info)
    v_en_res = res(v_en,D_M,w_xyz)
    sigma = sigma_1(v_en_res,weight_M)
    less_than_3s = norm(v_en_res,axis=1) < 3*sigma 
    data_i = sites_info[less_than_3s]

    # Second estimation
    w_lat_i,w_lon_i,w_i,w_lat_std_i,w_lon_std_i,w_std_i,w_xyz_i,w_xyz_std_i,v_en_i,weight_M_i = PMC(data_i)
    v_en_res_i = res(v_en[less_than_3s],D_M[less_than_3s],w_xyz_i)
    sigma_i = sigma_1(v_en_res_i,weight_M_i)
    v_en_res_ii = res(v_en,D_M,w_xyz_i)
    less_than_3s = norm(v_en_res_ii,axis=1) < 3*sigma_i
    data_ii = sites_info[less_than_3s]

    # loop
    while not np.array_equal(data_i,data_ii):
        data_i = data_ii
        w_lat_i,w_lon_i,w_i,w_lat_std_i,w_lon_std_i,w_std_i,w_xyz_i,w_xyz_std_i,v_en_i,weight_M_i = PMC(data_i) 
        v_en_res_i = res(v_en[less_than_3s],D_M[less_than_3s],w_xyz_i) 
        sigma_i = sigma_1(v_en_res_i,weight_M_i)
        v_en_res_ii = res(v_en,D_M,w_xyz_i) 
        less_than_3s = norm(v_en_res_ii,axis=1) < 3*sigma_i 
        data_ii = sites_info[less_than_3s] 
    rms = np.sqrt(np.mean(v_en_res_i**2,axis=0))*1e3*u.mm/u.yr
    
    df = data_i.reset_index(drop=True)
    df['ve'],df['vn'],df['ve_res'],df['vn_res'] = v_en_i[:,0],v_en_i[:,1],v_en_res_i[:,0],v_en_res_i[:,1]
    sites_retain_info = df
    
    n = len(df)
    epr = {'num_sites':n}
    epr['omega_cartesian'] = (w_xyz_i*u.rad/u.yr).to(u.deg/u.Ma)
    epr['omega_cartesian_std'] =  (w_xyz_std_i*u.rad/u.yr).to(u.deg/u.Ma)
    epr['omega_spherical'] = w_lat_i.to(u.deg), w_lon_i.to(u.deg), (w_i*u.rad/u.yr).to(u.deg/u.Myr),
    epr['omega_spherical_std'] = (w_lat_std_i*u.rad).to(u.deg), (w_lon_std_i*u.rad).to(u.deg), (w_std_i*u.rad/u.yr).to(u.deg/u.Myr)
    epr['rms'] = rms
    
    return epr,sites_retain_info    
