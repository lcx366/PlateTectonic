import numpy as np
from numpy.linalg import norm
from astropy.coordinates import spherical_to_cartesian,cartesian_to_spherical
from astropy.coordinates import EarthLocation,Longitude,Latitude
from astropy import units as u
from sphericalpolygon import Sphericalpolygon

import pandas as pd
from pkg_resources import resource_filename
from tqdm import trange

from os import path,makedirs
from pathlib import Path

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.font_manager import FontProperties

from ..utils import Const

def vel_unit_conver(vel):
    return vel.to(u.rad*u.mm/u.yr).value*(u.mm/u.yr)   

def H_unit_conver(H):
    return H.to(u.rad*u.km**2*u.kg/u.Ma).value*(u.km**2*u.kg/u.Ma)       

def plate_attribute(vertices,thickness,density):

    Re = Const.r0*u.km # volumetric radius of the Earth, [km] 
    rho = thickness*density # area density in km*g/Cm3
    
    vertices = np.array(vertices)
    if (vertices[0] != vertices[-1]).all():
        vertices = np.append(vertices,[vertices[0]],axis=0) # create a closed spherical polygon   
    polygon = Sphericalpolygon(vertices) 
    area = polygon.area(Re)
    perimeter = polygon.perimeter(Re)
    compactness = polygon.compactness()
    centroid = polygon.centroid(Re)
    inertia = polygon.inertia(Re,rho).to(u.km**2*u.kg)

    inertia_tensor = np.diag(inertia[:3])
    inertia_tensor[0,1] = inertia_tensor[1,0] = inertia[3]
    inertia_tensor[0,2] = inertia_tensor[2,0] = inertia[4]
    inertia_tensor[1,2] = inertia_tensor[2,1] = inertia[5]

    inertia_mag = norm(inertia_tensor)

    info = {'polygon':polygon,'area':area,'perimeter':perimeter,'compactness':compactness,'centroid':centroid,'inertia':inertia,\
    'inertia_tensor':inertia_tensor,'inertia_mag':inertia_mag,'thickness':thickness,'density':density}

    return info    

def read_bnds_pole(data_path,bnds_dir,pole_file,modelname):
   
    bnds_path = data_path + bnds_dir
    # prevent converting 'NA' to NaN by adding converters={ "column name":str} to read_csv
    data = pd.read_csv(data_path+pole_file,converters={ "Ab":str}) 

    info = {}
             
    df = data.dropna().reset_index(drop=True)
    n = len(df)
    plates_name,plates_Ab = df['Plate'],df['Ab']
    omega_lats = np.array(df['Lat. °N'])*u.deg
    omega_lons = np.array(df['Lon. °E'])*u.deg
    omega_mags = np.array(df['w °/Ma'])*u.deg/u.Ma  

    home = str(Path.home())
    direc = home + '/src/platemotion-data/'
    if not path.exists(direc): makedirs(direc)
    output_csv = open(direc+pole_file[:-4]+'_summary.csv','w')

    head = ','.join(('Plate','Ab','ωx','ωy','ωz','Lat.(°N)','Lon.(°E)','ω(°/Ma)','Hx','Hy','Hz','Lat.(°N)','Lon.(°E)','H(10^24 km^2 kg/Ma)',\
                'Q11','Q22','Q33','Q12','Q13','Q23','‖Q‖F(10^27 km^2 kg)','Area(km^2)','Perimeter(km)','Compactness','Lat.(°N)','Lon.(°E)','Depth(km)'))
    output_csv.write(head+'\n')

    for i in trange(n,ascii=True,desc='Loading plates from '+ modelname):

        plate = Plate.from_file(bnds_path + plates_Ab[i],skiprows=1) 
        plate_omega = (omega_lats[i],omega_lons[i],omega_mags[i])
        plate.set_omega(plate_omega,'spherical')
        plate.set_name(plates_name[i])

        info.update({plates_Ab[i]:plate})

        str_omega_cartesian = ','.join(['{:.3f}'.format(ele) for ele in plate.omega_cartesian.value])
        plate_omega_lat,plate_omega_lon,plate_omega_mag = plate.omega_spherical
        str_omega_spherical = '{:.2f},{:.2f},{:.3f}'.format(plate_omega_lat.value,plate_omega_lon.value,plate_omega_mag.value)
        str_H_cartesian = ','.join(['{:.3f}'.format(ele) for ele in plate.H_cartesian.value/1e24]) # 10^24 km2 kg / Ma
        plate_H_lat,plate_H_lon,plate_H_mag = plate.H_spherical
        str_H_spherical = '{:.2f},{:.2f},{:.3f}'.format(plate_H_lat.value,plate_H_lon.value,plate_H_mag.value/1e24)
        table1 = ','.join((plate.name,plates_Ab[i],str_omega_cartesian,str_omega_spherical,str_H_cartesian,str_H_spherical))

        str_inertia = ','.join(['{:.3f}'.format(ele) for ele in plate.inertia.value/1e27]) # 10^27 km2 kg
        str_inertia_mag = '{:.3f}'.format(plate.inertia_mag.value/1e27) # 10^27 km2 kg   
        str_area = '{:.2f}'.format(plate.area.value)
        str_perimeter = '{:.2f}'.format(plate.perimeter.value)
        str_compactness = '{:.2f}'.format(plate.compactness)
        str_centroid_latlon = ','.join(['{:.2f}'.format(ele.value) for ele in plate.centroid[:2]]) 
        str_centroid_depth = '{:.2f}'.format(plate.centroid[-1].value) 
        table2 = ','.join((str_inertia,str_inertia_mag,str_area,str_perimeter,str_compactness,str_centroid_latlon,str_centroid_depth))   
        output_csv.write(table1+','+table2+'\n')     
    output_csv.close()
    summary = pd.read_csv(direc+pole_file[:-4]+'_summary.csv') 
    info.update({'summary':summary})
    return PlateMotion(info)                    

def euler_matrix(lat,lon,r):
    # Design matrix of plate motion equations      

    A = np.empty((2,3))
    A[0] = -np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)                                                                                                    
    A[1] = np.sin(lon),-np.cos(lon),0
    return A*r   

def velocity_at_point(omega_cartesian, location, mode):
                
    if mode == 'geocentric':
        lat,lon,r = location
        pos_cartesian = spherical_to_cartesian(r,lat,lon)
    else:
        if mode == 'geodetic':
            lat_geode,lon,h = location
            pos_cartesian = EarthLocation.from_geodetic(lon,lat_geode,h)
            x,y,z = pos_cartesian.x,pos_cartesian.y,pos_cartesian.z
            pos_cartesian = (x,y,z)

        elif mode == 'cartesian':
            x,y,z = location
            pos_cartesian = location

        r, lat_rad, lon_rad = cartesian_to_spherical(x,y,z) 
        lat,lon = lat_rad.to(u.deg),lon_rad.to(u.deg)
                
    velocity_en = euler_matrix(lat,lon,r)@omega_cartesian
    velocity_cartesian = np.cross(omega_cartesian,pos_cartesian)

    velocity_xyz = vel_unit_conver(velocity_cartesian)
    velocity_en = vel_unit_conver(velocity_en)
    speed = norm(velocity_xyz)
    azimuth = np.arctan2(velocity_en[0],velocity_en[1]).to(u.deg)

    return velocity_xyz,velocity_en,speed,azimuth

def velocities_at_points(omega_cartesian, locations, mode):
    
    loca0,loca1,loca2 = locations
    n = len(loca0)
    __ = [[] for i in range(4)]
    for i in range(n):
        loca = loca0[i],loca1[i],loca2[i]
        _ = velocity_at_point(omega_cartesian,loca, mode)
        for j in range(4): __[j].append(_[j]) 
        
    velocity_xyz = u.Quantity(__[0])
    velocity_en = u.Quantity(__[1])
    speed = u.Quantity(__[2])
    azimuth = u.Quantity(__[3]) 
        
    return velocity_xyz,velocity_en,speed,azimuth             

class Velocity(object):
    '''
    class Velocity

    - attributes:
        - vertices: vertices of a closed spherical polygon in form of [[lat_0,lon_0],...,[lat_n,lon_n]]
        - lats: latitudes of the spherical polygon in degrees
        - lons: longitudes of the spherical polygon in degrees
        - orientation: vertices arrangement; it can be counterclockwise or clockwise

    - methods:
        - contains_points: determine if a single point or multiple points are inside a spherical polygon.
        - area: calculate the area or mass of a spherical polygon.
        - perimeter: calculate the perimeter of a spherical polygon.
        - centroid: identify the location of the centroid of a spherical polygon.
        - inertia: compute the geometrial or physical moment of inertia tensor of a spherical polygon.
        
    ''' 

    def __init__(self,info):

        self.info = info
        for key in info.keys():
            setattr(self, key, info[key])


    def __repr__(self):
    
        return 'instance of class Velocity'         

class Plate(object):
    '''
    class Plate

    - attributes:
        - vertices: vertices of a closed spherical polygon in form of [[lat_0,lon_0],...,[lat_n,lon_n]]
        - lats: latitudes of the spherical polygon in degrees
        - lons: longitudes of the spherical polygon in degrees
        - orientation: vertices arrangement; it can be counterclockwise or clockwise

    - methods:
        - contains_points: determine if a single point or multiple points are inside a spherical polygon.
        - area: calculate the area or mass of a spherical polygon.
        - perimeter: calculate the perimeter of a spherical polygon.
        - centroid: identify the location of the centroid of a spherical polygon.
        - inertia: compute the geometrial or physical moment of inertia tensor of a spherical polygon.
        
    ''' 

    def __init__(self,info):

        self.info = info
        for key in info.keys():
            setattr(self, key, info[key])


    def __repr__(self):
    
        return 'instance of class Plate'

    def from_array(vertices,thickness=100*u.km,density=3.1*u.g/u.cm**3):   
        '''
        Create an instance of class Sphericalpolygon from numpy array.
    
        Usage:
        polygon = Sphericalpolygon.from_array(vertices)
        thickness in km
        density in g/cm3

        Inputs:
        vertices -> [float 2d array] Vertices that make up the polygon in form of [[lat_0,lon_0],...,[lat_n,lon_n]] with unit of degrees. 
        If the first vertex is not equal to the last one, a point is automatically added to the end of the vertices sequence to form a closed polygon. 
        Vertices can be arranged either counterclockwise or clockwise.

        Outputs:
        polygon -> an instance of class Sphericalpolygon 

        Note: The spherical polygon has a latitude range of [-90°,90°] and a longitude range of [-180°,180°] or [0°,360°].
        '''
        info = plate_attribute(vertices,thickness,density)
        return Plate(info)
    
    def from_file(filename,thickness=100*u.km,density=3.1*u.g/u.cm**3,skiprows=0):
        '''
        Create an instance of class Sphericalpolygon from a file.
    
        Usage:
        polygon = Sphericalpolygon.from_file(filename,[skiprows])

        Inputs:
        filename -> [str] input file that lists vertices of a polygon in form of 

            # polygon info, such as name, and soure, etc.
            # comments
            #
            lat_0,lon_0
            ...
            lat_n,lon_n 

        with unit of degrees. If the first vertex is not equal to the last one, a point is automatically added to the end of the vertices sequence to form a closed polygon. 
        Vertices can be arranged either counterclockwise or clockwise.

        Parameters:
        skiprows -> [int, optional] skip the first `skiprows` lines, including comments; default: 0.

        Outputs:
        polygon -> an instance of class Sphericalpolygon 

        Note: The spherical polygon has a latitude range of [-90°,90°] and a longitude range of [-180°,180°] or [0°,360°].
        '''
        vertices = np.loadtxt(filename,skiprows=skiprows) 
        info = plate_attribute(vertices,thickness,density)
        return Plate(info)

    def set_name(self,name):    
        '''
        name: str
        '''
        info = self.info
        info['name'] = name
        setattr(self, 'name', name) 

    def set_epr(self,epr):    
        '''
        set Euler Pole parameters
        '''
        info = self.info
        info['epr'] = epr
        setattr(self, 'epr', epr)     

    def set_sites(self,sites):    
        '''
        set sites info
        '''
        info = self.info
        info['sites'] = sites
        setattr(self, 'sites', sites)     

    def set_omega(self,omega,mode):

        '''
        omega_spherical : lat, lon, omega_mag
        omega_cartesian : omega_x,omega_y,omega_z

        '''

        inertia_tensor = self.inertia_tensor
        info = self.info

        if mode == 'spherical':
            omega_lat, omega_lon, omega_mag = omega
            omega_spherical = Latitude(omega_lat), Longitude(omega_lon), omega_mag
            omega_cartesian = spherical_to_cartesian(omega_mag,omega_lat,omega_lon)
            omega_cartesian = u.Quantity(omega_cartesian)
        elif mode == 'cartesian':
            omega_x,omega_y,omega_z = omega  
            omega_cartesian = omega
            omega_mag, omega_lat_rad, omega_lon_rad = cartesian_to_spherical(omega_x,omega_y,omega_z) 
            omega_lat,omega_lon = omega_lat_rad.to(u.deg),omega_lon_rad.to(u.deg)
            omega_spherical = (omega_lat,omega_lon,omega_mag)
        else:
            raise Exception("Only 'spherical' and 'cartesian' are avaliable.")    

        # Angular Momentum 
        H_cartesian = np.dot(inertia_tensor,omega_cartesian) 
        H_cartesian = H_unit_conver(H_cartesian)
        H_x,H_y,H_z = H_cartesian
        H_mag, H_lat_rad, H_lon_rad = cartesian_to_spherical(H_x,H_y,H_z) 
        H_lat,H_lon = H_lat_rad.to(u.deg),H_lon_rad.to(u.deg)
        H_spherical = (H_lat,H_lon,H_mag)

        info.update({'omega_cartesian':omega_cartesian,'omega_spherical':omega_spherical,\
            'H_cartesian':H_cartesian,'H_spherical':H_spherical})

        for key in info.keys(): setattr(self, key, info[key])  

    def contains_points(self,points):
        '''
        Determine if a single point or multiple points are inside the given spherical polygon.

        Usage: 
        flag = polygon.contains_points([30,102])
        flags = polygon.contains_points([[30,102],[-75,33]])

        Inputs:
        points -> [float array with 2 elements or float 2d array] single point or multiple points to be determined in form of [lat,lon] or [[lat_0,lon_0],..,[lat_n,lon_n]] with unit of degrees.

        Outputs:
        flags -> [bool or bool array] If True, the point is inside the polygon, otherwise, it is outside.
        '''

        flags = self.polygon.contains_points(points)
        
        return flags   	     

    def velocity_at(self, location, mode):
        '''
        Calculate the geometrical or physical(if the area density is given) moment of inertia tensor of a specific spherical polygon over a sphere with a radius of R.

        Usage:
        inertia = polygon.inertia()
        inertia = polygon.inertia(6378.137,81)

        Parameters:
        R -> [optional, float, default = 1] sphere radius
        rho -> [optional, float, default = 1] area density of the spherical polygon

        Outputs:
        inertia -> [float array with 6 elements] symmetrical inertia tensor with six independent components.
        The first three components are located diagonally, corresponding to M_{11}, M_{22}, and M_{33}; the last three components correspond to M_{12}, M_{13}, and M_{23}.
        '''
                
        omega_cartesian = self.omega_cartesian

        try:
            n = len(location[0])
        except:
            n = 1
   
        if n > 1:
            velocity_xyz,velocity_en,speed,azimuth = velocities_at_points(omega_cartesian,location,mode)
        else:
            velocity_xyz,velocity_en,speed,azimuth = velocity_at_point(omega_cartesian,location,mode)
                
        v_info = {'xyz':velocity_xyz,'en':velocity_en,'speed':speed,'azimuth':azimuth}

        return Velocity(v_info)

    def plot(self,figname=None,zoom=1,scale=500,U=20):

        lats_bnds = self.polygon.lats
        lons_bnds = self.polygon.lons 
        centroid_lat,centroid_lon = self.centroid[:2]

        fig = plt.figure(dpi=200)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.NearsidePerspective(centroid_lon, centroid_lat,satellite_height=6.5e7/zoom))

        ax.set_global()
        # set the gridlines
        ax.gridlines(color='gray', linestyle = '--', xlocs = np.arange(-180,180,30), ylocs = np.linspace(-80,80,9))

        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.LAND,zorder=0)
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.RIVERS,zorder=1)
        ax.add_feature(cfeature.LAKES,zorder=2)

        ax.plot(lons_bnds, lats_bnds, color='m',transform=ccrs.Geodetic(),linewidth=1) 
        ax.fill(lons_bnds, lats_bnds, color='coral', transform=ccrs.Geodetic(),alpha=0.4) 

        if 'omega_cartesian' in self.info:

            lons_bnds_min,lons_bnds_max = lons_bnds.min(),lons_bnds.max()
            lats_bnds_min,lats_bnds_max = lats_bnds.min(),lats_bnds.max()

            if np.abs(lons_bnds_min - lons_bnds_max) > 355:
                if self.contains_points([89,0]): # North pole
                    lats_bnds_max = 85 
                elif self.contains_points([-89,0]): # South pole 
                    lats_bnds_min = -85

            lats_fill = np.linspace(lats_bnds_min,lats_bnds_max,10)
            lons_fill = np.linspace(lons_bnds_min,lons_bnds_max,10)    

            lons_mesh, lats_mesh = np.meshgrid(lons_fill,lats_fill)
            lons_flat = lons_mesh.flatten()
            lats_flat = lats_mesh.flatten()  

            inplate_flag = self.contains_points(np.array([lats_flat,lons_flat]).T)

            lons_inplate = lons_flat[inplate_flag]
            lats_inplate = lats_flat[inplate_flag]
            h_inplate = np.zeros_like(lats_inplate)*u.m  

            locations = [lats_inplate*u.deg,lons_inplate*u.deg,h_inplate]
            inplate_vel_en = self.velocity_at(locations,'geodetic').en
            inplate_ve,inplate_vn = inplate_vel_en[:,0],inplate_vel_en[:,1]

            font0 = FontProperties()
            font0.set_size('small')

            q = ax.quiver(lons_inplate, lats_inplate, inplate_ve, inplate_vn,transform = ccrs.PlateCarree(),scale=scale,width=0.004,color='g',regrid_shape=20)
            ax.quiverkey(q, X=0.65, Y=0.85, U=U,label=r'$20 \, mm/yr$', labelpos='E',coordinates='figure',fontproperties=font0)

        if figname is None:
            plt.show()
        else:
            plt.savefig(figname)    


class PlateMotion(object):
    '''
    class Plate

    - attributes:
        - vertices: vertices of a closed spherical polygon in form of [[lat_0,lon_0],...,[lat_n,lon_n]]
        - lats: latitudes of the spherical polygon in degrees
        - lons: longitudes of the spherical polygon in degrees
        - orientation: vertices arrangement; it can be counterclockwise or clockwise

    - methods:
        - contains_points: determine if a single point or multiple points are inside a spherical polygon.
        - area: calculate the area or mass of a spherical polygon.
        - perimeter: calculate the perimeter of a spherical polygon.
        - centroid: identify the location of the centroid of a spherical polygon.
        - inertia: compute the geometrial or physical moment of inertia tensor of a spherical polygon.
        
    ''' 

    def __init__(self,info):

        for key in info.keys():
            setattr(self, key, info[key])


    def __repr__(self):
    
        return 'instance of class PlateMotion'

    def loadmodel(modelname):
        data_path = resource_filename('platemotion', 'data/')

        if modelname == 'NNR-MORVEL56':      
            bnds_dir,pole_file = 'NnrMRVL_PltBndsLatLon/','NNR_MRV56.csv'

        elif modelname == 'GSRMv2.1':   
            bnds_dir,pole_file = 'NnrGSRMv2.1_PltBndsLatLon/','NNR_GSRMv2.1.csv' 
        res = read_bnds_pole(data_path,bnds_dir,pole_file,modelname)

        return res  

    def initialize(plates):
        info = {}
        try:
            n = len(plates)
            for i in range(n):
                info.update({plates[i].name:plates[i]})
        except:
            n = 1   
            info.update({plates.name:plates})     
        return PlateMotion(info)            

    def add_plate(self,plates):
        info = vars(self)
        try:
            n = len(plates)
            for i in range(n):
                info.update({plates[i].name:plates[i]})
        except:
            n = 1   
            info.update({plates.name:plates})     
        return PlateMotion(info)        