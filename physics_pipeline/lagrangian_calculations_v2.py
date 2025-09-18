# lagrangian_calculations.py
# 
# Calculate displacement, FTLE, and LAVD from backward-in-time Lagrangian trajectories initialized in a grid.

# LJK
#
# Date created: 11/09/22
# Last edited: 08/12/25

from glob import glob
from datetime import datetime
import numpy as np
import xarray as xr
import math
from config import *

def calc_ftle(ds,T,traj_lon_array,traj_lat_array):
    """ 
    Calculate the FTLE field for a time scale T.
    
    Input
        - ds: netCDF output from OceanParcels simulation
        - T : time of integration (in days)
        - traj_lon_array: longitude initialization array for Lagrangian simulation
        - traj_lat_array: latitude initialization array for Lagrangian simulation
    Output
        - ftle : 2d numpy array in units of days^-1
    """
    
    # Get the index for the end time of interest for the integration
    t0 = str(ds.time[0][0].values)
    t0_datetime = datetime.strptime(t0[:-3], '%Y-%m-%dT%H:%M:%S.%f')
    
    index = 0
    for t in ds.time[0]:
        t1 = str(t.values)
        t1_datetime = datetime.strptime(t1[:-3], '%Y-%m-%dT%H:%M:%S.%f')
        delta = (t0_datetime-t1_datetime)
        if str(delta) == '%s days, 0:00:00'%(T):
            break
        index += 1
        
    # Subset the data
    x_T_flat = ds.lon.isel(obs=index).values.squeeze()
    x_T = np.reshape(x_T_flat,(len(traj_lon_array),len(traj_lat_array)))                 
    y_T_flat = ds.lat.isel(obs=index).values.squeeze()
    y_T = np.reshape(y_T_flat,(len(traj_lon_array),len(traj_lat_array)))                  
    
    # Calculate backwards FTLE
    dxdy, dxdx = np.gradient(x_T, traj_lon_array, traj_lat_array)
    dydy, dydx = np.gradient(y_T, traj_lon_array, traj_lat_array)
    ny, nx = np.shape(dxdx)
    ftle = np.zeros([ny, nx])
    R = 6.3781*(10**6) #Earth's radius in meters
    for i in range(0, ny):
        for j in range(0, nx):
            J = np.array([[dxdx[i, j], dxdy[i, j]],[dydx[i, j], dydy[i, j]]])
            M = np.array([[R*R*np.cos(y_T[i, j]*np.pi/180.), 0], [0., R*R]])
            C = np.dot(np.dot(np.transpose(J), M), J)
            eig_lya, _ = np.linalg.eigh(C)
            ftle[i, j] = (1./np.abs(T))*np.log(np.sqrt(eig_lya.max()))
    return ftle,index

def distance_from_lat_lon(lat1,lon1,lat2,lon2):
    """
    Haversine distance in kilometers between two coordinate points. 
    Accepts negative (-180 to 180) or positive coordinate systems (0 to 360). 

    Input
        lat1,lon1: coordinates for point 1
        lat2,lon2: coordinates for point 2
        
    Output
        dist: distance in km between two points
    """

    R = 6371 # Radius of the earth in km
    delta_lat,delta_lon = math.radians(lat2-lat1),math.radians(lon2-lon1)
    lat1_radians,lat2_radians = math.radians(lat1),math.radians(lat2)
    
    a = math.sin(delta_lat/2) * math.sin(delta_lat/2) + math.cos(lat1_radians) * math.cos(lat2_radians) * ((math.sin(delta_lon/2))**2)
    dist = R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) # Distance in km
    
    return dist

def calc_displacement(ds,T,traj_lon_array,traj_lat_array):
    """ 
    Calculate 2D displacement field 
    
    Input
        - ds: netCDF output from OceanParcels simulation
        - T : time of integration (in days)
        - traj_lon_array: longitude initialization array for Lagrangian simulation
        - traj_lat_array: latitude initialization array for Lagrangian simulation
    Output
        - disp: 2d dimensional displacement array (km)
    """
    
    # Get the index for the end time of interest
    t0 = str(ds.time[0][0].values)
    t0_datetime = datetime.strptime(t0[:-3], '%Y-%m-%dT%H:%M:%S.%f')
    
    index = 0
    for t in ds.time[0]:
        t1 = str(t.values)
        t1_datetime = datetime.strptime(t1[:-3], '%Y-%m-%dT%H:%M:%S.%f')
        delta = (t0_datetime-t1_datetime)
        if str(delta) == '%s days, 0:00:00'%(T):
            break
        index += 1
        
    # Subset the data
    x0_flat = ds.lon.isel(obs=0).values.squeeze()
    y0_flat = ds.lat.isel(obs=0).values.squeeze()
    xT_flat = ds.lon.isel(obs=index).values.squeeze()
    yT_flat = ds.lat.isel(obs=index).values.squeeze()
  
    # Calculate backwards displacement
    disp_flat = []
    for l in np.arange(0,len(x0_flat)):
        disp_flat.append(distance_from_lat_lon(y0_flat[l],x0_flat[l],yT_flat[l],xT_flat[l]))
    
    return np.reshape(disp_flat,(len(traj_lon_array),len(traj_lat_array)))


def calc_LAVD_temporal_subset(T,output_freq,vort):
    """
    Input
        T: runtime of integration (days)
        output_freq: output frequency of OceanParcels simulation (hrs)
        vort: vorticity array from OceanParcels simulation
        
    Output
        LAVD: 2D Lagrangian averaged vorticity devation field (s^-1)
    """
    index = int((24/output_freq)*T)
    vort_avg_t = np.nanmean(np.array(vort),axis=0)[1:index+1] #Find average vorticity over the entire spatial domain at each time step
    LAVD = np.trapz(np.absolute(vort[:,1:index+1] - vort_avg_t), dx=output_freq*60*60, axis=1)/(T*24*60*60-output_freq*60*60)
    return LAVD

########################################################################################

# commenting out commands to use functions in other  scripts

#T = 28 # days to compute Lag. metrics

# OceanParcels initialization 
#output_freq = 6 # hours
#traj_lat_array = np.arange(10,35,0.03125)
#traj_lon_array = np.arange(215,245,0.03125)

#parcels_files = sorted(glob(data_dir+ parcels_trajs/*10.0_35.0_lon_215.0_245.0*.nc'))
#for file_path in parcels_files:
#    print(file_path)
#    filename = file_path.split('/')[-1]
#    ds = xr.open_dataset(file_path)
    
    # FTLE
#    ftle,index = calc_ftle(ds,T,traj_lon_array,traj_lat_array)
#    ftle_file_path = data_dir + 'FTLE/'' + filename[0:11] + '%s_day_FTLE'%(T) + filename[25:-3] 
#    np.save(ftle_file_path,ftle)
    
    # Displacement
#    disp = calc_displacement(ds,T,traj_lon_array,traj_lat_array)
#    disp_file_path = data_dir + 'disp/' + filename[0:11] + '%s_day_disp'%(T) + filename[25:-3] 
#    np.save(disp_file_path,disp)
    
    # Vorticity
#    vort_premask = ds.vort
#    vort = np.array(vort_premask.where(vort_premask != 0)) #filters out land values
#    LAVD = calc_LAVD_temporal_subset(T,output_freq,vort)
#    LAVD_file_path = data_dir + 'LAVD/' + filename[0:11] + '%s_day_LAVD'%(T) + filename[25:-3] 
#    np.save(LAVD_file_path,LAVD)
