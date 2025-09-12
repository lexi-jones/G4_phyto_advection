# Compute the distance travelled north or south of the G4 Gaussian clouds
# Note: v1 was incorrectly subsetting the data, not saving all 1000 particle distances; fixed in this version

#LJK
#Date created: 09/12/25
#Last edited: 09/12/25

import csv,math,itertools
import pandas as pd
import numpy as np
import xarray as xr
import scipy.stats as spy
import * from config

metadata_df = pd.read_csv(data_dir + 'AVISO_metadata_w_sst_sss.csv')
metadata_df = metadata_df.set_index('sample_id')

def angleFromCoordinate(lat1, lon1, lat2, lon2):
    
    """
    Lat/lon input expected to be in degrees
    0° is North, 90° is East, 180° is South,270° is West
    
    Angle output is in degrees
    """
    
    dLon = math.radians(lon2) - math.radians(lon1)
    y = math.sin(dLon) * math.cos(math.radians(lat2))
    x = ((math.cos(math.radians(lat1)) * math.sin(math.radians(lat2))) - 
        (math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(dLon)))
    brng = math.degrees(math.atan2(y, x))
    brng = (brng + 360) % 360

    return brng

def distance_from_lat_lon(lat1,lon1,lat2,lon2):
    """
    Lat/lon input expected to be in degrees
    """
    R = 6371 # Radius of the earth in km
    delta_lat,delta_lon = math.radians(lat2-lat1),math.radians(lon2-lon1)
    lat1_radians,lat2_radians = math.radians(lat1),math.radians(lat2)
    
    a = math.sin(delta_lat/2) * math.sin(delta_lat/2) + math.cos(lat1_radians) * math.cos(lat2_radians) * ((math.sin(delta_lon/2))**2)
    dist = R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) # Distance in km 
    return dist

def get_traj_angles(sample_id,num_days,dataset):
    """
    sample_id: sample id (string)
    num_days: int
    dataset: 'cmems' or 'oscar'
    """
    if dataset == 'cmems':
        save_dir = data_dir + '/parcels_trajs/gaussian_seeding_sample_sites/CMEMS_runs/'
        file_path = '%ssamplesite%s_200day_backward_runtime_20min_timestep_6hr_output_freq_1000_gaussian_particles_0.025_std_dist_cmems.nc'%(save_dir,sample_id)
        
    elif dataset == 'oscar':
        save_dir = data_dir + '/parcels_trajs/gaussian_seeding_sample_sites/OSCAR_final_runs/'
        file_path = '%ssamplesite%s_200day_backward_runtime_20min_timestep_6hr_output_freq_1000_gaussian_particles_0.025_std_dist_oscar_final_run3.nc'%(save_dir,sample_id)
        
    ds = xr.open_dataset(file_path)

    # Get trajectory data
    traj_lats,traj_lons = [],[]
    for p in np.arange(0, len(ds['traj'])):
        traj_lats_temp,traj_lons_temp = [],[]
        for t in [0,num_days]:
            traj_lats_temp.append(np.array(ds.lat[p,t*4])) 
            traj_lons_temp.append(np.array(ds.lon[p,t*4]))
        traj_lats.append(traj_lats_temp)
        traj_lons.append(traj_lons_temp)
        
    # Get the trajectory angle and distance between the initial and end point to calculate degrees North travelled
    dists_north = []
    for p in np.arange(0,np.shape(traj_lats)[0]): # iterate through each particle
        angle=angleFromCoordinate(traj_lats[p][0],traj_lons[p][0],traj_lats[p][-1],traj_lons[p][-1])
        dist=distance_from_lat_lon(traj_lats[p][0],traj_lons[p][0],traj_lats[p][-1],traj_lons[p][-1])
        dists_north.append(math.cos(math.radians(angle))*dist)
        
    return dists_north

def thirty_day_dists(dataset):
    """
    dataset: 'cmems' or 'oscar'
    """
    for x in np.arange(30,210,30): # number of days for advection 
        print(x)
        
        dayx_dists_north = []
        for s in metadata_df.index:
            dayx_dists_north.append(get_traj_angles(s,x,dataset))
        
        if dataset == 'cmems':
            save_dir = data_dir + '/parcels_trajs/gaussian_seeding_sample_sites/CMEMS_runs/'
        elif dataset == 'oscar':
            save_dir = data_dir + '/parcels_trajs/gaussian_seeding_sample_sites/OSCAR_final_runs/'
             
        np.save(save_dir+'%s_day%s_dists_north_v2'%(dataset,x),dayx_dists_north)
        
print('Computing cmems...')
thirty_day_dists('cmems')
print('Computing oscar...')
thirty_day_dists('oscar')
