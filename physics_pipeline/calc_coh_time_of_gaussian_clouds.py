# Calculate coherence time as the time when the mean particles spread from initial distance is 2x
# LJK
# Date created: 03/17/25
# Last edited: 08/12/25

import math
import numpy as np
import xarray as xr
import pandas as pd
from config import *


metadata_df = pd.read_csv(data_dir  + 'AVISO_eddies/AVISO_metadata.csv',index_col=0)
metadata_df = metadata_df[12:] # sample 16 & later
metadata_df

def distance_from_lat_lon(lat1,lon1,lat2,lon2):
    """
    Returns the distance between two geographic coordinate points. Accepts negative (-180 to 180) or positive coordinate systems (0 to 360). 

    Input
        lat1,lon1: geographic coordinates for point 1
        lat2,lon2: geographic coordinates for point 2
    Output
        dist: distance between the two point (units of kilometers)
    """
    R = 6371 # Radius of the earth in km
    delta_lat,delta_lon = math.radians(lat2-lat1),math.radians(lon2-lon1)
    lat1_radians,lat2_radians = math.radians(lat1),math.radians(lat2)
    a = math.sin(delta_lat/2) * math.sin(delta_lat/2) + math.cos(lat1_radians) * math.cos(lat2_radians) * ((math.sin(delta_lon/2))**2)
    dist = R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) # Distance in km
    return dist

def dispersal_calc(data_product,run_str):
    """
    data_product: 'cmems' or 'oscar_final'
    run_str: '','_run1','_run2','_run3'
    """

    if data_product == 'cmems':
        output_dir = data_dir + 'parcels_trajs/gaussian_seeding_sample_sites/CMEMS_runs/'
    elif data_product == 'oscar_final':
        output_dir = data_dir + 'parcels_trajs/gaussian_seeding_sample_sites/OSCAR_final_runs/'
    
    all_dist_tx_diff = []
    for m in list(metadata_df['sample_id']):
        
        output_file_path = '%ssamplesite%s_200day_backward_runtime_20min_timestep_6hr_output_freq_1000_gaussian_particles_0.025_std_dist_%s%s.nc'%(output_dir,m,data_product,run_str)
        ds = xr.open_dataset(output_file_path)
        
        traj_lats,traj_lons,traj_times = [],[],[]
        for i in np.arange(0,len(ds.lat)):
            traj_lats.append(np.array(ds.lat[i])) 
            traj_lons.append(np.array(ds.lon[i]))
            traj_times.append(np.array(ds.time[i]))

        # Get mean particle distance from center of mass in initial cloud
        lats_t0 = [i[0] for i in traj_lats]
        lons_t0 = [i[0] for i in traj_lons]
        center_lat_t0,center_lon_t0 = np.mean(lats_t0),np.mean(lons_t0)
        max_dist_t0 = np.max([distance_from_lat_lon(center_lat_t0,center_lon_t0,lats_t0[i],lons_t0[i]) for i in np.arange(0,len(lons_t0))])

        # Get mean particle distance from center of mass for each time step
        dist_tx = []
        for x in np.arange(0,len(traj_times[0])): # iterate through each time
            lats_tx = [t[x] for t in traj_lats]
            lons_tx = [t[x] for t in traj_lons]
            center_lat_tx,center_lon_tx = np.mean(lats_tx),np.mean(lons_tx)
            dists = [np.abs(distance_from_lat_lon(center_lat_tx,center_lon_tx,lats_tx[i],lons_tx[i])) for i in np.arange(0,len(lons_tx))]
            dist_tx.append(np.mean(dists))

        # fraction of distance compared to initialization size
        all_dist_tx_diff.append([i/max_dist_t0 for i in dist_tx])

    coherence_time = []
    for i in all_dist_tx_diff:
        dispersal_size = i[0]
        count = 1
        while (dispersal_size < 2) and (count < len(i)):
            dispersal_size = i[count]
            count += 1
        coherence_time.append(count)
        
    return coherence_time


print('Run0')
cmems_run0 = dispersal_calc('cmems','')

print('Run1')
cmems_run1 = dispersal_calc('cmems','_run1')

print('Run2')
cmems_run2 = dispersal_calc('cmems','_run2')

print('Run3')
cmems_run3 = dispersal_calc('cmems','_run3')

coh_df = pd.DataFrame({'sample_id': metadata_df.sample_id, 
                       'cmems_run0': cmems_run0, 
                       'cmems_run1': cmems_run1,
                       'cmems_run2': cmems_run2,
                       'cmems_run3': cmems_run3,})
coh_df = coh_df.set_index('sample_id')
coh_df.to_csv(data_dir + 'parcels_trajs/gaussian_seeding_sample_sites/CMEMS_coherence_time_by_run.csv')
