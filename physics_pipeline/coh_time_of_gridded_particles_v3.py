# Coherence Time of Gridded Particles

# Feeding a grid with eddy locations already removed, so there are "missing values" in the trajectory dataset.
# Original initialization was in an irregular shaped polygon in the NPSG

# LJK
# Date created: 06/09/25
# Last edited: 08/12/25

import time
import numpy as np
import pandas as pd
import xarray as xr
from numba import njit
from config import *

start_time = time.time()


@njit
def distance_from_lat_lon(lat1,lon1,lat2,lon2):
    """
    Returns the Haversine distance between two geographic coordinate points. Accepts negative (-180 to 180) or positive coordinate systems (0 to 360). 
    
    Using numpy instead of math package for numba efficiency

    Input
        lat1,lon1: geographic coordinates for point 1
        lat2,lon2: geographic coordinates for point 2
    Output
        dist: distance between the two point (units of kilometers)
    """
    R = 6371 # Radius of the earth in km
    delta_lat,delta_lon = np.radians(lat2-lat1),np.radians(lon2-lon1)
    lat1_radians,lat2_radians = np.radians(lat1),np.radians(lat2)
    a = np.sin(delta_lat/2) * np.sin(delta_lat/2) + np.cos(lat1_radians) * np.cos(lat2_radians) * ((np.sin(delta_lon/2))**2)
    dist = R * 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a)) # Distance in km
    return dist

# Open trajectory dataset (computed with `run_parcels_CMEMS_NPSG_outside_eddies.py`)
traj_ds = xr.open_dataset(data_dir + 'parcels_trajs/2021-11-23_200days_runtime_20min_timestep_0.02_spatial_res_24hr_output_freq_outside_15km_eddy_buff_NPSG_domain.nc')

# Read in lats and lons to help increase speed later when these are indexed frequently
traj_lats = traj_ds.lat.values
traj_lons = traj_ds.lon.values

# Assign indices to each particle based on what would be a regular initialization grid
lat_to_i = {lat: i for i, lat in enumerate(np.unique(traj_lats[:,0]))}
lon_to_j = {lon: j for j, lon in enumerate(np.unique(traj_lons[:,0]))}

# Map the actual coordinates to the gridded indeces
lat_index = np.array([lat_to_i[lat] for lat in traj_lats[:,0]])
lon_index = np.array([lon_to_j[lon] for lon in traj_lons[:,0]])

# Add the lat/lon indeces to a pandas dataframe, assigned to each traj
lookup_df = pd.DataFrame({'particle': np.arange(traj_ds.traj.size),
                    'lat_index': lat_index,
                    'lon_index': lon_index})

# MultiIndex for quick lookup (more than one variable can be set to index)
lookup_df.set_index(['lat_index', 'lon_index'], inplace=True) 

# Neighbor indeces; circle of points around center with 0.1 degree diameter (97 points total in this setup)

# V3 Neighbor indeces: circle of points around center; ~8.3 km radius in 0.02 res field
offsets = []
for y in range(-2,3):
    offsets.extend([(x,y) for x in range(-3,4)])
offsets.extend([(x,-3) for x in range(-2,3)])
offsets.extend([(x,3) for x in range(-2,3)])

# Calculate a weight based on the Euclidean distances for each point
weights = []
for point in offsets: 
    weights.append((4-np.linalg.norm(np.array(point)-np.array([0,0]))))
    
def weighted_mean(values, weights):
    weighted_sum = sum(v * w for v, w in zip(values, weights))
    total_weight = sum(weights)    
    return weighted_sum / total_weight

# Iterate through the particle trajectories
num_trajs = traj_ds.traj.size
coh_time = np.full(num_trajs,np.nan) # empty array to store coherence times

#successful_cloud_count = 0
for lat_i in np.unique(lat_index)[4::8]: # iterate through lats, skipping by 8 rows (cloud size + 1)
    elapsed_mins = round((time.time() - start_time)/60,1)
    print('On lat %s; %s mins elapsed ...'%(lat_i,elapsed_mins))
    #print('%s successful clouds'%(successful_cloud_count))
    print('%s coh times recorded'%(len(coh_time[~np.isnan(coh_time)])))
    
    lon_i = 4
    while lon_i < max(np.unique(lon_index))-4: # find viable clouds by lon
        center = (lat_i, lon_i) # center point of the cloud
        neighbors = [(center[0] + di, center[1] + dj) for di, dj in offsets] # all of the cloud neighbors
        
        try: # Check if there's all neighbor data available
            neighbor_ds_locs = np.array([lookup_df.loc[nbr][0] for nbr in neighbors])

        except KeyError:
            lon_i += 1
            continue # sends the code back up the while loop
            
        # Intial distance of neighbor particles from center 
        #(i == index of center particle, n == index of neighbor)
        # Converting to float for numba efficiency (it doesn't like xr data)

        center_traj = int(lookup_df.loc[center][0])
        max_dist_t0 = np.max([distance_from_lat_lon(traj_lats[center_traj,0],traj_lons[center_traj,0],
                                                traj_lats[n,0],traj_lons[n,0]) for n in neighbor_ds_locs])

        # Get mean particle distance from center of mass for each time step
        t = 0
        for t in np.arange(0,traj_ds.obs.size,10): # iterate every 10 days 

            # Get the distance of each neighbor from the central particle
            dists = [np.abs(distance_from_lat_lon(traj_lats[center_traj,t],traj_lons[center_traj,t],
                                            traj_lats[n,t],traj_lons[n,t])) for n in neighbor_ds_locs]
            
            # Weighted Mean = (Σ(value * weight)) / (Σ(weight))
            weighted_dist = weighted_mean(dists, weights)
            
            # Mean fraction of distance compared to initialization size
            dispersal_size = weighted_dist/max_dist_t0

            # Stop looking if dispersal spreads twice original distance
            if (dispersal_size > 2):
                coh_time[center_traj] = t
                break

            # In the case we get to the end, coh time is 200d
            elif t == traj_ds.obs.size-1:
                coh_time[center_traj] = t

        #successful_cloud_count += 1
        lon_i += 8 # success, go to the next possible cloud location
            
np.save(data_dir+'parcels_trajs/coh_time_of_outside_eddy_test_domain_0.02_res_weighted_mean_dist_7_by_7_cloud.npy',coh_time)
