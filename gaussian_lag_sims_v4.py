# Lagrangian simulations in a Gaussian distribution

# Simulate many particles from sample sites to get most probable pathways and dispersal
# v4: Include GLOBCURRENT option

# Lexi Jones-Kellett
# Date created: 10/21/24
# Last edited: 08/12/25

import csv,os,math,random,sys,cftime
import numpy as np
import xarray as xr
from glob import glob
from datetime import timedelta,datetime
from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, ErrorCode,Field
from config import *

#input_dataset = str(sys.argv[1]) # options: 'cmems','oscar','globcurrent0','globcurrent15'

print("Packages loaded...")
### Load Sample Site Locations ###

metadata = []
metadata_filename = data_dir + 'G4_metadata_no_physics.csv'
with open(metadata_dir + metadata_filename) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
#        if str(row[0]) == '60': # do this if you want to run a single site
        if (row[3] != 'nan') and (row[0] != '63') and (row[0] != '65'): # don't load in MQ blank sites, or closely spaced samples
            metadata.append(row)
metadata = np.array(metadata)
print("Sample locations loaded...")

### Generate Gaussian Cloud of Points ###
def generate_gaussian_circle_pts(std,num):
    """
    std = range; # of standard deviations
    num = number of points
    """
    points = []
    for i in np.arange(num):
        x = np.random.normal(loc=0, scale=std) # lons
        y = np.random.normal(loc=0, scale=std) # lats
        points.append([x,y])
    return np.array(points)

std, num = 0.025, 1000
points = generate_gaussian_circle_pts(std,num)
print("Gaussian cloud generated...")

### Advect Particles ### 

def DeleteParticle(particle, fieldset, time):
    """
    Stop a particle from running through the simulation and avoid an error that cancels the entire simulation.
    Used if a particle runs on land or out of bounds of the grid.
    """
    particle.delete()
    
runtime = 200 # days
timestep_mins = 20 # minutes
output_freq = 6 #hours

#for input_dataset in ['cmems','oscar_interim','globcurrent0']:
for input_dataset in ['cmems','oscar_final']:

    if input_dataset == 'cmems':
        parcels_input_files = sorted(glob(data_dir + 'CMEMS_data/dt_global_allsat_phy_l4_*.nc'))
        filenames = {'U': parcels_input_files,'V': parcels_input_files}
        variables = {'U': 'ugos','V': 'vgos'}
        dimensions = {'U': {'lon':'longitude','lat':'latitude','time':'time'},
                      'V': {'lon':'longitude','lat':'latitude','time':'time'}}
        output_dir = data_dir + 'parcels_trajs/gaussian_seeding_sample_sites/CMEMS_runs/'
        
    elif input_dataset == 'oscar_interim':
        parcels_input_files = sorted(glob(data_dir + 'OSCAR_2021/OSCAR_2021_transposed/oscar_currents_interim_*.nc'))
        filenames = {'U': parcels_input_files,'V': parcels_input_files}
        variables = {'U': 'u','V': 'v'}
        dimensions = {'U': {'lon':'lon','lat':'lat','time':'time'},
                      'V': {'lon':'lon','lat':'lat','time':'time'}}    
        output_dir = data_dir + 'parcels_trajs/gaussian_seeding_sample_sites/OSCAR_runs/'
        
    elif input_dataset == 'oscar_final':
        parcels_input_files = sorted(glob(data_dir + 'OSCAR_data/OSCAR_data_transposed/oscar_currents_final_*.nc'))
        filenames = {'U': parcels_input_files,'V': parcels_input_files}
        variables = {'U': 'u','V': 'v'}
        dimensions = {'U': {'lon':'lon','lat':'lat','time':'time'},
                      'V': {'lon':'lon','lat':'lat','time':'time'}}    
        output_dir = data_dir + 'parcels_trajs/gaussian_seeding_sample_sites/OSCAR_final_runs/'
        
    elif (input_dataset == 'globcurrent0') or (input_dataset == 'globcurrent15'):
        parcels_input_files = sorted(glob(data_dir + 'GLOBCURRENT/cmems_obs_mob_glo_phy-cur_my*.nc'))
        filenames = {'U': parcels_input_files,'V': parcels_input_files}
        variables = {'U': 'ugos','V': 'vgos'}
        dimensions = {'U': {'lon':'longitude','lat':'latitude','time':'time'},
                      'V': {'lon':'longitude','lat':'latitude','time':'time'}}
        output_dir = data_dir + 'parcels_trajs/gaussian_seeding_sample_sites/CMEMS_runs/'

    print("Advecting particles for sample site...")

    for m in metadata[1:]:
    #for m in metadata:
        print(m[0])

        fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, allow_time_extrapolation=True) 

        sample_date,sample_time = m[1],m[2]

        # CMEMS & GLORYS use np.datetime, and OSCAR uses cftime
        datetime_time = datetime(int(sample_date[0:4]),int(sample_date[4:6]),int(sample_date[6:8]),int(sample_time[0:2]),int(sample_time[3:5]))
        cftime_time = cftime.DatetimeJulian(int(sample_date[0:4]),int(sample_date[4:6]),int(sample_date[6:8]),int(sample_time[0:2]),int(sample_time[3:5]))

        lon_list = (float(m[4])+points[:,0]).tolist()
        lat_list = (float(m[3])+points[:,1]).tolist()

        # CMEMS & OSCAR are 0-360 lon, GLORYS is -180 to 180
        if input_dataset == 'glorys':
            lon_list = [i+360 for i in lon_list] #only works for negative lons

        if (input_dataset == 'oscar_interim') or (input_dataset == 'oscar_final'):
            pset = ParticleSet.from_list(fieldset = fieldset, pclass=JITParticle, 
                                         lon = lon_list, lat = lat_list, 
                                         depth = [0]*len(lon_list), time = [cftime_time]*len(lon_list))
        else:
            pset = ParticleSet.from_list(fieldset = fieldset, pclass=JITParticle, 
                                         lon = lon_list, lat = lat_list, 
                                         depth = [0]*len(lon_list), time = [datetime_time]*len(lon_list))

        output_file_path = '%ssamplesite%s_%sday_backward_runtime_%smin_timestep_%shr_output_freq_%s_gaussian_particles_%s_std_dist_%s_run3.nc'%(output_dir,m[0],runtime,timestep_mins,output_freq,num,std,input_dataset)

        ## Create output file
        if os.path.exists(output_file_path):
            f = open(output_file_path, 'w')
            f.close()
        output_file = pset.ParticleFile(name=output_file_path,outputdt=timedelta(hours=output_freq))

        ## Run the parcels simulation
        pset.execute(AdvectionRK4,
                    runtime=timedelta(days=runtime),
                    dt=-timedelta(minutes=timestep_mins),
                    output_file=output_file,
                    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

        #Write the trajectory information to the output file
        output_file.close() #deletes the Parcels output folder with npy files
