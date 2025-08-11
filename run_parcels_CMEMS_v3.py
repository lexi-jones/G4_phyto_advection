# Lexi Jones
# Date Created: 07/15/21
# Last Edited: 08/08/25

# Run an OceanParcels simulation with functions from 'functions_for_parcels'. 
# Format: python run_parcels_CMEMS_v3.py lat_start lat_stop lon_start lon_stop spatial_step date_input time_step output_freq runtime runtime_unit backwards
#	lat_start,lat_stop,lon_start,lon_stop: latitude and longitude bounds to intialize particles
#	spatial_step: lat/lon spacing between initial particles
#	date_input: start date with format YYYYMMDD
#	time_step: (minutes) the amount of time that passes in the fieldset before calculating new particle positions & attributes
#	output_freq: (hours) frequency to output information about the particle positions & attributes
#	runtime: how long to run the simulation
#	runtime_units: 'days' or 'hours', corresponds with runtime
#	backwards: 'y' if the particles should be run backwards in time, 'n' if the particles should be run forwards in time

from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile, DiffusionUniformKh,AdvectionRK4_3D,ErrorCode,Field
import numpy as np
import math,os,time,sys
from datetime import timedelta,datetime
import xarray as xr
from glob import glob
from config import *

############# FUNCTIONS ##############

from functions_for_parcels import CalcVort, DeleteParticle,particle_grid2d,simulate_particles2d,simulate_particles2d_backwards

class SpinnyParticle(JITParticle):
    u = Variable('u',dtype=np.float64)
    v = Variable('v',dtype=np.float64)
    vort = Variable('vort',dtype=np.float64)

########## RUN PARCELS #############

print(sys.argv)

#FROM USER INPUT
#Note: sys.argv[0] is the name of the script
lat_start,lat_stop = float(sys.argv[1]), float(sys.argv[2])
lon_start,lon_stop = float(sys.argv[3]), float(sys.argv[4])
spatial_step = float(sys.argv[5])

date_input = sys.argv[6]
start_year = int(date_input[0:4])
start_month = int(date_input[4:6])
start_day = int(date_input[6:8])

time_step = int(sys.argv[7])
output_freq = int(sys.argv[8])
runtime = int(sys.argv[9])
runtime_unit = sys.argv[10]
backwards = str(sys.argv[11])

start_date = datetime(start_year,start_month,start_day)

# Create fieldset
input_dir = data_dir + 'CMEMS_data/'
parcels_input_files = sorted(glob(input_dir+'dt_global_allsat_phy_l4_*.nc'))

filenames = {'U': parcels_input_files,'V': parcels_input_files}
#variables = {'U': 'uvel','V': 'vvel'} #name of the velocity variables in the netCDF file
variables = {'U': 'ugos','V': 'vgos'}

dimensions = {'U': {'lon':'longitude','lat':'latitude','time':'time'},
              'V': {'lon':'longitude','lat':'latitude','time':'time'}}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
print('Fieldset created.')

pset_dynamic,num_particles = particle_grid2d(fieldset,SpinnyParticle,[lat_start,lat_stop,spatial_step],[lon_start,lon_stop,spatial_step],start_date)
print('Particles initialized.')

output_dir = data_dir + 'parcels_trajs/'
output_file_path = '%s%s_%s%s_runtime_%smin_timestep_particle_start_lat_%s_%s_lon_%s_%s_spatial_step_%s_%shr_output_freq.nc'%(
output_dir,str(start_date)[0:10],runtime,runtime_unit,time_step,lat_start,lat_stop,lon_start,lon_stop,spatial_step,output_freq)

print('Output file: %s'%(output_file_path))

# Execute particle simulation
start_time = time.time()

if backwards == 'y':
    simulate_particles2d_backwards(pset_dynamic,output_file_path,runtime,runtime_unit,time_step,output_freq)
else:
    simulate_particles2d(pset_dynamic,output_file_path,runtime,runtime_unit,time_step,output_freq)
print("--- %s seconds ---" % (time.time() - start_time))
