# Set up for NPSG particles outside of eddies

# Lexi Jones-Kellett
# Date Created: 06/06/25
# Last Edited: 08/12/25

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

from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile, ErrorCode, Field
import numpy as np
import math,os,time,sys
from datetime import timedelta,datetime
import xarray as xr
from glob import glob

# Eddy stuff
from shapely.geometry import Point, Polygon
from shapely import geometry
import shapely
from shapely.ops import unary_union
from pyproj import Transformer, CRS
from shapely.ops import unary_union, transform
from shapely.strtree import STRtree

############# FUNCTIONS ##############

from functions_for_parcels import CalcVort, DeleteParticle,particle_grid2d,simulate_particles2d,simulate_particles2d_backwards

class SpinnyParticle(JITParticle):
    u = Variable('u',dtype=np.float64)
    v = Variable('v',dtype=np.float64)
    vort = Variable('vort',dtype=np.float64)
    
########## GET PARTICLE INIT COORDS #############
# Set up in `particle_map_of_coherence.ipynb`

### 1. Get NPSG particle grid
NPSG_poly = [[233, 33], [190, 33], [128, 24], [128, 9], [195, 9], [245, 21], [233,33]]
NPSG_poly_shapely = Polygon(NPSG_poly)
min_x, min_y, max_x, max_y = NPSG_poly_shapely.bounds

# Get grid of particles in the rectangle
res = 0.02 #0.01 #0.025
x_vals = np.arange(min_x, max_x, res)
y_vals = np.arange(min_y, max_y, res)
xx, yy = np.meshgrid(x_vals, y_vals)
grid_points = np.vstack((xx.ravel(), yy.ravel())).T

# Remove points outside the polygon
inside_NPSG_points = np.array([pt for pt in grid_points if NPSG_poly_shapely.contains(Point(pt))])
inside_NPSG_shapely_points = [Point(p) for p in inside_NPSG_points] # convert to shapely objects
print('%s NPSG points obtained!'%(len(inside_NPSG_shapely_points)))


### 2. Get NPSG Eddies & Edge (10km) buffer
# AVISO eddy data (already subsetted in `particle_map_of_coherence.ipynb`)
anti_ds = xr.open_dataset(AVISO_dir + 'META3.2_DT_allsat_Anticyclonic_long_20211123_NPSG_poly.nc')
cyc_ds = xr.open_dataset(AVISO_dir + 'META3.2_DT_allsat_Cyclonic_long_20211123_NPSG_poly.nc')

# Add 15km buffer around eddies to be consistent with the rest of the G4 analysis
def buffer_eddy_polygon(poly):
    """
    poly: a single eddy polygon in shapely format
    
    15km buffer around the eddy polygon using an Azimuthal Equidistant projection centered on the eddy
    """
    
    # Project the polygon from WGS84 to a Pacific coordinate projection
    center_lon, center_lat = poly.centroid.coords[0] # gets eddy center coords for the projection
    crs_pacific = CRS.from_proj4(f"+proj=aeqd +lat_0={center_lat} +lon_0={center_lon} +datum=WGS84 +units=m +no_defs")
    crs_wgs84 = CRS.from_epsg(4326) # EPSG code for the WGS84 coordinate system
    to_proj = Transformer.from_crs(crs_wgs84, crs_pacific, always_xy=True).transform
    poly_proj = transform(to_proj, poly)
        
    # Add a 15km buffer to the polgon using the new projection, then covert back to WGS84
    poly_buffered = poly_proj.buffer(15_000)  # 15km buffer (function expects meters)
    from_proj = Transformer.from_crs(crs_pacific, crs_wgs84, always_xy=True).transform
    poly_WGS84_wbuffer = transform(from_proj, poly_buffered)
    
    # Convert back to 360 again; necessary or polys across the dateline will be invalid 
    def shift_coords(coords):
        return [((x + 360 if x < 0 else x), y) for x, y in coords]

    exterior = shift_coords(poly_WGS84_wbuffer.exterior.coords)
    interiors = [shift_coords(ring.coords) for ring in poly_WGS84_wbuffer.interiors]
    final_poly = Polygon(exterior, interiors)
    return final_poly

eddy_polys = []
for eddy_ds in [anti_ds,cyc_ds]:
    for i in np.arange(0,len(eddy_ds.obs)):
        contour_lons = eddy_ds.effective_contour_longitude[i]
        contour_lats = eddy_ds.effective_contour_latitude[i]
        point_stack = np.column_stack((contour_lons.values.ravel(), contour_lats.values.ravel())).tolist()
        eddy_polys.append(Polygon(point_stack))
        
buff15km_eddy_polys = [buffer_eddy_polygon(poly) for poly in eddy_polys]
print('Eddy buffer computed...')


### 3. Remove particles in the eddies or edge
buffered_tree = STRtree(buff15km_eddy_polys) # trying STRtree to increase efficiency

def point_outside_all_buffers(point):
    possible_matches = buffered_tree.query(point) # Only check eddies close to the point
    return not any(poly.contains(point) for poly in possible_matches)

outside_eddy_wbuffer_points = [p for p in inside_NPSG_shapely_points if point_outside_all_buffers(p)]
lons = [p.x for p in outside_eddy_wbuffer_points]
lats = [p.y for p in outside_eddy_wbuffer_points]

print('%s in eddy particles removed!'%(len(inside_NPSG_shapely_points)-len(outside_eddy_wbuffer_points)))


########## RUN PARCELS #############

#time_step = 20 # minutes
#output_freq = 24 #hours
#runtime = 200 # days
#runtime_unit = 'days'
#start_date = datetime(2021,11,23)

# Create fieldset
#parcels_input_files = sorted(glob(data_dir + 'CMEMS_data/dt_global_allsat_phy_l4_*.nc'))

#filenames = {'U': parcels_input_files,'V': parcels_input_files}
#variables = {'U': 'ugos','V': 'vgos'}
#dimensions = {'U': {'lon':'longitude','lat':'latitude','time':'time'},
#              'V': {'lon':'longitude','lat':'latitude','time':'time'}}
#fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
#print('Fieldset created.')

#pset_dynamic = ParticleSet.from_list(fieldset=fieldset, pclass=SpinnyParticle, lon=lons, lat=lats, depth=[0]*len(lons), time=start_date)
#print('%s particles initialized.'%(len(lons)))

#output_dir = data_dir + 'parcels_trajs/'
#output_file_path = '%s%s_%s%s_runtime_%smin_timestep_%s_spatial_res_%shr_output_freq_outside_15km_eddy_buff_NPSG_domain.nc'%(
#output_dir,str(start_date)[0:10],runtime,runtime_unit,time_step,res,output_freq)

#print('Output file: %s'%(output_file_path))

# Execute particle simulation
#start_time = time.time()
#simulate_particles2d_backwards(pset_dynamic,output_file_path,runtime,runtime_unit,time_step,output_freq)
    
#print("--- %s seconds ---" % (time.time() - start_time))
