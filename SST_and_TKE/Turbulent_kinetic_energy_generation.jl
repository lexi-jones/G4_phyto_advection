# ...........................................................
#= 
   Mick Follows, 2025

  Evaluate TKE generation rate for NE Pacific

Flow:

(1) Read 1. ERA5 reanalysis surface heat fluxes and wind-stress
          (DAILY)
        2. ORAS5 MLD product (Regridded to ERA5 grid) 
          (MONTHLY)

          ERA5 is ECMWF climate reanalysis data. 
          Hourly means for midnight UTC over the region of interest were downloaded from
          Copernicus Climate Data Store (https://cds.climate.copernicus.eu/#!/home). 
          Input data fields are sensible heat flux, latent heat flux, net shortwave radiation,
          net longwave radiation, surface wind-stress (x- and y-).

          ORAS5 MLD is a related mixed layer depth product. Monthly mean reanalysis 
          obtained from Copernicus on original model grid (re-gridded to compatible grid). 

(2)  Evaluate TKE generation rate following Follows and Dutkiewicz (Deep-Sea Res, 2002),
     after Krauss and Turner (1967), Niiler and Krauss (1988) and other literature.

=#

# ...........................................................

using DataFrames
using NCDatasets
using DataStructures
using Dates
using CFTime
using GLMakie
 
# define some factors
secs_per_day = 86400.0 # seconds per day
secs_per_hour = 3600.0 # seconds per day
alpha = 0.15           # (kg m-3 oC-1) coeff for thermal contribution to 
                       #               equation of state for seawater
g = 9.80665            # (m s-2) acceleration due to gravity
rho = 1025.0           # (kg m-3) density of seawater (N.B fixed value)
Cp = 4000.0            # J kg-1 K-1 (varies 4000 at 0C to 4011 at 30C)
 
# Read the ERA5 surface fluxes etc................................
# 2020-2021 daily
# open the file
ds = Dataset("/INPUT_DIRECTORY_NAME/ERA5_INPUT_DATA_FILE.nc")
#
# read main surface flux data file and organize the data....
# latitude and longitude
lats = ds["latitude"][:];
lat = sort!(lats);  # lat array in ascending order for plotting
#lat_units = ds["latitude"].attrib["units"]
lons = ds["longitude"][:];
#lon_units = ds["longitude"].attrib["units"]
# date / time
# seconds since 1970-01-01
date = ds["valid_time"][:]
time_units = ds["valid_time"].attrib["units"]
# x,y,time fields
# surface upward sensible heat flux (J m-2)
sens_heat_flux_a = Array(ds["sshf"]) ;
sens_heat_flux_units = ds["sshf"].attrib["units"]
# surface upward latent heat flux (J m-2)
lat_heat_flux_a = Array(ds["slhf"]) ;
lat_heat_flux_units = ds["slhf"].attrib["units"]
# surface net downward shortwave radiation (J m-2)
short_wave_rad_a = Array(ds["ssr"]) ;
short_wave_rad_units = ds["ssr"].attrib["units"]
# surface net upward thermal radiation (J m-2)
long_wave_rad_a = Array(ds["str"]) ;
long_wave_rad_units = ds["str"].attrib["units"]
# surface downward eastward stress (N m-2 s)
tau_E_a = Array(ds["ewss"]) ;
tau_units = ds["ewss"].attrib["units"]
# surface downward northward stress (N m-2 s)
tau_N_a = Array(ds["nsss"]) ;
tau_units = ds["nsss"].attrib["units"]
# close the file
close(ds)
println("finished reading flux file ")           
# .................................................................

# Read land mask ..................................................
# 2D file!!
ds2 = Dataset("/INPUT_DIRECTORY_NAME/LAND_MASK_FILE.nc")
lsmask_a = ds2["lsm"];    # land surface mask (1.0, 0.0)
mask_a = 1.0 .- lsmask_a;   # mask out land, keep ocean
mask_a = floor.(mask_a) ;   # round down to 0.0 if < 1.0
close(ds2)                  # close the file
println("finished reading mask file ")           
# .................................................................

# Read MLD data ...................................................
# ORAS5 reanalysis MLD from Copernicus 
# 1990-2022 monthly means
# regridded by Jon Lauderdale to ERA5 standard grid
ds3 = Dataset("/INPUT_DIRECTORY_NAME/ORAS5_MLD_DATA_FILE_MONTHLYnc")
# lat, long not provided: same grid as other input data
MLD_b = Array(ds3["somxl010"]) ;    # MLD, criteria 10 (define?)
#MLD_b = Array(ds3["somxl030"]) ;   # MLD, criteria 30 (define?)
MLD_c = MLD_b ;
MLD_a = MLD_b ;
MLD_units = "m"
MLD_date = ds3["time"][:] ;     # time in MLD data set
#MLD_date = MLD_date_a ;
close(ds3) # close the MLD file
println("finished reading MLD file ")           
# .................................................................


# Read the global lat-long data for Mixed Layder Depth............
#   ... IF NECESSARY
ds4 = Dataset("/INPUT_DIRECTORY_NAME/LAT_LONG_DATA.nc")
# latitude
sMLD_lat = ds4["latitude"][:];
MLD_lat = sort!(sMLD_lat);  # lat array in ascending order for plotting
# longitude
MLD_lon = ds4["longitude"][:];
close(ds4)
println("read global MLD lat/long info")
# .................................................................


# reverse flux arrays in latitude dimension for plotting   
# (not necessary for MLD in my case )  
mask = reverse(mask_a, dims = 2) ;
sens_heat_flux = reverse(sens_heat_flux_a, dims = 2) ;
lat_heat_flux = reverse(lat_heat_flux_a, dims = 2) ;
short_wave_rad = reverse(short_wave_rad_a, dims = 2) ;
long_wave_rad = reverse(long_wave_rad_a, dims = 2) ;
tau_N = reverse(tau_N_a, dims = 2) ;
tau_E = reverse(tau_E_a, dims = 2) ;
#MLD_a = reverse(MLD_c, dims = 2) ;
println("finished reversing arrays...." )

# mask - excludes land : maskr - exludes ocean
maskr = 1.0 .- mask                            # make a reverse mask with NaNs
maskNaN = mask .* NaN                          # make NaN array of same size
mask[mask .< 0.9] .= maskNaN[mask .< 0.9]      # replace 0s in mask with NaNs
maskr[maskr .< 0.1] .= maskNaN[maskr .< 0.1]   # replace 0s in maskr with NaNs
println("finished mask operations...." )

# ......................................................................
# NOTE ON UNITS. From Copernicus ERA5 download page (latent heat flux entry)
# "... transfer of latent heat ... accumulated over a particular time period
# which depends on the data extracted. For the monthly averaged reanalysis and
# the monthly averaged ensemble members, the accumulation period is 1 day. For
# the monthly averaged reanalysis by hour of day, the accumulation period is 1
# hour and for the monthly averaged ensemble members by hour of day, the
# accumulation period is 3 hours. The units are J m-2. To convert to W m-2 the
# accumulated values should be divided by the accumulation period expressed in
# seconds. The ECMWF convention for vertical fluxes is positive downwards. "
# ......................................................................

# Evaluate net surface heat flux (J m-2)
# Note scaling by time for appropriate units
net_surf_heat_flux_A = sens_heat_flux .+ lat_heat_flux .+ 
			short_wave_rad .+ long_wave_rad ;
H_A = net_surf_heat_flux_A ./ secs_per_hour ; # scale units
net_surf_heat_flux_Wm2_units = "W m-2" ;

# Evaluate mod(wind-stress)
# Note scaling by time for appropriate units
mod_tau_b = ( tau_E.^2.0 + tau_N.^2.0).^0.5 ;
mod_tau = mod_tau_b ./ secs_per_hour ;
wind_stress_units = "N m-2" ;

# ...................................................................
# some details for cutting out region from global MLD field 
# really ugly code, but works...
# longitude is specified differently degrees W vs degrees E
# so need to figure out bounds of sub-array
        # longitude
        rownumber = collect(1:1440)
        aalon = [0.25*j for j in -1440:-1] ;  # "-ve" longitude going west
        y = [rownumber, MLD_lon, aalon] ;     # array of lon & -lon
        df = DataFrame(y, [:rownum, :mldlon, :alon]) ;  # dataframe
        df2 = df[df.alon .== -150.0, :]      # find western boundary row
        lonmin = Int(df2[1,1])                 # min index as integer 
        df2 = df[df.alon .== -110.0, :]            # find eastern boundary row
        lonmax = Int(df2[1,1])                 # max index as integer 
        # latitude
        rownumber = collect(1:721)
        aalat = [0.25*j for j in -360:360] ;  # latitude
        y = [rownumber, MLD_lat, aalat] ;     # array of lat
        dfa = DataFrame(y, [:rownum, :mldlat, :alat]) ;  # dataframe
        dfa2 = dfa[dfa.alat .== 10.0, :] ;    # southern boundary row
        latmin = Int(dfa2[1,1])          ;    # min index as integer 
        dfa2 = dfa[dfa.alat .== 40.0, :] ;    # find eastern boundary row
        latmax = Int(dfa2[1,1])               # max index as integer 
        println("cut out relevant area from global MLD")
        println("latmin = ",latmin," latmax = ",latmax," lonmin = ",lonmin,"   lonmax = ", lonmax) 
# ...................................................................

# select time slice 
itime_max = length(date)
MLD_itime_max = length(MLD_date)
println("itime_max = ", itime_max)           
println("MLD_itime_max = ", MLD_itime_max)           

# for every time slice in the sequence
    for itime in 1:itime_max 
    #for itime in 1:70
        # for each flux date in the sequence ... 
	dstring = Date(date[itime])          # current date for fluxes
	date_a = DateTime(date[itime])       # current date for fluxes
        # daily fluxes, monthly MLDs: Determine which MLD slice to use
        flux_year = Dates.year(date[itime]) ;
        flux_month = Dates.month(date[itime]) ;
        for itime_test in 1:MLD_itime_max 
           MLD_year = Dates.year(MLD_date[itime_test]) ;
           MLD_month = Dates.month(MLD_date[itime_test]) ;
           if MLD_year == flux_year && MLD_month == flux_month 
              MLD_itime = itime_test
	      global MLD_itime
           end 
        end
        # check matching time slices...
        MLD_year = Dates.year(MLD_date[MLD_itime]) ;
        MLD_month = Dates.month(MLD_date[MLD_itime]) ;
        println("flux_year, MLD_year ", flux_year,",", MLD_year,"     flux_month, MLD_month ", flux_month,",", MLD_month)        

# create single time-slice arrays for mask, net heat flux and wind-stress
	mask_t = mask[:,:]             ;    # masks
	maskr_t = maskr[:,:]           ;    # masks
	H_t = H_A[:,:,itime] .* mask[:,:] ; # net surface heat flux (J m-2)
        H_t[ismissing.(H_t)] .= NaN ;       # set missing to NaN
	mod_tau_t = mod_tau[:,:,itime] .* mask[:,:]; # wind stress

# NOTE: MLD is in global, monthly format ==================================
#       cut out appropriate local region (see above)
# mixed layer depth 
	MLD_d = MLD_a[:,:,MLD_itime]   ;    # MLD, approprite time slice
	MLD_t = MLD_d[:,:]' ;               # transpose grid for compatibility
        MLD_t[ismissing.(MLD_t)] .= NaN ; # set missing to NaN
	MLD = reverse(MLD_t, dims = 2) ;  # reverse lat axis for compatibility
# make new MLD array for NE Pacific region
        MLD_NEPac = MLD[lonmin:lonmax, latmin:latmax] ; 
      
# =========================================================================

	# evaluate TKE generation rate ....................................
	# Krauss and Turner model parameters m1 and m2 
	m1 = 1.25 ;             
	m2 = H_t .* 0.0 ;        # set up an array for m2
	m2[H_t .<= 0.0] .= 1.0;   # ocean losing heat
	m2[H_t .>= 0.0] .= 0.2;   # ocean gaining heat
	# eval TKE generation as Follows and Dutkiewicz (DSR, 2002), Eq. (1)
	# wind stress term
	fvsq = mod_tau_t ./ rho  ;          
	fv = fvsq.^0.5  ;                  
	TKEgen1_a = (m1 .* (fv.^3)) ;        
	# heat flux term: if H_t < 0, ocean losing heat and gaining TKE
        #      hence -1.0 factor
	TKEgen2_a = -1.0 .* H_t .* ((m2 .* g .*alpha .* MLD_NEPac) ./ (2.0 .* (rho^2) .* Cp)) ; 
	# total TKE generation rate
	TKEgen_sum_a = TKEgen1_a .+ TKEgen2_a ;
        # ......................................................................

        # Add time slice of data to the TKEgen NetCDF file

	if itime == 1
        # use first chunk to start the arrays
                TKEgen_wind = TKEgen1_a ;
                TKEgen_heat = TKEgen2_a ;
                TKEgen = TKEgen_sum_a ;
                MLDout = MLD_NEPac ;
                datex = date_a ;
                # make arrays and units available outside loop
                global TKEgen_wind
                global TKEgen_heat
                global TKEgen
                global MLDout
                global datex
                global lats
                global lons
                global TKEgen_units
        else
        # merge next chunk to main array
        # concatenate surface fluxes along time dimension
                TKEgen_wind = cat(TKEgen_wind, TKEgen1_a,dims=3) ;
                TKEgen_heat = cat(TKEgen_heat, TKEgen2_a,dims=3) ;
                TKEgen = cat(TKEgen, TKEgen_sum_a,dims=3) ;
                MLDout = cat(MLDout, MLD_NEPac,dims=3) ;
                datex = cat(datex, date_a, dims=1) ;
        end

   end
# finish stepping through monthly sequence

# write out TKEgen climatology to NetCDF 
flux_directory = "OUTPUT_DIRECTORY_NAME"
filename = "OUTPUT_FILE_NAME.nc"
data_set_name = flux_directory*filename
# Create NetCDF file. "c" = clobber
ds_out = Dataset(data_set_name,"c")

# determine dimensions of arrays
lon_size = size(TKEgen, 1)
lat_size = size(TKEgen, 2)
time_size = size(TKEgen, 3)
# set the dimesions for the file
defDim(ds_out,"lon",lon_size)
defDim(ds_out,"lat",lat_size)
defDim(ds_out,"times",time_size)

# set "missing" values to NaN for output to NetCDF ...
# NetCDF will take NaN but not "missing" (?) 
TKEgen[ismissing.(TKEgen)] .= NaN ;
TKEgen_wind[ismissing.(TKEgen_wind)] .= NaN ;
TKEgen_heat[ismissing.(TKEgen_heat)] .= NaN ;

# Define a global attribute
ds_out.attrib["title"] = "TKE generation rate"

# define and fill variables
#TKEgen
v1 = defVar(ds_out,"TKEgen",Float32,("lon","lat","times"))
TKEgen_units = "cm3 s-3"
v1.attrib["units"] = TKEgen_units            # add attributes
v1[:,:,:] = TKEgen * 1.0e6                   # fill the variable
                                             # *1.0e6 to scale units ...??
#TKEgen_wind
v2 = defVar(ds_out,"TKEgen_wind",Float32,("lon","lat","times"))
v2.attrib["units"] = TKEgen_units            # add attributes
v2[:,:,:] = TKEgen_wind  * 1.0e6             # fill the variable
                                             # *1.0e6 to scale units ...??

#TKEgen_heat
v3 = defVar(ds_out,"TKEgen_heat",Float32,("lon","lat","times"))
v3.attrib["units"] = TKEgen_units            # add attributes
v3[:,:,:] = TKEgen_heat * 1.0e6              # fill the variable
                                             # *1.0e6 to scale units ...CHECK
#MLD 
v7 = defVar(ds_out,"MLD",Float32,("lon","lat","times"))
v7.attrib["units"] = MLD_units               # add attributes
v7[:,:,:] = MLDout                           # fill the variable

# longitude -  need extra comma after "lon", ... why?
v4 = defVar(ds_out,"longitude",Float32,("lon",))
v4[:] = lons                                 # fill the variable
# latitude
v5 = defVar(ds_out,"latitude",Float32,("lat",))
v5[:] = lats                                 # fill the variable

# date / time
# note... made datex a "DateTime" variable above (not "Date")... ??
# in next line ,datex, is communicating format of variable ??
# ... anyway, seems to work...
v6 = defVar(ds_out,"date",datex,("times",))
v6[:] = datex                      # fill the variable

# close the NetCDF file
close(ds_out)



