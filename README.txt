`G4_phyto_advection` contains the scripts that accompany Jones-Kellett et al. 2026. Derived data used in the analysis is distributed at  https://doi.org/10.5281/zenodo.18867070 . 

~~~ ENVIRONMENTS ~~~
OceanParcels Lagrangian simulations (v2.2.2; https://doi.org/10.5194/gmd-12-3571-2019) were ran and stored on a computer cluster (environment `py3_parcels_v2.yaml`).

Genetic data was processed on a local computer (environment `G4.yaml`).


~~~ /physics_pipeline/ SCRIPTS (environment `py3_parcels_v2.yaml`) ~~~

`functions_for_parcels.py` : functions used to accompany OceanParcels Lagrangian simulations

GAUSSIAN CLOUDS
1. `gaussian_lag_sims_v4.py` : OceanParcels simulations of the Gaussian clouds
2. `calc_coh_time_of_gaussian_clouds.py` : Measure dispersal/coherence time of the clouds
3. `water_mass_traj_plots.ipynb` : Plot the Gaussian cloud trajectories
	- Fig 3 : Lagrangian coherency of the sampled water masses with examples
	- Fig S7 : CMEMS vs OSCAR trajectories
4. `animate_trajs.py` : Convert saved images into an animation/video
	- Video S1 : Gaussian cloud slideshow
5. `northern_dist_of_clouds_v2.py` : Compute the origin distance of the clouds, data outputted for `temp_lat_advection.ipynb`

GRIDDED LAGRANGIAN SIMS
1. `run_parcels_CMEMS_v3.py` : OceanParcels gridded simulations
2. `lagrangian_calculations_v2.py` : Calculates FTLE, LAVD from gridded particle simulation
3. `run_parcels_CMEMS_NPSG_outside_eddies.py` : High-res gridded simulation to compute outside-eddy dispersal rates
4. `coh_time_of_gridded_particles_v3.py` : Calculating a pseudo-coherence time of particle clusters across the gyre (outside of eddies)
5. `coh_grid_test.ipynb` : Compare coherence time from grid sim computed in `coh_time_of_gridded_particles_v3.py` to original Gaussian clouds
	- Fig S11a,b

MESOSCALE EDDIES & MAP FIGS 
1. `G4_AVISO_eddies_v3.ipnyb` : Identify intersected eddies along cruise track in AVISO dataset
2. `G4_paper_map_figs_v2.ipnyb` : Map visualizations
	- Fig 1 : Regimes of phytoplankton biomass and community structure in the North Pacific
	- Fig 2a : SLA tile plot with sample locations and eddy bounds
	- Fig S11c,d : map of outside-eddy water mass coherence times across the NPSG
3. `G4_eddy_structure_comparison_v2.ipynb` : Plot SLA and LAVD of the individual eddies
	- Fig S10 : LAVD coherent structures of each sampled eddy


~~~ /phytoplankton_pipeline/ SCRIPTS (environment `G4.yaml`) ~~~

Pre-prossessing of sequence data and internal standard correction described in Jones-Kellett et al. 2024 (https://doi.org/10.1093/ismeco/ycae115)
Internal standard correction scripts: https://github.com/lexi-jones/internal_std_correction

1. `lat_diel_running_means_v3.ipynb`: Calculates the spatiotemporal anomaly (STA) from abundance for each phytoplankton ASV 
	- Fig 2 : Stacked eukaryote groups	
	- Fig S1 : latitudinal trends
	- Fig S2 : STA methods
	- Fig S3 : Stacked cyanobacteria groups
	- Fig S4 : spatial and temporal distances between samples
2. `bootstrap_STA.py` : Bootstrap resampling (10,000x) w/ replacement to derive confidence intervals for STA 
3. `dispersal_vs_eco_metrics_v10.ipynb` : Compares coherence timescale with STA 
	- Fig 4 : Euk STA vs coherence time 
	- Fig 5 : Modeled exponential decay of euk major groups
	- Fig S8 : Cyan STA vs coherence time
	- Fig S93 : Euk STA vs coherence time (group level)
4. `sig_avg_STA_ASVs.py` : Calculcate average ASV STA by niche type and determine significance from bootstrap data
5. `compare_ASV_niche_avgs.ipynb` : Analyze average ASV STA by niche type
	- Fig 6 : Average ASV STA ranked

~~~ /SST_and_TKE/ SCRIPTS (environment `G4.yaml`) ~~~

Supplemental calculations

TEMPERATURE
1. `temp_lat_advection.ipynb` : Compare Gaussian cloud trajectory directions with temperature anomalies
	- Fig S5 : Underway SST anomalies versus water mass advection histories (Geostrophy only, Geostrophy + Ekman)
	- Fig S6 : steps to derive the spatiotemporal anomaly of SST
TKE
1. `Turbulent_kinetic_energy_generation.jl` : Calculate TKE; note that this Julia script was run outside the `G4.ymal` environment by MJF
2. `TKE_plot.ipynb` : Plot TKE
	- Fig S12 : TKE and MLD 2020-2022
