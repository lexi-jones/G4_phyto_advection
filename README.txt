`G4_phyto_advection` contains the scripts that accompany Jones-Kellett et al. 202X. 

~~~ ENVIRONMENTS ~~~
OceanParcels Lagrangian simulations (v2.2.2; https://doi.org/10.5194/gmd-12-3571-2019) were ran and stored on a computer cluster (environment `py3_parcels_v2.ymal`).

Genetic data was processed on a local computer (environment `G4.ymal`).


~~~ SCRIPTS (environment `py3_parcels_v2.ymal`) ~~~

`functions_for_parcels.py` : functions I use to accompany OceanParcels Lagrangian simulations

GAUSSIAN CLOUDS
1. `gaussian_lag_sims_v4.py` : OceanParcels simulations
2. `calc_coh_time_of_gaussian_clouds.py` : Measure dispersal of the clouds
3. `water_mass_traj_plots.ipynb` : Plot the Gaussian cloud trajectories
4. `animate_trajs.py` : Convert saved images into an animation/video

GRIDDED LAGRANGIAN SIMS
1. `run_parcels_CMEMS_v3.py` : OceanParcels gridded simulations
2. `lagrangian_calculations_v2.py` : Calculates FTLE, LAVD from gridded particle simulation
3. `run_parcels_CMEMS_NPSG_outside_eddies.py` : High-res gridded simulation to compute outside-eddy dispersal rates
4. `coh_time_of_gridded_particles_v3.py` : Calculating a pseudo-coherence time of particle clusters across the gyre (outside of eddies)
5. `coh_grid_test.ipynb` : Compare coherence time from grid sim computed in `coh_time_of_gridded_particles_v3.py` to original Gaussian clouds

MESOSCALE EDDIES & MAP FIGS
1. `G4_AVISO_eddies_v3.ipnyb` : Identify intersected eddies along cruise track (NOTE: NEED TO CLEAN UP, MOVE MAP FROM G4_paper_map_figs.ipynb TO HERE)
2. `G4_paper_map_figs_v2.ipnyb` : Map visualizations
3. `G4_eddy_structure_comparison_v2.ipynb` : Plot SLA, FTLE, LAVD of the individual eddies



~~~ SCRIPTS (environment `G4.ymal`) ~~~

PHYTOPLANKTON PIPELINE
1. `pp_group_annotation_v3.ipynb` : Annotate the phytoplankton taxa by group
2. `lat_diel_running_means_v3.ipynb`: Calculates the spatiotemporal anomaly (STA) from abundance for each phytoplankton ASV 
	- Fig S1 : latitudinal trends
	- Fig S4 : spatial and temporal distances between samples
	- Fig S2 : STA methods 
	- Fig 2 : Stacked eukaryote groups
	- Fig S3 : Stacked cyanobacteria groups 
3. `dispersal_vs_eco_metrics_v9.ipynb` : Compares coherence timescale with STA 
	- Fig 5 : Euk STA vs coherence time (total sum)
	- Fig S7 : Cyan STA vs coherence time
	- Fig S8 : Euk STA vs coherence time (group level)
4. `pcc_dist_by_physical_state_v2.ipynb` : Spearman Distance

TEMPERATURE
1. `temp_lat_advection.ipynb` : Compare Gaussian cloud trajectory directions with temperature anomalies

