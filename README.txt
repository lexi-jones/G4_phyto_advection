GAUSSIAN CLOUD PIPELINE
1. `` : OceanParcels simulations

GRIDDED LAGRANGIAN SIMS & EDDY ID PIPELINE
1. `run_parcels_CMEMS_v3.py` : OceanParcels gridded simulations
2. `gaussian_lag_sims_v4.py` : Calculates FTLE, LAVD from gridded particle simulation
3. `G4_AVISO_eddies_v3.ipnyb` : Identify intersected eddies along cruise track (NOTE: NEED TO CLEAN UP, MOVE MAP FROM G4_paper_map_figs.ipynb TO HERE)

BIOLOGY PIPELINE

1. `lat_diel_running_means_v3.ipynb`: Calculates the spatiotemporal anomaly (STA) from abundance for each phytoplankton ASV 
	- Fig S1 : latitudinal trends
	- Fig S4 : spatial and temporal distances between samples
	- Fig S2 : STA methods 
	- Fig 2 : Stacked eukaryote groups
	- Fig S3 : Stacked cyanobacteria groups 

RESULTS 

1. `dispersal_vs_eco_metrics_v9.ipynb` : Compares coherence timescale with STA 
	- Fig 5 : Euk STA vs coherence time (total sum)
	- Fig S7 : Cyan STA vs coherence time
	- Fig S8 : Euk STA vs coherence time (group level)
2. `pcc_dist_by_physical_state_v2.ipynb` : Spearman Distance

3. `G4_eddy_structure_comparison_v2.ipynb`