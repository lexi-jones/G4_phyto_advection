# Bootstrap STA 
#
# Bootstrap resampling (10,000x) w/ replacement the replicate data for each ASV & calculating an STA.
# Output: CSV for each ASV of 10,000x potential STAs at each sample site.
#
# LJK
# Date created : 12/22/25
# Last edited : 03/03/26

import json
from glob import glob
import numpy as np
import pandas as pd
from collections import OrderedDict
from config import * # directory file paths

#################################### DATA PREP #######################################

# Names of the technical replicates (all have 3 except for sample '20')
ordered_dict = data_dir + 'station_name_dict_NPSG_analysis.json'
with open(ordered_dict, 'r') as f:
    station_name_dict = json.load(f, object_pairs_hook=OrderedDict)
    
samples_in_analysis = []
for key,values in station_name_dict.items():
    samples_in_analysis.extend(values)
    
# Get the IDs of the 406 phytolankton ASVs included in the study
phyto_ASV_lat_temp_anom_sorted_df = pd.read_csv(highcov_dir + 'phyto_ASV_spatiotemporal_anoms_v3.csv', index_col=0)
ASVs = phyto_ASV_lat_temp_anom_sorted_df.columns

# Load count data from each internal standard estimate (BP, DR, TT)
def filter_replicate_files(ASV_file):
    ASV_df = pd.read_csv(highcov_dir + ASV_file,index_col=0).transpose().drop('taxonomy')
    ASV_df = ASV_df.loc[samples_in_analysis] # remove samples not in analysis
    ASV_df = ASV_df[ASVs] # remove ASVs not in analysis
    return ASV_df

BP_ASV_file = '221118-1309_LexiGradients-HighCov_2.09-fold-18S-correction_normalized_sequence_counts_abs_ASV_abundance_IS_BP.tsv'
DR_ASV_file = '221118-1309_LexiGradients-HighCov_2.09-fold-18S-correction_normalized_sequence_counts_abs_ASV_abundance_IS_DR.tsv'
TT_ASV_file = '221118-1309_LexiGradients-HighCov_2.09-fold-18S-correction_normalized_sequence_counts_abs_ASV_abundance_IS_TT.tsv'

BP_ASV_df = filter_replicate_files(BP_ASV_file)
DR_ASV_df = filter_replicate_files(DR_ASV_file)
TT_ASV_df = filter_replicate_files(TT_ASV_file)

# Load metadata for STA calculation
metadata_df = pd.read_csv(data_dir + 'AVISO_metadata_15km_near_eddy_w_sunrise_data.csv',index_col=0)
metadata_sorted_df = metadata_df.sort_values(by='time_since_sunrise')

#####################################################################################

def rolling_mean_looping_weighted(df, x_values, window):
    """
    Compute a weighted rolling mean with (looping) periodic boundary conditions,
    weighted by the distance between x-values. The window is centered on the center value.

    Params
    - df: pandas DataFrame to apply the rolling mean
    - x_values: x-values associated with the rows of df (e.g., time since sunrise)
    - window: window size (e.g. window=7 will be 3 on either side of center point)

    Returns
    - A df with the weighted rolling mean applied
    """
    # initialize a dataframe 
    result = pd.DataFrame(index=df.index, columns=df.columns, dtype=float)

    # Iterate through each ASV
    for col in df.columns:
        rolling_mean_column = []

        # Iterate through the samples to compute looping rolling mean
        for i in range(len(df)):

            # Define current the window
            inds = [(i + j - window // 2) % len(df) for j in range(window)] # window indeces with looping
            window_vals = df.iloc[inds, df.columns.get_loc(col)].values
            x_window = x_values.iloc[inds].values # window x's
            x_current = x_values.iloc[i] # center x
            distances = np.abs(x_window - x_current) # distance between the current point and other vals in the window
            normalized_distances = distances/np.max(distances)
            
            # Exponential decay for weights based on distance
            weights = np.exp(-normalized_distances)
            normalized_weights = weights/weights.sum() 

            # Compute the weighted mean
            rolling_mean_column.append(np.sum(window_vals*normalized_weights))

        result[col] = rolling_mean_column

    return result


boot_array = np.arange(0,10000) # number of times to resample

# Iterate through each ASV
count = 0
for a in ASVs:
    if count % 10 == 0:
        print(count)
    
    #################### Bootstrap Abundance ####################
    
    # Start an empty pandas dataframe to store data from each bootstrap iteration
    boot_abund_df = pd.DataFrame(index=list(station_name_dict.keys()), columns=boot_array)
    
    # Iterate through each sample to get the replicates 
    for key, values in station_name_dict.items():
        sample_reps = [] # measured replicates 
        for v in values: 
            sample_reps.append(BP_ASV_df.loc[v,a])
            sample_reps.append(DR_ASV_df.loc[v,a])
            sample_reps.append(TT_ASV_df.loc[v,a])
        sample_reps = np.array(sample_reps)
        sample_reps = sample_reps[~np.isnan(sample_reps)] # remove any nan values 

        # Resample number of times in boot array
        for n in boot_array:
        
            # Take mean of randomly resampled replicates with replacement
            mean_draw_w_replace = np.mean(np.random.choice(sample_reps,size=len(sample_reps),replace=True))
            boot_abund_df.loc[key,n] = mean_draw_w_replace
            
    #################### Calculate STA (see also `lat_diel_running_means_v3.ipynb`) ####################
    # 1. Rolling mean
    boot_abund_df = boot_abund_df+1 # add buffer to avoid error with 0 abundance
    boot_lat_rolling_df = boot_abund_df.rolling(window=9,center=True,min_periods=4).mean() # lat rolling mean

    # 2. Normalized latitudinal anomaly
    boot_lat_anom_df = np.subtract(boot_abund_df,boot_lat_rolling_df) # subtract lat rolling mean
    boot_lat_norm_anom_df = boot_lat_anom_df / boot_lat_rolling_df # normalize by lat rolling mean

    # 3. Temporal anomaly
    boot_lat_norm_anom_sorted = boot_lat_norm_anom_df.reindex(metadata_sorted_df.index)
    boot_temp_rolling_df = rolling_mean_looping_weighted(boot_lat_norm_anom_sorted, metadata_sorted_df['time_to_sunrise'],7)
    boot_lat_temp_anom_df = np.subtract(boot_lat_norm_anom_sorted,boot_temp_rolling_df)

    # 4. Sort back to latitudal ordering and rescale by abundances from lat rolling mean 
    boot_lat_temp_anom_sorted_df = np.multiply(boot_lat_temp_anom_df.reindex(metadata_df.index), (boot_lat_rolling_df))
    
    # 5. Save output (one file per ASV)
    boot_lat_temp_anom_sorted_df.to_csv(highcov_dir + '/bootstrapped_STAs/%s_bootstrapped_STA.csv'%(a)) 
    
    count += 1
