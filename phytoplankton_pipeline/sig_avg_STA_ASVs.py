# Significant Average STA Eukaryotic Phytoplankton ASVs by Niche Type
# 
# LJK
#
# Date created: 02/27/26
# Last edited: 05/05/26
#

import numpy as np
import pandas as pd
from config import * # directory file paths

# Load metadata
metadata_df = pd.read_csv(data_dir + 'AVISO_metadata_15km_near_eddy_w_coh_time.csv',index_col=0)

# Sample niche categories
eddy_samples = list(metadata_df[(metadata_df.eddy_indicator != 'O')].index)
outside_mixed_samples = list(metadata_df[(metadata_df.eddy_indicator == 'O') & (metadata_df.mean_coh < 90)].index)
outside_coh_samples = list(metadata_df[(metadata_df.eddy_indicator == 'O') & (metadata_df.mean_coh >= 90)].index)

# Load euk STAs 
euk_ASV_STA = pd.read_csv(highcov_dir + 'euk_phyto_ASV_spatiotemporal_anoms_v3.csv',index_col=0)

# Init a data table to hold significance of avg STAs for each euk ASV within each niche;
# Significance defined as 95% CI having same sign (positive or negative)
sig_ASV_STA = pd.DataFrame(index=euk_ASV_STA.columns,columns=['Eddy','Outside Mixed','Outside Coherent'])
CI_ASV_STA = pd.DataFrame(index=euk_ASV_STA.columns,columns=['Eddy','Outside Mixed','Outside Coherent'])

def boot_avg_STA(sample_set, ASV_boot_df):
    """
    Check for 95% CI significance 
    """
    
    sample_set_df = ASV_boot_df[ASV_boot_df['sample_id'].isin(sample_set)].iloc[:,1:] # crop out first column of sample names
    boot_means = sample_set_df.mean(axis=0) # get mean STA for each boot iteration
    boot_mean_lower = np.quantile(boot_means,0.025)
    boot_mean_upper = np.quantile(boot_means,0.975)

    if np.sign(boot_mean_lower) == np.sign(boot_mean_upper): # significant
        ASV_sig = 1
        ASV_CI = '[%s, %s]$^*$'%(round(boot_mean_lower,1),round(boot_mean_upper,1))
    else: # not significant
        ASV_sig = 0
        ASV_CI = '[%s, %s]'%(round(boot_mean_lower,1),round(boot_mean_upper,1))
    return ASV_sig, ASV_CI

# Iterate through each euk ASV
for ASV in list(sig_ASV_STA.index):
    print(ASV)
    # Open bootstrap data for the specific ASV 
    ASV_boot_df = pd.read_csv(highcov_dir + 'bootstrapped_STAs/%s_bootstrapped_STA.csv'%(ASV))
    
    # Niche 1: Eddy
    eddy_ASV_sig, eddy_ASV_CI = boot_avg_STA(eddy_samples, ASV_boot_df)
    sig_ASV_STA.loc[ASV]['Eddy'] = eddy_ASV_sig
    CI_ASV_STA.loc[ASV]['Eddy'] = eddy_ASV_CI

    # Niche 2: Outside mixed
    mixed_ASV_sig, mixed_ASV_CI = boot_avg_STA(outside_mixed_samples, ASV_boot_df)
    sig_ASV_STA.loc[ASV]['Outside Mixed'] = mixed_ASV_sig
    CI_ASV_STA.loc[ASV]['Outside Mixed'] = mixed_ASV_CI
    
    # Niche 3: Outside coherent
    coh_ASV_sig, coh_ASV_CI = boot_avg_STA(outside_coh_samples, ASV_boot_df)
    sig_ASV_STA.loc[ASV]['Outside Coherent'] = coh_ASV_sig
    CI_ASV_STA.loc[ASV]['Outside Coherent'] = coh_ASV_CI

#sig_ASV_STA.to_csv(highcov_dir + 'sig_avg_STA_ASVs_by_niche.csv')
CI_ASV_STA.to_csv(highcov_dir + 'CI_avg_STA_ASVs_by_niche.csv')




