# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import pandas as pd


# %%
def trtmb_calc(aachange_file):
    """Truncal-bTMB calculation
    Reference:
    Li, S., Noor, Z.S., Zeng, W. et al. Sensitive detection of tumor mutations from blood and its application to immunotherapy prognosis.
    Nat Commun 12, 4172 (2021). https://doi.org/10.1038/s41467-021-24457-2

    """
    mut_df = pd.read_csv(aachange_file, sep='\t')
    mut_df.sort_values(by='VAF', inplace=True)
    vaf_cutoff = mut_df['VAF'][:5].mean() * 0.6
    truncal_mut_df = mut_df[mut_df['VAF'] >= vaf_cutoff]
    trtmb = truncal_mut_df['VAF'].sum() / vaf_cutoff
    return trtmb
