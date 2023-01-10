# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import numpy as np
import pandas as pd


def sampleinfo_stat(sample_info, cols):
    stat_df = []
    for col in cols:
        if sample_info[col].dtype == 'object':
            for key, value in sample_info[col].value_counts().to_dict().items():
                stat_df.append([col, key, value])
        elif sample_info[col].dtype == np.int or sample_info[col].dtype == np.float32 or sample_info[
            col].dtype == np.float64:
            stat_df.extend([[col, 'mean', sample_info[col].mean()], [col, 'median', sample_info[col].median()],
                            [col, 'std', sample_info[col].std()]])
    return pd.DataFrame(stat_df)