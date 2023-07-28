# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import numpy as np
import pandas as pd


def sampleinfo_stat(sample_info, cols=None):
    """统计样本信息"""
    if not cols:
        cols = sample_info.columns
    stat_df = []
    for col in cols:
        stat_df.append([col, ''])
        if sample_info[col].dtype == 'object':
            for key, value in sample_info[col].value_counts().to_dict().items():
                stat_df.append([key, value])
        elif sample_info[col].dtype == np.int or sample_info[col].dtype == np.float32 or sample_info[
            col].dtype == np.float64:
            stat_df.extend([['mean', sample_info[col].mean()], ['quantile25', sample_info[col].quantile(0.25)],
                            ['median', sample_info[col].median()], ['quantile75', sample_info[col].quantile(0.75)],
                            ['std', sample_info[col].std()]])
    return pd.DataFrame(stat_df)
