# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import numpy as np
import pandas as pd


def sampleinfo_stat(sample_info, save=True, cols=None):
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
            stat_df.extend([['interval', f'{sample_info[col].min()}~{sample_info[col].max()}'],
                            ['mean', sample_info[col].mean()], ['median', sample_info[col].median()]])

    return pd.DataFrame(stat_df)
