# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import numpy as np
import pandas as pd


def get_co_occurrence_df(adata, groupby):
    n_rows = adata.uns[f'{groupby}_co_occurrence']['occ'].shape[0]
    n_cols = adata.uns[f'{groupby}_co_occurrence']['occ'].shape[1]
    n_distances = adata.uns[f'{groupby}_co_occurrence']['occ'].shape[2]
    row_indices, col_indices = np.meshgrid(np.arange(n_rows), np.arange(n_cols), indexing='ij')
    row_indices = row_indices.reshape(-1, 1)  # (100, 1)
    col_indices = col_indices.reshape(-1, 1)  # (100, 1)
    arr_reshaped = adata.uns[f'{groupby}_co_occurrence']['occ'].reshape(n_rows * n_cols, n_distances)
    arr_final = np.hstack((row_indices, col_indices, arr_reshaped))
    co_occurrence_df = pd.DataFrame(arr_final, columns=['group1', 'group2'] + adata.uns[f'{groupby}_co_occurrence'][
                                                                                  'interval'].tolist()[:-1])
    co_occurrence_df[['group1', 'group2']] = co_occurrence_df[['group1', 'group2']].astype(int)
    categories = adata.obs['cell_type'].cat.categories
    co_occurrence_df['group1'] = co_occurrence_df['group1'].map(lambda x: categories[x])
    co_occurrence_df['group2'] = co_occurrence_df['group2'].map(lambda x: categories[x])
    return co_occurrence_df
