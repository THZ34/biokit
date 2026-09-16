# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.neighbors import radius_neighbors_graph


def calculate_radius_neighbor_counts(adata, radius=0.05, cell_type_col='cell_type', sample_col='sample',
                                     coord_cols=('x_slide_mm', 'y_slide_mm'), cell_types=None, n_jobs=-1,
                                     store_key=None, ):
    required_columns = [cell_type_col, sample_col, *coord_cols]
    missing_columns = [column for column in required_columns if column not in adata.obs.columns]
    if missing_columns:
        raise KeyError(f'adata.obs 缺少必要字段：{missing_columns}')
    if not adata.obs_names.is_unique:
        raise ValueError('adata.obs_names 不唯一，无法保证邻居计数结果与细胞一一对应。')
    if adata.obs[cell_type_col].isna().any():
        raise ValueError(f'{cell_type_col} 存在缺失值，不能计算邻居细胞类型计数。')
    if adata.obs[sample_col].isna().any():
        raise ValueError(f'{sample_col} 存在缺失值，不能按样本分别计算邻域。')
    if adata.obs[list(coord_cols)].isna().any().any():
        raise ValueError(f'空间坐标存在缺失值：{adata.obs[list(coord_cols)].isna().sum().to_dict()}')
    cell_type_values = adata.obs[cell_type_col].astype(str)
    sample_values = adata.obs[sample_col].astype(str)
    if cell_types is None:
        cell_types = sorted(cell_type_values.unique())
    else:
        cell_types = list(cell_types)
        unexpected_types = sorted(set(cell_type_values.unique()) - set(cell_types))
        if unexpected_types:
            raise ValueError(f'cell_types 未包含以下实际细胞类型：{unexpected_types}')
    if len(cell_types) == 0:
        raise ValueError('未检测到可用于计算的细胞类型。')
    radius_um = int(round(radius * 1000))
    if store_key is None:
        store_key = f'neighbor_count_{radius_um}um'
    count_matrix = np.zeros((adata.n_obs, len(cell_types)), dtype=np.int32)
    print(f'邻域半径：{radius_um} μm；样本数：{sample_values.nunique()}；细胞类型数：{len(cell_types)}')
    for sample in sample_values.unique():
        sample_positions = np.flatnonzero(sample_values.to_numpy() == sample)
        sample_obs = adata.obs.iloc[sample_positions]
        coords = sample_obs.loc[:, coord_cols].to_numpy(dtype=np.float64)
        if not np.isfinite(coords).all():
            raise ValueError(f'样本 {sample} 存在非有限空间坐标。')
        type_codes = pd.Categorical(sample_obs[cell_type_col].astype(str), categories=cell_types, ordered=True, ).codes
        if (type_codes < 0).any():
            invalid_types = sample_obs.loc[type_codes < 0, cell_type_col].astype(str).unique().tolist()
            raise ValueError(f'样本 {sample} 存在未成功编码的细胞类型：{invalid_types}')
        adjacency = radius_neighbors_graph(coords, radius=radius, mode='connectivity', metric='euclidean',
                                           include_self=False, n_jobs=n_jobs, )
        type_indicator = csr_matrix(
            (np.ones(len(sample_positions), dtype=np.uint8), (np.arange(len(sample_positions)), type_codes),),
            shape=(len(sample_positions), len(cell_types)), dtype=np.uint8, )
        sample_count_matrix = (adjacency @ type_indicator).toarray().astype(np.int32)
        if not np.array_equal(sample_count_matrix.sum(axis=1), adjacency.getnnz(axis=1)):
            raise AssertionError(f'样本 {sample} 的邻居总数与分类邻居数不一致。')
        count_matrix[sample_positions] = sample_count_matrix
        print(f'{sample}: {len(sample_positions)} cells, '
              f'{adjacency.nnz:,} neighbor pairs, '
              f'mean neighbors = {adjacency.getnnz(axis=1).mean():.2f}')
    neighbor_count_df = pd.DataFrame(count_matrix, index=adata.obs_names, columns=cell_types, )
    adata.obsm[store_key] = count_matrix
    adata.obs[f'{store_key}_total'] = count_matrix.sum(axis=1)
    adata.uns[store_key] = {'radius_mm': radius, 'radius_um': radius_um, 'coordinate_columns': list(coord_cols),
                            'sample_column': sample_col, 'cell_type_column': cell_type_col, 'cell_types': cell_types,
                            'include_self': False, }
    return neighbor_count_df
