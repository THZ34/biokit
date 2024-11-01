# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def marker_positive_contour(adata, marker, pos_cutoff=0, coord_key='X_umap', ax=None):
    """marker阳性细胞的密度图

    :param adata: AnnData
    """
    temp_adata = adata[adata[:, marker].X > pos_cutoff]
    if ax is None:
        fig, ax = plt.subplots()
    scatter_df = pd.DataFrame(temp_adata.obsm[coord_key], columns=['x', 'y'], index=temp_adata.obs_names)
    scatter_df['color'] = temp_adata[:, marker].to_df()
    sns.scatterplot(data=scatter_df, x='x', y='y', hue='color', ax=ax)
    sns.kdeplot(data=scatter_df, x='x', y='y', hue='color', ax=ax, )
    return ax
