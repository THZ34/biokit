# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import math

import numpy as np
import pandas as pd
import seaborn as sns
from palantir.plot import _validate_gene_trend_input


def degrees_to_vector(angle_degrees):
    angle_radians = math.radians(angle_degrees)
    dx = math.cos(angle_radians)
    dy = math.sin(angle_radians)
    return dx, dy


def distance_in_direction(point1, point2, angle_degrees):
    direction = degrees_to_vector(angle_degrees)
    # 计算两个点的坐标向量
    vector_point1 = [point1[0], point1[1]]
    vector_point2 = [point2[0], point2[1]]
    # 计算两个点的坐标向量在指定方向上的距离
    distance = 0
    for i in range(len(direction)):
        distance += (vector_point2[i] - vector_point1[i]) * direction[i]

    return distance


def farthest_point_at_angle(points, angle_degrees):
    # 计算中心点
    center = [sum([point[0] for point in points]) / len(points), sum([point[1] for point in points]) / len(points)]

    if angle_degrees is not None:  # 有角度则计算距离中心点最远的点
        max_distance = 0
        farthest_point = None
        for point in points:
            # 计算点与中心点的距离
            distance = distance_in_direction(center, point, angle_degrees)
            if distance > max_distance:
                max_distance = distance
                farthest_point = point
        return farthest_point
    else:  # 否则选择离中心点最近的点
        max_distance = np.inf
        farthest_point = None
        for point in points:
            point1 = center
            point2 = point
            # 计算点与中心点的距离
            distance = np.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)
            if distance < max_distance:
                max_distance = distance
                farthest_point = point
        return farthest_point


def set_terminal_state(adata, groupby, group, angle):
    x, y = farthest_point_at_angle(adata[adata.obs[groupby] == group].obsm['X_umap'], angle_degrees=angle)
    cellid = adata[adata.obs[groupby] == group].obs_names[
        np.argmin(np.sum((adata[adata.obs[groupby] == group].obsm['X_umap'] - [x, y]) ** 2, axis=1))]
    return cellid


def set_terminal_states(adata, groupby, groups=None, angles=None):
    if not groups:
        groups = list(adata.obs[groupby].unique())
        groups.sort()
    if not angles:
        angles = [None] * len(groups)
    terminal_state_dict = {}
    for cluster, angle in zip(groups, angles):
        terminal_state_dict[set_terminal_state(adata, groupby, cluster, angle)] = cluster
    return pd.Series(terminal_state_dict)


def plot_gene_trends(data, gene=None, branches=None, gene_trend_key="gene_trends", branch_names="branch_masks",
                     ax=None, color_dict=None):
    """
    Plot the gene trends for each gene across different panels.

    Parameters
    ----------
    data : Union[Dict, sc.AnnData]
        An AnnData object or a dictionary that contains the gene trends.
    genes : Optional[List[str]], optional
        A list of genes to plot. If not provided, all genes will be plotted. Default is None.
    gene_trend_key : str, optional
        The key to access gene trends in the varm of the AnnData object. Default is 'gene_trends'.
    branch_names : Union[str, List[str]], optional
        Key to retrieve branch names from the AnnData object, or a list of branch names. If a string is provided,
        it will be treated as a key in either AnnData.uns or AnnData.obsm. For AnnData.obsm, the column names will
        be treated as branch names. If it cannot be found in AnnData.obsm and AnnData.uns then branch_names + "_columns"
        will be looked up in AnnData.uns.
        Default is 'branch_masks'.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object containing the plot.

    Raises
    ------
    KeyError
        If 'branch_names' as a string is not found in either .uns or .obsm, or if 'gene_trend_key + "_" + branch_name'
        is not found in .varm.
    ValueError
        If 'data' is not an AnnData object and not a dictionary.
    """
    gene_trends = _validate_gene_trend_input(data, gene_trend_key, branch_names)
    # Branches and genes
    if not branches:
        branches = list(gene_trends.keys())
    if not color_dict:
        color_dict = dict(zip(branches, sns.color_palette('tab10', len(branches))))
        color_dict = pd.Series(color_dict)
    for branch in branches:
        trends = gene_trends[branch]["trends"]
        ax.plot(trends.columns.astype(float), trends.loc[gene, :], color=color_dict[branch], label=branch, )
        ax.set_xticks([0, 1])
        ax.set_title(gene)
