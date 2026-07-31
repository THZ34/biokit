# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


"""
适用于单细胞和空转数据的探索性分析

需要提供的内容包括:
1.adata: 单细胞或者空转的表达量数据
2.sample_col: 样本信息列名, 需要在adata.obs中
3.groupby
4.cell_type
5.tissue

"""

import os
import sys
from biokit.plot import create_fig
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import matplotlib.patches as mpatches_legend
import matplotlib.patches as mpatches
from colorcet import LinearSegmentedColormap


# %%
def _prepare_groups(adata, groupby, groups):
    if not groups:
        groups = adata.obs[groupby].unique()
    return groups


def _prepare_color_dict(groups, color_dict):
    if not color_dict:
        color_dict = dict(zip(groups, sns.color_palette('rainbow_r', len(groups))))
    return color_dict


def _prepare_cmap(cmap):
    if not cmap:
        cmap = LinearSegmentedColormap.from_list('default_cmap', ['#2E86C1', '#FFFFFF', '#E66B55'], N=256)
    return cmap


# %%
def group_heatmap(sample_info, groupby, groups, heatmap_df, color_dict=None, cmap=None, zscore=False):

    groups = _prepare_groups(sample_info, groupby, groups)
    color_dict = _prepare_color_dict(groups, color_dict)
    cmap = _prepare_cmap(cmap)

    # 数据预处理
    heatmap_df = heatmap_df.copy()
    if zscore:
        _z = heatmap_df.sub(heatmap_df.mean(axis=0), axis=1).div(heatmap_df.std(axis=0).replace(0, 1),
            axis=1)
        heatmap_df = _z.T  # 14 program × 12 patient
    vmax = float(np.abs(heatmap_df.values).max())

    """
    设计规则:
    以ax_heatmap为中心, ax的长宽与program_score_df.shape / 2相等
    """
    # === Step 1: sns.heatmap square=True ===
    ax_width = heatmap_df.shape[1] / 2
    ax_height = heatmap_df.shape[0] / 2
    fig, ax_heatmap = create_fig(ax_width=ax_width, ax_height=ax_height, left=0.13, bottom=0.08, right=0.95, top=0.95)
    sns.heatmap(heatmap_df, cmap=cmap, center=0, vmin=-vmax, vmax=vmax, linewidths=0.3, linecolor='white',
                square=True, cbar_kws={'label': 'Program z-score'}, ax=ax_heatmap)
    # === Step 2: 读 ax_heatmap 实际位置, 算 line_height ===
    actual = ax_heatmap.get_position()
    left, bottom, width, height = actual.x0, actual.y0, actual.width, actual.height
    n_rows = heatmap_df.shape[0]  # 14 program
    line_height = height / n_rows

    # === Step 3: ax_group 顶部色块 (12 patient 按 Dual/Triple) ===
    ax_group = fig.add_axes([left, bottom + height + line_height * 0.2, width, line_height * 0.5])
    # 画色块 + label
    group_series = sample_info.loc[heatmap_df.columns,groupby]
    for i, p in enumerate(heatmap_df.columns):
        c = group_palette[group_series.loc[p]]
        ax_group.add_patch(
            mpatches.Rectangle((i, 0), 1, 1, facecolor=c, edgecolor='white', linewidth=0.2, label=group_series.loc[p]))
    ax_group.set_xlim(0, len(heatmap_df.columns))
    ax_group.set_ylim(0, 1)
    ax_group.axis('off')

    # === Step 4: cax 高度 = height*0.5, bottom = bottom ===
    cax = ax_heatmap.collections[0].colorbar.ax
    cax.set_position([left + width + 0.02, bottom, 0.025, height * 0.5])
    cax.tick_params(labelsize=7.5)
    cax.set_ylabel('Program z-score', fontsize=8.5, labelpad=6)

    # === Step 5: ticks + labels (含 TLS 副标 vermillion 加粗) ===
    ax_heatmap.set_xticks(np.arange(len(heatmap_df.columns)) + 0.5)
    ax_heatmap.set_xticklabels(heatmap_df.columns, rotation=45, ha='right', fontsize=8.5)
    # y 轴: program + 副标 (in cancer cells / in immune cells)
    program_subtitle_map = {'TLS_signature': 'immune'}
    y_labels = []
    y_colors = []
    for p in heatmap_df.index:
        label = program_label_dict[p]
        if program_subtitle_map.get(p, 'cancer') == 'immune':
            y_labels.append(f'{label}\nin immune cells')
            y_colors.append('#D55E00')  # vermillion
        else:
            y_labels.append(f'{label}\nin cancer cells')
            y_colors.append('black')
    ax_heatmap.set_yticks(np.arange(len(y_labels)) + 0.5)
    ax_heatmap.set_yticklabels(y_labels, rotation=0, fontsize=7.5)
    for tick_label, color in zip(ax_heatmap.get_yticklabels(), y_colors):
        tick_label.set_color(color)
        tick_label.set_fontweight('bold')
    ax_heatmap.tick_params(axis='both', length=3, width=0.6, direction='in')

    # === Step 6: title y = bottom + height + line_height ===
    title_y = bottom + height + line_height
    fig.suptitle('14 functional-state programs (patient-level z-score)', fontsize=11, fontweight='bold', y=title_y,
                 va='bottom')

    # === Step 7: group legend (ax_group 内右上 + bbox 推到图外) ===
    legend_handles = [mpatches_legend.Patch(color=group_palette[k], label=k) for k in
                      ['Dual therapy', 'Triple therapy']]
    ax_group.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=8, frameon=True)
