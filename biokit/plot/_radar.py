# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import numpy as np
from matplotlib import pyplot as plt


def radarplot(radar_df, groupby, value, groups=None, color_dict=None, hue=None, hue_color_dict=None, ax=None):
    radar_df = radar_df.copy()
    if not groups:
        groups = radar_df[groupby].unique()
    if not color_dict:
        import seaborn as sns
        color_dict = dict(zip(groups, sns.color_palette('Set1', n_colors=len(groups))))
    if hue:
        if not hue_color_dict:
            hue_color_dict = dict(
                zip(radar_df[hue].unique(), sns.color_palette('Set1', n_colors=len(radar_df[hue].unique()))))
    if not ax:
        fig, ax = plt.subplots(subplot_kw=dict(polar=True))

    if not hue:
        radarplot_base(radar_df, groupby, value, groups=groups, color_dict=color_dict, ax=ax)
    else:
        hue_groups = radar_df[hue].unique()
        print(hue_groups)
        for i, hue_group in enumerate(hue_groups):
            radarplot_base(radar_df[radar_df[hue] == hue_group], groupby, value, groups=groups,
                           color_dict=color_dict, ax=ax, label=hue_group, facecolor=hue_color_dict[hue_group])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    return ax


def radarplot_base(radar_df, groupby, value, groups=None, figsize=None, color_dict=None, facecolor='deepskyblue',
                   ax=None, label=None):
    if not figsize:
        figsize = (6, 6)
    if groups is None:
        groups = radar_df[groupby].unique()
    if not color_dict:
        import seaborn as sns
        color_dict = dict(zip(groups, sns.color_palette('Set1', n_colors=len(groups))))
    if not ax:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))

    radar_df = radar_df[radar_df[groupby].isin(groups)]
    labels = radar_df[groupby].to_numpy()
    stats = radar_df[value].to_numpy()
    # 计算角度
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    # 闭合图形
    stats = np.concatenate((stats, [stats[0]]))
    angles += angles[:1]
    # 绘图
    ax.fill(angles, stats, color=facecolor, alpha=0.25)
    ax.plot(angles, stats, color=facecolor, linewidth=2, label=label)
    ax.scatter(angles[:-1], stats[:-1], color=[color_dict[i] for i in labels], s=50, zorder=2)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)
    ax.grid(color='lightgrey', linestyle='-', linewidth=2)
    ax.spines['polar'].set_visible(False)
    return ax
