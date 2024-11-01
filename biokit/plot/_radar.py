# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import numpy as np
from matplotlib import pyplot as plt


def radarplot(radar_df, groupby, value, groups=None, figsize=None, color_dict=None, facecolor='deepskyblue', ax=None,
              sns=None, label=None):
    if not figsize:
        figsize = (6, 6)
    if not groups:
        groups = radar_df[groupby].unique()
    if not color_dict:
        color_dict = dict(zip(groups, sns.color_palette('Set1', n_colors=len(groups))))

    radar_df = radar_df[radar_df[groupby].isin(groups)]
    if not ax:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))
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
