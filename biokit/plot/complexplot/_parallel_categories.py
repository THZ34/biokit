# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from biokit.data_convert import bezier_curve_S


def parallel_categories(df, x_groups=None, y_groups=None, width=0.2, x_color_dict=None, y_color_dict=None, ax=None,
                        normalize=True, interval_pct=3, edgewidth=1, n_smooth=1000, annot=True, annot_fmt='{:.2f}',
                        annot_loc=-0.2, alpha=1):
    if not x_groups:
        x_groups = df.columns.to_list()
    if not y_groups:
        y_groups = df.index.to_list()
    if not ax:
        fig, ax = plt.subplots(figsize=(len(x_groups), 4))
    if not x_color_dict:
        x_color_dict = dict.fromkeys(x_groups, 'none')
    if not y_color_dict:
        y_color_dict = dict(zip(y_groups, sns.color_palette(n_colors=df.shape[0])))
    if normalize:
        df = df.div(df.sum(axis=0), axis=1)
    df = df.mul(100).round(2)
    if annot is True:
        annot_fmt = '{:.2f}%'

    df = df.loc[y_groups,x_groups]
    assert len(x_groups) > 1, 'x_groups must be more than 1'

    interval = df.sum().max() * interval_pct / 100
    # 累计柱
    bottoms = np.array([0] * df.shape[1])
    for y_group in y_groups:
        heights = df.loc[y_group]
        ax.bar(x=np.arange(df.shape[1]), height=heights, width=width, bottom=bottoms, color=y_color_dict[y_group],
               label=y_group, alpha=alpha, edgecolor=[x_color_dict[x_group] for x_group in x_groups],
               linewidth=edgewidth,
               zorder=1)
        for x, y, height in zip(np.arange(df.shape[1]), bottoms, heights):
            ax.text(x=x + annot_loc, y=y + height / 2, s=annot_fmt.format(height), ha='center', va='center',
                    color='black', fontsize=8, fontweight='bold')
        bottoms += heights
        bottoms += interval
    # 组间变化曲线
    bottoms = np.array([0] * df.shape[1])
    for y_group in y_groups:
        heights = df.loc[y_group]
        for i in range(1, df.shape[1]):
            x1, y1 = i - 1 + width / 2, bottoms[i - 1]
            x2, y2 = i - width / 2, bottoms[i]
            h1, h2 = heights[i - 1], heights[i]
            x3, y3 = x2, y2 + h2
            x4, y4 = x1, y1 + h1
            control1 = (x1 + x2) / 2, y1
            control2 = (x1 + x2) / 2, y2
            control3 = (x3 + x4) / 2, y3
            control4 = (x3 + x4) / 2, y4
            curve1 = bezier_curve_S(np.array([[x1, y1], control1, control2, [x2, y2]]), n_smooth)
            curve2 = bezier_curve_S(np.array([[x3, y3], control3, control4, [x4, y4]]), n_smooth)
            curve_x = np.concatenate([[x1], curve1[::-1, 0], [x2, x2], curve2[::-1, 0], [x1, x1]])
            curve_y = np.concatenate([[y1], curve1[::-1, 1], [y2, y3], curve2[::-1, 1], [y4, y1]])
            ax.fill(curve_x, curve_y, color=y_color_dict[y_group], linewidth=0.5, alpha=alpha * 0.75, zorder=0)
        bottoms += heights
        bottoms += interval

    ax.set_xticks(np.arange(df.shape[1]))
    ax.set_xticklabels(x_groups)
    ylength = (df.sum().max() - df.min().min()) * (1 + interval_pct / 100 * (df.shape[0] - 1))
    # ax.set_ylim(df.min().min() + ylength, df.min().min())
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.tight_layout()
    return ax
