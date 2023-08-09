# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


# %%
def sigmoid(y1, y2, n_smooth=100):
    s_curve_x = np.linspace(0, 10, n_smooth) - 5
    s_curve_y = 1 / (1 + np.exp(-s_curve_x))
    curve_y = s_curve_y * (y2 - y1) + y1
    return curve_y


def mutation_timescape(df, timepoints=None, color_dict=None, figsize=None, n_smooth=100, alpha=0.5):
    if not color_dict:
        color_dict = dict(zip(df.index, sns.hls_palette(df.shape[0])))
    if not figsize:
        fig_length = df.shape[1] * 2 + 1
        fig_height = df.shape[0] + 2.5
        figsize = (fig_length, fig_height)
    if not timepoints:
        timepoints = list(range(df.shape[1]))

    raw_df = df.copy()
    df = raw_df.copy()
    df.insert(0, 'day0', 0)
    base_df = pd.DataFrame([[0] * df.shape[1]], index=['base'], columns=df.columns)
    df = pd.concat([base_df, df], axis=0)
    cum_df = df.cumsum(axis=0)
    center_df = cum_df - cum_df.mean()

    fig = plt.figure(figsize=figsize)
    grid_rows = 3 + df.shape[0]
    grid_cols = 1
    grid = plt.GridSpec(grid_rows, grid_cols)
    ax_dict = dict()

    # timescape
    ax = fig.add_subplot(grid[:3, :])
    ax_dict['timescape'] = ax
    timescape_ax_xlim = \
        timescape(raw_df, timepoints=timepoints, color_dict=color_dict, n_smooth=n_smooth, alpha=alpha, ax=ax)[1]
    ax.set_xlim(timescape_ax_xlim)
    # 突变频率
    for i, mut in enumerate(df.index[1:]):
        ax = fig.add_subplot(grid[4 + i:5 + i, :])
        ax_dict[mut] = ax
        ylim = (-0.2, 1.2)
        ax.bar(x=-0.25, height=ylim[1] - ylim[0], bottom=ylim[0], width=0.25, color=color_dict[mut], align='edge')
        ax.text(x=0.05, y=0.5, s=mut, fontsize=12, ha='left', va='center')
        points_x = np.linspace(2, timescape_ax_xlim[1] - 0.25, raw_df.shape[1])
        points_y = df.loc[mut][1:].to_numpy()
        points_y = (points_y - points_y.min()) / (points_y.max() - points_y.min())
        ax.plot(points_x, points_y, color=color_dict[mut], linewidth=2)
        for point_x, point_y, vaf in zip(points_x, points_y, df.loc[mut][1:]):
            ax.scatter(point_x, point_y, s=10, color=color_dict[mut])
            ax.text(x=point_x, y=point_y + 0.05 if point_y < 0.5 else point_y - 0.05,
                    s=f'{vaf * 100:0.2f}' if vaf > 0 else 'ND', ha='center',
                    va='bottom' if point_y < 0.5 else 'top')
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylim(ylim)
        ax.set_xlim(timescape_ax_xlim)

    plt.subplots_adjust(hspace=0, bottom=0.1, left=0.1, right=0.9)


def timescape(df, timepoints=None, color_dict=None, n_smooth=100, alpha=0.5, ax=None):
    """timescape图
    """
    if not color_dict:
        color_dict = dict(zip(df.index, sns.hls_palette(df.shape[0])))
    if not ax:
        fig, ax = plt.subplots(figsize=(8, 4))
    if not timepoints:
        timepoints = list(range(df.shape[1]))
    # 标准化vaf
    raw_df = df.copy()
    df = raw_df.copy()
    df.insert(0, 'day0', 0)
    base_df = pd.DataFrame([[0] * df.shape[1]], index=['base'], columns=df.columns)
    df = pd.concat([base_df, df], axis=0)
    cum_df = df.cumsum(axis=0)
    center_df = cum_df - df.sum() / 2
    # 标准化时间点横坐标
    timepoints = np.array(timepoints)
    timepoints = (timepoints - timepoints.min()) / (timepoints.max() - timepoints.min())
    timepoints = timepoints * raw_df.shape[1]
    timepoints = timepoints + 1
    timepoints = [0] + timepoints.tolist()
    # 生成曲线
    curves = []
    curvex = []
    for i in range(df.shape[1] - 1):
        curvex.extend(np.linspace(timepoints[i], timepoints[i + 1], n_smooth))
    for mut in df.index:
        curve_y = []
        for i in range(df.shape[1] - 1):
            x1, x2 = i, i + 1
            x1_name = df.columns[x1]
            x2_name = df.columns[x2]
            y1, y2 = center_df.loc[mut, x1_name], center_df.loc[mut, x2_name]
            temp_curve_y = sigmoid(y1, y2, n_smooth=n_smooth)
            curve_y.extend(temp_curve_y)
        curves.append(curve_y)
    # 填充曲线
    for i in range(raw_df.shape[0]):
        mut = df.index[i + 1]
        curvey1 = curves[i]
        curvey2 = curves[i + 1]
        ax.fill_between(curvex, curvey1, curvey2, color=color_dict[mut], alpha=alpha)
    timescape_ax_xlim = (-0.25, max(timepoints) + 0.25)
    ax.set_xticks(timepoints)
    xticklabels_bottom = ['timepoint'] + list(df.columns[1:])
    ax.set_xticklabels(xticklabels_bottom)
    ax.set_yticks([])
    ax.set_xlim(timescape_ax_xlim)
    ax.grid(axis='x', linestyle='--', alpha=1, linewidth=1)  # 网格线
    ax2 = ax.twiny()
    ax2.set_xticks(timepoints)
    xticklabels_upper = ['Max variant \nallele fraction (%)'] + [f'{vaf * 100:0.2f}' for vaf in df.max()[1:]]
    ax2.set_xticklabels(xticklabels_upper)
    ax2.set_xlim(timescape_ax_xlim)

    return (ax, ax2), timescape_ax_xlim, center_df
