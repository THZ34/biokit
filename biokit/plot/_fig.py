# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import warnings
import matplotlib.pyplot as plt
import numpy as np


def create_fig(ax_width=4, ax_height=4, n_rows=1, n_cols=1, bottom=0.1, top=0.85, left=0.05, right=0.95, wspace=0.2,
               hspace=0.3):
    """
    基于期望单轴尺寸(ax_width/ax_height, 单位英寸)计算figure大小，并用subplots生成Axes；
    - constrained_layout=True 时走官方自动布局（忽略 left/right/top/bottom/wspace/hspace）；
    - 否则用 subplots_adjust 精确应用边距和间距。
    """
    # 总轴区尺寸（含轴间距，不含页边距）
    all_ax_w = n_cols * ax_width + (n_cols - 1) * wspace * ax_width
    all_ax_h = n_rows * ax_height + (n_rows - 1) * hspace * ax_height
    # 依据边距反推整页尺寸
    fig_w = all_ax_w / max(right - left, 1e-6)
    fig_h = all_ax_h / max(top - bottom, 1e-6)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
    if n_rows == 1 and n_cols == 1:
        axes = axes if isinstance(axes, plt.Axes) else axes[0]
        axes = axes if isinstance(axes, plt.Axes) else axes[0]
        return fig, axes
    axes = np.array(axes).reshape(n_rows, n_cols)
    return fig, axes


def create_fig_old(ax_width=4, ax_height=4, n_rows=1, n_cols=1, bottom=0.1, top=0.85, left=0.05, right=0.95, wspace=0.2,
                   hspace=0.3, constrained_layout=True):
    warnings.warn("create_fig_old 已过时，请使用 create_fig_new", DeprecationWarning, stacklevel=2)
    """Create a figure and axes"""
    all_ax_height = (n_rows + hspace * (n_rows - 1)) * ax_height
    fig_height = all_ax_height / (top - bottom)
    all_ax_width = (n_cols + wspace * (n_cols - 1)) * ax_width
    fig_width = all_ax_width / (right - left)
    ax_width_ratio = ax_width / fig_width
    ax_height_ratio = ax_height / fig_height
    axes = []
    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=constrained_layout)
    for row in range(n_rows):
        ax_line = []
        for col in range(n_cols):
            ax_line.append(fig.add_axes([left + (1 + wspace) * ax_width_ratio * col,
                                         bottom + (1 + hspace) * ax_height_ratio * (n_rows - row - 1), ax_width_ratio,
                                         ax_height_ratio]))
        axes.append(ax_line)
    axes = np.array(axes)
    if n_rows == 1 and n_cols == 1:
        axes = axes[0, 0]
    return fig, axes
