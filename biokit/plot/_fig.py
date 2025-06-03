# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import matplotlib.pyplot as plt
import numpy as np

def create_fig(ax_width=4, ax_height=4, n_rows=1, n_cols=1, bottom=0.1, top=0.85, left=0.05, right=0.75, wspace=0.2,
               hspace=0.3):
    """Create a figure and axes"""
    all_ax_height = (n_rows + wspace * (n_rows - 1)) * ax_height
    fig_height = all_ax_height / (top - bottom)
    all_ax_width = (n_cols + hspace * (n_cols - 1)) * ax_width
    fig_width = all_ax_width / (right - left)
    ax_width_ratio = ax_width / fig_width
    ax_height_ratio = ax_height / fig_height

    axes = []
    fig = plt.figure(figsize=(fig_width, fig_height))
    for row in range(n_rows):
        ax_line = []
        for col in range(n_cols):
            ax_line.append(fig.add_axes([left + (1 + wspace) * ax_width_ratio * col,
                                         bottom + (1 + hspace) * ax_height_ratio * (n_rows - row - 1), ax_width_ratio,
                                         ax_height_ratio]))
        axes.append(ax_line)
    if n_rows == 1 and n_cols == 1:
        axes = axes[0][0]
    return fig, np.array(axes)
