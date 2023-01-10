# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %% 
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.container import Container
from matplotlib.patches import Wedge
from scipy.sparse import coo_matrix

config = {"font.family": 'Microsoft YaHei', 'pdf.fonttype': 42}
rcParams.update(config)


# %%
def heatmap(mutations, ax, color_dict, style='simple', sep=None):
    heatmapfunc_dict = {'simple': heatmap_simple, 'circle': heatmap_circledot, 'box': heatmap_cumulativebox}
    ax.invert_yaxis()
    ax.set_xlim(-0.5, mutations.shape[1] - 0.5)
    ax.set_ylim(-0.5, mutations.shape[0] - 0.5)
    ax.axis('off')
    return heatmapfunc_dict[style](mutations, ax, color_dict, sep)


def heatmap_simple(mutations, ax, color_dict=None, sep=None):
    """

    :param mutations:
    :param ax:
    :param color_dict:
    :param sep:
    :return:
    """
    sparse_df = coo_matrix(mutations)
    heatmap_sparse_df = pd.DataFrame([sparse_df.row, sparse_df.col, sparse_df.data]).T
    heatmap_sparse_df.columns = ['y', 'x', 'variant']
    for variant in color_dict.keys():
        temp_df = sparse_df[sparse_df['variant'] == variant]
        ax.bar(x=temp_df['x'], bottom=temp_df['y'] - 0.4, height=0.8, width=0.8, edgecolor='grey',
               color=color_dict[variant], label=variant, align='center')
    return ax


def heatmap_circledot(mutations, ax, color_dict=None, sep=','):
    """

    :param mutations:
    :param ax:
    :param color_dict:
    :param sep:
    :return:
    """
    sparse_df = coo_matrix(mutations)
    new_sparse_df = []
    for row, col, vars in zip(sparse_df.row, sparse_df.col, sparse_df.data):
        if sep in vars:
            vars = vars.split(sep)
            interval = 360 / len(vars)
            for t1, t2, var in zip(np.arange(-interval, 360 - interval, interval), np.arange(0, 360, interval),
                                   vars):
                var = var.capitalize()
                new_sparse_df.append([row, col, t1, t2, var])
        else:
            var = vars
            new_sparse_df.append([row, col, 0, 360, var])
    sparse_df = pd.DataFrame(new_sparse_df, columns=['y', 'x', 'theta1', 'theta2', 'var'])

    for var in color_dict.keys():
        patches = []
        temp_sparse_df = sparse_df[sparse_df['var'] == var]
        for y, x, t1, t2, var in temp_sparse_df.to_numpy():
            patch = Wedge(center=(x, y), r=0.4, theta1=t1, theta2=t2, fc=color_dict[var], label='_nolegend_', ec='grey')
            patches.append(patch)
            ax.add_patch(patch)
        wedge_container = Container(patches, label=var)
        ax.add_container(wedge_container)
    return ax


def heatmap_cumulativebox(mutations, ax, color_dict, sep=','):
    """

    :param mutations:
    :param ax:
    :param color_dict:
    :param sep:
    :return:
    """
    sparse_df = coo_matrix(mutations)
    new_sparse_df = []
    for row, col, vars in zip(sparse_df.row, sparse_df.col, sparse_df.data):
        height = 0.8
        if sep in vars:
            vars = vars.split(sep)
            interval = height / len(vars)
            for bottom, var in zip(np.arange(-height / 2, height / 2, interval), vars):
                new_sparse_df.append([col, row + bottom, interval, var])
        else:
            var = vars
            new_sparse_df.append([col, row - height / 2, height, var])
    sparse_df = pd.DataFrame(new_sparse_df, columns=['x', 'bottom', 'height', 'var'])
    for var in color_dict:
        temp_sparse_df = sparse_df[sparse_df['var'] == var]
        ax.bar(x=temp_sparse_df['x'], bottom=temp_sparse_df['bottom'],
               height=temp_sparse_df['height'], label=var, color=color_dict[var], ec='grey', align='center')
    return ax
