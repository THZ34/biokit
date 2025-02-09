# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %% 
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
import numpy as np

config = {"font.family": 'Microsoft YaHei', 'pdf.fonttype': 42}
rcParams.update(config)


def cumulative_bar(df, color_dict=None, ax=None, normalize=False, proportion=False, counts=False, alpha=1):
    """累加柱状图"""
    df = df.copy()
    if not ax:
        fig, ax = plt.subplots()
    if not color_dict:
        color_dict = dict(zip(df.index, sns.hls_palette(df.shape[0])))
    if normalize:
        df = df.div(df.sum(axis=0), axis=1)
    if proportion:
        proportion_df = df.div(df.sum(axis=1), axis=0)
        proportion_df = proportion_df.mul(100).round(2)

    bottom = [0] * df.shape[1]
    for sample in df.index:
        ax.bar(x=range(df.shape[1]), height=df.loc[sample], bottom=bottom, color=color_dict[sample],
               label=sample, edgecolor='black', alpha=alpha)
        bottom = [i + j for i, j in zip(bottom, df.loc[sample])]
    ylim = 1 if normalize else df.sum().max() * 1.2
    ax.set_ylim(0, ylim)
    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns)
    handles, labels = ax.get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))
    if normalize:
        ax.set_yticklabels([f'{int(i * 100)}%' for i in ax.get_yticks()])

    return ax


def cumulative_barh(df, color_dict=None, ax=None, normalize=False):
    """累加柱状图"""
    df = df.copy()
    if not ax:
        fig, ax = plt.subplots()
    if not color_dict:
        color_dict = dict(zip(df.index, sns.hls_palette(df.shape[0])))
    if normalize:
        df = df.div(df.sum(axis=0), axis=1)

    left = [0] * df.shape[1]
    for sample in df.index:
        ax.barh(y=range(df.shape[1]), width=df.loc[sample], left=left, color=color_dict[sample],
                label=sample, edgecolor='black')
        left = [i + j for i, j in zip(left, df.loc[sample])]
    xlim = 1 if normalize else df.sum().max() * 1.2
    ax.set_xlim(0, xlim)
    ax.set_yticks(range(df.shape[1]))
    ax.set_yticklabels(df.columns)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    if normalize:
        ax.set_xticklabels([f'{int(i * 100)}%' for i in ax.get_xticks()])
    return ax
