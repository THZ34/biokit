# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %% 
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

config = {"font.family": 'Microsoft YaHei', 'pdf.fonttype': 42}
rcParams.update(config)


def cumulative_bar(df, color_dict=None, ax=None, normalize=False):
    """累加柱状图"""
    df = df.copy()
    if not ax:
        fig, ax = plt.subplots()
    if not color_dict:
        color_dict = dict(zip(df.index, sns.hls_palette(df.shape[0])))
    if normalize:
        df = df.div(df.sum(axis=0), axis=1)

    bottom = [0] * df.shape[1]
    for sample in df.index:
        ax.bar(x=range(df.shape[1]), height=df.loc[sample], bottom=bottom, color=color_dict[sample],
               label=sample, edgecolor='black')
        bottom = [i + j for i, j in zip(bottom, df.loc[sample])]
    ylim = 1 if normalize else df.sum().max() * 1.2
    ax.set_ylim(0, ylim)
    ax.set_xticks([])
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
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
    ax.set_yticklabels([])
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    if normalize:
        ax.set_xticklabels([f'{int(i * 100)}%' for i in ax.get_xticks()])
    return ax
