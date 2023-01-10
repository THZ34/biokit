# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %% 
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

config = {"font.family": 'Microsoft YaHei', 'pdf.fonttype': 42}
rcParams.update(config)


def cumulative_bar(df, color_dict=None, ax=None):
    """累加柱状图"""
    if not ax:
        fig, ax = plt.subplots()
    if not color_dict:
        color_dict = dict(zip(df.index, sns.hls_palette(df.shape[0])))
    bottom = [0] * df.shape[1]
    for sample in df.index:
        ax.bar(x=range(df.shape[1]), height=df.loc[sample], bottom=bottom, color=color_dict[sample],
               label=sample, edgecolor='black')
        bottom = [i + j for i, j in zip(bottom, df.loc[sample])]
    ax.set_ylim(0, df.sum().max() * 1.2)
    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    return ax
