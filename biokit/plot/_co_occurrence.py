# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import seaborn as sns
from matplotlib import pyplot as plt

from data._co_occurrence import get_co_occurrence_df
from plot import create_fig


def cooccurrence_lineplot(adata, groupby, groups1=None, groups2=None, color_dict=None):
    if groups1 is None:
        groups1 = adata.obs[groupby].cat.categories
    if groups2 is None:
        groups2 = adata.obs[groupby].cat.categories
    if color_dict is None:
        groups = list(set(groups1) | set(groups2))
        color_dict = dict(zip(groups, sns.color_palette('Set2', n_colors=len(groups))))

    co_occurrence_df = get_co_occurrence_df(adata, groupby)
    fig, ax = create_fig(12, 4, len(groups1), 1)
    for group1 in groups1:
        for group2 in groups2:
            temp_cooccurrence_df = co_occurrence_df[(co_occurrence_df['group1'] == group1) & (
                    co_occurrence_df['group2'] == group2)]
            ax.plot(temp_cooccurrence_df.columns[2:],
                    temp_cooccurrence_df.values[0][2:],
                    color=color_dict[group2], label=group2)
            ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
            plt.tight_layout()
            plt.show()
