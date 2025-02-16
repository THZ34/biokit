# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import matplotlib.pyplot as plt
import seaborn as sns


def contour_map(df, value, x='x', y='y', levels1=None, levels2=None, cmap1='Blues_r', cmap2='Reds', alpha=1, ax=None):
    if levels1 is None:
        levels1 = [0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2]
    if levels2 is None:
        levels2 = [0.2, 0.4, 0.6, 0.8, 1]
    if ax is None:
        fig, ax = plt.subplots()

    for levels, cmap in zip([levels1, levels2], [cmap1, cmap2]):
        sns.kdeplot(data=df, x=x, y=y, weights=value, fill=True, levels=levels, ax=ax, bw_adjust=0.2,
                    cmap=cmap, alpha=alpha)

    return ax
