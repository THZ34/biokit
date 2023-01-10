import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text


def volcano_plot(df, x='log2fc', y='-log10p', color=None, color_dict=None, anno=None, ax=None):
    """

    :param anno:
    :param df:
    :param x:
    :param y:
    :param color:
    :param color_dict:
    :param ax:
    :return:
    """
    df = df.copy()
    if not color:
        color = 'color'
        df[color] = 'other'
        df.loc[(df[x] >= 1) & (df[y] >= -np.log10(0.05)), color] = 'up'
        df.loc[(df[x] <= -1) & (df[y] >= -np.log10(0.05)), color] = 'down'

    if not color_dict:
        color_dict = {'up': 'orangered', 'down': 'deepskyblue', 'other': 'grey'}

    if not ax:
        fig, ax = plt.subplots(figsize=(6, 4))

    if not anno:
        anno = 20
    if type(anno) == float:
        anno = df[(df[x] < -1) | (df[x] > 1)].sort_values(by=y, ascending=False).index[:int(df.shape[0] * anno)]
    elif type(anno) == int:
        anno = df[(df[x] < -1) | (df[x] > 1)].sort_values(by=y, ascending=False).index[:anno]
    elif type(anno) == list:
        anno = sorted(list(set(anno) & set(df.index)))

    # 画点
    for key in color_dict:
        temp_df = df[df[color] == key]
        ax.scatter(x=temp_df[x], y=temp_df[y], s=5, c=color_dict[key])

    # 注释基因
    texts = []
    for gene in anno:
        texts.append(ax.text(x=df.loc[gene][x], y=df.loc[gene][y], s=gene, fontsize=5,
                             bbox={'facecolor': color_dict[df.loc[gene][color]], 'alpha': 0.3, 'pad': 2,
                                   'linewidth': 0}))
        adjust_text(texts)

    ax.set_xlabel(x)
    ax.set_ylabel(y)
    return ax
