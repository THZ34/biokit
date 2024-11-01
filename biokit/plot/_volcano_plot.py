import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text


def volcano_plot(df, x='log2fc', y='-log10p', color=None, color_dict=None, anno=None, ax=None, groupnames=None,
                 textadjust=True):
    """

    :param textadjust:
    :param groupnames:
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
        fig, ax = plt.subplots(figsize=(6, 5))

    if anno is None:
        anno = 20
    if type(anno) == float:
        anno = df[(df[x] < -1) | (df[x] > 1)].sort_values(by=y, ascending=False).index[:int(df.shape[0] * anno)]
    elif type(anno) == int:
        anno = df[(df[x] < -1) | (df[x] > 1)].sort_values(by=y, ascending=False).index[:anno]
    elif type(anno) == list:
        anno = sorted(list(set(anno) & set(df.index)))
    if not groupnames:
        groupnames = ('control', 'case')
    controlname, casename = groupnames

    # 画点
    for key in color_dict:
        temp_df = df[df[color] == key]
        ax.scatter(x=temp_df[x], y=temp_df[y], s=5, c=color_dict[key])

    # 注释基因
    texts = []
    for gene in anno:
        texts.append(ax.text(x=df.loc[gene][x], y=df.loc[gene][y], s=gene, fontsize=8,
                             bbox={'facecolor': color_dict[df.loc[gene][color]], 'alpha': 0.3, 'pad': 2,
                                   'linewidth': 0}))
    if textadjust:
        adjust_text(texts, only_move={'points': 'y', 'texts': 'y'})

    # case control 箭头
    xmax = df[x].abs().max()
    ymax = df[y].max()
    # width = 0.3
    width = ymax * 0.04
    head_width = width * 1.5
    head_length = 0.15 * xmax
    shape = 'full'
    ax.arrow(x=-0.1 * xmax, y=ymax * 1.1, dx=-0.8 * xmax, dy=0, color=color_dict['down'], length_includes_head=True,
             width=width, head_width=head_width, shape=shape, head_length=head_length)
    ax.arrow(x=0.1 * xmax, y=ymax * 1.1, dx=0.8 * xmax, dy=0, color=color_dict['up'], length_includes_head=True,
             width=width, head_width=head_width, shape=shape, head_length=head_length)
    ax.text(x=0.12 * xmax, y=ymax * 1.15, s=f'high-exp in {casename}', fontsize=10, fontweight='bold',
            ha='left')
    ax.text(x=-0.12 * xmax, y=ymax * 1.15, s=f'high-exp in {controlname}', fontsize=10, fontweight='bold',
            ha='right')
    ax.text(x=0.5 * xmax, y=ymax * 1.1, s=f'n_genes = {df[color].value_counts().to_dict().get("up", 0)}', fontsize=10,
            fontweight='bold', ha='center', va='center')
    ax.text(x=-0.5 * xmax, y=ymax * 1.1, s=f'n_genes = {df[color].value_counts().to_dict().get("down", 0)}',
            fontsize=10,
            fontweight='bold', ha='center', va='center')
    #
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_xlim(xmax * -1.05, xmax * 1.05)
    ax.set_ylim(0, ymax * 1.23)
    # ax.set_ylim(0, ymax)
    return ax
