from matplotlib.patches import FancyBboxPatch
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import Iterable


def crosstab_plot(crosstab=None, value_x=None, value_y=None, xlabel=None, ylabel=None, pad=-0.05,
                  rounding_size=0.15, x_order=None, y_order=None, table_cmap='seagreen', row_cmap='orangered',
                  col_cmap='deepskyblue', ax=None, background_param=None):
    """

    :param crosstab:
    :param ax:
    :param background_param:
    :param value_x:
    :param value_y:
    :param xlabel:
    :param ylabel:
    :param pad:
    :param rounding_size:
    :param x_order:
    :param y_order:
    :param table_cmap:
    :param row_cmap:
    :param col_cmap:
    :return:
    """

    if not crosstab and (isinstance(value_x, Iterable) and isinstance(value_y, Iterable)):
        if not xlabel:
            xlabel = 'X label'
        if not ylabel:
            ylabel = 'Y label'
        if not x_order:
            x_order = sorted(list(set(value_x)))
        if not y_order:
            y_order = sorted(list(set(value_y)))
        value_x = list(value_x)
        value_y = list(value_y)
        crosstab_input = pd.DataFrame([value_x, value_y], index=[xlabel, ylabel]).T
        crosstab = pd.crosstab(crosstab_input[ylabel], crosstab_input[xlabel])
        crosstab = crosstab[x_order].loc[y_order]

    elif isinstance(crosstab, pd.DataFrame):
        xlabel = xlabel or crosstab.columns.name
        ylabel = ylabel or crosstab.index.name

    else:
        raise ValueError('Require crosstab or (value_x and value_y)')

    if not ax:
        fig, ax = plt.subplots(figsize=(crosstab.shape[1] + 1, crosstab.shape[0] + 1))

    if not background_param:
        background_param = {'facecolor': 'white', 'linestyle': '--', 'edgecolor': 'black', 'alpha': 1}

    fig_width = crosstab.shape[1] + 1
    fig_height = crosstab.shape[0] + 1
    table_width = crosstab.shape[1]
    table_height = crosstab.shape[0]
    ax.axis(False)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    patchs = []

    # ������
    ax.table(crosstab.values, cellLoc='center', edges='open',
             bbox=((1 / fig_width), 0, table_width / fig_width, table_height / fig_height), zorder=2)

    # ������������
    ax.table([[i] for i in crosstab.index], cellLoc='center', edges='open', zorder=2,
             bbox=(0.5 / fig_width, 0, 0.5 / fig_width, table_height / fig_height))
    ax.text(x=0.25, y=table_height / 2, s=xlabel, rotation=90, ha='center', zorder=2, va='center')

    # ������������
    ax.table([crosstab.columns.to_list()], cellLoc='center', zorder=2,
             bbox=((1 / fig_width), table_height / fig_height, table_width / fig_width, 0.5 / fig_height), edges='open')
    ax.text(x=1 + table_width / 2, y=table_height + 0.75, s=ylabel, ha='center', zorder=2, va='center')

    # Բ�Ǿ���ɫ��
    # �з�����
    patch = FancyBboxPatch(xy=(0, 0), width=0.5, height=table_height, color='peru', zorder=1,
                           boxstyle=f'round,pad={pad},rounding_size={rounding_size}')
    patchs.append(patch)
    # �з���
    row_colors = sns.light_palette(row_cmap, table_height)
    for y, color in zip(range(table_height), row_colors):
        patch = FancyBboxPatch(xy=(0.5, y), width=0.5, height=1, color=color, zorder=1,
                               boxstyle=f'round,pad={pad},rounding_size={rounding_size}')
        patchs.append(patch)

    # �з�����
    patch = FancyBboxPatch(xy=(1, table_height + 0.5), width=table_width, height=0.5, color='steelblue', zorder=1,
                           boxstyle=f'round,pad={pad},rounding_size={rounding_size}')
    patchs.append(patch)
    # �з���
    col_colors = sns.light_palette(col_cmap, table_width)
    for x, color in zip(reversed(range(table_width)), col_colors):
        patch = FancyBboxPatch(xy=(x + 1, table_height), width=1, height=0.5, color=color, zorder=1,
                               boxstyle=f'round,pad={pad},rounding_size={rounding_size}')
        patchs.append(patch)

    # ������
    table_max = crosstab.max().max() + 1
    table_colors = sns.light_palette(table_cmap, 256)
    table_colors.reverse()
    for x in reversed(range(table_width)):
        for y in range(table_height):
            value = crosstab.iloc[-y, x]
            color = table_colors[int(np.floor(256 * value / table_max))]
            patch = FancyBboxPatch(xy=(x + 1, y), width=1, height=1, color=color, zorder=1,
                                   boxstyle=f'round,pad={pad},rounding_size={rounding_size}', ec='gainsboro')
            patchs.append(patch)

    # ������
    patch = FancyBboxPatch(xy=(0, 0), width=fig_width, height=fig_height, facecolor=background_param['facecolor'],
                           boxstyle=f'round,pad=0,rounding_size={rounding_size}', zorder=0, linewidth=2,
                           edgecolor=background_param['edgecolor'], linestyle=background_param['linestyle'])
    patchs.append(patch)

    # ��ax��������patch
    for patch in patchs:
        ax.add_patch(patch)

    return ax