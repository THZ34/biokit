import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import f_oneway, ttest_ind
import scipy


def p2text(p, cutoff):
    for cutoff_value in sorted(list(cutoff.keys())):
        if p <= cutoff_value:
            return cutoff[cutoff_value]
    return 'ns'


def testbox(data, y, x=0, ylim=None, groupby=None, groups=None, testfunc=ttest_ind, kind='box', cutoff=None,
            width=0.8, ax=None, colors=None, cutoff_color=None):
    """

    :param data:
    :param y:
    :param x:
    :param groupby:
    :param groups:
    :param testfunc:
    :param kind:
    :param cutoff:
    :param width:
    :param ax:
    :param colors:
    :param cutoff_color:
    :return:
    """
    from seaborn import color_palette
    import math
    # 生成默认参数

    if not groups:
        groups = data[groupby].unique().tolist()
    n_groups = len(groups)
    if not ax:
        fig, ax = plt.subplots(figsize=(2 * n_groups, 8))
    if not ylim:
        ymax = data[data[groupby].isin(groups)][y].max()
        ymin = data[data[groupby].isin(groups)][y].min()
        ylength = ymax - ymin
        ptext_y_interval = ylength * 0.05
    else:
        ymin, ymax = ylim
        ylength = ymax - ymin
        ptext_y_interval = ylength * 0.05
    if not set(groups).issubset(set(data[groupby])):
        raise ValueError('Groups must be subset of data[groupby].values')
    if not colors:
        colors = color_palette(n_colors=len(groups))
    if not cutoff:
        cutoff = {0.05: '*', 0.01: '**', 0.001: '***', 0.0001: '****'}
    if not cutoff_color:
        cutoff_color = dict(
            zip(['*', '**', '***', '****', 'ns'], ['orange', 'darkorange', 'orangered', 'red', 'deepskyblue']))

    # 确定box坐标
    data = data.copy()
    box_width = width / n_groups
    position_start = x - width / 2 + box_width / 2
    positions = [position_start + i * box_width for i in range(n_groups)]
    # 画图box
    for group, position, color in zip(groups, positions, colors):
        ax.boxplot(data[data[groupby] == group][y], positions=[position],
                   widths=box_width, patch_artist=True, boxprops={'facecolor': color})
    # 依次计算组间p-value并在图中标注
    ptext_y_bottom = ymax + ptext_y_interval
    n_text = 0
    for interval in range(1, n_groups):
        for i in range(n_groups - interval):
            group1 = groups[i]
            group2 = groups[i + interval]
            x1 = positions[i]
            x2 = positions[i + interval]
            pvalue = testfunc(data[data[groupby] == group1][y], data[data[groupby] == group2][y]).pvalue
            ptext = p2text(pvalue, cutoff)
            ptext_y = ptext_y_bottom + n_text * ptext_y_interval
            ax.text((x1 + x2) / 2, ptext_y, ptext, ha='center', va='bottom',
                    color=cutoff_color[ptext])
            ax.plot([x1, x1, x2, x2], [ptext_y - 0.6 * ptext_y_interval, ptext_y,
                                       ptext_y, ptext_y - 0.6 * ptext_y_interval], color=cutoff_color[ptext])
    return ax

