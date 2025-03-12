import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import seaborn as sns
from biokit.data_convert import p2text as p2text_func


def testbox(data, y, x=0, ylim=None, groupby=None, groups=None, testfunc=ttest_ind, kind='box', cutoff=None,
            width=0.8, interval='auto', ax=None, colors=None, cutoff_color=None, p2text=True, hide_flier=False,
            **kwargs):
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

    if not groups:
        groups = data[groupby].unique().tolist()
    n_groups = len(groups)
    if not ax:
        fig, ax = plt.subplots(figsize=(2 * n_groups, 8))

    sym = None
    if not ylim:
        if not hide_flier:
            ymax = data[data[groupby].isin(groups)][y].max()
            ymin = data[data[groupby].isin(groups)][y].min()
            ylength = ymax - ymin
            ptext_y_interval = ylength * 0.05
        elif hide_flier:
            iqr = data[data[groupby].isin(groups)][y].quantile(0.75) - data[data[groupby].isin(groups)][y].quantile(
                0.25)
            ymax = iqr * 1.5 + data[data[groupby].isin(groups)][y].quantile(0.75)
            ymin = data[data[groupby].isin(groups)][y].min()
            ylength = ymax - ymin
            ptext_y_interval = ylength * 0.05
            sym = ''
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
    data = data.copy()
    if interval == 'auto':
        interval = 1 - width
    box_width = width / n_groups
    interval_width = box_width * interval / width
    position_start = x - (width + interval) / 2 + box_width / 2
    positions = [position_start + i * box_width + i * interval_width for i in range(n_groups)]
    if kind == 'box':
        fig_objects = ax.boxplot([data[data[groupby] == group][y] for group in groups], positions=positions,
                                 widths=box_width, patch_artist=True, sym=sym, **kwargs)
        for patch, color in zip(fig_objects['boxes'], colors):
            patch.set_facecolor(color)
    elif kind == 'violin':
        for group, position in zip(groups, positions):
            temp_data = data[data[groupby] == group]
            temp_data['x'] = position
            sns.violinplot(data=temp_data, x=position, y=y, ax=ax, color=colors[groups.index(group)], **kwargs)
    elif kind == 'bar':
        data['x'] = x
        data = data[[groupby, 'x', y]]
        sns.barplot(data=data, x='x', y=y, hue=groupby, ax=ax, palette=colors, width=width, **kwargs)

    ptext_y_bottom = ymax + ptext_y_interval
    n_text = 0
    for interval in range(1, n_groups):
        for i in range(n_groups - interval):
            group1 = groups[i]
            group2 = groups[i + interval]
            x1 = positions[i]
            x2 = positions[i + interval]
            pvalue = testfunc(data[data[groupby] == group1][y].to_list(),
                              data[data[groupby] == group2][y].to_list()).pvalue
            if p2text:
                ptext = p2text_func(pvalue, cutoff)
                color = cutoff_color[ptext]
            else:
                ptext = f'{pvalue:0.4f}' if pvalue > 0.0001 else f'{pvalue:0.2e}'
                color = 'red' if pvalue < 0.05 else 'black'
            ptext_y = ptext_y_bottom + n_text * ptext_y_interval
            ax.text((x1 + x2) / 2, ptext_y, ptext, ha='center', va='bottom', color=color)
            if interval == 1:
                ax.plot([x1, x1, x2, x2],
                        [ptext_y - 0.6 * ptext_y_interval, ptext_y, ptext_y, ptext_y - 0.6 * ptext_y_interval],
                        color=color)
            elif interval > 1:
                ax.plot([x1, x1, x2, x2],
                        [ptext_y - 0.3 * ptext_y_interval, ptext_y, ptext_y, ptext_y - 0.3 * ptext_y_interval],
                        color=color)
            n_text += 1
    ax.set_xticks(positions)
    if len(ax.get_xticks()) == len(groups):
        ax.set_xticklabels(groups)

    ax.set_ylim(ymin - 0.05 * ylength, ymax + 0.15 * ylength)
    return ax
