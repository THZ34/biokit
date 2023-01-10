import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact, ttest_ind


def fisher_bar(fisher_df, x, y, y_target, x_order=None, ax=None, y_max=100):
    if not ax:
        fig, ax = plt.subplots()
    if not x_order:
        x_order = sorted(list(fisher_df[x].unique()))

    fisher_table = pd.crosstab(fisher_df[x], fisher_df[y])
    p_value = fisher_exact(fisher_table)[1]
    for x_target in x_order:
        percent = fisher_table.loc[x_target][y_target] / fisher_table.loc[x_target].sum() * 100
        ax.bar(x=x_target, height=percent)
    ax.plot([0, 0, 1, 1], [1.02 * y_max, 1.05 * y_max, 1.05 * y_max, 1.02 * y_max], color='black', linewidth=1)
    ax.text(x=0.5, y=1.05 * y_max, s=f'{p_value:0.3f}' if p_value > 0.001 else f'{p_value:0.3e}', ha='center',
            va='bottom', color='red' if p_value < 0.05 else 'black')

    ax.set_ylim(0, y_max * 1.15)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(x_order)
    ax.set_title(x)
    ax.set_ylabel(f'{y}={y_target}(%)')
    return ax, p_value


def ttest_box(ttest_df, x, y, x_order=None, ax=None, kind='boxen', y_min=None, y_max=None, colors=None):
    ttest_df = ttest_df.copy()
    if not ax:
        fig, ax = plt.subplots(figsize=(3, 6))
    if not x_order:
        x_order = sorted(list(ttest_df[x].unique()))
    if len(x_order) != 2:
        raise ValueError
    if not y_max:
        y_max = ttest_df[y].max()
    if not y_min:
        y_min = ttest_df[y].min()

    if not kind in ['box', 'boxen', 'violin']:
        raise ValueError('kind allowed box, boxen and violin')
    if not colors:
        colors = ['deepskyblue', 'orangered']

    replace_dict = {x_order[0]: 0, x_order[1]: 1}
    ttest_df[x].replace(replace_dict, inplace=True)

    p_value = ttest_ind(ttest_df[ttest_df[x] == 0][y], ttest_df[ttest_df[x] == 1][y])[1]
    if kind == 'box':
        sns.boxplot(data=ttest_df, x=x, y=y, ax=ax, palette=colors, width=0.6)
    elif kind == 'boxen':
        sns.boxenplot(data=ttest_df, x=x, y=y, ax=ax, outlier_prop=0.01, k_depth="proportion",
                      scale='exponential', palette=colors)
    elif kind == 'violin':
        sns.violinplot(data=ttest_df, x=x, y=y, ax=ax, outlier_prop=0.01, k_depth="proportion", cut=0,
                       scale='width', palette=colors)
    ax.plot([0, 0, 1, 1], [y_max * 1.02 - y_min * 0.02, y_max * 1.05 - y_min * 0.05, y_max * 1.05 - y_min * 0.05,
                           y_max * 1.02 - y_min * 0.02], color='black', linewidth=1)
    ax.text(x=0.5, y=y_max * 1.05 - y_min * 0.05, s=f'{p_value:0.3f}' if p_value > 0.001 else f'{p_value:0.3e}',
            ha='center',
            va='bottom', color='red' if p_value < 0.05 else 'black')

    ax.set_ylim(y_min, y_max * 1.15 - y_min * 0.15)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(x_order)
    ax.set_title(y)
    return ax, p_value
