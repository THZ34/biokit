from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import pandas as pd


def fisher_bar(fisher_df, x, y, y_target, x_order=None, ax=None):
    if not ax:
        fig, ax = plt.subplots()
    if not x_order:
        x_order = sorted(list(fisher_df[x].unique()))

    fisher_table = pd.crosstab(fisher_df[x], fisher_df[y])
    p_value = fisher_exact(fisher_table)[1]
    for x_target in x_order:
        percent = fisher_table.loc[x_target][y_target] / fisher_table.loc[x_target].sum() * 100
        ax.bar(x=x_target, height=percent)
    ax.plot([0, 0, 1, 1], [102, 105, 105, 102], color='black')
    ax.text(x=0.5, y=105, s=f'{p_value:0.3f}', ha='center', va='bottom', color='red' if p_value < 0.05 else 'black')

    ax.set_ylim(0, 115)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(x_order)
    ax.set_title(x)
    ax.set_ylabel(f'{y}={y_target}(%)')
    return ax
