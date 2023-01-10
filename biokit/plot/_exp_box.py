import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import f_oneway
from sympy.physics.quantum.identitysearch import scipy


def exp_box(data, var, groupby, order=None, cutoff=None, test='t', no_ns=False, kind='box', ax=None,
            color_dict=None, meanline=False, cutoff_color=None):
    """

    :param meanline:
    :param color_dict:
    :param data: pandas DataFrame
    :param var: variation name
    :param groupby: x axis
    :param figsize: tuple e.p. (4,8)
    :param order: sort x
    :param cutoff: dictory of cutoff and label
    :return: axes
    """
    if not cutoff:
        cutoff = {0.05: '*', 0.01: '**', 0.001: '***', 0.0001: '****'}
    if not cutoff_color:
        cutoff_color = dict(
            zip(['*', '**', '***', '****', 'ns'], ['orange', 'darkorange', 'orangered', 'red', 'deepskyblue']))
    if not order:
        order = sorted(list(set(data[groupby])))

    # box and strip (violin and strip)
    if kind == 'box':
        sns.boxplot(x=groupby, y=var, data=data, width=0.7, order=order, ax=ax)
    elif kind == 'violin':
        sns.violinplot(x=groupby, y=var, data=data, width=0.7, order=order, ax=ax)
    elif kind == 'boxen':
        sns.boxenplot(x=groupby, y=var, data=data, width=0.7, order=order, ax=ax)
    else:
        print('Support only box, violin and boxen')

    sns.stripplot(x=groupby, y=var, data=data, color='black', size=1, jitter=0.4, order=order, ax=ax)
    # mean line
    if meanline:
        for i in range(len(order)):
            group = order[i]
            mean = data[data[groupby] == group][var].mean()
            ax.plot([i + 0.3, i - 0.3], [mean, mean], c='red', label='mean value', linewidth=3)

    ax.set_xlim(-0.5, len(order) - 0.5)
    # plt.xticks(rotation=30)
    ax.set_title(var, fontsize=30)
    ax.set_ylabel('')
    test_df = pd.DataFrame(columns=['pval', 'stat'])
    # define Statistical function
    if test == 't':
        def tmp_func(value_i, value_j):
            f = scipy.stats.levene(list(values_i), list(values_j))  # test for homogeneity of variance
            if f.pvalue <= 0.05:
                ttest = scipy.stats.ttest_ind(list(values_i), list(values_j), equal_var=False)
            else:
                ttest = scipy.stats.ttest_ind(list(values_i), list(values_j), equal_var=True)
            return ttest

        test_func = tmp_func
    elif test == 'anova':
        test_func = f_oneway
    else:
        raise ValueError('t or anova')
    # Statistical significance
    cutoff_list = list(cutoff.keys())
    cutoff_list.sort()
    iter_list = []
    print(cutoff_list)
    for span in range(1, len(order)):
        for i in range(len(order) - span):
            iter_list.append((i, i + span))
    height = 1.2  # height of test label
    text_list = []
    for i, j in iter_list:
        obs_i = order[i]
        obs_j = order[j]
        values_i = data[data[groupby] == obs_i][var]
        values_j = data[data[groupby] == obs_j][var]
        test = test_func(values_i, values_j)
        print(test)
        # significant
        max = data[var].max()
        if set(values_i) == set(values_j):
            text = 'ns'
        else:
            for c in cutoff_list:
                if test.pvalue <= c:
                    text = cutoff[c]
                    break
                else:
                    text = 'ns'
        text_list.append(text)
        if text == 'ns' and no_ns == True:
            pass
        else:
            ax.plot([obs_i, obs_j], [height * max, height * max], c='black')
            ax.plot([obs_i, obs_i], [(height - 0.02) * max, height * max], c='black')
            ax.plot([obs_j, obs_j], [(height - 0.02) * max, height * max], c='black')
            ax.text(x=(i + j) / 2, y=(height + 0.01) * max, s=text, horizontalalignment='center',
                    fontweight='bold', fontsize=16, c=cutoff_color[text])
            height += 0.1
        if no_ns:
            ax.set_ylim(0, data[var].max() * (
                    1.2 + 0.1 * len(order) * (len(order) - 1) / 2 - 0.1 * pd.value_counts(text).get('ns', 0)))
        else:
            ax.set_ylim(0, data[var].max() * (1.2 + 0.1 * len(order) * (len(order) - 1) / 2))
    plt.tight_layout()
    return plt.gcf()
