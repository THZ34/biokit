# coding

# 作者：唐宏振
# 邮箱：tanghongzhen34@gmail.com
# 这个脚本的主要目的是打包一些常用的函数和画图功能，简化操作

import copy
import math
from math import pi

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scipy
import scipy.stats
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test
from scipy.stats import f_oneway


def annular(df, groupby, order=None, colors=None):
    """环状统计图

    :param df:
    :param groupby:
    :param order:
    :param colors:
    :return:
    """
    if not order:
        order = list(pd.value_counts(df[groupby]).index)
    if not colors:
        colors = [None] * len(order)
    # matplotlib 极坐标系
    ax = plt.subplot(projection='polar')
    ax.grid(False)
    ax.set_thetagrids([])
    ax.set_rgrids([])
    ax.spines['polar'].set_visible(False)
    # 统计
    percent = [(df[df[groupby] == group].shape[0] / df.shape[0]) for group in order]
    width = [i * 2 * pi for i in percent]
    left = []
    i = 0
    for length in width:
        left.append(i)
        i += length
    # 环 + 百分比
    for group, color, i in zip(order, colors, range(len(order))):
        plt.barh(y=2.5, width=width[i], left=left[i], height=1, label=group,
                 color=color)
        plt.text(x=width[i] / 2 + left[i], y=2.5, s='%0.2f%%' % (percent[i] * 100), ha='center', va='center',
                 fontdict={'fontsize': 20, 'weight': 'bold', 'family': 'SimHei'})
    plt.ylim(0, 3)
    plt.legend(bbox_to_anchor=(0.2, 0.2, 1.1, 0.9))
    plt.tight_layout()
    return ax


# box/violin plot + t-test
# coding='utf-8'
def basic_analysis_pipeline(adata, name):
    """Basic single cell analysis and plot

    :param file:
    :param name:
    :param highest_n:
    :return:
    """
    # Basic statistics
    exp_df = adata.to_df()
    mito_genes = [i for i in exp_df.columns if 'MT-' in i]
    mito_df = exp_df[mito_genes]
    adata.obs['n_counts'] = exp_df.sum(1)
    adata.obs['mito_percent'] = mito_df[mito_genes].sum(1) / adata.obs['n_counts']
    adata.obs['n_genes'] = [(len(exp_df.columns) - pd.value_counts(exp_df.loc[sample]).get(0, 0)) for sample in
                            exp_df.index]
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'mito_percent'], jitter=0.4, multi_panel=True, save=' ' + name + '.pdf',
                 show=False)
    sc.pl.scatter(adata, x='n_counts', y='mito_percent', save=' count-mito_percent ' + name + '.pdf', show=False)
    sc.pl.scatter(adata, x='n_counts', y='n_genes', save=' count-gene ' + name + '.pdf', show=False)

    # quality control
    adata.qc = sc.pp.calculate_qc_metrics(adata)
    sc.pp.filter_cells(adata, 50)
    sc.pp.filter_genes(adata, 3)
    sc.pp.log1p(adata)
    adata = adata[adata.obs['n_genes'] < 4000, :]
    adata = adata[adata.obs['mito_percent'] < 0.3, :]

    # normalize and filter high variable gene
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.pl.highly_variable_genes(adata, save=' ' + name + '.pdf', show=False)
    adata = adata[:, adata.var['highly_variable']]

    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, show=False, save=' ' + name + '.pdf')
    sc.pl.pca_variance_ratio(adata, show=False, save=' ' + name + '.pdf')

    # reduce dimension
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
    # sc.tl.tsne(adata, random_state=10)
    sc.tl.umap(adata, random_state=10)

    # cluster
    sc.tl.leiden(adata)
    sc.tl.louvain(adata)
    adata.write(name + ' result.h5ad')


def hypergeometric_distribution_subtype(adata, cell_markers, groupby, top=20, p=0.05, other_cells='other cells'):
    """define cell subtype with cell markers

    :param adata: anndata class
    :param cell_markers:
    :param groupby:
    :param top: number of gene use as cluster marker
    :return:
    """
    marker_dict = {}
    with open(cell_markers) as f:
        for line in f:
            line = line.strip().split('\t')
            marker_dict[line[0]] = line[1:]
    df_dict = rank_reconstruct(adata, groupby, top)
    del df_dict['params']
    cluster_list = list(df_dict.keys())

    M = adata.shape[1]
    N = top
    subtype_matrix = []
    subtype_list = list(marker_dict.keys())
    for subtype in subtype_list:
        subtype_marker = marker_dict[subtype]
        n = len(subtype_marker)
        x = []
        rv = scipy.stats.hypergeom(M, n, N)
        for cluster in cluster_list:
            df = df_dict[cluster]
            cluster_marker = df.index
            k = len(set(subtype_marker) & set(cluster_marker))
            x.append(k)
        subtype_matrix.append(rv.pmf(x))
    subtype_df = pd.DataFrame(subtype_matrix, columns=cluster_list, index=subtype_list)
    adata.uns['hyge_subtype'] = subtype_df
    adata.obs['hyge_subtype'] = adata.obs[groupby]
    subtype_df.dropna(axis=1, how='all', inplace=True)
    for cluster in cluster_list:
        a = subtype_df.nsmallest(1, cluster)[cluster][0]
        if (cluster in subtype_df.columns) and (a < p):
            adata.obs['hyge_subtype'].replace(cluster, subtype_df.nsmallest(1, cluster).index[0], inplace=True)
        else:
            adata.obs['hyge_subtype'].replace(cluster, other_cells, inplace=True)


def expression_subtype(adata, cell_markers, groupby):
    marker_dict = {}
    with open(cell_markers) as f:
        for line in f:
            line = line.strip().split('\t')
            marker_dict[line[0]] = line[1:]

    subtype_list = list(marker_dict.keys())

    cluster_list = list(set(adata.obs[groupby]))
    cluster_list.sort()
    exp_df = adata.to_df()
    subexpression_df = pd.DataFrame(index=subtype_list, columns=cluster_list)
    for cluster in cluster_list:
        tmp_df_1 = exp_df.loc[adata.obs[groupby] == cluster]
        for subtype in subtype_list:
            genelist = marker_dict[subtype]
            genelist = list(set(genelist) & set(tmp_df_1.columns))
            tmp_df_2 = tmp_df_1[genelist]
            if tmp_df_2.empty:
                subexpression_df[cluster].loc[subtype] = 0
            else:
                subexpression_df[cluster].loc[subtype] = tmp_df_2.mean().mean()

    subexpression_df = pd.DataFrame(subexpression_df, dtype='float')
    subexpression_df.dropna(axis=1, how='all', inplace=True)
    adata.uns['expression_subtype'] = subexpression_df
    adata.obs['expression_subtype'] = adata.obs[groupby]

    for cluster in cluster_list:
        if (cluster in subexpression_df.columns) and subexpression_df[cluster].sum() != 0:
            adata.obs['expression_subtype'].replace(cluster, subexpression_df.nlargest(1, cluster).index[0],
                                                    inplace=True)
        else:
            adata.obs['expression_subtype'].replace(cluster, 'unrecognized cells', inplace=True)


def overlap_subtype(adata, cell_markers, groupby, top=20, rate=False, other_cells='other cells'):
    """define cell subtype with cell markers

    :param adata: anndata class
    :param cell_markers:
    :param groupby:
    :param top: number of gene use as cluster marker
    :return:
    """

    marker_dict = {}
    all_marker = []
    with open(cell_markers) as f:
        for line in f:
            line = line.strip().split('\t')
            marker_dict[line[0]] = line[1:]
            all_marker.extend(line[1:])

    df_dict = rank_reconstruct(adata, groupby, top)
    del df_dict['params']
    cluster_list = list(df_dict.keys())
    cell_list = list(marker_dict.keys())
    intersection_matrix = []

    for cluster in cluster_list:
        intersection_line = []
        cluster_marker = set(df_dict[cluster].index)
        for cell in cell_list:
            cell_marker = set(marker_dict[cell])
            intersection_line.append(len(cluster_marker & cell_marker))
        intersection_matrix.append(intersection_line)
    intersection_df = pd.DataFrame(intersection_matrix, index=cluster_list, columns=cell_list)
    adata.uns['overlap_subtype'] = intersection_df
    adata.obs['overlap_subtype'] = adata.obs[groupby]
    if not rate:
        for cluster in cluster_list:
            line = intersection_df.loc[cluster]
            if line.sum() != 0:
                cell_type = []
                for cell in cell_list:
                    if line[cell] != 0:
                        cell_type.append(cell)

                cell_type = ' + '.join(cell_type)
                adata.obs['overlap_subtype'].replace(cluster, cell_type, inplace=True)
            else:
                adata.obs['overlap_subtype'].replace(cluster, other_cells, inplace=True)

    else:
        for cluster in cluster_list:
            line = intersection_df.loc[cluster]
            if line.sum() != 0:
                total = line.sum()
                cell_type = []
                for cell in cell_list:
                    if line[cell] != 0:
                        cell_type.append(str(line[cell] / total) + ' ' + cell)

                cell_type = ' + '.join(cell_type)
                adata.obs['overlap_subtype'].replace(cluster, cell_type, inplace=True)
            else:
                adata.obs['overlap_subtype'].replace(cluster, other_cells, inplace=True)


def color_piecewise(adata, var_names, cmap, lower=0, upper=None, coor_name='X_umap', s=10, alpha=1,
                    ground=(0.898438, 0.898438, 0.898438, 1), surface=(0, 0, 0, 1), show_limit=False, figsize=(6, 4)):
    # init
    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    else:
        pass
    plt.figure(figsize=figsize)
    coor_df = pd.DataFrame(adata.obsm[coor_name], index=adata.obs_names, columns=['x', 'y'])
    # var filter
    not_in_var_obs = list(set(var_names) - (set(adata.var_names) ^ set(adata.obs.columns)))
    not_in_var_obs.sort()
    print(str(not_in_var_obs) + ' not in var_names or obs.columns')
    genes = list(set(var_names) & set(adata.var_names))
    # create plot data frame
    var_names = list(set(var_names) & (set(adata.var_names) ^ set(adata.obs.columns)))
    plot_shape = math.ceil(math.sqrt(len(var_names)))
    exp_df = pd.concat([adata.obs, adata.to_df()[genes]], axis=1, sort=False)
    if plot_shape * (plot_shape - 1) > len(var_names):
        plot_shape_x = plot_shape
        plot_shape_y = plot_shape - 1
    else:
        plot_shape_x = plot_shape
        plot_shape_y = plot_shape

    for i in range(len(var_names)):
        # subplot
        plt.subplot(plot_shape_x, plot_shape_y, i + 1)
        gene = var_names[i]
        color_list = cmap(exp_df[gene])
        for j in range(len(exp_df[gene])):
            if exp_df[gene][j] <= lower:
                color_list[j] = ground
        if upper:
            for j in range(len(exp_df[gene])):
                if exp_df[gene][j] >= upper:
                    color_list[j] = surface
        color_df = pd.DataFrame(color_list, index=exp_df.index, columns=['r', 'g', 'b', 'a'])
        plot_df = pd.concat([exp_df[genes], adata.obs, color_df, coor_df], axis=1, sort=False)
        plot_df.sort_values(by=gene, inplace=True)
        plt.scatter(plot_df['x'], plot_df['y'], c=plot_df[['r', 'g', 'b', 'a']].to_numpy(), s=s, alpha=alpha)
        if show_limit:
            plt.title('%s  limit:%s~%s' % (gene, lower, upper))
        else:
            plt.title(gene)
        plt.xticks([])
        plt.yticks([])


def volcano(rank_df, key_genes=(), s=1, p=0.05, lfc=0.5, log_p=False, lfc_limit=None, bar_len=1, angle=45,
            key_gene_filter=False):
    """volcano plot

    :param rank_df: gene rank dataframe, at least contains p-value and log fold change
    :param key_genes: importent genes need a label in plot
    :param s: size of point
    :param p: cutoff of p-values
    :param lfc: cutoff of log fold change
    :param log_p: None or a real number
    :return: matplotlib axes
    """
    from math import log, sin, cos, pi
    df = rank_df.copy()
    if lfc_limit:
        df = df[df['logfoldchanges'] > lfc_limit[0]]
        df = df[df['logfoldchanges'] < lfc_limit[1]]
    df[['pvals', 'logfoldchanges']] = pd.DataFrame(df[['pvals', 'logfoldchanges']], dtype='float')
    # classification
    upregulate_genes = df[(df['pvals'] < p) & (df['logfoldchanges'] > lfc)].index
    downregulate_genes = df[(df['pvals'] < p) & (df['logfoldchanges'] < -lfc)].index
    other_genes = df[~df.index.isin(upregulate_genes ^ downregulate_genes)].index
    other_genes = list(set(df.index) - set(upregulate_genes) - set(downregulate_genes))
    p_min = df[df['pvals'] != 0]['pvals'].min()
    # log
    if log_p:
        df['pvals'] = [-log(i, log_p) if i != 0 else -log(p_min, log_p) for i in df['pvals']]
        p = -log(p, log_p)
    # plot
    plt.plot([df['logfoldchanges'].min(), df['logfoldchanges'].max()], [p, p], linestyle='--', c='black')
    plt.plot([lfc, lfc], [p, df['pvals'].max()], linestyle='--', c='black')
    plt.plot([-lfc, -lfc], [p, df['pvals'].max()], linestyle='--', c='black')
    plt.scatter(x=df.loc[upregulate_genes]['logfoldchanges'], y=df.loc[upregulate_genes]['pvals'], s=s, c='red',
                label='up-regulated genes')
    plt.scatter(x=df.loc[downregulate_genes]['logfoldchanges'], y=df.loc[downregulate_genes]['pvals'], s=s, c='blue',
                label='down-regulated genes')
    plt.scatter(x=df.loc[other_genes]['logfoldchanges'], y=df.loc[other_genes]['pvals'], s=s / 10, c='black')
    plt.legend(bbox_to_anchor=(0.5, 0.5, 0.5, 0.5))
    plt.xlabel('log2 (fold change)')
    plt.ylabel('log' + str(log_p) + ' (p-value)')
    # plt.grid()
    # key
    if key_genes:
        if key_gene_filter:
            key_genes = [i for i in key_genes if i in upregulate_genes]
            key_genes.extend([i for i in key_genes if i in downregulate_genes])
        for gene in key_genes:
            if gene in df.index:
                if df.loc[gene]['logfoldchanges'] > 0:
                    plt.text(x=df.loc[gene]['logfoldchanges'] + bar_len * cos(angle * pi / 360),
                             y=df.loc[gene]['pvals'] + bar_len * sin(angle * pi / 360), s=gene,
                             horizontalalignment='left')
                    plt.plot([df.loc[gene]['logfoldchanges'],
                              df.loc[gene]['logfoldchanges'] + bar_len * cos(angle * pi / 360)],
                             [df.loc[gene]['pvals'], df.loc[gene]['pvals'] + bar_len * sin(angle * pi / 360)],
                             linestyle='solid', c='black')
                if df.loc[gene]['logfoldchanges'] < 0:
                    plt.text(x=df.loc[gene]['logfoldchanges'] + bar_len * cos(angle * pi / 360),
                             y=df.loc[gene]['pvals'] + bar_len * sin(angle * pi / 360), s=gene,
                             horizontalalignment='right')
                    plt.plot([df.loc[gene]['logfoldchanges'],
                              df.loc[gene]['logfoldchanges'] - bar_len * cos(angle * pi / 360)],
                             [df.loc[gene]['pvals'], df.loc[gene]['pvals'] + bar_len * sin(angle * pi / 360)],
                             linestyle='solid', c='black')


def hyper(exp_df, gene_1, gene_2):
    M = exp_df.shape[0]
    N = exp_df[exp_df[gene_1] > 0].shape[0]
    n = exp_df[exp_df[gene_2] > 0].shape[0]
    k = exp_df[exp_df[gene_1] > 0][exp_df[gene_2] > 0].shape[0]
    rv = scipy.stats.hypergeom(M, N, n)
    p = rv.pmf(k)
    return p, N - k, n - k, k


def exp_box(data, var, groupby, figsize=(4, 8), order=None,
            cutoff={0.05: '*', 0.01: '**', 0.001: '***', 0.0001: '****'}, test='t', no_ns=False, kind='box', ax=None):
    """

    :param data: pandas DataFrame
    :param var: variation name
    :param groupby: x axis
    :param figsize: tuple e.p. (4,8)
    :param order: sort x
    :param cutoff: dictory of cutoff and label
    :return: axes
    """
    if not order:
        order = sorted(list(set(data[groupby])))

    # box and strip (violin and strip)
    plt.figure(figsize=figsize)
    if kind == 'box':
        sns.boxplot(x=groupby, y=var, data=data, color='grey', width=0.7, order=order, whis=100)
    elif kind == 'violin':
        sns.violinplot(x=groupby, y=var, data=data, color='grey', width=0.7, order=order, whis=100)
    else:
        print('Support only box or violin')
    sns.stripplot(x=groupby, y=var, data=data, color='black', size=1, jitter=0.4, order=order)
    # mean line
    for i in range(len(order)):
        group = order[i]
        mean = data[data[groupby] == group][var].mean()
        plt.plot([i + 0.3, i - 0.3], [mean, mean], c='red', label='mean value', linewidth=3)
    plt.xlim(-0.5, len(order) - 0.5)
    # plt.xticks(rotation=30)
    plt.title(var, fontsize=30)
    plt.ylabel('')
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
        f = scipy.stats.levene(list(values_i), list(values_j))  # test for homogeneity of variance
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
            plt.plot([obs_i, obs_j], [height * max, height * max], c='black')
            plt.plot([obs_i, obs_i], [(height - 0.02) * max, height * max], c='black')
            plt.plot([obs_j, obs_j], [(height - 0.02) * max, height * max], c='black')
            plt.text(x=(i + j) / 2, y=(height + 0.01) * max, s=text, horizontalalignment='center',
                     fontweight='bold',
                     fontsize=16)
            height += 0.1
        if no_ns == True:
            plt.ylim(0, data[var].max() * (
                    1.2 + 0.1 * len(order) * (len(order) - 1) / 2 - 0.1 * pd.value_counts(text).get('ns', 0)))
        else:
            plt.ylim(0, data[var].max() * (1.2 + 0.1 * len(order) * (len(order) - 1) / 2))
    plt.tight_layout()
    return plt.gcf()


def kaplan_meier(grouped_df, groupby, ax=None, order=None):
    kmf = KaplanMeierFitter()
    te_list = []
    if not ax:
        ax = plt.subplot()
    if not order:
        group_list = list(set(grouped_df[groupby]))
        group_list.sort()
    else:
        group_list = order

    for group in group_list:
        T, E = grouped_df[grouped_df[groupby] == group][['time', 'status']].T.to_numpy()
        kmf.fit(T, E, label=group)
        kmf.plot(ax=ax)

    if len(te_list) == 2:
        statistical_parameter = logrank_test(te_list[0][0], te_list[1][0], te_list[0][1], te_list[1][1], alpha=.99)
        ax.text(s='p-value:%s' % round(statistical_parameter.p_value, 4), x=grouped_df['time'].max() * 0.75,
                y=0.8)
    return ax


def kaplan_meier(grouped_df, groupby, time='time', status='status', groups=None, ax=None, cox=False,
                 table=False, pathway=None):
    """

    :param grouped_df: pandas.DataFrame, contain at least three columns(time,status,group)
    :param groupby: column name of group, necessary
    :param time: column name of time, default: 'time'
    :param status: column name of status, default: 'status'
    :param groups: groups and order
    :param ax: matplotlib.axes
    :param cox: cox-regression, default: False
    :param cox_event: group of events occurred
    :param table: plot number of risk
    :return: ax
    """
    kmf = KaplanMeierFitter()
    if not groups:
        groups = list(set(grouped_df[groupby]))
        groups.sort()
        if len(groups) <= 1:
            raise ValueError('Only one group')
    if not ax:
        ax = plt.subplot()

    matrix = []
    kmf_dict = {}
    # kaplan-meier curve
    for group in groups:
        T, E = grouped_df[grouped_df[groupby] == group][[time, status]].to_numpy().T
        kmf.fit(T, E)
        kmf_dict[group] = copy.copy(kmf)
        kmf.plot(ax=ax, ci_show=False)
        os_median = kmf._median
        upper = os_median + 1.96 * T.std() / math.sqrt(len(T))
        lower = os_median - 1.96 * T.std() / math.sqrt(len(T))
        median_ci95 = f'{os_median:4.1f}({lower:4.1f}-{upper:4.1f})'
        matrix.append([groupby, group, len(E), E[E == 1].shape[0], median_ci95])

    str_event_to_number = dict(zip(groups, range(len(groups))))
    cox_df = cox(grouped_df.replace(str_event_to_number), [groups[0]])
    statistical_parameter = multivariate_logrank_test(grouped_df['time'],
                                                      [str_event_to_number[i] for i in grouped_df[groupby]],
                                                      grouped_df['status'])
    # set plot
    ax.set_xlabel('')
    ax.set_xlim(0, ax.get_xticks()[-1])

    ax.set_ylabel('Survival Probability (%)')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels([str(int(i * 100)) for i in ax.get_yticks()])
    ax.set_ylim(0, 1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().set_visible(False)
    # cox
    raw_xticks = ax.get_xticks()
    xticks = ax.get_xticks()[1:-1]
    cox_table = [['', groups[0], groups[1]],
                 ['OS median', matrix[0][-1].split('(')[0], matrix[1][-1].split('(')[0]],
                 ['HR(95CI)', f"{cox_df['p-value'][pathway]:0.3f}", cox_df['HR(CI95)'][pathway]]]
    scale = xticks[-1] * 0.4 / xticks[-1]
    start = xticks[-1] * 0.6 / ax.get_xticks()[-1]
    cox_table_plot = ax.table(cox_table, cellLoc='center', bbox=(start, 0.80, scale, 0.21), edges='open')
    cox_table_plot.set_fontsize(10)
    end = ax.get_xticks()[-1] * 6 / 7
    ax.plot([xticks[-1] * 0.6, end], [0.8, 0.8], color='black', linewidth=1)
    ax.plot([xticks[-1] * 0.6, end], [0.94, 0.94], color='black', linewidth=1)
    # table number of risk
    time_number_table = []
    for group in groups:
        kmf = kmf_dict[group]
        number_of_sample = len(kmf.event_observed)
        time_number_table.append(
            (round(number_of_sample * (1 - kmf.cumulative_density_at_times(xticks)))).astype(dtype='int').to_list())
        plt.text(s=group, x=xticks[-1], y=(1 - kmf.cumulative_density_at_times(xticks[-1])), va='bottom', ha='right')

    xticks = ax.get_xticks()[1:-1]
    width = (xticks[-1] - xticks[0]) / raw_xticks[-1] * (len(xticks)) / (len(xticks) - 1)
    start = xticks[0] / ax.get_xticks()[-1] - (width / (len(xticks) * 2))
    ax.text(x=ax.get_xticks().mean(), y=-0.2, s='Number at Risk', fontsize=10, ha='center')
    bbox = (start, -0.25 - 0.1 * len(groups), width, 0.1 * len(groups))
    ax.table(time_number_table, cellLoc='center', bbox=bbox,
             edges='open', rowLabels=['mutated', 'nonmutation'])
    # table_list = [
    #     ax.table([cox_table[0]], cellLoc='center', bbox=(start, 0.90, scale, 0.05)),
    #     ax.table([cox_table[1]], cellLoc='center', bbox=(start, 0.85, scale, 0.05)),
    #     ax.table([cox_table[2]], cellLoc='center', bbox=(start, 0.80, scale, 0.05))]
    # table
    # title color
    if statistical_parameter.p_value < 0.06:
        title_color = 'red'
    else:
        title_color = 'black'
    ax.set_title(f'{groupby} p-value:{statistical_parameter.p_value:3.5f}', fontsize=10, color=title_color)

    # if cox and len(groups) == 2:
    #     cox_df = cox(grouped_df, [pathway])
    #     cox_table =
    ax.set_xticks(xticks)
    ax.set_xlim(raw_xticks[0], raw_xticks[-1])
    plt.tight_layout()
    return matrix, ax
