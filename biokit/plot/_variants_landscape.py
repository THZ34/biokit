import math

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy.sparse import coo_matrix
from scipy.stats import pearsonr


def variants_landscape_mutation(mutation_file, clinical_file, save=True, n_genes=None, sort_gene=True):
    """展示突变和临床信息的复合热图

    :param sort_gene: 按突变数量排序
    :param n_genes: 热图基因数量
    :param mutation_file: 突变矩阵，列：基因，行：样本，值：突变类型
    :param clinical_file: 临床信息
    :param save: 是否保存pdf，默认保存为landscape.pdf,如果不保存，则返回plt.figure
    :return: None
    """
    # 修改字体,子图间距
    plt.subplots_adjust(left=0.1, bottom=0, top=1, right=0.8, hspace=3, wspace=0)
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    # 读入数据
    if isinstance(mutation_file, pd.DataFrame):
        mut_df = mutation_file
    else:
        mut_df = pd.read_csv(mutation_file, sep='\t', index_col=0)
        mut_df.replace({np.nan: None}, inplace=True)
    if isinstance(clinical_file, pd.DataFrame):
        sample_info = clinical_file
    else:
        sample_info = pd.read_csv(clinical_file, sep='\t', index_col=0)

    # 对样本和基因进行排序,过滤
    mut_df['sum'] = mut_df.astype(bool).sum(1)
    mut_df.loc['sum'] = mut_df.astype(bool).sum()
    if sort_gene:
        mut_df.sort_values('sum', axis=0, inplace=True, ascending=True)
    # mut_df.sort_values('sum', axis=1, inplace=True, ascending=False)
    mut_df.drop('sum', axis=0, inplace=True)
    mut_df.drop('sum', axis=1, inplace=True)
    # sample_info = sample_info.loc[mut_df.columns]
    if n_genes:
        mut_df = mut_df.iloc[-n_genes:, :]
    mut_df = mut_df[sample_info.index]
    mut_df.fillna('None', inplace=True)
    # 区分连续数据和离散数据
    continuous_columns = sample_info.dtypes[~(sample_info.dtypes == 'object')].index
    discrete_columns = sample_info.dtypes[sample_info.dtypes == 'object'].index
    # 定义连续、离散变量和突变的colormap
    continuous_colors = dict(zip(continuous_columns, sns.husl_palette(len(continuous_columns))))
    discrete_colors = {}
    for column in discrete_columns:
        values = sorted(list(set(sample_info[column])))
        discrete_colors[column] = dict(zip(values, sns.hls_palette(len(values))))
    mut_types = []
    for column in mut_df:
        for value in mut_df[column]:
            if value:
                mut_types.extend(value.split(';'))
    mut_types = sorted(list(set(mut_types)))
    mut_colors = dict(zip(mut_types, sns.hls_palette(len(mut_types))))
    mut_colors['None'] = 'grey'
    # 原始数据备份
    raw_data = {'mut': mut_df.copy(), 'info': sample_info.copy()}
    # 画图 - 定义画板和子图
    figsize = (sample_info.shape[0] / 4, sample_info.shape[1] + mut_df.shape[0] / 4)
    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(int(sample_info.shape[1] + mut_df.shape[0]), 1)
    # 画图 - 离散变量
    series = 0
    sample_info.index = range(sample_info.shape[0])
    for column in discrete_columns:
        ax = plt.subplot(grid[series, 0])
        for subtype in discrete_colors[column]:
            tmp_df = sample_info[sample_info[column] == subtype]
            ax.bar(x=tmp_df.index, height=[1] * tmp_df.shape[0], color=discrete_colors[column][subtype], width=1,
                   label=subtype)
            # ax.axis('off')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_ylabel(column, rotation=0, ha='right', fontsize=10)
            ax.set_xlim(-0.8, sample_info.shape[0])
        ax.legend(title=column, bbox_to_anchor=(1, 1))
        series += 1
    # 画图 - 连续变量
    for column in continuous_columns:
        ax = plt.subplot(grid[series:series + 2, 0])
        ax.bar(x=sample_info.index, height=sample_info[column], color=continuous_colors[column])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_ylabel(column, rotation=0, ha='right', fontsize=10)
        ax.set_xlim(-0.8, sample_info.shape[0])
        series += 2
    # 画图 - 变异
    mut_df.index = range(mut_df.shape[0])
    mut_df.columns = range(mut_df.shape[1])
    ax = plt.subplot(grid[series:, 0])
    heatmap_df = []
    for sample in mut_df.columns:
        for gene in mut_df.index:
            muts = mut_df[sample][gene]
            x = sample
            bottom = gene

            if type(muts) == str:
                muts = muts.split(';')
                height = 0.7 / len(muts)
                for mut in muts:
                    heatmap_df.append([x, height, bottom, mut])
                    bottom += height
            else:
                height = 0.7
                heatmap_df.append([x, height, bottom, 'None'])

    heatmap_df = pd.DataFrame(heatmap_df)
    for mut in mut_colors:
        color = mut_colors[mut]
        tmp_df = heatmap_df[heatmap_df[3] == mut]
        ax.bar(x=tmp_df[0], height=tmp_df[1], bottom=tmp_df[2], color=color, label=mut, edgecolor='grey')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_ylabel('Gene', rotation=0, ha='right', fontsize=10)
    ax.set_xlim(-0.8, sample_info.shape[0])
    ax.set_ylim(0, mut_df.shape[0])
    ax.set_yticks(mut_df.index)
    ax.set_yticklabels(raw_data['mut'].index)
    ax.legend(title='Mutation Type', bbox_to_anchor=(1, 0.98))
    # 保存
    if save:
        if type(save) == str:
            pdf = PdfPages(f'{save}.pdf')
        else:
            pdf = PdfPages(f'landscape.pdf')
        pdf.savefig()
        pdf.close()
    else:
        return fig


def variants_landscape_compare(mut_df1: [pd.DataFrame, str], mut_df2: [pd.DataFrame, str],
                               sample_info: [pd.DataFrame, str], save=True, n_genes=None, marker=('1', '2'),
                               sort=False):
    """展示突变和临床信息的复合热图

    :param sort: 基因排序方式 : 'sum':突变数量,'corr':组间相关性,'score':相关性*数量
    :param sample_info:
    :param marker: 区分同一个基因的标记
    :param mut_df1: 突变矩阵1
    :param mut_df2: 突变矩阵2
    :param n_genes: 基因数量
    :param save: 是否保存pdf，默认保存为landscape.pdf,如果不保存，则返回plt.figure
    :return: None
    """
    if sort in ['sum', 'corr', 'score', False]:
        pass
    elif sort:
        sort = 'score'
    else:
        raise ValueError('unexpected value of sort')

    # 修改字体
    plt.subplots_adjust(left=0.1, bottom=0, top=1, right=0.8, hspace=3, wspace=0)
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    # 读入数据
    if type(sample_info) == str:
        sample_info = pd.read_csv(sample_info, sep='\t', index_col=0)
    if type(mut_df1) == str:
        mut_df1 = pd.read_csv(mut_df1, sep='\t', index_col=0)
    if type(mut_df2) == str:
        mut_df2 = pd.read_csv(mut_df2, sep='\t', index_col=0)

    mut_df1.replace({np.nan: None}, inplace=True)
    mut_df2.replace({np.nan: None}, inplace=True)
    # 对样本和基因进行排序,过滤
    genes = mut_df1.index.intersection(mut_df2.index)
    mut_df1 = mut_df1.loc[genes]
    mut_df2 = mut_df2.loc[genes]
    corr_df = pd.DataFrame(
        [pearsonr(mut_df1.loc[gene].astype(bool).astype(int), mut_df2.loc[gene].astype(bool).astype(int))[0] *
         (1 - pearsonr(mut_df1.loc[gene].astype(bool).astype(int), mut_df2.loc[gene].astype(bool).astype(int))[0])
         for gene in genes], index=genes)

    mut_df1['sum'] = mut_df1.astype(bool).sum(1)
    mut_df2['sum'] = mut_df2.astype(bool).sum(1)
    mut_df1['corr'] = corr_df
    mut_df2['corr'] = corr_df
    mut_df1['score'] = mut_df1['sum'] * mut_df1['corr'] + mut_df2['sum'] * mut_df2['corr']
    mut_df2['score'] = mut_df1['sum'] * mut_df1['corr'] + mut_df2['sum'] * mut_df2['corr']
    if sort:
        mut_df1.sort_values(by=sort, inplace=True)
        mut_df2.sort_values(by=sort, inplace=True)
    if n_genes:
        mut_df1 = mut_df1.iloc[-n_genes:, :]
        mut_df2 = mut_df2.iloc[-n_genes:, :]
    genes = list(mut_df1.index)
    mut_df1.drop(['sum', 'corr', 'score'], axis=1, inplace=True)
    mut_df2.drop(['sum', 'corr', 'score'], axis=1, inplace=True)
    mut_df1 = mut_df1[sample_info.index]
    mut_df2 = mut_df2[sample_info.index]
    mut_df1.index = mut_df1.index + '_' + marker[0]
    mut_df2.index = mut_df2.index + '_' + marker[1]
    # 区分连续数据和离散数据
    continuous_columns = sample_info.dtypes[~(sample_info.dtypes == 'object')].index
    discrete_columns = sample_info.dtypes[sample_info.dtypes == 'object'].index
    # 定义连续、离散变量和突基因的colormap
    continuous_colors = dict(zip(continuous_columns, sns.husl_palette(len(continuous_columns))))
    discrete_colors = {}
    for column in discrete_columns:
        values = sorted(list(set(sample_info[column])))
        discrete_colors[column] = dict(zip(values, sns.hls_palette(len(values))))

    gene_colors = {0: 'deepskyblue', 1: 'hotpink'}
    # 画图 -  定义画板和子图
    figsize = (sample_info.shape[0] / 4, sample_info.shape[1] + len(genes) / 2)
    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(int(sample_info.shape[1] * 2 + len(genes)), 1)
    # 画图 - 离散变量
    series = 0
    sample_info.index = range(sample_info.shape[0])
    for column in discrete_columns:
        ax = plt.subplot(grid[series, 0])
        for subtype in discrete_colors[column]:
            tmp_df = sample_info[sample_info[column] == subtype]
            ax.bar(x=tmp_df.index, height=[1] * tmp_df.shape[0], color=discrete_colors[column][subtype], width=1,
                   label=subtype)
            # ax.axis('off')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_ylabel(column, rotation=0, ha='right', fontsize=10)
            ax.set_xlim(-0.6, sample_info.shape[0] - 0.4)
        ax.legend(title=column, bbox_to_anchor=(1, 1))
        series += 1
    # 画图 - 连续变量
    for column in continuous_columns:
        ax = plt.subplot(grid[series:series + 2, 0])
        ax.bar(x=sample_info.index, height=sample_info[column], color=continuous_colors[column])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_ylabel(column, rotation=0, ha='right', fontsize=10)
        ax.set_xlim(-0.6, sample_info.shape[0] - 0.4)
        series += 2
    # 画图 - 变异
    mut_df1.index = range(mut_df1.shape[0])
    mut_df1.columns = range(mut_df1.shape[1])
    mut_df2.index = range(mut_df2.shape[0])
    mut_df2.columns = range(mut_df2.shape[1])
    ax = plt.subplot(grid[series:, 0])
    heatmap_df = []
    for sample in sample_info.index:
        for gene_index in mut_df1.index:
            muts_1 = mut_df1[sample][gene_index]
            muts_2 = mut_df2[sample][gene_index]
            color = gene_colors[gene_index % 2]
            x = sample
            height = 0.6
            bottom1 = gene_index * 2 + 0.8
            bottom2 = gene_index * 2
            color1 = 'grey' if muts_1 == 0 else color
            color2 = 'grey' if muts_2 == 0 else color
            label1 = 'None' if muts_1 == 0 else 'Mutated'
            label2 = 'None' if muts_2 == 0 else 'Mutated'
            heatmap_df.extend([[x, height, bottom1, color1, label1], [x, height, bottom2, color2, label2]])

    heatmap_df = pd.DataFrame(heatmap_df)
    for label in ['None', 'Mutated']:
        tmp_df = heatmap_df[heatmap_df[4] == label]
        ax.bar(x=tmp_df[0], height=tmp_df[1], bottom=tmp_df[2], color=tmp_df[3], label=label)
    # for bottom in range(len(genes) - 1):
    #     bottom = bottom * 2 + 1.7
    #     plt.plot([-0.6, sample_info.shape[0] - 0.4], [bottom, bottom], color='grey')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_ylabel('Gene', rotation=0, ha='right', fontsize=10)
    ax.set_xlim(-0.6, sample_info.shape[0] - 0.4)
    ax.set_ylim(0, len(genes) * 2)
    yticks = [i * 2 + 0.7 for i in range(len(genes))]
    # yticklabels = sum([[gene + '_' + marker[1], gene + '_' + marker[0]] for gene in genes], [])
    ax.set_yticks(yticks)
    ax.set_yticklabels(genes)
    ax.legend(title='Mutation Type', bbox_to_anchor=(1, 0.98))
    # 保存
    if save:
        if type(save) == str:
            pdf = PdfPages(f'{save}.pdf')
        else:
            pdf = PdfPages(f'landscape.pdf')
        pdf.savefig()
        pdf.close()
    else:
        return fig


# %% oncoplot
def oncoplot(mutations, sample_info, figsize=None, color_dict=None, discrete_colors=None, variants=None,
             fraction_lim=None, info_loc=None, fraction_annot=False, legend_cc=1):
    """oncoplot 瀑布图


    :param info_loc: 临床信息位置
    :param variants:
    :param fraction_lim:
    :param legend_cc:
    :param fraction_annot:
    :param discrete_colors:
    :param mutations: 突变矩阵，使用biokit.preprocessing.read_aachange()读取
    :param sample_info:样本临床信息
    :param figsize:
    :param color_dict:
    :return:
    """
    # 矩阵预处理，生成heatmap稀疏矩阵，提取样本名，基因名和突变类型
    mutations = mutations.copy()
    sample_info = sample_info.copy()

    xticklabels = mutations.columns
    yticklabels = mutations.index
    mutations.index = range(mutations.shape[0])
    mutations.columns = range(mutations.shape[1])
    heatmap_sparse_df = coo_matrix(mutations)
    heatmap_sparse_df = pd.DataFrame([heatmap_sparse_df.row, heatmap_sparse_df.col, heatmap_sparse_df.data]).T
    heatmap_sparse_df.columns = ['y', 'x', 'variant']

    if not variants:
        variants = list(heatmap_sparse_df['variant'].unique())
        # 重新排序,将multi_hits和no mutate排到最后
        for variant in ['multi_hits', 'no mutate']:
            if variant in variants:
                variants.remove(variant)
                variants.append(variant)

    sample_info = sample_info.loc[xticklabels]
    continuous_columns = sample_info.dtypes[~(sample_info.dtypes == 'object')].index
    discrete_columns = sample_info.dtypes[sample_info.dtypes == 'object'].index

    # 临床信息的位置
    if not info_loc:
        info_loc = {'upper': [], 'bottom': []}

    # 如果列名不在info_loc中，则加入info_loc['upper']
    for col in sample_info.columns:
        if col not in info_loc['upper'] and col not in info_loc['bottom']:
            info_loc['upper'].append(col)

    # 画图前准备
    if not figsize:
        figsize = (mutations.shape[1] + 3, len(discrete_columns) + len(continuous_columns) * 2 + len(yticklabels) + 3)
        figsize = (figsize[0] / 4, figsize[1] / 4)
    fig = plt.figure(figsize=figsize)
    grid_rows = len(discrete_columns) * 2 + len(continuous_columns) * 3 + len(yticklabels) + 3 + 1
    grid_cols = mutations.shape[1] + 3
    grid = plt.GridSpec(grid_rows, grid_cols)

    # 确定各个子图部分的位置，保存到字典
    top_of_ax = 0
    ax_dict = {'upper': {}, 'bottom': {}}
    # 热图上方的临床信息
    for col in info_loc['upper']:
        if col in discrete_columns:
            ax = plt.subplot(grid[top_of_ax, :mutations.shape[1]])
            top_of_ax += 2
        elif col in continuous_columns:
            ax = plt.subplot(grid[top_of_ax:top_of_ax + 2, :mutations.shape[1]])
            top_of_ax += 3
        ax_dict['upper'][col] = ax

    # 热图及突变比例
    ax_dict['heatmap'] = plt.subplot(grid[top_of_ax + 3:, :mutations.shape[1]])
    if fraction_annot:
        ax_dict['mut_stat_sample'] = plt.subplot(grid[top_of_ax:top_of_ax + 2, :mutations.shape[1]])
        ax_dict['mut_stat_gene'] = plt.subplot(grid[top_of_ax + 3:, mutations.shape[1] + 1:])
    else:
        ax_dict['mut_stat_sample'] = plt.subplot(grid[top_of_ax:top_of_ax + 3, :mutations.shape[1]])
        ax_dict['mut_stat_gene'] = plt.subplot(grid[top_of_ax + 3:, mutations.shape[1]:])

    # 热图下方的临床信息
    for col in info_loc['bottom']:
        if col in discrete_columns:
            ax = plt.subplot(grid[top_of_ax, :mutations.shape[1]])
            top_of_ax += 2
        elif col in continuous_columns:
            ax = plt.subplot(grid[top_of_ax:top_of_ax + 2, :mutations.shape[1]])
            top_of_ax += 3
        ax_dict['bottom'][col] = ax

    # 突变颜色， multi_hits标记为黑色
    if not color_dict:
        colors = sns.hls_palette(len(variants))
        color_dict = dict(zip(variants, colors))
        color_dict['multi_hits'] = (0, 0, 0)
        color_dict['no mutate'] = (0.9, 0.9, 0.9)

    # 临床信息（画图）
    # 定义连续、离散变量和突变类型的colormap
    continuous_colors = dict(zip(continuous_columns, sns.husl_palette(len(continuous_columns))))
    if discrete_colors is None:
        discrete_colors = {}

    for column in discrete_columns:
        values = sorted(list(set(sample_info[column])))
        discrete_colors.setdefault(column, dict(zip(values, sns.hls_palette(len(values)))))
    sample_info.index = range(sample_info.shape[0])

    # 画图 - 临床信息
    ax_info = {**ax_dict['upper'], **ax_dict['bottom']}
    for column in info_loc['upper'] + info_loc['bottom']:
        ax = ax_info[column]
        if column in discrete_columns:
            for subtype in discrete_colors[column]:
                tmp_df = sample_info[sample_info[column] == subtype]
                ax.bar(x=tmp_df.index, height=[1] * tmp_df.shape[0], color=discrete_colors[column][subtype], width=1,
                       label=subtype)
            ax.spines['left'].set_visible(False)
            ax.set_yticks([])
            ax.set_ylabel(column, rotation=0, ha='right', va='center', fontsize=10)

        if column in continuous_columns:
            ax.bar(x=sample_info.index, height=sample_info[column], color=continuous_colors[column])
            ax.set_ylabel(column, rotation=0, ha='right', fontsize=10)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_xlim(-0.6, sample_info.shape[0] - 0.4)

    # 突变统计
    patient_mut_stat = pd.concat([pd.value_counts(mutations[i]) for i in range(len(xticklabels))], axis=1)
    patient_mut_stat.fillna(0, inplace=True)
    patient_mut_stat = patient_mut_stat * 100 / mutations.shape[0]

    gene_mut_stat = pd.concat([pd.value_counts(mutations.loc[i]) for i in range(len(yticklabels))], axis=1)
    gene_mut_stat.fillna(0, inplace=True)
    gene_mut_stat = gene_mut_stat * 100 / mutations.shape[1]

    # 样本突变统计(柱状图)
    ax = ax_dict['mut_stat_sample']
    ploted_variants = []
    for variant in variants:
        if len(ploted_variants) == 0:
            bottom = [0] * patient_mut_stat.shape[1]
        else:
            bottom = patient_mut_stat.loc[ploted_variants].sum(0)
        ax.bar(x=patient_mut_stat.columns, height=patient_mut_stat.loc[variant], bottom=bottom,
               color=color_dict[variant])
        ploted_variants.append(variant)
    if fraction_lim:
        ax.set_ylim(0, fraction_lim[0])
    if fraction_annot:
        ax.set_xticks(range(patient_mut_stat.shape[1]))
        ax.set_xticklabels([f'{i:0.0f}%' for i in patient_mut_stat.iloc[1:, :].sum(0)])
    else:
        ax.set_xticks([])
    for axis in ['bottom', 'top', 'right']:
        ax.spines[axis].set_visible(False)

    # 基因突变统计(柱状图 )
    ax = ax_dict['mut_stat_gene']
    ploted_variants = []
    for variant in variants:
        if len(ploted_variants) == 0:
            bottom = [0] * gene_mut_stat.shape[1]
        else:
            bottom = gene_mut_stat.loc[ploted_variants].sum(0)
        ax.barh(y=gene_mut_stat.columns, width=gene_mut_stat.loc[variant], left=bottom, color=color_dict[variant],
                align='center', height=0.8)
        ploted_variants.append(variant)

    ax.invert_yaxis()
    if fraction_lim:
        ax.set_xlim(0, fraction_lim[1])
    ax.set_ylim(len(yticklabels) - 1.5, -0.5)
    if fraction_annot:
        ax.set_yticks(range(gene_mut_stat.shape[1]))
        ax.set_yticklabels([f'{i:0.0f}%' for i in gene_mut_stat.iloc[1:, :].sum(0)])
    else:
        ax.set_yticks([])

    ax.xaxis.set_ticks_position('top')
    for axis in ['bottom', 'left', 'right']:
        ax.spines[axis].set_visible(False)
    # 突变热图
    ax = ax_dict['heatmap']
    for variant in color_dict.keys():
        temp_df = heatmap_sparse_df[heatmap_sparse_df['variant'] == variant]
        ax.bar(x=temp_df['x'], bottom=temp_df['y'] - 0.4, height=0.8, width=0.8, edgecolor='grey',
               color=color_dict[variant],
               label=variant, align='center')
    ax.invert_yaxis()
    ax.set_yticks(range(len(yticklabels)))
    ax.set_yticklabels(yticklabels)
    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels)
    for axis in ['bottom', 'left', 'right', 'top']:
        ax.spines[axis].set_visible(False)
    ax.set_ylim(len(yticklabels) - 1.5, -0.5)
    ax.set_xlim(-0.5, sample_info.shape[0])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # 调整突变比例的坐标系
    ax_dict['mut_stat_gene'].set_ylim(ax_dict['heatmap'].get_ylim())  # 基因突变统计的纵坐标与热图保持一致
    for ax in sum(
            [[ax_dict['mut_stat_sample']], list(ax_dict['upper'].values()), list(ax_dict['upper'].values())],
            []):
        ax.set_xlim(-0.5, sample_info.shape[0])

    # 转移xtickslabel
    if len(info_loc['bottom']) > 0:
        ax_heatmap = ax_dict['heatmap']
        ax = ax_dict['bottom'][-1]
        ax.set_xticks(ax_heatmap.get_xticks())
        ax.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=90)
        ax_heatmap.set_xticks([])
        ax_heatmap.set_xticklabels([])

    # 设置图例
    top_of_legend = 0
    axes_require_legend = {'Variants Type': ax_dict['heatmap']}
    axes_require_legend.update(ax_dict['upper'])
    axes_require_legend.update(ax_dict['bottom'])
    titles = ['Variants Type']
    titles.extend([col for col in info_loc['upper'] + info_loc['bottom'] if col in discrete_columns])

    ax_heatmap = ax_dict['heatmap']
    for title in titles:
        ax = axes_require_legend[title]
        handles, labels = ax.get_legend_handles_labels()
        legend = ax_heatmap.legend(handles=handles, labels=labels, title=title, bbox_to_anchor=(
            (mutations.shape[1] + 3) / mutations.shape[1], 1 - (top_of_legend * legend_cc / 50)), ncol=1,
                                   loc='upper left',
                                   fontsize=15, title_fontsize=15)
        ax_heatmap.add_artist(legend)
        top_of_legend += (len(labels) + 2)

    return fig, discrete_columns, continuous_columns, figsize, ax_dict
