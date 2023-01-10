import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix


def oncoplot(mutations, sample_info, figsize=None, color_dict=None, discrete_colors=None, fraction_lim=None,
             info_loc=None, fraction_annot=False, legend_cc=1, allow_multi_hits=True, heatmap_kind='box'):
    """oncoplot 瀑布图


    :param heatmap_kind: 热图风格，可选box和circle
    :param allow_multi_hits: 如果不允许基因多位点突变，将所有多位点突变替换成multi_hits
    :param fraction_percent:
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

    # 如果不允许基因多位点突变，将所有多位点突变替换成multi_hits
    if not allow_multi_hits:
        np.where(pd.concat([mutations[col].str.contains(',') for col in mutations.columns], axis=1), 'multi_hits',
                 mutations)

    # 转换成稀疏矩阵
    heatmap_sparse_df = coo_matrix(mutations)
    heatmap_sparse_df = pd.DataFrame([heatmap_sparse_df.row, heatmap_sparse_df.col, heatmap_sparse_df.data]).T
    heatmap_sparse_df.columns = ['y', 'x', 'variant']

    # 获取所有突变类型
    variants = heatmap_sparse_df['variant'].unique()
    new_variants = []
    for mut_type in variants:
        if ',' in mut_type:
            new_variants.extend(mut_type.split(','))
        else:
            new_variants.append(mut_type)
    variants = sorted(list(set(new_variants)))
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
        figsize = (figsize[0] * 1.5, figsize[1] * 1.5)
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

    # 突变比例
    if fraction_annot:
        ax_dict['mut_stat_sample'] = plt.subplot(grid[top_of_ax:top_of_ax + 2, :mutations.shape[1]])
        top_of_ax += 3
        ax_dict['mut_stat_gene'] = plt.subplot(grid[top_of_ax:top_of_ax + mutations.shape[0], mutations.shape[1] + 1:])
    else:
        ax_dict['mut_stat_sample'] = plt.subplot(grid[top_of_ax:top_of_ax + 3, :mutations.shape[1]])
        top_of_ax += 3
        ax_dict['mut_stat_gene'] = plt.subplot(grid[top_of_ax:top_of_ax + mutations.shape[0], mutations.shape[1]:])

    # 热图
    ax_dict['heatmap'] = plt.subplot(grid[top_of_ax:top_of_ax + mutations.shape[0], :mutations.shape[1]])
    top_of_ax += mutations.shape[0] + 1

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
        if 'no mutate' in color_dict:
            color_dict['no mutate'] = (0.9, 0.9, 0.9)
        if 'multi_hits' in color_dict:
            color_dict['multi_hits'] = (0, 0, 0)

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

        for axis in ['bottom', 'right', 'top']:
            ax.spines[axis].set_visible(False)
        ax.set_xticks([])
        ax.set_xlim(-0.6, sample_info.shape[0] - 0.4)

    # 突变统计
    patient_mut_stat = pd.DataFrame()
    for patient_num in mutations.columns:
        patient_muts = mutations[patient_num]
        single_mut = patient_muts[~patient_muts.str.contains(',')].to_list()
        multi_muts = patient_muts[patient_muts.str.contains(',')]
        for multi_mut in multi_muts:
            single_mut.extend(multi_mut.split(','))
        patient_mut_stat = pd.concat([patient_mut_stat, pd.value_counts(single_mut)], axis=1)
    patient_mut_stat.columns = mutations.columns
    patient_mut_stat.fillna(0, inplace=True)
    patient_mut_stat = patient_mut_stat.loc[variants]

    gene_mut_stat = pd.DataFrame()
    for gene_num in mutations.index:
        gene_muts = mutations.loc[gene_num]
        single_mut = gene_muts[~gene_muts.str.contains(',')].to_list()
        multi_muts = gene_muts[gene_muts.str.contains(',')]
        for multi_mut in multi_muts:
            single_mut.extend(multi_mut.split(','))
        gene_mut_stat = pd.concat([gene_mut_stat, pd.value_counts(single_mut)], axis=1)

    gene_mut_stat.columns = mutations.index
    gene_mut_stat.fillna(0, inplace=True)
    gene_mut_stat = gene_mut_stat.loc[variants]

    if not allow_multi_hits:
        patient_mut_stat = patient_mut_stat * 100 / mutations.shape[0]
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
    else:
        ax.set_ylim(0, patient_mut_stat.sum(0).max())

    if fraction_annot:
        ax.set_xticks(range(patient_mut_stat.shape[1]))
        ax.set_xticklabels([f'{i:0.0f}%' for i in (mutations != 'no mutate').sum(0) * 100 / mutations.shape[0]])
        ax.tick_params(bottom=False)
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
    else:
        ax.set_xlim(0, gene_mut_stat.sum(0).max())
    ax.set_ylim(len(yticklabels) - 1.5, -0.5)
    if fraction_annot:
        ax.set_yticks(range(gene_mut_stat.shape[1]))
        ax.set_yticklabels([f'{i:0.0f}%' for i in (mutations != 'no mutate').sum(1) * 100 / mutations.shape[1]],
                           fontsize=10, ha='center', x=-0.1)
        ax.tick_params(left=False)
    else:
        ax.set_yticks([])

    ax.xaxis.set_ticks_position('top')
    for axis in ['bottom', 'left', 'right']:
        ax.spines[axis].set_visible(False)
    # 突变热图
    ax = ax_dict['heatmap']
    if heatmap_kind == 'box':
        from biokit.plot._heatmap import heatmap_cumulativebox
        heatmap_cumulativebox(mutations, ax, color_dict, sep=',')
    elif heatmap_kind == 'circle':
        from biokit.plot._heatmap import heatmap_circledot
        heatmap_circledot(mutations, ax, color_dict, sep=',')

    ax.invert_yaxis()
    ax.set_yticks(range(len(yticklabels)))
    ax.set_yticklabels(yticklabels)
    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels)
    for axis in ['bottom', 'left', 'right', 'top']:
        ax.spines[axis].set_visible(False)
    ax.set_ylim(len(yticklabels) - 0.5, -0.5)
    ax.set_xlim(-0.6, sample_info.shape[0] - 0.4)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # 调整突变比例的坐标系
    ax_dict['mut_stat_gene'].set_ylim(ax_dict['heatmap'].get_ylim())  # 基因突变统计的纵坐标与热图保持一致
    for ax in sum([[ax_dict['mut_stat_sample']], list(ax_dict['upper'].values()), list(ax_dict['upper'].values())], []):
        ax.set_xlim(-0.6, sample_info.shape[0] - 0.4)

    # 转移xtickslabel
    if len(info_loc['bottom']) > 0:
        ax_heatmap = ax_dict['heatmap']
        ax = ax_dict['bottom'][info_loc['bottom'][-1]]
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
