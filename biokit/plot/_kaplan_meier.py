import copy
import math
import warnings

import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.exceptions import ConvergenceError
from lifelines.statistics import multivariate_logrank_test
from matplotlib import pyplot as plt
from matplotlib import rcParams

from biokit.analysis import cox


# %%

def kaplan_meier(grouped_df, groupby, time='time', status='status', groups=None, cox_analysis=True, figsize=None,
                 color_dict=None, cox_table_loc=(0.45, 1), dropna=True, cox_ref=None, show_censors=False, ci_show=False):
    """绘制km曲线
    :param cox_ref:
    :param dropna: 去掉缺失值
    :param cox_table_loc: cox表格 左上角的坐标
    :param figsize:
    :param color_dict:
    :param cox_analysis:
    :param grouped_df: pandas.DataFrame, contain at least three columns(time,status,group)
    :param groupby: column name of group, necessary
    :param time: column name of time, default: 'time'
    :param status: column name of status, default: 'status'
    :param groups: groups and order
    :return: matrix,fig
    """
    # 去掉缺失值
    if dropna:
        # 如果有缺失值，去除并警告
        if grouped_df[[time, status, groupby]].isna().sum().sum() > 0:
            warnings.warn(f'DataFrame contains missing values\n{grouped_df[[time, status, groupby]].isna().sum()}')
            grouped_df = grouped_df.copy()
            grouped_df.dropna(subset=[time, status, groupby], axis=0, inplace=True)

    if not groups:
        groups = list(set(grouped_df[groupby]))
        groups.sort()

    # 筛除无用的组
    grouped_df = grouped_df[[(i in groups) for i in grouped_df[groupby]]]

    # 创建画板
    if not figsize:
        figsize = (6, 3 + len(groups) / 2)
    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(nrows=7 + len(groups), ncols=1)
    ax_dict = {'km': plt.subplot(grid[:5, :]), 'table': plt.subplot(grid[7:, :])}

    # 设置颜色
    if not color_dict:
        colors = sns.husl_palette(len(groups))
        color_dict = dict(zip(groups, colors))

    # kmf拟合，存储，绘制km曲线
    ax_km = ax_dict['km']
    matrix = []
    kmf_dict = {}
    # log-rank 检验
    log_rank = multivariate_logrank_test(grouped_df[time].to_list(), grouped_df[groupby].to_list(),
                                         grouped_df[status].to_list())
    for group in groups:
        T, E = grouped_df[grouped_df[groupby] == group][[time, status]].to_numpy().T
        kmf = KaplanMeierFitter()
        kmf.fit(T, E)
        kmf_dict[group] = copy.copy(kmf)
        kmf.plot_survival_function(ax=ax_km, ci_show=ci_show, color=color_dict[group], label=group,
                                   show_censors=show_censors,)
        os_median = kmf.median_survival_time_
        upper = os_median + 1.96 * T.std() / math.sqrt(len(T))
        lower = os_median - 1.96 * T.std() / math.sqrt(len(T))
        median_ci95 = f'{os_median:4.1f}({lower:4.1f}-{upper:4.1f})'
        matrix.append([groupby, group, len(E), E[E == 1].shape[0], median_ci95, log_rank.p_value])

    # get xticks
    raw_xticks = ax_km.get_xticks()
    xticks = ax_km.get_xticks()

    # 表格在图的右上角显示，三线表格式，以全图长宽为参考，表格长度=0.5，高度= (组数+1) x 0.08
    if cox_analysis:
        # 整理表格
        try:
            cox_df = cox(grouped_df[[time, status, groupby]], time=time, status=status, mod='multiple',
                         ref_dict={groupby: cox_ref or groups[0]})
            cox_table = [['group', f'median {time}', 'cox HR(95CI)', 'cox p-value']]
            for group in groups:
                hr, hr_l, hr_h, cox_p = cox_df.loc[(groupby, group), ['HR', 'HR(95CI-Low)', 'HR(95CI-High)', 'p-value']]
                cox_table.append(
                    [group, f'{kmf_dict[group].median_survival_time_:0.1f}', f'{hr:0.1f}({hr_l:0.1f}~{hr_h:0.1f})',
                     f'{cox_p:0.2e}' if cox_p < 0.0001 else f'{cox_p:0.4f}'])
            length = 0.5
            height = (len(groups) + 1) * 0.08
            left = cox_table_loc[0]
            bottom = cox_table_loc[1] - height
            cox_table_plot = ax_km.table(cox_table, cellLoc='center', bbox=(left, bottom, length, height), edges='open')
            cox_table_plot.auto_set_font_size(False)
            cox_table_plot.set_fontsize(7)

            # 三线表，线
            ax_length = raw_xticks[-1] - (xticks[0] - xticks[1])
            line_start = ax_length * left + (xticks[0] - xticks[1]) / 4
            line_end = ax_length * 0.95
            line_heights = [bottom + height, bottom + height - 0.08, bottom]
            for line_height in line_heights:
                ax_km.plot([line_start, line_end], [line_height, line_height], color='black', linewidth=1)
        except ConvergenceError:
            print(f'{groupby}, {time}, {status} ConvergenceError')

    # 设置图的坐标标签，范围等参数
    ax_km.set_xlabel('')
    ax_km.set_xticks(raw_xticks[1:])
    ax_km.set_xlim((xticks[0] - xticks[1]) / 4, raw_xticks[-1] - 3 * (xticks[0] - xticks[1]) / 4)
    ax_km.set_ylabel('Survival Probability (%)')
    ax_km.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax_km.set_yticklabels([str(int(i * 100)) for i in ax_km.get_yticks()])
    ax_km.set_ylim(0, 1)

    ax_km.spines['right'].set_visible(False)
    ax_km.spines['top'].set_visible(False)
    ax_km.get_legend().set_visible(False)

    # 标题颜色，低于0.05设置为红色，否则为黑色
    ax_km.set_title(
        f'{groupby} p-value:{log_rank.p_value:0.2e}' if log_rank.p_value < 0.0001 else f'{groupby} p-value:{log_rank.p_value:0.4f}',
        fontsize=rcParams['font.size'] + 2, color='red' if log_rank.p_value < 0.05 else 'black', y=1.1)

    # 获取 number in risk
    time_number_table = []
    for group in groups:
        number_at_risk = [grouped_df[(grouped_df[groupby] == group) & (grouped_df[time] >= time_stamp)][status].shape[0]
                          for time_stamp in raw_xticks[1:]]
        time_number_table.append(number_at_risk)

    # table number of risk, 在km曲线下方显示
    ax_table = ax_dict['table']
    ax_table.set_xlim(ax_km.get_xlim()[0], ax_km.get_xlim()[1])
    interval = (ax_km.get_xticks()[1] - ax_km.get_xticks()[0]) / (ax_km.get_xlim()[1] - ax_km.get_xlim()[0])
    table_length = interval * len(ax_km.get_xticks())
    table_start = -interval / 4
    # table 标签
    for group, i in zip(reversed(groups), range(len(groups))):
        ax_table.text(x=-(ax_km.get_xticks()[1] - ax_km.get_xticks()[0]) / 2, y=i + 0.5, s=group,
                      color=color_dict[group], ha='right', va='center')

    ax_table.table(time_number_table, bbox=(table_start, 0, table_length, 1), cellLoc='center', edges='open')
    ax_table.axis('off')
    ax_table.set_title('Number at Risk', fontsize=rcParams['font.size'] + 2, ha='center')
    ax_table.set_ylim(0, len(groups))
    # 返回生存分析矩阵和ax
    # plt.subplots_adjust(top=0.85, bottom=0, left=0.15, right=0.95, hspace=-0.1)
    plt.subplots_adjust(top=0.85, bottom=0, left=0.15, right=0.95, hspace=-0.2)
    return matrix, log_rank.p_value, fig
