import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import rcParams


# %%
def forest_plot(cox_result, ax=None):
    """cox回归森林图

    :param cox_result: biokit.analysis.cox 函数返回值
    :param ax:
    :return: ax
    """
    if not ax:
        figsize = (10, (cox_result.shape[0] + 3) / 2)
        fig, ax = plt.subplots(figsize=figsize)

    config = {"font.size": 10}
    rcParams.update(config)

    # 检查无穷大值,找出Inf之外的最大值
    max_except_inf = cox_result['HR(95CI-High)'][~np.isinf(cox_result['HR(95CI-High)'])].max()

    # 整理分组表格
    index_df = cox_result.index.to_frame()
    index_dict = {}
    table = []
    for groupby in index_df[0].unique():
        # print(groupby)
        # print(index_df[0])
        groups = index_df.loc[groupby].index.to_list()
        groups_cox_table = []
        group_ref_cox_table = []
        if len(groups) > 1:
            index_dict[groupby] = groups
            for group in groups:
                hr, hr_l, hr_h, n_sample = cox_result[['HR', 'HR(95CI-Low)', 'HR(95CI-High)', 'n_sample']].loc[groupby].loc[group]
                if (groupby, group) in cox_result.multi_ref_variables:
                    group_ref_cox_table.append([groupby, f'{group}\n(n={n_sample})', 'reference'])
                else:
                    if hr_h < 1000:
                        groups_cox_table.append(
                            ['', f'{group}\n(n={n_sample})', f'{hr:0.2f}\n({hr_l:0.2f} ~ {hr_h:0.2f})'])
                    else:
                        groups_cox_table.append(
                            ['', f'{group}\n(n={n_sample})', f'{hr:0.2f}\n({hr_l:0.2f} ~ {hr_h:0.2e})'])
                        ax.ticklabel_format(style='sci', axis='y')
        else:
            group = groupby
            n_sample = cox_result.df.shape[0]
            hr, hr_l, hr_h = cox_result[['HR', 'HR(95CI-Low)', 'HR(95CI-High)']].loc[groupby, group]
            if hr_h < 1000:
                groups_cox_table.append(
                    [groupby, f'{group}\n(n={n_sample})', f'{hr:0.2f}\n({hr_l:0.2f} ~ {hr_h:0.2f})'])
            else:
                groups_cox_table.append(
                    [groupby, f'{group}\n(n={n_sample})', f'{hr:0.2f}\n({hr_l:0.2f} ~ {hr_h:0.2e})'])
                ax.ticklabel_format(style='sci', axis='y')

        table.extend(group_ref_cox_table)
        table.extend(groups_cox_table)

    # 画图
    ax.set_ylim(0, index_df.shape[0])
    ax.set_yticks([])
    ax.set_xlim(-max_except_inf * 1.2, max_except_inf * 1.2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # 分区条
    width = max_except_inf * 2.4
    for y in range(0, index_df.shape[0], 2):
        ax.barh(y=y, width=width, left=-max_except_inf * 1.2, color='grey', alpha=0.3, align='edge',
                height=1, zorder=0)

    # 分组表
    ax.table(table, cellLoc='center', bbox=(0, 0, 0.45, 1), edges='open', zorder=1)

    # HR(CI95)
    # 6段参考线
    ax.plot([1, 1], [0, index_df.shape[0]], linestyle='--', color='red', alpha=0.7, zorder=1)
    lines = [0]
    interval = max_except_inf // 4
    if interval > 0:
        lines.extend(list(
            np.arange(1 + interval, max_except_inf, interval)))
        for x in lines:
            ax.plot([x, x], [0, index_df.shape[0]], color='grey', alpha=0.7, zorder=1)
        # 参考线坐标
        xticklabels = [0, 1]
        xticklabels.extend(list(np.arange(1 + interval, max_except_inf, interval)))
        xticks = [i for i in xticklabels]

        ax.set_xticks(xticks)
        if ax.get_xlim()[1] < 1000:
            ax.set_xticklabels([str(int(i)) for i in xticklabels])
        else:
            ax.set_xticklabels([f'{i:0.2e}' for i in xticklabels])

    else:
        ax.set_xticks([i for i in ax.get_xticks() if 0 <= i < max_except_inf * 1.2])
        ax.set_xticklabels([f'{i:0.2f}' for i in ax.get_xticks()])

    # HR + error bar
    cox_table = cox_result.to_numpy().tolist()
    cox_table.reverse()
    significant_p = cox_result['p-value'][cox_result['p-value'] <= 0.05]
    # 如果有p值显著，按p值标HR方块颜色，如果没有，全部标黑色
    if len(significant_p) > 0:
        cmap = sns.blend_palette(['grey', 'red'], as_cmap=True)
        p_interval = -np.log10(min(significant_p)) - -np.log10(0.05)
        color_dict = dict(
            zip(significant_p, cmap([256 * (-np.log10(i) - -np.log10(0.05)) / p_interval for i in significant_p])))
    else:
        color_dict = {}

    # 自上而下绘制HR和p值
    y = len(table) - 1
    text_bottom = max_except_inf + max_except_inf * 0.02
    for line in table:
        if line[0] != '':
            groupby = line[0]
        group = line[1].split('\n')[0]
        hr, hr_l, hr_h, p = cox_result.loc[groupby, group]

        # 如果HR high是无穷大，画箭头
        if hr_l == 0:
            ax.arrow(x=0, y=y + 0.5, dx=max_except_inf, dy=0, color=color_dict.get(p, 'black'), width=0.03,
                     length_includes_head=True, zorder=2)
            ax.plot([0, 0], [y + 0.4, y + 0.6], color='black', zorder=2)
        else:
            ax.plot([hr_l, hr_h], [y + 0.5, y + 0.5], color='black', zorder=2)
            if hr != 1:
                ax.plot([hr_l, hr_l], [y + 0.4, y + 0.6], color='black', zorder=2)
                ax.plot([hr_h, hr_h], [y + 0.4, y + 0.6], color='black', zorder=2)
        ax.scatter(x=hr, y=y + 0.5, s=50, color=color_dict.get(p, 'black'), zorder=3)
        # p值
        ax.text(y=y + 0.5, x=text_bottom, s=f'{p:0.4f}' if p > 0.0001 else f'{p:0.2e}',
                color=color_dict.get(p, 'black'), ha='left', va='center')
        y -= 1

    # 在上方各留出两行的空间，在下方各留出一行的空间
    edge_width = 1 / (cox_result.shape[0] + 3)
    title_y = 1 - edge_width
    top = 1 - 2 * edge_width
    bottom = edge_width
    plt.suptitle('Hazard Ratio', fontsize=20, y=title_y, va='bottom')
    plt.subplots_adjust(left=0.05, right=0.95, bottom=bottom, top=top)
    return ax
