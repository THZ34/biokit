import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import rcParams


# %%
def cox(cox_df, df, ax=None):
    """cox回归森林图

    :param ax:
    :param cox_df: biokit.analysis.cox 函数返回值
    :param df: cox分析的原始数据
    :return: ax
    """
    config = {
        "font.family": 'Microsoft YaHei',
        "font.size": 10,
        "mathtext.fontset": 'stix',
        "font.serif": ['SimSun'],
    }
    rcParams.update(config)

    # 整理分组表格
    index_df = cox_df.index.to_frame()
    index_dict = {}
    table = []
    for groupby in index_df[0].unique():
        groups = index_df.loc[groupby].index.to_list()
        index_dict[groupby] = groups
        if len(groups) == 1:
            hr, hr_l, hr_h = cox_df[['HR', 'HR(95CI-Low)', 'HR(95CI-High)']].loc[groupby].loc[groups[0]]
            table.append([groupby, groups[0], f'{hr:0.2f}\n({hr_l:0.2f} ~ {hr_h:0.2f})'])
        else:
            group = groups[0]
            n_samples = df[df[groupby] == group].shape[0]
            hr, hr_l, hr_h = cox_df[['HR', 'HR(95CI-Low)', 'HR(95CI-High)']].loc[groupby].loc[group]
            if hr == hr_l == hr_h == 1:
                table.append([groupby, f'{group}\n(n={n_samples})', 'reference'])
            else:
                table.append([groupby, f'{group}\n(n={n_samples})', f'{hr:0.2f}\n({hr_l:0.2f} ~ {hr_h:0.2f})'])

            for group in groups[1:]:
                n_samples = df[df[groupby] == group].shape[0]
                hr, hr_l, hr_h = cox_df[['HR', 'HR(95CI-Low)', 'HR(95CI-High)']].loc[groupby].loc[group]
                if hr == hr_l == hr_h == 1:
                    table.append(['', f'{group}\n(n={n_samples})', 'reference'])
                else:
                    table.append(['', f'{group}\n(n={n_samples})', f'{hr:0.2f}\n({hr_l:0.2f} ~ {hr_h:0.2f})'])
    # 画图
    if not ax:
        figsize = (10, index_df.shape[0] * 1.1 / 2)
        fig, ax = plt.subplots(figsize=figsize)
    ax.set_ylim(0, index_df.shape[0])
    ax.set_yticks([])
    ax.set_xlim(-cox_df['HR(95CI-High)'].max() * 1.2, cox_df['HR(95CI-High)'].max() * 1.2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # 分区条
    width = cox_df['HR(95CI-High)'].max() * 2.4
    for y in range(0, index_df.shape[0], 2):
        plt.barh(y=y, width=width, left=-cox_df['HR(95CI-High)'].max() * 1.2, color='grey', alpha=0.3, align='edge',
                 height=1)

    # 分组表
    ax.table(table, cellLoc='center', bbox=(0, 0, 0.45, 1), edges='open')

    # HR(CI95)
    # 6段参考线
    lines = [0]
    interval = cox_df['HR(95CI-High)'].max() // 4
    lines.extend(list(np.arange(1 + interval, cox_df['HR(95CI-High)'].max(), interval)))
    for x in lines:
        ax.plot([x, x], [0, index_df.shape[0]], color='grey', alpha=0.7)
    ax.plot([1, 1], [0, index_df.shape[0]], linestyle='--', color='red', alpha=0.7)
    # 参考线坐标
    xticklabels = [0, 1]
    xticklabels.extend(list(np.arange(1 + interval, cox_df['HR(95CI-High)'].max(), interval)))
    xticks = [i for i in xticklabels]
    ax.set_xticks(xticks)
    ax.set_xticklabels(map(lambda x: str(int(x)), xticklabels))
    # HR + error bar
    cox_table = cox_df.to_numpy().tolist()
    cox_table.reverse()
    significant_p = cox_df['p-value'][cox_df['p-value'] <= 0.05]
    cmap = sns.blend_palette(['grey', 'red'], as_cmap=True)
    p_interval = -np.log10(min(significant_p)) - -np.log10(0.05)
    color_dict = dict(
        zip(significant_p, cmap([256 * (-np.log10(i) - -np.log10(0.05)) / p_interval for i in significant_p])))
    y = 0
    for hr, hr_l, hr_h, p in cox_table:
        # plt.barh(y=y + 0.5, width=box_width, height=box_height, left=hr, color=color_dict.get(p, 'black'),
        # align='center')
        plt.errorbar(y=y + 0.5, x=hr, xerr=hr_h - hr, fmt='o', color=color_dict.get(p, 'black'), ecolor='black')
        y += 1
    # p值
    y = 0
    text_bottom = max(cox_df['HR(95CI-High)']) + cox_df['HR(95CI-High)'].max() * 0.02
    for hr, hr_l, hr_h, p in cox_table:
        plt.text(y=y + 0.5, x=text_bottom, s=f'{p:0.5f}', color=color_dict.get(p, 'black'))
        y += 1

    plt.suptitle('Hazard Ratio', fontsize=20)
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=cox_df.shape[0] / (cox_df.shape[0] + 2))
    return ax
