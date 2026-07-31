from collections.abc import Iterable
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
from scipy.stats import fisher_exact, chi2_contingency


def crosstab_plot(crosstab=None, value_x=None, value_y=None, xlabel=None, ylabel=None, pad=-0.05,
                  rounding_size=0.15, x_order=None, y_order=None, table_cmap='seagreen', row_cmap='orangered',
                  col_cmap='deepskyblue', ax=None, background_param=None, show_pvalue=False,
                  standard_scale=None, adaptive_font=True, xlabel_color='black', ylabel_color='black'):
    """
    绘制列联表（带色块背景、比例归一化和 p 值标注）
    --------------------------------------------------
    参数：
        crosstab: DataFrame 或 None
            已经生成的列联表。
        value_x, value_y: Iterable
            若未提供 crosstab，则从两个分类变量构造列联表。
        xlabel, ylabel: str
            行列标签名称。
        standard_scale: {'group', 'var', None}
            - 'group': 按行最大值归一化 (横向比较)
            - 'var'  : 按列最大值归一化 (纵向比较)
            - None   : 全表归一化
        show_pvalue: bool
            显示统计学显著性（自动 Fisher/χ²）。
        adaptive_font: bool
            行列标签字号随表格大小自动调整。
    """

    # ==== 数据预处理 ====
    if isinstance(crosstab, pd.DataFrame):
        xlabel = xlabel or crosstab.columns.name
        ylabel = ylabel or crosstab.index.name

    elif crosstab is None and (isinstance(value_x, Iterable) and isinstance(value_y, Iterable)):
        xlabel = xlabel or 'X label'
        ylabel = ylabel or 'Y label'
        x_order = x_order or sorted(list(set(value_x)))
        y_order = y_order or sorted(list(set(value_y)))
        crosstab_input = pd.DataFrame([value_x, value_y], index=[xlabel, ylabel]).T
        crosstab = pd.crosstab(crosstab_input[ylabel], crosstab_input[xlabel])
        crosstab = crosstab.loc[y_order, x_order]
    else:
        raise ValueError('Require crosstab or (value_x and value_y)')

    # ==== 图形初始化 ====
    if not ax:
        fig, ax = plt.subplots(figsize=(crosstab.shape[1] + 1, crosstab.shape[0] + 1))
    else:
        fig = ax.figure

    if not background_param:
        background_param = {'facecolor': 'white', 'linestyle': '--', 'edgecolor': 'black', 'alpha': 1}

    fig_width = crosstab.shape[1] + 1
    fig_height = crosstab.shape[0] + 1
    table_width, table_height = crosstab.shape[1], crosstab.shape[0]
    ax.axis(False)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    patchs = []
    # ==== 主体表 ====
    ax.table(crosstab.values, cellLoc='center', edges='open',
             bbox=((1 / fig_width), 0, table_width / fig_width, table_height / fig_height), zorder=2)

    # ==== 行标签 ====
    row_table = ax.table([[i] for i in crosstab.index], cellLoc='center', edges='open', zorder=2,
                         bbox=(0.5 / fig_width, 0, 0.5 / fig_width, table_height / fig_height))
    for cell in row_table.get_celld().values():
        cell.get_text().set_rotation('vertical')

    # ==== 列标签 ====
    col_table = ax.table([crosstab.columns.to_list()], cellLoc='center', zorder=2,
                         bbox=((1 / fig_width), table_height / fig_height, table_width / fig_width, 0.5 / fig_height),
                         edges='open')

    # ==== 字号自适应 ====
    if adaptive_font:
        for cell in row_table.get_celld().values():
            text_obj = cell.get_text()
            bbox_height = cell.get_height() * fig_height * fig.dpi / 72
            fontsize = max(6, min(12, bbox_height * 0.55))
            text_obj.set_fontsize(fontsize)

        for cell in col_table.get_celld().values():
            text_obj = cell.get_text()
            bbox_width = cell.get_width() * fig_width * fig.dpi / 72
            fontsize = max(6, min(12, bbox_width * 0.3))
            text_obj.set_fontsize(fontsize)

    ax.text(x=0.25, y=table_height / 2, s=ylabel, rotation=90, ha='center', va='center',
            fontsize=10, zorder=2, fontweight='bold')
    ax.text(x=1 + table_width / 2, y=table_height + 0.75, s=xlabel, ha='center', va='center',
            fontsize=10, zorder=2, fontweight='bold')

    # ==== 背景色块 ====
    patchs.append(FancyBboxPatch(xy=(0, 0), width=0.5, height=table_height, color=row_cmap, zorder=1,
                                 boxstyle=f'round,pad={pad},rounding_size={rounding_size}'))
    row_colors = sns.light_palette(row_cmap, table_height)
    for y, color in zip(range(table_height), row_colors):
        patchs.append(FancyBboxPatch(xy=(0.5, y), width=0.5, height=1, color=color, zorder=1,
                                     boxstyle=f'round,pad={pad},rounding_size={rounding_size}'))

    patchs.append(FancyBboxPatch(xy=(1, table_height + 0.5), width=table_width, height=0.5, color=col_cmap, zorder=1,
                                 boxstyle=f'round,pad={pad},rounding_size={rounding_size}'))
    col_colors = sns.light_palette(col_cmap, table_width)
    for x, color in zip(reversed(range(table_width)), col_colors):
        patchs.append(FancyBboxPatch(xy=(x + 1, table_height), width=1, height=0.5, color=color, zorder=1,
                                     boxstyle=f'round,pad={pad},rounding_size={rounding_size}'))

    # ==== 表格主色块：比例标准化 ====
    color_table = crosstab.copy().astype(float)
    if standard_scale == 'group':  # 行归一化
        color_table = color_table.div(color_table.max(axis=1).replace(0, np.nan), axis=0)
    elif standard_scale == 'var':  # 列归一化
        color_table = color_table.div(color_table.max(axis=0).replace(0, np.nan), axis=1)
    else:
        standard_scale = None  # 全表归一化

    norm = plt.Normalize(vmin=color_table.min().min(), vmax=color_table.max().max())
    if isinstance(table_cmap, str):
        cmap = sns.light_palette(table_cmap, as_cmap=True)
    else:
        cmap = table_cmap

    for x in range(table_width):
        for y in range(table_height):
            val = color_table.iloc[y, x]
            color = cmap(norm(val))
            patchs.append(FancyBboxPatch(xy=(x + 1, table_height - y - 1), width=1, height=1, color=color, zorder=1,
                                         boxstyle=f'round,pad={pad},rounding_size={rounding_size}', ec='gainsboro'))

    # ==== 外框 ====
    patchs.append(FancyBboxPatch(xy=(0, 0), width=fig_width, height=fig_height,
                                 facecolor=background_param['facecolor'],
                                 boxstyle=f'round,pad=0,rounding_size={rounding_size}', zorder=0, linewidth=2,
                                 edgecolor=background_param['edgecolor'], linestyle=background_param['linestyle']))

    for patch in patchs:
        ax.add_patch(patch)

    # ==== p 值标注 ====
    if show_pvalue:
        arr = crosstab.to_numpy()
        if arr.shape == (2, 2):
            pvalue = fisher_exact(arr)[1]
            method = "Fisher's exact"
        else:
            _, pvalue, _, _ = chi2_contingency(arr)
            method = "Chi-square"

        if pvalue == 0:
            p_text = '<0.0001'
        elif pvalue < 0.001:
            p_text = '<0.001'
        elif pvalue < 0.01:
            p_text = '<0.01'
        elif pvalue < 0.05:
            p_text = '<0.05'
        else:
            p_text = f'= {pvalue:.3f}'
        text = f"{method}\nP{p_text}"
        textx = 0.5 / (crosstab.shape[1] + 1)
        texty = (0.5 + crosstab.shape[1]) / (crosstab.shape[0] + 1)
        ax.text(textx, texty, text, transform=ax.transAxes, ha='center', va='center',
                fontsize=8, fontweight='bold', color='red' if pvalue < 0.05 else 'black')

    return ax
