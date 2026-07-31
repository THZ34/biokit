# -*- coding: utf-8 -*-
"""Nature 风格的 Kaplan-Meier 生存曲线绘图。

布局（自上而下，高度随分组数自适应）：
    ┌─────────────────────────────┐
    │  统计/Cox 表 (cox_ax)        │  行数 = 分组数 + 表头
    ├─────────────────────────────┤
    │  KM 曲线 (km_ax)             │  固定高度，不被多组压缩
    ├─────────────────────────────┤
    │  Number at risk (risk_ax)    │  行数 = 分组数；底部为时间轴
    └─────────────────────────────┘
"""
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.exceptions import ConvergenceError
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import median_survival_times

# ggsci npg 配色（Nature 系列常用的定性配色）
NPG_PALETTE = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#7E6148', '#DC0000',
               '#B09C85']


def _fmt_logrank_p(p):
    if p is None or np.isnan(p):
        return 'Log-rank P = NA'
    return 'Log-rank P < 0.001' if p < 0.001 else f'Log-rank P = {p:.3f}'


def _fmt_cox_p(p):
    if p is None or np.isnan(p):
        return 'NA'
    return f'{p:.1e}' if p < 0.001 else f'{p:.3f}'


def _fmt_median_ci(kmf):
    """中位生存期 (95% CI)；用 lifelines 正确估计 CI，未达中位记为 NR。"""
    m = kmf.median_survival_time_
    if m is None or np.isinf(m) or np.isnan(m):
        return 'NR'
    try:
        ci = median_survival_times(kmf.confidence_interval_)
        lo, hi = ci.iloc[0, 0], ci.iloc[0, 1]
    except Exception:
        lo = hi = np.nan
    s = lambda x: 'NR' if (x is None or np.isinf(x) or np.isnan(x)) else f'{x:.1f}'
    return f'{m:.1f} ({s(lo)}\u2013{s(hi)})'


def _draw_stats_table(ax, header, cell_text, groups, color_dict, fontsize):
    ax.axis('off')
    nrows = len(cell_text) + 1
    ncols = len(header)
    # 按各列最长文本自适应列宽（matplotlib 表格不会自动缩字，需手动分配）
    cols = list(zip(*([header] + cell_text)))
    maxlens = [max(len(str(x)) for x in col) for col in cols]
    total = float(sum(maxlens))
    col_widths = [ml / total for ml in maxlens]
    tbl = ax.table(cellText=[header] + cell_text, cellLoc='center', colWidths=col_widths, bbox=[0, 0, 1, 1],
                   edges='open')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(max(fontsize - 1, 6))
    for c in range(ncols):
        tbl[0, c].set_text_props(weight='bold')
    for r, g in enumerate(groups, start=1):
        tbl[r, 0].set_text_props(color=color_dict.get(g, 'black'))
    # 三线表：顶线、表头下线、底线
    for y in (1.0, 1.0 - 1.0 / nrows, 0.0):
        ax.plot([0, 1], [y, y], transform=ax.transAxes, color='black', lw=0.8, clip_on=False)


def kaplan_meier(grouped_df, groupby, time='time', status='status', groups=None, cox_analysis=True, cox_ref=None,
                 cox_fn=None, color_dict=None, width=5.8, km_height=3.2, row_height=0.30, figsize=None, base_fontsize=8,
                 show_censors=True, ci_show=True, ci_alpha=0.15, xlabel=None, ylabel='Survival probability (%)',
                 title=None, p_loc=(0.04, 0.06), dropna=True):
    """绘制 Nature 风格 KM 曲线。

    参数（仅列出与原版不同/新增的）
    ----
    title : str | None
        图标题，置于最上方并预留空间。默认 None 时用 groupby（分组依据的
        指标名，如基因名），让读者一眼看出按什么分的高低组；传 '' 可去掉。
    cox_fn : callable | None
        Cox 分析函数；None 时延迟从 ``biokit.analysis`` 导入。签名需为
        ``cox_fn(df, time, status, mod, ref_dict) -> (cox_df, _, _)``。
    width / km_height / row_height : float
        画布宽、KM 曲线区高度、每行表高度（单位英寸）；figsize=None 时据此
        与分组数自动算出总高度，保证多组不被压缩。
    p_loc : tuple
        log-rank p 值在 KM 区内的位置（axes 比例坐标），默认左下。

    返回 (stats_df, logrank_p, fig)。stats_df 为含每组 N/事件/中位/HR 的表。
    """
    rc = {'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': base_fontsize, 'axes.linewidth': 0.8, 'axes.unicode_minus': False, 'pdf.fonttype': 42,
        # 矢量导出时文字可在 AI 中编辑
        'ps.fonttype': 42, 'svg.fonttype': 'none', }

    with plt.rc_context(rc):
        # ---- 缺失值 ----
        if dropna:
            sub = grouped_df[[time, status, groupby]]
            if sub.isna().any().any():
                warnings.warn(f'DataFrame 含缺失值，已删除：\n{sub.isna().sum()}')
                grouped_df = grouped_df.dropna(subset=[time, status, groupby]).copy()

        # ---- 分组 ----
        if not groups:
            groups = sorted(grouped_df[groupby].dropna().unique().tolist())
        grouped_df = grouped_df[grouped_df[groupby].isin(groups)]
        n = len(groups)

        # ---- 配色 ----
        if color_dict is None:
            palette = NPG_PALETTE[:n] if n <= len(NPG_PALETTE) else sns.husl_palette(n)
            color_dict = dict(zip(groups, palette))

        # ---- 标题（默认用分组指标名）----
        title = str(groupby) if title is None else title
        # title = None
        title_h = 0.42 if title else 0.0

        # ---- 自适应画布与三段布局 ----
        cox_h = (n + 1) * row_height + 0.15 if cox_analysis else 0.0
        risk_h = n * row_height + 0.45
        if figsize is None:
            fig_h = km_height + cox_h + risk_h + 0.5 + title_h
            figsize = (width, fig_h)
        else:
            fig_h = figsize[1]
        fig = plt.figure(figsize=figsize)

        ratios, has_cox = [], cox_analysis
        if has_cox:
            ratios.append(cox_h)
        ratios += [km_height, risk_h]
        gs = fig.add_gridspec(len(ratios), 1, height_ratios=ratios, hspace=0.12)
        i = 0
        cox_ax = fig.add_subplot(gs[i]) if has_cox else None
        i += 1 if has_cox else 0
        km_ax = fig.add_subplot(gs[i]);
        i += 1
        risk_ax = fig.add_subplot(gs[i])

        # ---- log-rank ----
        log_rank = multivariate_logrank_test(grouped_df[time].to_list(), grouped_df[groupby].to_list(),
            grouped_df[status].to_list())
        p = log_rank.p_value

        # ---- 拟合 + 画曲线 ----
        kmf_dict = {}
        for g in groups:
            sg = grouped_df[grouped_df[groupby] == g]
            T, E = sg[time].to_numpy(), sg[status].to_numpy()
            kmf = KaplanMeierFitter()
            kmf.fit(T, E, label=str(g))
            kmf_dict[g] = kmf
            kmf.plot_survival_function(ax=km_ax, color=color_dict[g], linewidth=1.6, ci_show=ci_show, ci_alpha=ci_alpha,
                show_censors=show_censors, censor_styles={'ms': 4, 'marker': '|'})

        # ---- KM 区样式 ----
        t_max = float(grouped_df[time].max())
        x_lo = -0.02 * t_max  # 左侧留白，避免风险表首列贴轴
        ticks = [t for t in mticker.MaxNLocator(nbins=6).tick_values(0, t_max) if 0 <= t <= t_max]
        km_ax.set_xlim(x_lo, t_max)
        km_ax.set_ylim(0, 1.02)
        km_ax.set_xticks(ticks)
        km_ax.set_xticklabels([])  # 时间轴统一放到最底部
        km_ax.set_xlabel('')  # 清掉 lifelines 默认的 'timeline'
        km_ax.set_yticks([0, .2, .4, .6, .8, 1])
        km_ax.set_yticklabels([0, 20, 40, 60, 80, 100])
        km_ax.set_ylabel(ylabel)
        km_ax.spines['top'].set_visible(False)
        km_ax.spines['right'].set_visible(False)
        km_ax.tick_params(direction='out', length=3, width=0.8)
        leg = km_ax.get_legend()
        if leg:
            leg.remove()
        km_ax.text(*p_loc, _fmt_logrank_p(p), transform=km_ax.transAxes, ha='left', va='bottom', color='black')

        # ---- Cox / 统计表 ----
        stats_rows = []
        if has_cox:
            include_cox, cox_df = False, None
            try:
                if cox_fn is None:
                    from biokit.analysis import cox as cox_fn
                cox_df, _, _ = cox_fn(grouped_df[[time, status, groupby]], time=time, status=status, mod='multiple',
                                      ref_dict={groupby: cox_ref or groups[0]})
                include_cox = True
            except ConvergenceError:
                warnings.warn(f'{groupby}/{time}/{status} Cox 不收敛，仅显示中位生存')
            except Exception as e:
                warnings.warn(f'Cox 失败（{e}），仅显示中位生存')

            header = ['Group', 'N', 'Event', 'Median (95% CI)']
            if include_cox:
                header += ['HR (95% CI)', 'Cox P']
            ref = cox_ref or groups[0]
            cell_text = []
            for g in groups:
                mask = grouped_df[groupby] == g
                n_g = int(mask.sum())
                e_g = int((mask & (grouped_df[status] == 1)).sum())
                med = _fmt_median_ci(kmf_dict[g])
                row = [str(g), str(n_g), str(e_g), med]
                hr_str, cp = '', np.nan
                if include_cox:
                    if g == ref:
                        row += ['Ref', '\u2014']
                    else:
                        try:
                            hr, lo, hi, cp = cox_df.loc[
                                (groupby, g), ['HR', 'HR(95CI-Low)', 'HR(95CI-High)', 'p-value']]
                            hr_str = f'{hr:.2f} ({lo:.2f}\u2013{hi:.2f})'
                            row += [hr_str, _fmt_cox_p(cp)]
                        except Exception:
                            row += ['NA', 'NA']
                cell_text.append(row)
                stats_rows.append(
                    {'group': g, 'n': n_g, 'events': e_g, 'median_ci': med, 'hr': hr_str, 'cox_p': cp, 'logrank_p': p})
            _draw_stats_table(cox_ax, header, cell_text, groups, color_dict, base_fontsize)
        else:
            for g in groups:
                mask = grouped_df[groupby] == g
                stats_rows.append(
                    {'group': g, 'n': int(mask.sum()), 'events': int((mask & (grouped_df[status] == 1)).sum()),
                        'median_ci': _fmt_median_ci(kmf_dict[g]), 'logrank_p': p})

        # ---- Number at risk + 底部时间轴 ----
        risk_ax.set_xlim(km_ax.get_xlim())
        risk_ax.set_ylim(-0.5, n - 0.5)
        risk_ax.set_xticks(ticks)
        risk_ax.set_xlabel(xlabel or time)
        risk_ax.set_yticks([])
        for sp in ('top', 'left', 'right'):
            risk_ax.spines[sp].set_visible(False)
        risk_ax.spines['bottom'].set_linewidth(0.8)
        risk_ax.tick_params(direction='out', length=3, width=0.8)
        # risk_ax.text(0.0, 1.0, 'Number at risk', transform=risk_ax.transAxes,
        #              ha='left', va='bottom', fontweight='bold')
        for idx, g in enumerate(groups):
            y = n - 1 - idx
            sg = grouped_df[grouped_df[groupby] == g]
            risk_ax.text(-0.05, y, str(g), transform=risk_ax.get_yaxis_transform(), ha='right', va='center',
                         color=color_dict[g], clip_on=False)
            for t in ticks:
                risk_ax.text(t, y, str(int((sg[time] >= t).sum())), ha='center', va='center')

        if title:
            band = title_h / fig_h
            fig.subplots_adjust(left=0.20, right=0.96, top=1 - band, bottom=0.10)
            fig.suptitle(title, y=1 - band * 0.45, fontsize=base_fontsize + 3, fontweight='bold')
        else:
            fig.subplots_adjust(left=0.20, right=0.96, top=0.97, bottom=0.10)

    import pandas as pd
    return pd.DataFrame(stats_rows), p, fig
