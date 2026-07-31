import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from biokit.data_convert import p2text as p2text_func


def testbox(data, y, x=None, groupby=None, groups=None, kind='box', testfunc=ttest_ind, cutoff=None, colors=None,
            cutoff_color=None, p2text=True, showfliers=False, ax=None, width=0.8, ylim=None, **kwargs):
    """
    绘制分组比较图 (box / violin / bar)，并在 hue 组之间添加显著性标注

    修正点：
    1. x 顺序由 order 显式控制（柱子/xtick/显著性永远一致）
    2. bar 模式下 ymax / base_y 由 errorbar 上限决定（而不是原始值最大）
    """

    data = data.copy()

    # ---------- x ----------
    if x is None:
        data['_x'] = 'all'
        x = '_x'

    order = data[x].unique().tolist()

    # ---------- hue ----------
    if groups is None:
        groups = data[groupby].unique().tolist()

    k = len(groups)

    # ---------- ax ----------
    if ax is None:
        fig, ax = plt.subplots(figsize=(2.5 * len(order), 6))

    # ---------- colors ----------
    if colors is None:
        colors = sns.color_palette("Set2", n_colors=k)
    palette = dict(zip(groups, colors))

    # ---------- cutoff ----------
    if cutoff is None:
        cutoff = {0.05: '*', 0.01: '**', 0.001: '***', 0.0001: '****'}

    if cutoff_color is None:
        cutoff_color = {'*': 'orange', '**': 'darkorange', '***': 'orangered', '****': 'red', 'ns': 'deepskyblue'}

    # ---------- plot ----------
    plot_kws = dict(data=data, x=x, y=y, order=order,  # ⭐ 顺序锁死
                    hue=groupby, hue_order=groups, palette=palette, ax=ax, legend=False, )

    if kind == 'box':
        sns.boxplot(**plot_kws, width=width, showfliers=showfliers, **kwargs)

    elif kind == 'violin':
        sns.violinplot(**plot_kws, width=width, cut=0, **kwargs)

    elif kind == 'bar':
        sns.barplot(capsize=0.2, gap=0.2, **plot_kws, **kwargs)

    # ---------- x positions ----------
    x_levels = np.arange(len(order))
    offsets = (np.arange(k) + 0.5) * (width / k) - width / 2
    positions_per_x = [xj + offsets for xj in x_levels]

    # ---------- y base / ymax ----------
    # ---------- y base / ymax ----------
    if kind == 'bar':
        base_y_per_x = [0.0] * len(order)

        # 1) bar 本体
        for container in ax.containers:
            if not hasattr(container, 'patches'):
                continue
            for patch in container.patches:
                try:
                    x_mid = patch.get_x() + patch.get_width() / 2
                    y_top = patch.get_y() + patch.get_height()
                except Exception:
                    continue
                x_idx = int(round(x_mid))
                if 0 <= x_idx < len(base_y_per_x):
                    base_y_per_x[x_idx] = max(base_y_per_x[x_idx], float(y_top))

        # 2) errorbar from collections
        for coll in ax.collections:
            if not hasattr(coll, 'get_segments'):
                continue
            try:
                segs = coll.get_segments()
            except Exception:
                continue
            for seg in segs:
                if seg is None or len(seg) == 0:
                    continue
                y_upper = float(np.max(seg[:, 1]))
                x_mid = float(np.mean(seg[:, 0]))
                x_idx = int(round(x_mid))
                if 0 <= x_idx < len(base_y_per_x):
                    base_y_per_x[x_idx] = max(base_y_per_x[x_idx], y_upper)

        # 3) errorbar from lines
        for line in ax.lines:
            xdata = np.asarray(line.get_xdata(), dtype=float)
            ydata = np.asarray(line.get_ydata(), dtype=float)
            if xdata.size == 0 or ydata.size == 0:
                continue
            if not np.isfinite(xdata).all() or not np.isfinite(ydata).all():
                continue
            x_mid = float(np.mean(xdata))
            y_upper = float(np.max(ydata))
            x_idx = int(round(x_mid))
            if 0 <= x_idx < len(base_y_per_x):
                base_y_per_x[x_idx] = max(base_y_per_x[x_idx], y_upper)

        ymin = 0
        ymax = max(base_y_per_x) if len(base_y_per_x) > 0 else 1.0

    else:
        ymin = data[y].min()
        ymax = data[y].max()
        base_y_per_x = [ymax] * len(order)

    ylength = ymax - ymin if ymax > ymin else max(abs(ymax), 1.0)
    pstep = ylength * 0.05

    # ---------- significance ----------
    if testfunc is not None:
        for j, xval in enumerate(order):
            subset = data[data[x] == xval]
            cur_top = base_y_per_x[j] * 1.05
            for i in range(k - 1):
                for m in range(i + 1, k):
                    g1, g2 = groups[i], groups[m]
                    v1 = subset[subset[groupby] == g1][y].dropna()
                    v2 = subset[subset[groupby] == g2][y].dropna()
                    if len(v1) == 0 or len(v2) == 0:
                        continue

                    pval = testfunc(v2, v1).pvalue
                    if p2text:
                        ptxt = p2text_func(pval, cutoff)
                        color = cutoff_color.get(ptxt, 'black')
                    else:
                        ptxt = f'{pval:.3g}'
                        color = 'red' if pval < 0.05 else 'black'

                    cur_top += pstep
                    x1, x2 = positions_per_x[j][i], positions_per_x[j][m]

                    ax.plot([x1, x1, x2, x2], [cur_top - 0.3 * pstep, cur_top, cur_top, cur_top - 0.3 * pstep],
                            color=color)
                    ax.text((x1 + x2) / 2, cur_top, ptxt, ha='center', va='bottom', color=color)

                    ymax = max(ymax, cur_top)

    # ---------- ticks ----------
    if x == '_x' and groupby is not None:
        # 单一 x，但用 hue 名作为 xticklabels
        ax.set_xticks(positions_per_x[0])
        ax.set_xticklabels(groups)
    else:
        ax.set_xticks(x_levels)
        ax.set_xticklabels(order)

    # ---------- ylim ----------
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim(ymin - 0.05 * ylength, ymax + 2 * pstep)

    return ax


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from biokit.data_convert import p2text as p2text_func


def testbox(data, y, x=None, x_order=None, groupby=None, groups=None, kind='box', testfunc=ttest_ind, cutoff=None,
            colors=None, cutoff_color=None, p2text=True, showfliers=False, ax=None, width=0.8, ylim=None, **kwargs):
    data = data.copy()
    if x is None:
        data['_x'] = 'all'
        x = '_x'
    if not x_order:
        order = data[x].unique().tolist()
    else:
        order = x_order
    if groups is None:
        groups = data[groupby].unique().tolist()
    k = len(groups)
    if ax is None:
        fig, ax = plt.subplots(figsize=(2.5 * len(order), 6))
    if colors is None:
        colors = sns.color_palette("Set2", n_colors=k)
    palette = dict(zip(groups, colors))
    if cutoff is None:
        cutoff = {0.05: '*', 0.01: '**', 0.001: '***', 0.0001: '****'}
    if cutoff_color is None:
        cutoff_color = {'*': 'orange', '**': 'darkorange', '***': 'orangered', '****': 'red', 'ns': 'deepskyblue'}

    plot_kws = dict(data=data, x=x, y=y, order=order, hue=groupby, hue_order=groups, palette=palette, ax=ax,
                    legend=False)
    if kind == 'box':
        sns.boxplot(width=width, showfliers=showfliers, gap=0.2, **kwargs, **plot_kws, )
    elif kind == 'violin':
        sns.violinplot(cut=0, linewidth=0, gap=0.2, **plot_kws, )
        sns.boxplot(width=0.8, showfliers=False, color=None, gap=0.7, boxprops={"facecolor": "none"},
                    capprops={"linewidth": 0}, **plot_kws, )
    elif kind == 'bar':
        sns.barplot(capsize=0.2, gap=0.2, **kwargs, **plot_kws)
    elif kind == 'strip':
        sns.stripplot(size=5, jitter=True, **kwargs, **plot_kws, )
    elif kind == 'swarm':
        sns.swarmplot(size=5, **kwargs, **plot_kws, )

    x_levels = np.arange(len(order))
    offsets = (np.arange(k) + 0.5) * (width / k) - width / 2
    positions_per_x = [xj + offsets for xj in x_levels]

    def nearest_x_idx(xpos):
        return int(np.argmin(np.abs(x_levels - xpos)))

    if kind == 'bar':
        data_ymin = float(np.nanmin(data[y])) if np.isfinite(np.nanmin(data[y])) else 0.0
        ymin = min(0.0, data_ymin)
        base_y_per_x = np.full(len(order), ymin, dtype=float)

        for patch in ax.patches:
            try:
                x_mid = patch.get_x() + patch.get_width() / 2
                y0 = patch.get_y()
                y1 = y0 + patch.get_height()
                y_top = max(y0, y1)
            except Exception:
                continue
            if not np.isfinite(x_mid) or not np.isfinite(y_top):
                continue
            x_idx = nearest_x_idx(x_mid)
            base_y_per_x[x_idx] = max(base_y_per_x[x_idx], float(y_top))

        for coll in ax.collections:
            if not hasattr(coll, 'get_segments'):
                continue
            try:
                segs = coll.get_segments()
            except Exception:
                continue
            for seg in segs:
                if seg is None or len(seg) == 0:
                    continue
                seg = np.asarray(seg, dtype=float)
                if seg.ndim != 2 or seg.shape[1] != 2:
                    continue
                if not np.isfinite(seg).all():
                    continue
                x_mid = float(np.mean(seg[:, 0]))
                y_top = float(np.max(seg[:, 1]))
                x_idx = nearest_x_idx(x_mid)
                base_y_per_x[x_idx] = max(base_y_per_x[x_idx], y_top)

        for line in ax.lines:
            xdata = np.asarray(line.get_xdata(), dtype=float)
            ydata = np.asarray(line.get_ydata(), dtype=float)
            if xdata.size == 0 or ydata.size == 0:
                continue
            if not np.isfinite(xdata).all() or not np.isfinite(ydata).all():
                continue
            x_mid = float(np.mean(xdata))
            y_top = float(np.max(ydata))
            x_idx = nearest_x_idx(x_mid)
            base_y_per_x[x_idx] = max(base_y_per_x[x_idx], y_top)

        ymax = float(np.max(base_y_per_x)) if len(base_y_per_x) > 0 else 1.0
    else:
        ymin = float(data[y].min())
        ymax = float(data[y].max())
        base_y_per_x = np.full(len(order), ymax, dtype=float)

    ylength = ymax - ymin
    pstep = ylength * 0.05
    raw_cur_top = ylim[1] if ylim is not None else ymax + pstep
    if testfunc is not None:
        for j, xval in enumerate(order):
            cur_top = raw_cur_top
            subset = data[data[x] == xval]
            if not np.isfinite(ylength) or ylength <= 0:
                ylength = max(abs(ymax), 1.0)
            for i in range(k - 1):
                for m in range(i + 1, k):
                    g1, g2 = groups[i], groups[m]
                    v1 = subset[subset[groupby] == g1][y].dropna()
                    v2 = subset[subset[groupby] == g2][y].dropna()
                    if len(v1) == 0 or len(v2) == 0:
                        continue
                    pval = testfunc(v2, v1).pvalue
                    if p2text:
                        ptxt = p2text_func(pval, cutoff)
                        color = cutoff_color.get(ptxt, 'black')
                    else:
                        ptxt = f'{pval:.3g}'
                        color = 'red' if pval < 0.05 else 'black'
                    color = 'black'
                    x1, x2 = positions_per_x[j][i], positions_per_x[j][m]
                    y_bracket = cur_top + pstep
                    ax.plot([x1, x1, x2, x2], [y_bracket - 0.3 * pstep, y_bracket, y_bracket, y_bracket - 0.3 * pstep],
                            color=color)
                    ax.text((x1 + x2) / 2, y_bracket, ptxt, ha='center', va='bottom', color=color)
                    cur_top = y_bracket
                    ymax = max(ymax, y_bracket)

    if x == '_x' and groupby is not None:
        ax.set_xticks(positions_per_x[0])
        ax.set_xticklabels(groups)
        ax.set_xlabel(groupby)
    else:
        print(x)
        ax.set_xticks(x_levels)
        ax.set_xticklabels(order)

    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1] + 4 * pstep)
    else:
        cur_top = raw_cur_top
        ax.set_ylim(ymin - 0.05 * ylength, cur_top + 2 * pstep)

    # 去掉右上边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax
