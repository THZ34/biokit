# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import sys

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_template import FigureCanvas
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from numpy import unique


def gsea_pathwayplot(df, column='Adjusted P-value', title='', color='-log10 Padj', cutoff=0.05, top_term=10,
                     sizes=None, norm=None, legend=True, figsize=(6, 5.5),
                     cmap='RdBu_r', ofname=None, **kwargs):
    """Visualize enrichr results.

    :param df: GSEApy DataFrame results.
    :param column: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title
    :param cutoff: p-adjust cut-off.
    :param top_term: number of enriched terms to show.
    :param ascending: bool, the order of y axis.
    :param sizes: tuple, (min, max) scatter size. Not functional for now
    :param norm: maplotlib.colors.Normalize object.
    :param legend: bool, whether to show legend.
    :param figsize: tuple, figure size.
    :param cmap: matplotlib colormap
    :param ofname: output file name. If None, don't save figure

    """

    colname = column
    # sorting the dataframe for better visualization
    if colname in ['Adjusted P-value', 'P-value']:
        df = df[df[colname] <= cutoff]
        if len(df) < 1:
            msg = "Warning: No enrich terms when cutoff = %s" % cutoff
            return msg
        df = df.assign(logAP=lambda x: - x[colname].apply(np.log10))
        colname = 'logAP'
    df = df.sort_values(by=colname).iloc[-top_term:, :]
    #
    temp = df['Overlap'].str.split("/", expand=True).astype(int)
    df = df.assign(Hits=temp.iloc[:, 0], Background=temp.iloc[:, 1])
    df = df.assign(Hits_ratio=lambda x: x.Hits / x.Background)
    # x axis values
    x = df.loc[:, colname].values
    # combined_score = df['Combined Score'].round().astype('int')
    combined_score = df.loc[:, color].values
    # y axis index and values
    y = [i for i in range(0, len(df))]
    ylabels = df['Term'].values
    # Normalise to [0,1]
    # b = (df['Count']  - df['Count'].min())/ np.ptp(df['Count'])
    # area = 100 * b

    # control the size of scatter and legend marker
    levels = numbers = np.sort(df.Hits.unique())
    if norm is None:
        norm = Normalize()
    elif isinstance(norm, tuple):
        norm = Normalize(*norm)
    elif not isinstance(norm, Normalize):
        err = ("``size_norm`` must be None, tuple, "
               "or Normalize object.")
        raise ValueError(err)
    min_width, max_width = np.r_[20, 100] * plt.rcParams["lines.linewidth"]
    norm.clip = True
    if not norm.scaled():
        norm(np.asarray(numbers))
    size_limits = norm.vmin, norm.vmax
    scl = norm(numbers)
    widths = np.asarray(min_width + scl * (max_width - min_width))
    if scl.mask.any():
        widths[scl.mask] = 0
    sizes = dict(zip(levels, widths))
    df['sizes'] = df.Hits.map(sizes)
    area = df['sizes'].values

    # creat scatter plot
    if hasattr(sys, 'ps1') and (ofname is None):
        # working inside python console, show figure
        fig, ax = plt.subplots(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
    vmin = np.percentile(combined_score.min(), 2)
    vmax = np.percentile(combined_score.max(), 98)
    sc = ax.scatter(x=x, y=y, s=area, edgecolors='face', c=combined_score,
                    cmap=cmap, vmin=vmin, vmax=vmax)

    if column in ['Adjusted P-value', 'P-value']:
        xlabel = "-log$_{10}$(%s)" % column
    else:
        xlabel = column
    ax.set_xlabel(xlabel, fontsize=14, fontweight='bold')
    ax.yaxis.set_major_locator(plt.FixedLocator(y))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(ylabels))
    ax.set_yticklabels(ylabels, fontsize=16)

    # ax.set_ylim([-1, len(df)])
    ax.grid()
    # colorbar
    cax = fig.add_axes([0.9, 0.20, 0.03, 0.22])
    cbar = fig.colorbar(sc, cax=cax, )
    cbar.ax.tick_params(right=True)
    cbar.ax.set_title(color, loc='left', fontsize=12)

    # for terms less than 3
    if len(df) >= 3:
        # find the index of the closest value to the median
        idx = [area.argmax(), np.abs(area - area.mean()).argmin(), area.argmin()]
        idx = unique(idx)
    else:
        idx = df.index.values
    label = df.iloc[idx, df.columns.get_loc('Hits')]

    if legend:
        handles, _ = ax.get_legend_handles_labels()
        legend_markers = []
        for ix in idx:
            legend_markers.append(ax.scatter([], [], s=area[ix], c='k'))
        # artist = ax.scatter([], [], s=size_levels,)
        ax.legend(legend_markers, label, title='Gene Counts')
    ax.set_title(title, fontsize=20, fontweight='bold')

    if ofname is not None:
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches='tight', dpi=300)
        return
    return ax


def kobas_dotplot(df, column='Adjusted P-value', title='', color='-log10 Padj', cutoff=0.05, top_term=10,
                  sizes=None, norm=None, legend=True, figsize=(6, 5.5),
                  cmap='RdBu_r', ofname=None, **kwargs):
    """Visualize enrichr results.

    :param df: GSEApy DataFrame results.
    :param column: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title
    :param cutoff: p-adjust cut-off.
    :param top_term: number of enriched terms to show.
    :param ascending: bool, the order of y axis.
    :param sizes: tuple, (min, max) scatter size. Not functional for now
    :param norm: maplotlib.colors.Normalize object.
    :param legend: bool, whether to show legend.
    :param figsize: tuple, figure size.
    :param cmap: matplotlib colormap
    :param ofname: output file name. If None, don't save figure

    """

    colname = column
    # sorting the dataframe for better visualization
    if colname in ['Correct P-value', 'P-value']:
        df = df[df[colname] <= cutoff]
        if len(df) < 1:
            msg = "Warning: No enrich terms when cutoff = %s" % cutoff
            return msg
        df = df.assign(logAP=lambda x: - x[colname].apply(np.log10))
        colname = 'logAP'
    df = df.sort_values(by=colname).iloc[-top_term:, :]
    #
    temp = df[['Input number', 'Background number']]
    df = df.assign(Hits=temp.iloc[:, 0], Background=temp.iloc[:, 1])
    df = df.assign(Hits_ratio=lambda x: x.Hits / x.Background)
    # x axis values
    x = df.loc[:, colname].values
    combined_score = df.loc[:, color].values
    # y axis index and values
    y = [i for i in range(0, len(df))]
    ylabels = df['Term'].values
    # Normalise to [0,1]
    # b = (df['Count']  - df['Count'].min())/ np.ptp(df['Count'])
    # area = 100 * b

    # control the size of scatter and legend marker
    levels = numbers = np.sort(df.Hits.unique())
    if norm is None:
        norm = Normalize()
    elif isinstance(norm, tuple):
        norm = Normalize(*norm)
    elif not isinstance(norm, Normalize):
        err = ("``size_norm`` must be None, tuple, "
               "or Normalize object.")
        raise ValueError(err)
    min_width, max_width = np.r_[20, 100] * plt.rcParams["lines.linewidth"]
    norm.clip = True
    if not norm.scaled():
        norm(np.asarray(numbers))
    size_limits = norm.vmin, norm.vmax
    scl = norm(numbers)
    widths = np.asarray(min_width + scl * (max_width - min_width))
    if scl.mask.any():
        widths[scl.mask] = 0
    sizes = dict(zip(levels, widths))
    df['sizes'] = df.Hits.map(sizes)
    area = df['sizes'].values

    # creat scatter plot
    if hasattr(sys, 'ps1') and (ofname is None):
        # working inside python console, show figure
        fig, ax = plt.subplots(figsize=figsize)
    else:
        # If working on commandline, don't show figure
        fig = Figure(figsize=figsize)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
    vmin = np.percentile(combined_score.min(), 2)
    vmax = np.percentile(combined_score.max(), 98)
    sc = ax.scatter(x=x, y=y, s=area, edgecolors='face', c=combined_score,
                    cmap=cmap, vmin=vmin, vmax=vmax)

    if column in ['Adjusted P-value', 'P-value']:
        xlabel = "-log$_{10}$(%s)" % column
    else:
        xlabel = column
    ax.set_xlabel(xlabel, fontsize=14, fontweight='bold')
    ax.yaxis.set_major_locator(plt.FixedLocator(y))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(ylabels))
    ax.set_yticklabels(ylabels, fontsize=16)

    # ax.set_ylim([-1, len(df)])
    ax.grid()
    # colorbar
    cax = fig.add_axes([0.9, 0.20, 0.03, 0.22])
    cbar = fig.colorbar(sc, cax=cax, )
    cbar.ax.tick_params(right=True)
    cbar.ax.set_title(color, loc='left', fontsize=12)

    # for terms less than 3
    if len(df) >= 3:
        # find the index of the closest value to the median
        idx = [area.argmax(), np.abs(area - area.mean()).argmin(), area.argmin()]
        idx = unique(idx)
    else:
        idx = df.index.values
    label = df.iloc[idx, df.columns.get_loc('Hits')]

    if legend:
        handles, _ = ax.get_legend_handles_labels()
        legend_markers = []
        for ix in idx:
            legend_markers.append(ax.scatter([], [], s=area[ix], c='k'))
        # artist = ax.scatter([], [], s=size_levels,)
        ax.legend(legend_markers, label, title='Gene Counts')
    ax.set_title(title, fontsize=20, fontweight='bold')

    if ofname is not None:
        # canvas.print_figure(ofname, bbox_inches='tight', dpi=300)
        fig.savefig(ofname, bbox_inches='tight', dpi=300)
        return
    return ax


def metascape_dotplot(df, y='Adjusted P-value', title='', color='-log10 Padj', cutoff=0.05, top_term=10,
                      sizes=None, norm=None, legend=True, figsize=(6, 5.5),
                      cmap='RdBu_r', ofname=None, **kwargs):
    """Visualize enrichr results.

    :param df: GSEApy DataFrame results.
    :param y: which column of DataFrame to show. Default: Adjusted P-value
    :param title: figure title
    :param cutoff: p-adjust cut-off.
    :param top_term: number of enriched terms to show.
    :param ascending: bool, the order of y axis.
    :param sizes: tuple, (min, max) scatter size. Not functional for now
    :param norm: maplotlib.colors.Normalize object.
    :param legend: bool, whether to show legend.
    :param figsize: tuple, figure size.
    :param cmap: matplotlib colormap
    :param ofname: output file name. If None, don't save figure

    """

    colname = y
    # sorting the dataframe for better visualization
    df = df.sort_values(by=colname).iloc[-top_term:, :]
    #
    temp = df[['Input number', 'Background number']]
    df = df.assign(Hits=temp.iloc[:, 0], Background=temp.iloc[:, 1])
    df = df.assign(Hits_ratio=lambda x: x.Hits / x.Background)
    # x axis values
    x = df.loc[:, colname].values
    combined_score = df.loc[:, color].values
    # y axis index and values
    y = [i for i in range(0, len(df))]
    ylabels = df['Description'].values

    levels = numbers = np.sort(df.Hits.unique())
    if norm is None:
        norm = Normalize()
    elif isinstance(norm, tuple):
        norm = Normalize(*norm)
    elif not isinstance(norm, Normalize):
        err = ("``size_norm`` must be None, tuple, "
               "or Normalize object.")
        raise ValueError(err)
    min_width, max_width = np.r_[20, 100] * plt.rcParams["lines.linewidth"]
    norm.clip = True
    if not norm.scaled():
        norm(np.asarray(numbers))
    size_limits = norm.vmin, norm.vmax
    scl = norm(numbers)
    widths = np.asarray(min_width + scl * (max_width - min_width))
    if scl.mask.any():
        widths[scl.mask] = 0
    sizes = dict(zip(levels, widths))
    df['sizes'] = df.Hits.map(sizes)
    area = df['sizes'].values

    # creat scatter plot
    fig, ax = plt.subplots(figsize=figsize)

    vmin = np.percentile(combined_score.min(), 2)
    vmax = np.percentile(combined_score.max(), 98)
    sc = ax.scatter(x=x, y=y, s=area, edgecolors='face', c=combined_score,
                    cmap=cmap, vmin=vmin, vmax=vmax)

    if y in ['Adjusted P-value', 'P-value']:
        xlabel = "-log$_{10}$(%s)" % y
    else:
        xlabel = y
    ax.set_xlabel(xlabel, fontsize=14, fontweight='bold')
    ax.yaxis.set_major_locator(plt.FixedLocator(y))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(ylabels))
    ax.set_yticklabels(ylabels, fontsize=16)

    ax.grid()
    # colorbar
    cax = fig.add_axes([0.91, 0.20, 0.02, 0.22])
    cbar = fig.colorbar(sc, cax=cax, )
    cbar.ax.tick_params(right=True)
    cbar.ax.set_title(color, loc='left', fontsize=12)

    # for terms less than 3
    if len(df) >= 3:
        # find the index of the closest value to the median
        idx = [area.argmax(), np.abs(area - area.mean()).argmin(), area.argmin()]
        idx = unique(idx)
    else:
        idx = df.index.values
    label = df.iloc[idx, df.columns.get_loc('Hits')]

    if legend:
        handles, _ = ax.get_legend_handles_labels()
        legend_markers = []
        for ix in idx:
            legend_markers.append(ax.scatter([], [], s=area[ix], c='k'))
        ax.legend(legend_markers, label, title='Gene Counts')
    ax.set_title(title, fontsize=20, fontweight='bold')

    if ofname is not None:
        fig.savefig(ofname, bbox_inches='tight', dpi=300)
        return
    return ax



def gseapy_dotplot(df, x, color, size, cmap, size_scale, size_log):
    xlabel = x
    cbar_title = color

    df.sort_values(by=x, ascending=False, inplace=True)
    norm = Normalize(df[color].min(), df[color].max())
    df['color'] = df[color].apply(lambda x: cmap(norm(x)))
    df['size'] = df[size] / df[size].max()

    if size_log:
        df['size'] = np.log10(df['size'] + 1)
        df['size'] = df['size'] / df['size'].max()
        df['size'] = df['size'] * size_scale

    fig, ax = plt.subplots(figsize=(8, 4))
    ymax = df.shape[0]
    for y, pathway in enumerate(df.index):
        x, size, color = df.loc[pathway, ['Odds Ratio', 'size', 'color']]
        ax.scatter(x, ymax - y, s=size, color=color, edgecolor='black', linewidth=0.5,
                   label=pathway)

    plt.subplots_adjust(left=0.5, right=0.78)
    # 设置colorbar
    cax = fig.add_axes([0.8, 0.15, 0.02, 0.3])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, cax=cax)
    cbar.set_label(cbar_title)
    # 设置yticklabels,xlabel
    ax.set_yticks(range(1, df.shape[0] + 1))
    ax.set_yticklabels(df.index[::-1])
    ax.set_xlabel(xlabel)
    # 设置网格
    ax.grid()
    return ax

def single_dotplot(adata, var, groupby, cmap=None, vmax=None, vmin=None, dot_max=1, dot_min=0, color_scale=None,
                   size_scale=None, ax=None):
    if var.isinstance(dict):
        genes = sum(list(var.values()), [])
    elif var.isinstance(str):
        genes = [var]
    else:
        genes = var
    temp_exp_df = adata[:, genes].to_df()
    mean_exp_df = temp_exp_df.groupby(adata.obs[groupby]).mean()
    size_df = (temp_exp_df > 0).groupby(adata.obs[groupby]).sum() / temp_exp_df.groupby(adata.obs[groupby]).count()
