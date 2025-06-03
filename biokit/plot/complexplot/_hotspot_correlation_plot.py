# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list


def hotspot_correlation_plot(local_correlation_z, modules, linkage, mod_cmap='Set2', vmin=-8, vmax=8, z_cmap='RdBu_r',
                           yticklabels=False, hull=False):
    colors = sns.color_palette(mod_cmap, n_colors=10)
    module_colors = {i: colors[(i - 1) % len(colors)] for i in modules.unique()}
    module_colors[-1] = '#ffffff'

    row_colors1 = pd.Series([module_colors[i] for i in modules], index=local_correlation_z.index, )

    row_colors = pd.DataFrame({"Modules": row_colors1, })

    cm = sns.clustermap(local_correlation_z, row_linkage=linkage, col_linkage=linkage, vmin=vmin, vmax=vmax,
                        cmap=z_cmap, xticklabels=False, yticklabels=yticklabels, row_colors=row_colors,
                        rasterized=True, )

    fig = plt.gcf()
    plt.sca(cm.ax_heatmap)
    plt.ylabel("")
    plt.xlabel("")

    cm.ax_row_dendrogram.remove()

    # Add 'module X' annotations
    ii = leaves_list(linkage)
    mod_reordered = modules.iloc[ii]
    mod_map = {}
    y = np.arange(modules.size)

    for x in mod_reordered.unique():
        if x == -1:
            continue
        mod_map[x] = y[mod_reordered == x].mean()

    plt.sca(cm.ax_row_colors)
    for mod, mod_y in mod_map.items():
        plt.text(-.5, y=mod_y, s="Module {}".format(mod), horizontalalignment='right', verticalalignment='center')
    plt.xticks([])

    if hull:

        for mod, mod_y in mod_map.items():
            print(mod)
            edge_length = modules.value_counts()[mod] / 2
            center = mod_y
            x1, y1 = center - edge_length, center - edge_length
            x2, y2 = center - edge_length, center + edge_length
            x3, y3 = center + edge_length, center + edge_length
            x4, y4 = center + edge_length, center - edge_length
            cm.ax_heatmap.plot([x1, x2, x3, x4, x1], [y1, y2, y3, y4, y1], color=module_colors[mod], lw=2)

    # Find the colorbar 'child' and modify
    min_delta = 1e99
    min_aa = None
    for aa in fig.get_children():
        try:
            bbox = aa.get_position()
            delta = (0 - bbox.xmin) ** 2 + (1 - bbox.ymax) ** 2
            if delta < min_delta:
                delta = min_delta
                min_aa = aa
        except AttributeError:
            pass

    min_aa.set_ylabel('Z-Scores')
    min_aa.yaxis.set_label_position("left")
