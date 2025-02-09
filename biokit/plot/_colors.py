# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from itertools import combinations

# %% ∑÷¿‡—’…´
color_style = {'style1': ['#000000', '#13213c', '#fca311', '#e5e5e5'],
               'style2': ['#f4f1de', '#df7a5e', '#3c405b', '#82b29a', '#f2cc8e'],
               'style3': ['#264653', '#2a9d8e', '#e9c46b', '#f3a261', '#e66f51'],
               'style4': ['#f66f69', '#feb3ae', '#fff4f2', '#1597a5', '#0e606b', '#ffc24b'],
               'style5': ['#90c9e7', '#219ebc', '#136783', '#02304a', '#feb705', '#ff9e02', '#fa8600'],
               'style6': ['#73bad6', '#0d4c6d', '#033250', '#02263e', '#ef4143', '#bf1e2e', '#c4323f'],
               'style7': ['#e73847', '#f0faef', '#a8dadb', '#457b9d', '#1d3557'],
               'style8': ['#b7b5a0', '#44757a', '#452a3d', '#d44c3c', '#dd6c4c', '#e5855d', '#eed5b7'],
               'style9': ['#4eab90', '#8eb69c', '#edddc3', '#eebf6d', '#d94f33', '#834026'],
               'style10': ['#db3124', '#fc8c5a', '#ffdf92', '#e6f1f3', '#90bee0', '#4b10b2'],
               'style11': ['#cb997e', '#ddbea9', '#fde8d5', '#b8b7a3', '#a5a58d', '#6b705c'],
               'style12': ['#fbf0c3', '#54686f', '#e57b7f', '#9e3150', '#87bba4'],
               'style13': ['#90c9e6', '#219ebc', '#023047', '#ffb703', '#fb8402'],
               'style14': ['#010713', '#01244c', '#013565', '#ffc300', '#fb8402'],
               'style15': ['#780001', '#c11221', '#fef0d5', '#002f49', '#669bbb'],
               'style16': ['#264653', '#287271', '#2a9d8c', '#8ab07d', '#e9c46b', '#f3a261', '#e66f51'],
               'style17': ['#44045a', '#413e85', '#30688d', '#1f928b', '#35b777', '#91d542', '#f8e620'],
               'style18': ['#000000', '#3a0064', '#7a1b6d', '#bd3752', '#ed6825', '#fbb41a'],
               'style19': ['#000000', '#410e73', '#8c2a81', '#df4a68', '#fc9a6b', '#fcf8bb'],
               'style20': ['#16068a', '#6300a9', '#9e189d', '#cc4975', '#ec7853', '#fdb32e'],
               'style21': ['#a40545', '#f46f44', '#fdd985', '#e9fea1', '#7fcba4', '#4b65af'],
               'biomega1': ['#a2d2e7', '#67a8cd', '#ffc17f', '#cf9f88', '#6fb3a8', '#b3e19b', '#50aa4b', '#ff9d9f',
                            '#f36569', '#3581b7', '#cdb6da', '#704ba3', '#9a7fbd', '#dba9a8', '#e43030', '#e99b78',
                            '#ff8831'],
               'biomega2': ['#9ec6dc', '#f4acac', '#fdcbaa', '#d7c1e8', '#b2d9d9', '#f2a3ba', '#f4da69', '#6ebdbd',
                            '#f9b17d', '#ddda91']
               }


def scicolors(n_colors, style='style1'):
    if n_colors <= len(color_style[style]):
        fold = n_colors // len(color_style[style])
        colors = color_style[style] * (fold + 1)
    elif n_colors == len(color_style[style]):
        colors = color_style[style]
    else:
        step = len(color_style[style]) // n_colors
        colors = color_style[style][:step:]
    return colors[:n_colors]


def show_scicolors():
    max_n_colors = max([len(i) for i in color_style.values()])
    fig, ax = plt.subplots(figsize=(max_n_colors, len(color_style)))
    for i, key in enumerate(color_style):
        for j, color in enumerate(color_style[key]):
            ax.barh(y=i, width=0.8, height=0.8, left=j, color=color)
    ax.set_xticks([])
    ax.set_yticks(range(len(color_style)))
    ax.set_yticklabels(color_style.keys(), fontsize=16)
    ax.spines.clear()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()


# %% cmap
cmap_colors_dict = {'cmap1': ['#27787F', 'white', '#D26156'],
                    'cmap2': ['#277F62', 'white', '#D2B356'],
                    'cmap3': ['#AF98D5', 'white', '#E6C35E'],
                    }
low_colors = ['#013565', '#13213c', '#219ebc', '#264653', '#27787F', '#35b777', '#4b65af', '#a8dadb']
high_colors = ['#D26156', '#df4a68', '#df7a5e', '#e66f51', '#e9c46b', '#f46f44', '#fca311', '#ffc300']

i = 4
for low_color in low_colors:
    for high_color in high_colors:
        cmap_colors_dict[f'cmap{i}'] = [low_color, 'white', high_color]
        i += 1
#

twocolor_cmap_colors_dict = {'cmap1': ['#27787F', 'white', '#D26156'], 'cmap2': ['#277F62', 'white', '#D2B356'],
                             'cmap3': ['#AF98D5', 'white', '#E6C35E'], 'cmap4': ['#013565', 'white', '#D26156'],
                             'cmap5': ['#013565', 'white', '#df4a68'], 'cmap6': ['#013565', 'white', '#df7a5e'],
                             'cmap7': ['#013565', 'white', '#e66f51'], 'cmap8': ['#013565', 'white', '#e9c46b'],
                             'cmap9': ['#013565', 'white', '#f46f44'], 'cmap10': ['#013565', 'white', '#fca311'],
                             'cmap11': ['#013565', 'white', '#ffc300'], 'cmap12': ['#13213c', 'white', '#D26156'],
                             'cmap13': ['#13213c', 'white', '#df4a68'], 'cmap14': ['#13213c', 'white', '#df7a5e'],
                             'cmap15': ['#13213c', 'white', '#e66f51'], 'cmap16': ['#13213c', 'white', '#e9c46b'],
                             'cmap17': ['#13213c', 'white', '#f46f44'], 'cmap18': ['#13213c', 'white', '#fca311'],
                             'cmap19': ['#13213c', 'white', '#ffc300'], 'cmap20': ['#219ebc', 'white', '#D26156'],
                             'cmap21': ['#219ebc', 'white', '#df4a68'], 'cmap22': ['#219ebc', 'white', '#df7a5e'],
                             'cmap23': ['#219ebc', 'white', '#e66f51'], 'cmap24': ['#219ebc', 'white', '#e9c46b'],
                             'cmap25': ['#219ebc', 'white', '#f46f44'], 'cmap26': ['#219ebc', 'white', '#fca311'],
                             'cmap27': ['#219ebc', 'white', '#ffc300'], 'cmap28': ['#264653', 'white', '#D26156'],
                             'cmap29': ['#264653', 'white', '#df4a68'], 'cmap30': ['#264653', 'white', '#df7a5e'],
                             'cmap31': ['#264653', 'white', '#e66f51'], 'cmap32': ['#264653', 'white', '#e9c46b'],
                             'cmap33': ['#264653', 'white', '#f46f44'], 'cmap34': ['#264653', 'white', '#fca311'],
                             'cmap35': ['#264653', 'white', '#ffc300'], 'cmap36': ['#27787F', 'white', '#D26156'],
                             'cmap37': ['#27787F', 'white', '#df4a68'], 'cmap38': ['#27787F', 'white', '#df7a5e'],
                             'cmap39': ['#27787F', 'white', '#e66f51'], 'cmap40': ['#27787F', 'white', '#e9c46b'],
                             'cmap41': ['#27787F', 'white', '#f46f44'], 'cmap42': ['#27787F', 'white', '#fca311'],
                             'cmap43': ['#27787F', 'white', '#ffc300'], 'cmap44': ['#35b777', 'white', '#D26156'],
                             'cmap45': ['#35b777', 'white', '#df4a68'], 'cmap46': ['#35b777', 'white', '#df7a5e'],
                             'cmap47': ['#35b777', 'white', '#e66f51'], 'cmap48': ['#35b777', 'white', '#e9c46b'],
                             'cmap49': ['#35b777', 'white', '#f46f44'], 'cmap50': ['#35b777', 'white', '#fca311'],
                             'cmap51': ['#35b777', 'white', '#ffc300'], 'cmap52': ['#4b65af', 'white', '#D26156'],
                             'cmap53': ['#4b65af', 'white', '#df4a68'], 'cmap54': ['#4b65af', 'white', '#df7a5e'],
                             'cmap55': ['#4b65af', 'white', '#e66f51'], 'cmap56': ['#4b65af', 'white', '#e9c46b'],
                             'cmap57': ['#4b65af', 'white', '#f46f44'], 'cmap58': ['#4b65af', 'white', '#fca311'],
                             'cmap59': ['#4b65af', 'white', '#ffc300'], 'cmap60': ['#a8dadb', 'white', '#D26156'],
                             'cmap61': ['#a8dadb', 'white', '#df4a68'], 'cmap62': ['#a8dadb', 'white', '#df7a5e'],
                             'cmap63': ['#a8dadb', 'white', '#e66f51'], 'cmap64': ['#a8dadb', 'white', '#e9c46b'],
                             'cmap65': ['#a8dadb', 'white', '#f46f44'], 'cmap66': ['#a8dadb', 'white', '#fca311'],
                             'cmap67': ['#a8dadb', 'white', '#ffc300']}

cmap_colors_dict = {**color_style, **twocolor_cmap_colors_dict}


def scicmap(style='cmap1'):
    colors = cmap_colors_dict[style]
    cmap = LinearSegmentedColormap.from_list(style, colors)
    return cmap


def show_scicmap(figsize=None):
    n_cmaps = len(cmap_colors_dict)
    nax_cols = math.ceil(math.sqrt(n_cmaps / 8))
    nax_rows = math.ceil(n_cmaps / nax_cols)

    if not figsize:
        figsize = (nax_cols * 4, nax_rows)

    fig, axes = plt.subplots(nax_rows, nax_cols, figsize=figsize)
    for i, (key, value) in enumerate(cmap_colors_dict.items()):
        ax = axes.flatten()[i]
        cmap = LinearSegmentedColormap.from_list(key, value)
        ax.imshow([np.linspace(0, 1, 100)] * 10, cmap=cmap)
        ax.axis('off')
        ax.set_title(key)

    for ax in axes.flatten()[n_cmaps:]:
        ax.axis('off')
    plt.tight_layout()
    plt.show()
