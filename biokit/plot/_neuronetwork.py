# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import seaborn as sns
import matplotlib.pyplot as plt


def draw_neural_net_double(n_neuro_layers, featurenames=None, yname='y', x_interval=2, neurosize=300, colors=None,
                           arrowcolors=None, arrowwidth=0.3, ax=None):
    """
    绘制神经网络图
    :param n_neuro_layers:
    :param featurenames:
    :param yname:
    :param x_interval:
    :param neurosize:
    :param colors:
    :param arrowcolors:
    :param arrowwidth:
    :param ax:
    :return:
    """
    # 计算参数
    ymax = max(n_neuro_layers) * 2 + 3
    left = 0
    ycenter = ymax / 2
    if colors is None:
        colors = sns.hls_palette(len(n_neuro_layers))
    if arrowcolors is None:
        arrowcolors = colors
    arrowcolors = [''] + arrowcolors
    arrow_length = 0.6
    if featurenames is None:
        featurenames = [f'feature{i + 1}' for i in range(n_neuro_layers[0] * 2)]
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    # 特征层 取n个特征
    arrow_left = left + x_interval * 0.1
    n_neuro = n_neuro_layers[0]
    bottom = ycenter + 1.5
    for gene in featurenames[:n_neuro][::-1]:  # 上半部分基因
        y = bottom + 0.5
        ax.text(left, y, gene, ha='right', va='center')
        ax.arrow(arrow_left, y, arrow_length * x_interval, 0, linewidth=arrowwidth, head_width=0.3, head_length=0.3,
                 fc='black', ec='black', length_includes_head=True)
        bottom += 1
    top = ycenter - 1.5
    for gene in featurenames[-n_neuro:]:  # 下半部分基因
        y = top - 0.5
        ax.text(left, y, gene, ha='right', va='center')
        ax.arrow(arrow_left, y, arrow_length * x_interval, 0, linewidth=arrowwidth, head_width=0.3, head_length=0.3,
                 fc='black', ec='black', length_includes_head=True)
        top -= 1
    for i in [ycenter - 1, ycenter, ycenter + 1]:  # 省略号
        ax.scatter(left - 0.5, i, s=neurosize / 20, c='black', edgecolors='black', zorder=2)

    # 神经网络
    previous_y_list = []
    for n_neuro, color, arrowcolor in zip(n_neuro_layers, colors, arrowcolors):
        left += x_interval
        bottom = ycenter + 1.5
        y_list = []
        for i in range(n_neuro):  # 上半部分神经元
            ax.scatter(left, bottom + 0.5, s=neurosize, c=color, edgecolors='black', zorder=2)
            y_list.append(bottom + 0.5)
            bottom += 1
        top = ycenter - 1.5
        for i in range(n_neuro):  # 下半部分神经元
            ax.scatter(left, top - 0.5, s=neurosize, c=color, edgecolors='black', zorder=2)
            y_list.append(top - 0.5)
            top -= 1
        for y in y_list:  # 连接上一层的箭头
            for previous_y in previous_y_list:
                ax.arrow(left - x_interval, previous_y, x_interval, y - previous_y, linewidth=arrowwidth,
                         head_width=0.1,
                         head_length=0, fc=arrowcolor, ec=arrowcolor, length_includes_head=True)
        for i in [ycenter - 1, ycenter, ycenter + 1]:  # 省略号
            ax.scatter(left, i, s=neurosize / 20, c='black', edgecolors='black', zorder=2)
        previous_y_list = y_list

    # 输出层
    left += x_interval
    y = ycenter
    color = colors[len(n_neuro_layers)]
    arrowcolor = arrowcolors[len(n_neuro_layers)]
    ax.scatter(left, y, s=neurosize, c=color, edgecolors='black', zorder=2)
    for previous_y in previous_y_list:
        ax.arrow(left - x_interval, previous_y, x_interval, y - previous_y, linewidth=arrowwidth, head_width=0.1,
                 head_length=0, fc=arrowcolor, ec=arrowcolor, length_includes_head=True)

    # 输出层标签
    left += x_interval
    arrow_left = left + x_interval * 0.1
    ax.text(left, y, yname, ha='left', va='center')
    ax.arrow(arrow_left - x_interval, ycenter, arrow_length * x_interval, 0, linewidth=arrowwidth, head_width=0.3,
             head_length=0.3, fc='black', ec='black', length_includes_head=True)

    ax.set_xlim(-4, (len(n_neuro_layers) + 3) * x_interval + 1)
    ax.set_ylim(0, ymax)
    ax.axis('off')
    return ax


def draw_neural_net_single(n_neuro_layers, featurenames=None, yname='y', x_interval=2, neurosize=300, colors=None,
                           arrowcolors=None, arrowwidth=0.3, ax=None):
    # 计算参数
    ymax = max(n_neuro_layers)
    left = 0
    ycenter = ymax / 2
    if colors is None:
        colors = sns.hls_palette(len(n_neuro_layers))
    if arrowcolors is None:
        arrowcolors = colors
    arrowcolors = [''] + arrowcolors
    arrow_length = 0.6
    if featurenames is None:
        featurenames = [f'feature{i + 1}' for i in range(n_neuro_layers[0])]

    # 特征层 取n个特征
    arrow_left = left + x_interval * 0.1
    n_neuro = n_neuro_layers[0]
    bottom = ycenter - n_neuro / 2
    for gene in featurenames[:n_neuro][::-1]:
        y = bottom + 0.5
        ax.text(left, y, gene, ha='right', va='center')
        ax.arrow(arrow_left, y, arrow_length * x_interval, 0, linewidth=arrowwidth, head_width=0.3, head_length=0.3,
                 fc='black', ec='black', length_includes_head=True)
        bottom += 1
    # 神经网络
    previous_y_list = []
    for n_neuro, color, arrowcolor in zip(n_neuro_layers, colors, arrowcolors):
        left += x_interval
        bottom = ycenter - n_neuro / 2
        y_list = []
        for i in range(n_neuro):  # 上半部分神经元
            plt.scatter(left, bottom + 0.5, s=neurosize, c=color, edgecolors='black', zorder=2)
            y_list.append(bottom + 0.5)
            bottom += 1
        for y in y_list:  # 连接上一层的箭头
            for previous_y in previous_y_list:
                plt.arrow(left - x_interval, previous_y, x_interval, y - previous_y, linewidth=arrowwidth,
                          head_width=0.1,
                          head_length=0, fc=arrowcolor, ec=arrowcolor, length_includes_head=True)
        previous_y_list = y_list

    # 输出层
    left += x_interval
    y = ycenter
    color = colors[len(n_neuro_layers)]
    arrowcolor = arrowcolors[len(n_neuro_layers)]
    plt.scatter(left, y, s=neurosize, c=color, edgecolors='black', zorder=2)
    for previous_y in previous_y_list:
        plt.arrow(left - x_interval, previous_y, x_interval, y - previous_y, linewidth=arrowwidth, head_width=0.1,
                  head_length=0, fc=arrowcolor, ec=arrowcolor, length_includes_head=True)

    # 输出层标签
    left += x_interval
    arrow_left = left + x_interval * 0.1
    ax.text(left, y, yname, ha='left', va='center')
    plt.arrow(arrow_left - x_interval, ycenter, arrow_length * x_interval, 0, linewidth=arrowwidth, head_width=0.3,
              head_length=0.3, fc='black', ec='black', length_includes_head=True)

    ax.set_xlim(-4, (len(n_neuro_layers) + 3) * x_interval + 1)
    ax.set_ylim(0, ymax)
    ax.axis('off')
    return ax
