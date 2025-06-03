# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import math
import matplotlib.pyplot as plt

base_color = {'A': '#91d542', 'T': '#db3124', 'C': '#457b9d', 'G': '#fa8600', ' ': 'white'}


def reverse_complement(seq):
    """
    反向互补序列
    :param seq:
    :return:
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[i] for i in seq[::-1]])


def plot_blast_result(blast_result, query_df, subject_seq, subject_name, full_length=40, ax=None):
    if not ax:
        fig, axes = plt.subplots(blast_result.shape[0], 1, figsize=(12, 3 * blast_result.shape[0]))
        if blast_result.shape[0] == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
    else:
        axes = [ax]
        fig = ax.figure

    # 在两端加上full_length个空格, 并且修正sub_start和sub_end, 避免left_pad和right_pad序列超出范围
    subject_seq = ' ' * full_length + subject_seq + ' ' * full_length

    for i, ax in zip(blast_result.index, axes):
        query_name = blast_result.loc[i, 'Query ID']
        query_seq = query_df.loc[query_df['miRNA'] == query_name, 'seq'].values[0]
        sub_start, sub_end, query_start, query_end = \
            blast_result.loc[i, ['Subject Start', 'Subject End', 'Query Start', 'Query End']].values

        pos_strand = True
        # 判断比对的正负链
        if sub_start < sub_end:
            sub_start, sub_end = sub_end, sub_start
            pos_strand = False

        sub_end += full_length
        sub_start += full_length

        mapped_sub_seq = subject_seq[sub_end - 1:sub_start]
        mapped_query_seq = query_seq[query_start - 1:query_end]
        # 计算两端需要补齐的长度
        left_pad = math.ceil((full_length - len(mapped_sub_seq)) / 2)
        right_pad = math.floor((full_length - len(mapped_sub_seq)) / 2)
        plot_sub_seq = subject_seq[sub_end - 1 - left_pad:sub_start + right_pad]

        # 画参考序列
        y = 0.6
        for x, base in enumerate(plot_sub_seq):
            x_in_mapped_region = left_pad - 1 < x < left_pad + len(mapped_sub_seq)
            if x_in_mapped_region:
                color = base_color[base]
                fontweight = 'bold'
            else:
                color = 'k'
                fontweight = 'normal'
            ax.text(x, y, base, color=color, fontsize=18, ha='center', va='center', fontweight=fontweight)

        # 画查询序列
        y = 0
        for x, base in enumerate(query_seq):
            x_in_mapped_region = query_start - 1 <= x <= query_end - 1
            if x_in_mapped_region:
                color = base_color[base]
                fontweight = 'bold'
            else:
                color = 'k'
                fontweight = 'normal'
            if pos_strand:
                plot_x = 38 - right_pad - x + query_start
                ax.text(plot_x, y, base, color=color, fontsize=18, ha='center', va='center',
                        fontweight=fontweight)
                if reverse_complement(base) == plot_sub_seq[plot_x] and x_in_mapped_region:
                    ax.plot([plot_x, plot_x], [0.25, 0.45], color='k', lw=2)
            elif not pos_strand and x_in_mapped_region:
                plot_x = left_pad - query_start + x + 1
                ax.text(plot_x, y, base, color=color, fontsize=18, ha='center', va='center',
                        fontweight=fontweight)
                if base == plot_sub_seq[plot_x]:
                    ax.plot([plot_x, plot_x], [0.25, 0.45], color='k', lw=2)

        # 补充头尾以及序列名
        # 3'5'
        ax.text(-1, 0.6, "5'-", fontsize=18, ha='right', va='center')
        ax.text(40, 0.6, "-3'", fontsize=18, ha='left', va='center')

        if pos_strand:
            ax.text(-1, 0, "3-'", fontsize=18, ha='right', va='center')
            ax.text(40, 0, "-5'", fontsize=18, ha='left', va='center')
        else:
            ax.text(-1, 0, "5-'", fontsize=18, ha='right', va='center')
            ax.text(40, 0, "-3'", fontsize=18, ha='left', va='center')

        # 参考基因名及范围
        ax.text(20, 1, subject_name, fontsize=20, ha='center', va='center',
                bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
        ax.text(-0.5, 1, sub_end - full_length, fontsize=18, ha='left', va='center',
                bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
        ax.text(39.5, 1, sub_start - full_length, fontsize=18, ha='right', va='center',
                bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
        ax.plot([0, 39], [1.05, 1.05], color='k', lw=2)

        # 查询序列名
        ax.text(20, -0.4, query_name, fontsize=20, ha='center', va='center')

        ax.set_xlim(-3, 43)
        ax.set_ylim(-1, 1.5)
        ax.set_xticks([])
        ax.set_yticks([])
        # plt.show()
    plt.subplots_adjust(0, 0, 1, 1, 0, 0)
    return fig, ax
