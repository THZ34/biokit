# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %% 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch

config = {"font.family": 'Microsoft YaHei', 'pdf.fonttype': 42}
rcParams.update(config)


def flow_chart(df, datasets, fold_result, rank_df):
    """

    :param df:
    :param datasets:
    :param fold_result:
    :param rank_df:
    :return:
    """
    train_df = df[df['dataset'].isin(datasets)]

    best_n_feature = {}
    for model_name in fold_result:
        train_scores = []
        test_scores = []
        scores_list = []
        x = []
        for n_feature in sorted([int(i) for i in list(fold_result[model_name].keys())]):
            train_score = np.mean(fold_result[model_name][str(n_feature)]['train'])
            test_score = np.mean(fold_result[model_name][str(n_feature)]['test'])
            train_scores.append(train_score)
            test_scores.append(test_score)
            scores_list.append(train_score + test_score)
            x.append(n_feature)

        best_n_feature[model_name] = (scores_list.index(max(scores_list)), max(scores_list))
    best_n_feature = pd.DataFrame(best_n_feature, index=['n_features', 'max score']).T
    best_n_feature.sort_values(by='max score', ascending=False, inplace=True)
    snsax = sns.clustermap(rank_df)
    plt.close()

    plt.figure(figsize=(24, 7))
    gridspec = plt.GridSpec(14, 24)

    #  第一步：整理数据集
    ax = plt.subplot(gridspec[2:14, 0:6])
    data_integrate(datasets, train_df, ax)
    plt.subplots_adjust(0, 0, 1, 1)

    #  衔接箭头
    ax = plt.subplot(gridspec[2:14, 6:8])
    ax.axis('off')
    ax.set(xlim=(0, 2), ylim=(0, 6))
    ax.arrow(x=0, y=3, dx=2, dy=0, width=0.1, fc='black', length_includes_head=True, shape='full', zorder=1)

    #  第二步，RFE获得基因rank,画热图
    ax = plt.subplot(gridspec[2:14, 10:16])
    sns.heatmap(rank_df.iloc[snsax.dendrogram_row.reordered_ind, snsax.dendrogram_col.reordered_ind].T,
                cmap='viridis_r',
                cbar=False, ax=ax)
    ax.set(ylabel=None, xlabel=None, xticks=[])
    ax.set_yticklabels(labels=ax.get_yticklabels(), fontsize=10, weight='bold')
    ax.set_ylim(ax.get_ylim()[0] * 1.125, ax.get_ylim()[0] * -0.125)

    #  衔接箭头
    ax = plt.subplot(gridspec[2:14, 16:18])
    ax.axis('off')
    ax.set(xlim=(0, 2), ylim=(0, 6))
    ax.arrow(x=0, y=3, dx=2, dy=0, width=0.1, fc='black', length_includes_head=True, shape='full', zorder=1)

    # 特征数量筛选
    axes = [plt.subplot(gridspec[(i * 4 + 4):(i * 4 + 7), 18:24]) for i in range(2)]
    for ax, model_name in zip(axes[:2], best_n_feature.index[:2]):
        x = range(1, len(fold_result[model_name]) + 1)
        train_y = [np.mean(fold_result[model_name][str(i)]['train']) for i in x]
        test_y = [np.mean(fold_result[model_name][str(i)]['test']) for i in x]
        sum_y = [i + j for i, j in zip(train_y, test_y)]
        ax.plot(x, train_y, label='train')
        ax.plot(x, test_y, label='test')
        ax.plot(x, sum_y, label='train+test')
        ax.scatter(x=best_n_feature.loc[model_name]['n_features'], y=max(sum_y), color='red', s=20)
        ax.text(x=best_n_feature.loc[model_name]['n_features'], y=max(sum_y) - 0.05,
                s=f"{best_n_feature.loc[model_name]['n_features']:0.0f}", va='top', ha='center')
        ax.set(ylim=(0, 2.2), yticks=[0, 1, 2], ylabel='f1 score', title=model_name)

    ax = plt.subplot(gridspec[11:13, 18:24])
    ax.axis('off')
    ax.set(xlim=(0, 4), ylim=(0, 5.5))
    ax.text(x=0, y=1, s='x: Number of features\nRed point: Number of features with max score', va='top')
    ax.scatter(x=[2, 2, 2], y=[4, 3, 2], c='black')
    #  步骤标题
    for x1, x2, title in [(0, 6, 'Data Integrate'), (8, 16, 'Gene Importance Rank'), (18, 24, 'Feature Selection')]:
        ax = plt.subplot(gridspec[0:2, x1:x2])
        ax.axis('off')
        ax.set(xlim=(0, 2), ylim=(0, 2))
        ax.text(x=1, y=1, s=title, fontsize=30, ha='center', va='center')

    pdf = PdfPages('The flow chart.pdf')
    pdf.savefig()
    plt.close()
    pdf.close()


def data_integrate(datasets, train_df, ax):
    """数据集整合示意图

    :param datasets:
    :param train_df:
    :param ax:
    :return:
    """
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-0.5, 4.5)
    ax.axis('off')
    inte_x, inte_y, inte_width, inte_height = 0, 0, 4, 4
    # 边框
    patch = FancyBboxPatch(xy=(inte_x, inte_y), width=inte_width, height=inte_height, facecolor='white', linestyle='--',
                           linewidth=1, boxstyle='round,pad=0,rounding_size=0.05', zorder=3)
    ax.add_patch(patch)
    # ax.text(inte_x + inte_width / 2, inte_y + inte_height + 0.25, s='Datasets Integrate', va='bottom', ha='center',
    #         fontsize=15)
    # 数据集
    facecolors = sns.hls_palette(4)
    bottom = 0
    for facecolor, dataset in zip(facecolors, datasets):
        sub_x = inte_x + 2
        sub_y = inte_y + 0.25 + bottom
        sub_width = 1.9
        sub_height = (inte_height - 0.5) * train_df[train_df['dataset'] == dataset].shape[0] / train_df.shape[0]
        patch = FancyBboxPatch(xy=(sub_x, sub_y), width=sub_width, height=sub_height, facecolor=facecolor,
                               linewidth=1, boxstyle='round,pad=-0.05,rounding_size=0.05', zorder=4)
        ax.add_patch(patch)
        ax.text(inte_x + 0.1, sub_y + sub_height / 2, dataset, ha='left', va='center', zorder=4, fontsize=15)
        bottom += sub_height
    # 标注整合数据集维度
    ax.text(x=inte_x + inte_width + 0.15, y=inte_y + inte_height / 2, s=f'{train_df.shape[0]} samples', va='center',
            ha='center', rotation=-90, bbox={'facecolor': 'white', 'linewidth': 0}, zorder=2, fontsize=15)
    ax.plot([inte_x + inte_width, inte_x + inte_width + 0.3], [inte_y + inte_height, inte_y + inte_height],
            color='black')
    ax.plot([inte_x + inte_width, inte_x + inte_width + 0.3], [inte_y, inte_y], color='black')

    ax.text(x=inte_x + inte_width / 2, y=inte_y - 0.15, s=f'{train_df.shape[1]} genes', va='center', ha='center',
            bbox={'facecolor': 'white', 'linewidth': 0}, zorder=2, fontsize=15)
    ax.plot([inte_x, inte_x], [inte_y - 0.3, inte_y], color='black')
    ax.plot([inte_x + inte_width, inte_x + inte_width], [inte_y - 0.3, inte_y], color='black')
    # 维度标注箭头
    ax.arrow(x=inte_x + 0.05, y=inte_y - 0.15, dx=inte_width - 0.1, dy=0, width=0.03, shape='full',
             length_includes_head=True, zorder=1)
    ax.arrow(x=inte_x + inte_width + 0.15, y=inte_y + 0.05, dx=0, dy=inte_height - 0.1, width=0.03,
             length_includes_head=True, shape='full', zorder=1)
    return ax
