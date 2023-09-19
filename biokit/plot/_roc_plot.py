import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve


# ROC曲线

def rocplots(states=None, values=None, colors=None, plot=True, ax=None, linestyles=None, labels=None, annot_aucs=False,
             annot_loc=(1, 0)):
    """分析ROC并画图，返回cutoff,AUC和ROC曲线图

    :param annot_loc:
    :param annot_aucs:
    :param labels:
    :param linestyles:
    :param ax:
    :param colors:
    :param states:标记状态的列名，二分类，1代表事件发生
    :param values:用于预测状态的值的列名
    :param plot:False不画图；True画图并返回plt.figure对象；输入文本则保存1200dpi图片
    :return:cutoff,auc,fig
    """

    # 检查画板
    if not ax:
        fig, ax = plt.subplots(figsize=(7.142857, 5.8823))
        plt.subplots_adjust(0.05, 0.05, 0.75, 0.9)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Rositive Rate')
        ax.set_title('ROC(Receiver Operating Characteristic)')

    if not colors:
        colors = sns.hls_palette(len(states))
    if not linestyles:
        linestyles = ['--'] * len(states)

    cutoffs, aucs = [], []
    for state, value, color, linestyle, label in zip(states, values, colors, linestyles, labels):
        cutoff, auc = rocplot(value=value, status=state, ax=ax, label=label, color=color, linestyle=linestyle)
        cutoffs.append(cutoff)
        aucs.append(auc)
    ax.legend(bbox_to_anchor=(1, 1))
    if annot_aucs:
        max_string_length = max([len(label) for label in labels])
        auc_table = [f'{"AUC":<{max_string_length}}'] + [
            f'{label}{" " * ((max_string_length - len(label)) * 2)} : {auc:0.2f}' for label, auc in zip(labels, aucs)]
        auc_strings = '\n'.join(auc_table)
        ax.text(x=annot_loc[0], y=annot_loc[1], s=auc_strings, backgroundcolor='lightcyan', ha='right', va='bottom')
    return cutoffs, aucs


def rocplot(value, status, ax=None, label='', color='deepskyblue', linestyle=(0, (2, 2, 1, 1)), fill=False,
            linewidth=10):
    fpr, tpr, thresholds = roc_curve(status, value)
    ax.plot(fpr, tpr, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    if fill:
        fill_x = fpr
        fill_x = np.append(fill_x, [[1]])
        fill_y = tpr
        fill_y = np.append(fill_y, [[0]])
        ax.fill(fill_x, fill_y, color='deepskyblue')
    cutoff = thresholds[np.where((tpr - fpr) == (tpr - fpr).max())]
    auc = roc_auc_score(status, value)
    return cutoff, auc
