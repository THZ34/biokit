# ROC曲线
def roc_plot(dataframe, status, value, color=None, plot=True, ax=None):
    """分析ROC并画图，返回cutoff,AUC和ROC曲线图

    :param ax:
    :param color:
    :param dataframe: pandas.DataFrame对象
    :param status:标记状态的列名，二分类，1代表事件发生
    :param value:用于预测状态的值的列名
    :param plot:False不画图；True画图并返回plt.figure对象；输入文本则保存1200dpi图片
    :return:cutoff,auc,fig
    """
    from sklearn.metrics import roc_curve, roc_auc_score
    import numpy as np
    import matplotlib.pyplot as plt

    # 检查画板
    if not ax:
        fig, ax = plt.subplots()
    fpr, tpr, thresholds = roc_curve(dataframe[status], dataframe[value])
    # 计算最佳cutoff
    cutoff = thresholds[np.where((tpr - fpr) == (tpr - fpr).max())]
    # 计算AUC
    auc = roc_auc_score(dataframe[status], dataframe[value])
    # 画图
    if plot:
        ax.plot(fpr, tpr, color=color, label=f'AUC={auc:.3f}')
        ax.plot([0, 1], [0, 1], color='grey', linestyle='--')
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(f'{value} ~ {status}')
        ax.legend()
    return {'cutoff': cutoff, 'auc': auc}
