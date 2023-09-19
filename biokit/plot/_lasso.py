import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator
from sklearn.linear_model import LassoCV


def lassocv(X=None, y=None, ax=None, alphas=None, cv=10, random_state=0, linewidth=1):
    """

    :param X:
    :param y:
    :param ax:
    :param alphas:
    :param cv:
    :param random_state:
    :return:
    """
    if not ax:
        fig, ax = plt.subplots(figsize=(12, 4))
    if alphas is None:
        alphas = 10 ** np.linspace(-3, 1, 100)

    lasso_result = []
    mse_path = []
    for alpha in alphas:
        model = LassoCV(alphas=[alpha], cv=cv, random_state=random_state)
        model.fit(X, y)
        lasso_result.append([alpha, (model.coef_ != 0).sum(), model.mse_path_.mean()])
        mse_path.append(model.mse_path_)
    lasso_result = pd.DataFrame(lasso_result, columns=['alpha', 'n_feature', 'MSE'])
    lasso_result.index = lasso_result['alpha']
    mse_path = pd.DataFrame(mse_path).T

    model = LassoCV(alphas=alphas, cv=cv, random_state=random_state)
    model.fit(X, y)
    width = 0.6
    for mse_row, alpha, i in zip(model.mse_path_, model.alphas, range(len(model.alphas))):
        mse_max = np.max(mse_row)
        mse_min = np.min(mse_row)
        ax.plot([i, i], [mse_min, mse_max], color='k', linewidth=linewidth)
        ax.plot([i - width, i + width], [mse_max, mse_max], color='k', linewidth=linewidth)
        ax.plot([i - width, i + width], [mse_min, mse_min], color='k', linewidth=linewidth)
    ax.scatter(x=lasso_result.index.to_list().index(model.alpha_) + 1, y=lasso_result.loc[model.alpha_]['MSE'], s=5,
               marker='D')
    ax.set_xlabel('log10 alpha')
    ax.set_ylabel('Mean Square Error')
    ax.set_xticks(range(len(model.alphas)))
    ax.set_xticklabels(model.alphas)

    ax2 = ax.twiny()
    xlim = ax.get_xlim()
    ax2.set_xlim(xlim[0], xlim[1])
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(lasso_result.iloc[ax.get_xticks(), 'n_feature'])
    ax2.minorticks_off()

    print(lasso_result)
    ax.set_xticklabels(np.round(np.log10(alphas), 2).tolist())
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax2.xaxis.set_major_locator(MultipleLocator(2))
    return model, lasso_result, mse_path, ax
