import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator
from seaborn import color_palette
from sklearn.linear_model import LassoCV


def lassocv(X=None, y=None, ax=None, alphas=None, cv=10, random_state=0):
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
    colors = color_palette("Oranges", as_cmap=True)(lasso_result['MSE'])
    bplot1 = ax.boxplot(mse_path)
    ax.scatter(x=lasso_result.index.to_list().index(model.alpha_) + 1, y=lasso_result.loc[model.alpha_]['MSE'], s=5,
               marker='D')
    # for patch, color in zip(bplot1['boxes'], colors):
    #     patch.set_facecolor(color)

    ax.set_xlabel('log10 alpha')
    ax.set_ylabel('Mean Square Error')

    ax2 = ax.twiny()
    xlim = ax.get_xlim()
    ax2.set_xlim(xlim[0], xlim[1])
    ax2.set_xticks(ax.get_xticks())
    ax2.minorticks_off()

    ax.set_xticklabels(np.round(np.log10(alphas), 2).tolist())
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax2.set_xticklabels(lasso_result['n_feature'])
    ax2.xaxis.set_major_locator(MultipleLocator(2))
    return model, lasso_result, mse_path, ax
