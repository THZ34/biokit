# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import matplotlib.pyplot as plt
import math


def pca_variance_ratio(pca, ax=None):
    """ªÊ÷∆PCA∑Ω≤Óπ±œ◊¬ Õº"""
    if ax is None:
        fig, ax = plt.subplots()
    fig, ax = plt.subplots()
    for i, variance in enumerate(pca.explained_variance_ratio_):
        ax.text(i, variance + 0.02, f'{variance:.2f}', ha='center', va='bottom')
    ax.bar(range(pca.n_components), pca.explained_variance_ratio_, color='skyblue', width=0.6)
    ax.set_xlim(-0.5, pca.n_components - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Principal component')
    ax.set_ylabel('Variance ratio')
    ax.set_xticks(range(pca.n_components))
    ax.set_xticklabels([f'PC{i + 1}' for i in range(pca.n_components)])
    ax.set_title('PCA Variance Ratio')
    return ax


def pca_loadings(pca, featurenames=None, n_pcs=2, n_features=None, axes=None):
    """ªÊ÷∆PCA‘ÿ∫…Õº"""
    if not axes:
        n_row = math.floor(math.sqrt(n_pcs * 2))
        n_col = math.ceil(n_pcs / n_row)
        fig, axes = plt.subplots(n_row, n_col, figsize=(6 * n_col, 4 * n_row))
    if featurenames is None:
        featurenames = [f'feature{i + 1}' for i in range(n_features)]
    if n_features is None:
        n_features = min(10, pca.n_components)

    ymin = pd.DataFrame(pca.components_[:n_features, :]).min().min()
    ymax = pd.DataFrame(pca.components_[:n_features, :]).max().max()
    axes = axes.flatten()
    for i in range(n_pcs):
        ax = axes[i]
        loading = pca.components_[i]
        loading_df = pd.DataFrame({'feature': featurenames, 'loading': loading})
        loading_df.sort_values('loading', ascending=False, inplace=True)

        for j, (feature, loading) in enumerate(
                zip(loading_df['feature'][:n_features], loading_df['loading'][:n_features])):
            ax.text(j, loading, feature, ha='center', va='center', rotation=90)

        ymean = (loading_df.iloc[n_features]['loading'] + loading_df.iloc[-n_features]['loading']) / 2
        for j in range(n_features, n_features + 3):
            ax.scatter(j, ymean, color='black', s=30)
        for j, (feature, loading) in enumerate(
                zip(loading_df['feature'][-n_features:], loading_df['loading'][-n_features:])):
            ax.text(j + 3 + n_features, loading, feature, ha='center', va='center', rotation=90)

        ax.set_xlim(-0.5, n_features * 2 + 2.5)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks([])
        ax.set_title(f'PC{i + 1} Loading')
    return axes
