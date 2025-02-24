# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import numpy as np
import scanpy as sc


def get_rank_df(adata, group, key='rank_genes_groups', logfc_cutoff=1, pval_cutoff=0.05, fix_pval=True):
    """
    从adata中提取差异基因排序结果
    :param adata: AnnData object
    :param group: group name
    :param key: key in adata.uns
    :return: rank_df
    """
    rank_df = sc.get.rank_genes_groups_df(adata, group=group, key=key)
    rank_df.index = rank_df['names']
    if fix_pval:
        rank_df.loc[rank_df['pvals'] == 0, 'pvals'] = rank_df.loc[rank_df['pvals'] != 0, 'pvals'].min()
        rank_df.loc[rank_df['pvals_adj'] == 0, 'pvals_adj'] = rank_df.loc[rank_df['pvals_adj'] != 0, 'pvals_adj'].min()

    rank_df['-log10(padj)'] = -np.log10(rank_df['pvals_adj'])
    rank_df['-log10(p)'] = -np.log10(rank_df['pvals'])
    rank_df['change'] = 'other'
    rank_df.loc[(rank_df['logfoldchanges'] > logfc_cutoff) & (rank_df['pvals_adj'] < pval_cutoff), 'change'] = 'up'
    rank_df.loc[(rank_df['logfoldchanges'] < -logfc_cutoff) & (rank_df['pvals_adj'] < pval_cutoff), 'change'] = 'down'

    return rank_df
