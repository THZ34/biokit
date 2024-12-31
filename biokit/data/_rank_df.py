# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import numpy as np
import scanpy as sc


def get_rank_df(adata, group, key='rank_genes_groups'):
    """
    从adata中提取差异基因排序结果
    :param adata: AnnData object
    :param group: group name
    :param key: key in adata.uns
    :return: rank_df
    """
    rank_df = sc.get.rank_genes_groups_df(adata, group=group, key=key)
    rank_df['-log10(padj)'] = -np.log10(rank_df['pvals_adj'])
    rank_df['-log10(p)'] = -np.log10(rank_df['pvals'])

    return rank_df
