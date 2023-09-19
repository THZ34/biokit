# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import numpy as np
import pandas as pd


def calculate_average(data, method='arithmetic'):
    if method == 'arithmetic':
        return np.mean(data)
    elif method == 'geometric':
        return np.prod(data) ** (1 / len(data))
    elif method == 'harmonic':
        return len(data) / np.sum(1 / data)
    else:
        raise ValueError("Unknown method")


def signature_score(exp_df, method='geometric', signature_dict=None):
    """

    :param method: 评分方法
    :param exp_df: pd.DataFrame, index为基因名，columns为样本名
    :param signature_dict:
    :return:
    """
    if not signature_dict:
        from biokit.data import load_rna_signature
        signature_dict = load_rna_signature()

    signature_df = pd.DataFrame(index=exp_df.index, columns=signature_dict.keys())
    for sig_name in signature_dict:
        genes = signature_dict[sig_name]
        genes = [x for x in genes if x in exp_df.columns]
        if len(genes) > 0:
            signature_df[sig_name] = calculate_average(exp_df[genes].T, method=method)
        else:
            signature_df[sig_name] = np.nan
    return pd.DataFrame(signature_df)
