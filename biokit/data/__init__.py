# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %% 
import pandas as pd
import os

data_path = os.path.dirname(os.path.abspath(__file__))


# %%
def load_hg19_karyo():
    """读取数据库中的hg19 karyoband"""
    karyoband_file = os.path.join(data_path, 'ref/hg19/hg19_karyo_band.txt')
    return pd.read_csv(karyoband_file, sep='\t')


def load_hg19_ref():
    """读取数据库中的hg19 karyoband"""
    karyoband_file = os.path.join(data_path, 'ref/hg19/HG19 refGene.txt')
    return pd.read_csv(karyoband_file, sep='\t')


def load_hg38_karyo():
    """读取数据库中的hg38 karyoband"""
    karyoband_file = os.path.join(data_path, 'ref/hg38/hg38_karyo_band.txt')
    return pd.read_csv(karyoband_file, sep='\t')


def load_hg38_ref():
    """读取数据库中的hg38 karyoband"""
    karyoband_file = os.path.join(data_path, 'ref/hg38/HG38 refGene.txt')
    return pd.read_csv(karyoband_file, sep='\t')


def load_xcell_signature():
    import json
    return json.load(open(os.path.join(data_path, 'xCell/xcell_sig.json')))
