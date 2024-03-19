# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import os
def load_go_pathway_geneset():
    """读取GO pathway geneset"""
    go_file = os.path.join(data_path, 'pathway_geneset/GO_geneset.txt')
    return pd.read_csv(go_file, sep='\t')