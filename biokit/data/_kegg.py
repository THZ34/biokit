# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import os
import json
import requests

def download_kegg_pathway(hsaid,datapath):
    """Download KEGG pathway data by hsa id"""
    if not os.path.exists(os.path.join(datapath, f'{hsaid}.json')):
        os.makedirs(datapath)
        url = f'http://togows.dbcls.jp/entry/pathway/{hsaid}/genes.json'
        r = requests.get(url)
        with open(os.path.join(datapath, f'{hsaid}.json'), 'w') as f:
            f.write(r.text)
        pathway_gene_dict = json.loads(r.text)
        pathway_gene_dict = {k: v.split(';')[0] for k, v in pathway_gene_dict.items() if v}
    else:
        with open(os.path.join(datapath, f'{hsaid}.json')) as f:
            pathway_gene_dict = json.load(f)
            pathway_gene_dict = {k: v.split(';')[0] for k, v in pathway_gene_dict.items() if v}
    return pathway_gene_dict
