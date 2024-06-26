# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import os
import json
import requests

datapath = os.path.join(os.path.dirname(__file__), 'kegg')


def download_kegg_pathway(hsaid, datapath=datapath):
    """Download KEGG pathway data by hsa id"""
    if not os.path.exists(os.path.join(datapath, f'{hsaid}.json')):
        os.makedirs(datapath, exist_ok=True)
        url = f'http://togows.dbcls.jp/entry/pathway/{hsaid}/genes.json'
        r = requests.get(url)
        with open(os.path.join(datapath, f'{hsaid}.json'), 'w') as f:
            f.write(r.text)
        pathway_gene_dict = json.loads(r.text)
        pathway_gene_dict = {k: v.split(';')[0] for k, v in pathway_gene_dict[0].items() if v}
    else:
        with open(os.path.join(datapath, f'{hsaid}.json')) as f:
            pathway_gene_dict = json.load(f)
            pathway_gene_dict = {k: v.split(';')[0] for k, v in pathway_gene_dict[0].items() if v}
    return pathway_gene_dict
