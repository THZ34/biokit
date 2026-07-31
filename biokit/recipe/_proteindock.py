# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import requests
import os
import pandas as pd
import numpy as np


def download_alphafold_pdb(gene_uniprotid_dict, outdir=''):
    os.makedirs(outdir, exist_ok=True)

    for gene, uni in gene_uniprotid_dict.items():
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uni}-F1-model_v4.pdb"
        dst = os.path.join(outdir, f"{gene}.pdb")
        resp = requests.get(url, stream=True)
        if resp.status_code == 200:
            with open(dst, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Downloaded: {dst}")
        else:
            print(f"ERROR: Failed to download for {gene} (UniProt:{uni}), status code {resp.status_code}")
