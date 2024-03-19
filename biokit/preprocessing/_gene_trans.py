# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import gzip
import os
import re
import pandas as pd


# %%
def download_gencode(versions=None, outdir='ref/gencode'):
    """下载gencode GRCh38 """
    os.makedirs(outdir, exist_ok=True)
    if not versions:
        versions = range(34, 45)
    with open('download_gencode.sh', 'w') as f:
        for release_version in versions:
            url = f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release_version}/gencode.v{release_version}.annotation.gtf.gz'
            filename = f'gencode.v{release_version}.annotation.gtf.gz'
            unzip_filename = f'gencode.v{release_version}.annotation.gtf'
            if not os.path.exists(f'{outdir}/{filename}') and not os.path.exists(f'{outdir}/{unzip_filename}'):
                f.write(f'wget {url} -O {outdir}/{filename}\n')
            else:
                print(f'{filename}已存在')


def download_ensembl(versions=None, outdir='ref/ensembl'):
    """下载ensembl GRCh38 """
    os.makedirs(outdir, exist_ok=True)
    if not versions:
        versions = range(76, 111)
    with open('download_ensembl.sh', 'w') as f:
        for release_version in versions:
            url = f'https://ftp.ensembl.org/pub/release-{release_version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release_version}.gtf.gz'
            filename = f'Homo_sapiens.GRCh38.{release_version}.gtf.gz'
            unzip_filename = f'Homo_sapiens.GRCh38.{release_version}.gtf'
            if not os.path.exists(f'{outdir}/{filename}') and not os.path.exists(f'{outdir}/{unzip_filename}'):
                f.write(f'wget {url} -O {outdir}/{filename}\n')
            else:
                print(f'{filename}已存在')


def get_gencode_genename_df(versions, outdir):
    genename_df = {}
    genetype_df = {}
    geneid_df = {}
    for release_version in versions:
        version_genename = {}
        version_genetype = {}
        version_geneid = {}
        # 提取genename
        gencode_gz_file = f'ref/gencode/gencode.v{release_version}.annotation.gtf.gz'
        with gzip.open(gencode_gz_file, 'rt') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                try:
                    genename = re.findall('gene_name "(.*?)";', line)[0]
                except IndexError:
                    continue
                genetype = re.findall('gene_type "(.*?)";', line)[0]
                geneid = re.findall('gene_id "(.*?)";', line)[0]
                version_genename[genename] = 1
                version_genetype[genename] = genetype
                version_geneid[genename] = geneid
        genename_df[release_version] = pd.DataFrame.from_dict(version_genename, orient='index',
                                                              columns=[f'gencode_v{release_version}'])
        genetype_df[release_version] = pd.DataFrame.from_dict(version_genetype, orient='index',
                                                              columns=[f'gencode_v{release_version}'])
        geneid_df[release_version] = pd.DataFrame.from_dict(version_geneid, orient='index',
                                                            columns=[f'gencode_v{release_version}'])
    genename_df = pd.concat(genename_df.values(), axis=1)
    genetype_df = pd.concat(genetype_df.values(), axis=1)
    geneid_df = pd.concat(geneid_df.values(), axis=1)
    genename_df.fillna(0, inplace=True)
    genename_df.to_csv('ref/gencode/genename.csv')
    genetype_df.to_csv('ref/gencode/genetype.csv')
    geneid_df.to_csv('ref/gencode/geneid.csv')


def get_ensembl_genename_df(versions):
    genename_df = {}
    genetype_df = {}
    geneid_df = {}
    for release_version in versions:
        version_genename = {}
        version_genetype = {}
        version_geneid = {}
        # 提取genename
        ensembl_gz_file = f'ref/ensembl/Homo_sapiens.GRCh38.{release_version}.gtf.gz'
        with gzip.open(ensembl_gz_file, 'rt') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                try:
                    genename = re.findall('gene_name "(.*?)";', line)[0]
                except IndexError:
                    continue
                genetype = re.findall('gene_biotype "(.*?)";', line)[0]
                geneid = re.findall('gene_id "(.*?)";', line)[0]
                version_genename[genename] = 1
                version_genetype[genename] = genetype
                version_geneid[genename] = geneid
        genename_df[release_version] = pd.DataFrame.from_dict(version_genename, orient='index',
                                                              columns=[f'ensembl_v{release_version}'])
        genetype_df[release_version] = pd.DataFrame.from_dict(version_genetype, orient='index',
                                                              columns=[f'ensembl_v{release_version}'])
        geneid_df[release_version] = pd.DataFrame.from_dict(version_geneid, orient='index',
                                                            columns=[f'ensembl_v{release_version}'])
    genename_df = pd.concat(genename_df.values(), axis=1)
    genetype_df = pd.concat(genetype_df.values(), axis=1)
    geneid_df = pd.concat(geneid_df.values(), axis=1)
    genename_df.fillna(0, inplace=True)
    genename_df.to_csv('ref/ensembl/genename.csv')
    genetype_df.to_csv('ref/ensembl/genetype.csv')
    geneid_df.to_csv('ref/ensembl/geneid.csv')


def detect_version(genes, genename_df):
    genesum = genename_df.loc[genes].sum()
    version = genesum[genesum == genesum.max()].index[0]
    return version


def genename_version_convert(genes, geneid_df, from_version, to_version):
    genes = [i for i in genes if i in geneid_df.index]
    from_geneid_df = geneid_df.loc[genes, from_version].dropna()
    genes = from_geneid_df.index
    from_version_dict = dict(zip(from_geneid_df.index, from_geneid_df))
    to_geneid_df = geneid_df.loc[geneid_df[to_version].isin(from_geneid_df.to_list()), to_version]
    to_version_dict = dict(zip(to_geneid_df, to_geneid_df.index))
    gene_trans_dict = {gene: to_version_dict.get(from_version_dict[gene], gene) for gene in genes}
    gene_trans_dict = {key: value for key, value in gene_trans_dict.items() if (key != value)}
    return gene_trans_dict
