# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import pandas as pd
import os
import numpy as np
import math


# %%
def genecnr(infile):
    gcnv = {}
    with open(infile, "r") as fd:
        for line in fd.readlines():
            arr = line.strip().split("\t")
            if arr[0] == "chromosome":
                continue
            num = 2 ** float(arr[4])
            region = (int(arr[2]) - int(arr[1]) + 1) * float(arr[7])
            if arr[3] == "Background" or arr[3] == "-" or float(arr[7]) < 0.25:
                continue
            if not arr[3] in gcnv:
                gcnv[arr[3]] = {}
                gcnv[arr[3]]["c"] = arr[0]
                gcnv[arr[3]]["s"] = int(arr[1])
                gcnv[arr[3]]["l"] = 0
                gcnv[arr[3]]["n"] = 0
                gcnv[arr[3]]["b"] = 0
            gcnv[arr[3]]["e"] = int(arr[2])
            gcnv[arr[3]]["l"] += region
            gcnv[arr[3]]["n"] += num * region
            gcnv[arr[3]]["b"] += 1
    gcnv = pd.DataFrame(gcnv)
    return gcnv


def check_aachange_file(record_df):
    # 检查vardict aachange file是否存在
    aachange_files = []
    for workdir, pairid in record_df[['workdir', 'pairid']].to_numpy():
        cloudid = workdir.split('/')[-2]
        aachange_file_path = None
        for file_path in [f'{workdir}/somatic/vardict/{pairid}/{pairid}.snv.sindel.AAchange.VAF.xls'
                          f'{workdir}/call-vardict_anno_process/{cloudid}_{pairid}.snv.sindel.AAchange.VAF.xls',
                          f'{workdir}/call-vardict_anno_process/execution/{cloudid}_{pairid}.snv.sindel.AAchange.VAF.xls',
                          f'{workdir}/call-vardict_anno_process/cacheCopy/execution/{cloudid}_{pairid}.snv.sindel.AAchange.VAF.xls']:
            if os.path.exists(file_path):
                aachange_file_path = file_path
                break
        aachange_files.append(aachange_file_path)
    record_df['aachange_file'] = aachange_files
    # 没有找到aachange文件的样本，解析缓存命中
    for tumorid in record_df[record_df['aachange_file'].isna()]['tumorid']:
        workdir, pairid = record_df.loc[tumorid, ['workdir', 'pairid']]
        cloudid = workdir.split('/')[-2]
        project_path = '/'.join(workdir.split('/')[:-1])
        caching_file = f'{workdir}/call-vardict_anno_process/call_caching_placeholder.txt'
        previous_cormwellid = open(caching_file).read().split()[-1].split('/')[-2]
        file_path1 = f'{project_path}/{previous_cormwellid}/call-vardict_anno_process/{cloudid}_{pairid}.snv.sindel.AAchange.VAF.xls'
        file_path2 = f'{project_path}/{previous_cormwellid}/call-vardict_anno_process/execution/{cloudid}_{pairid}.snv.sindel.AAchange.VAF.xls'
        if os.path.exists(file_path1):
            record_df.loc[tumorid, 'aachange_file'] = file_path1
        elif os.path.exists(file_path2):
            record_df.loc[tumorid, 'aachange_file'] = file_path2
        else:
            record_df.loc[tumorid, 'aachange_file'] = None


def check_cnr_file(record_df):
    # CNV文件
    cnv_files = []
    for workdir, pairid in record_df[['workdir', 'pairid']].to_numpy():
        cloudid = workdir.split('/')[-2]
        cnv_file = None
        for file in [f'{workdir}/call-cnvkit/{cloudid}_{pairid}.call.cnr',
                     f'{workdir}/call-cnvkit/execution/{cloudid}_{pairid}.call.cnr',
                     f'{workdir}/call-cnvkit/cacheCopy/execution/{cloudid}_{pairid}.call.cnr']:
            if os.path.exists(file):
                cnv_file = file
                break
        cnv_files.append(cnv_file)
    record_df['cnv_file'] = cnv_files

    # 没有找到cnr文件的样本，解析缓存命中
    for tumorid in record_df[record_df['cnv_file'].isna()]['tumorid']:
        workdir, pairid = record_df.loc[tumorid, ['workdir', 'pairid']]
        cloudid = workdir.split('/')[-2]
        project_path = '/'.join(workdir.split('/')[:-1])
        caching_file = f'{workdir}/call-cnvkit/call_caching_placeholder.txt'
        previous_cormwellid = open(caching_file).read().split()[-1].split('/')[-2]
        file_path1 = f'{project_path}/{previous_cormwellid}/call-cnvkit/{cloudid}_{pairid}.call.cnr'
        file_path2 = f'{project_path}/{previous_cormwellid}/call-cnvkit/execution/{cloudid}_{pairid}.call.cnr'
        if os.path.exists(file_path1):
            record_df.loc[tumorid, 'cnv_file'] = file_path1
        elif os.path.exists(file_path2):
            record_df.loc[tumorid, 'cnv_file'] = file_path2
        else:
            record_df.loc[tumorid, 'cnv_file'] = None


def read_cnv(record_df):
    cnv_df = pd.DataFrame()
    for file, tumorid in record_df[['cnv_file', 'tumorid']].to_numpy():
        temp_df = genecnr(file).T
        temp_df['sample'] = tumorid
        cnv_df = pd.concat([cnv_df, temp_df], axis=0)

    cnv_df['gene'] = cnv_df.index
    cnv_df = cnv_df[cnv_df['gene'] != 'Background']
    cnv_df['ratio'] = cnv_df['n'] / cnv_df['l']
    cnv_df['ratio'] = cnv_df['ratio'].astype(float)
    cnv_df['logr'] = np.log(cnv_df['ratio']) / math.log(2)
    cnv_df['cn'] = 2 ** cnv_df['logr'] * 2

    return cnv_df


def confirm_gainloss(cnv_df, cnvlimit=None):
    if not cnvlimit:
        cnvlimit = (4, 1)

    cnv_up, cnv_low = cnvlimit
    cn_df = cnv_df.copy()
    cn_df['effect'] = 'no mutate'
    cn_df.loc[cn_df['cn'] > cnv_up, 'effect'] = 'CNV_amp'
    cn_df.loc[cn_df['cn'] < cnv_low, 'effect'] = 'CNV_del'


def float_qc_df(qc_df):
    percent_cols = ['Target_uniformity', 'Mapping_rate', 'Duplication_rate', 'Capture_rate', 'Target_coverage',
                    'Target_1000X', 'Target_500X', 'Target_300X', 'Target_100X', 'Target_10X', 'Q20', 'Q30', 'GC',
                    'Clean_base_ratio']
    bp_cols = ['Insert_size', 'Clean_Base', 'Average_read_length']
    depth_cols = ['Depth_in_target']

    for col in percent_cols:
        qc_df[col] = qc_df[col].str.rstrip('%').astype(float) / 100
    for col in bp_cols:
        qc_df[col] = qc_df[col].str.rstrip(' bp').astype(float)
    for col in depth_cols:
        qc_df[col] = qc_df[col].str.rstrip('X').astype(float)
    qc_df['Coding_size'] = qc_df['Coding_size'].apply(lambda x: x.replace(',', ''))
    qc_df = qc_df.astype(float)
    return qc_df
