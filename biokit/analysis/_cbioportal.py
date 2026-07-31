# -*- coding: utf-8 -*-
"""
跨 cBioPortal 数据集的「基因–生存关联」粗筛 (rough screen)

对给定基因列表，在指定目录下的多个 cBioPortal 数据集中逐一做 KM best-cutoff
生存分析，汇总每个 (数据集, 基因, 生存终点) 的最佳 cutoff 与最优 logrank p 值，
用于挑选样本量更大、可用于验证生存结论的数据集。

与生存分析在 *sample* 层面对齐（与你更正后的逻辑一致）：先把病人级生存信息
广播到样本表，再按 SAMPLE_ID 与表达矩阵 join。
"""

import os
import tarfile

import numpy as np
import pandas as pd

# 默认优先选用的 RNA 表达文件名（按优先级从高到低）
DEFAULT_RNA_FILENAMES = ('data_mrna_seq_v2_rsem.txt', 'data_mrna_seq_tpm.txt',)


# --------------------------------------------------------------------------- #
# 辅助函数
# --------------------------------------------------------------------------- #
def _find_package_by_label(cbioportal_data_dir, dataset_label):
    """在各子目录里找到 split('.')[0] == dataset_label 的 .tar.gz，返回完整路径。"""
    for d in _list_datasets(cbioportal_data_dir):
        ddir = os.path.join(cbioportal_data_dir, d)
        for f in os.listdir(ddir):
            if f.endswith('.tar.gz') and f.split('.')[0] == dataset_label:
                return os.path.join(ddir, f)
    return None


def _list_datasets(cbioportal_data_dir):
    return [d for d in os.listdir(cbioportal_data_dir) if os.path.isdir(os.path.join(cbioportal_data_dir, d))]


def _pick_rna_file(member_names, rna_filenames):
    """挑选 RNA 表达文件：排除 z-score 文件，优先匹配 rna_filenames。"""
    candidates = [f for f in member_names if
        'data_mrna' in f and ('tpm' in f or 'rsem' in f) and 'zscore' not in f.lower()]
    if not candidates:
        return None
    preferred = [f for f in candidates if f.split('/')[-1] in rna_filenames]
    return preferred[0] if preferred else candidates[0]


def _load_sample2patient(handle, dataset_label, tumor_only=True):
    """读取 data_clinical_sample，返回 SAMPLE_ID -> PATIENT_ID 映射 (Series)。"""
    sample_info_file = f'{dataset_label}/data_clinical_sample.txt'
    sample_info_df = pd.read_csv(handle.extractfile(sample_info_file), sep='\t', skiprows=4, index_col='SAMPLE_ID')
    if tumor_only and 'SAMPLE_TYPE' in sample_info_df.columns:
        mask = ~sample_info_df['SAMPLE_TYPE'].astype(str).str.contains('normal', case=False, na=False)
        sample_info_df = sample_info_df[mask]
    return sample_info_df['PATIENT_ID']


def _load_expression_sample_level(handle, rna_exp_file):
    """读取表达矩阵，返回 DataFrame(index=SAMPLE_ID, columns=gene)。"""
    expr = pd.read_csv(handle.extractfile(rna_exp_file), sep='\t')
    gene_col = expr.columns[0]  # 一般是 Hugo_Symbol
    expr = expr.drop(columns=[c for c in ('Entrez_Gene_Id',) if c in expr.columns])
    expr = expr.set_index(gene_col)
    expr = expr.apply(pd.to_numeric, errors='coerce')
    expr = expr.groupby(level=0).mean()  # 同名基因取均值
    return expr.T  # -> sample x gene


def _build_survival_df(patient_info_df, endpoints=None):
    """构建病人级生存表，返回 (sur_df, [(endpoint, time_col, status_col), ...])。"""
    time_cols = patient_info_df.columns[patient_info_df.columns.str.contains('MONTHS')]
    triples = []
    for tcol in time_cols:
        endpoint = tcol.split('_')[0]  # OS / DFS / PFS / DSS ...
        scol = endpoint + '_STATUS'
        if scol not in patient_info_df.columns:
            continue
        if endpoints is not None and endpoint not in endpoints:
            continue
        triples.append((endpoint, tcol, scol))
    if not triples:
        return None, []

    cols = sorted({c for _, t, s in triples for c in (t, s)})
    sur_df = patient_info_df[cols].copy()
    sur_df = sur_df.replace({'[Not Available]': np.nan, '[Not Applicable]': np.nan})
    for _, tcol, scol in triples:
        sur_df[tcol] = pd.to_numeric(sur_df[tcol], errors='coerce')
        # 期望形如 "1:DECEASED"/"0:LIVING"；取冒号前数字，非数值 -> NaN
        sur_df[scol] = pd.to_numeric(sur_df[scol].astype(str).str.split(':').str[0], errors='coerce')
    sur_df = sur_df.dropna(how='all')
    return sur_df, triples


def _median_survival(durations, events):
    """KM 估计的中位生存期；事件不足无法达到中位时返回 inf（lifelines 行为）。"""
    from lifelines import KaplanMeierFitter
    if len(durations) == 0:
        return np.nan
    kmf = KaplanMeierFitter()
    kmf.fit(durations, event_observed=events)
    return kmf.median_survival_time_


def _bh_fdr(pvals):
    """Benjamini-Hochberg FDR，忽略 NaN。"""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


# --------------------------------------------------------------------------- #
# 主函数
# --------------------------------------------------------------------------- #
def screen_gene_survival_across_datasets(cbioportal_data_dir, genes, km_best_cutoff, datasets=None,
        rna_filenames=DEFAULT_RNA_FILENAMES, endpoints=None,  # 例如 ['OS']；None 则使用全部 *_MONTHS 终点
        tumor_only=True,  # 仅保留肿瘤样本（依据 data_clinical_sample 的 SAMPLE_TYPE）
        min_group_size=10,  # 高/低分组的最小样本数，过小则跳过该 (基因, 终点)
        add_fdr=True,  # 是否对汇总 p 值做 BH-FDR 校正
        verbose=True, ):
    """
    在多个 cBioPortal 数据集中粗筛基因对生存的影响。

    km_best_cutoff : callable
        签名 ``km_best_cutoff(df, value, time, status) -> (cutoff, result_df)``，
        其中 result_df 含一列 ``'pvalue'``（每个候选 cutoff 的 logrank p）。
        最优 p 取 ``result_df['pvalue'].min()``。

    返回 DataFrame：每行一个 (dataset, gene, endpoint)，含 cutoff、pvalue、
    样本/事件数；若 add_fdr 则附 BH 校正 q 值并按 p 排序。
    """
    if datasets is None:
        datasets = _list_datasets(cbioportal_data_dir)
        if verbose:
            print(f"Available datasets: {datasets}")

    genes = list(genes)
    records = []

    for dataset in datasets:
        dataset_dir = os.path.join(cbioportal_data_dir, dataset)
        if not os.path.isdir(dataset_dir):
            continue
        tarballs = [f for f in os.listdir(dataset_dir) if f.endswith('.tar.gz')]
        if not tarballs:
            continue
        pkg = os.path.join(dataset_dir, tarballs[0])
        dataset_label = tarballs[0].split('.')[0]

        try:
            with tarfile.open(pkg, 'r:gz') as handle:
                patient_info_df = pd.read_csv(handle.extractfile(f'{dataset_label}/data_clinical_patient.txt'),
                    sep='\t', skiprows=4, index_col='PATIENT_ID')
                sur_df, triples = _build_survival_df(patient_info_df, endpoints)
                if sur_df is None:
                    if verbose:
                        print(f"  [skip] 无生存终点: {dataset_label}")
                    continue

                rna_exp_file = _pick_rna_file(handle.getnames(), rna_filenames)
                if rna_exp_file is None:
                    if verbose:
                        print(f"  [skip] 无 RNA 表达文件: {dataset_label}")
                    continue
                if verbose:
                    print(f"{dataset_label} -> 终点 {[e for e, _, _ in triples]}; "
                          f"表达 {rna_exp_file}")

                sample2patient = _load_sample2patient(handle, dataset_label, tumor_only)
                expr_sample = _load_expression_sample_level(handle, rna_exp_file)
        except Exception as e:
            if verbose:
                print(f"  [error] 读取 {dataset_label} 失败: {e}")
            continue

        # 把病人级生存广播到样本级（reindex 容忍缺失病人，不会 KeyError）
        sample_surv = sur_df.reindex(sample2patient.values)
        sample_surv.index = sample2patient.index  # 回到 SAMPLE_ID

        for gene in genes:
            if gene not in expr_sample.columns:
                continue
            gene_expr = expr_sample[gene]
            for endpoint, tcol, scol in triples:
                try:
                    km_df = sample_surv[[tcol, scol]].copy()
                    km_df[gene] = gene_expr.reindex(km_df.index)
                    km_df = km_df.dropna(how='any')
                    if len(km_df) < 2 * min_group_size:
                        continue

                    cutoff, cut_df = km_best_cutoff(km_df, value=gene, time=tcol, status=scol)
                    best_p = float(cut_df['pvalue'].min())

                    # km_best_cutoff 用 (value > cutoff) 划高组，这里保持一致
                    vals = km_df[gene]
                    high = vals > cutoff
                    n_high = int(high.sum())
                    n_low = int((~high).sum())
                    if min(n_high, n_low) < min_group_size:
                        continue

                    med_high = _median_survival(km_df.loc[high, tcol], km_df.loc[high, scol])
                    med_low = _median_survival(km_df.loc[~high, tcol], km_df.loc[~high, scol])
                    # 方向：高表达组中位生存更短 -> 高表达预后更差
                    if np.isnan(med_high) or np.isnan(med_low) or med_high == med_low:
                        direction = 'unclear'
                    elif med_high < med_low:
                        direction = 'high_worse'
                    else:
                        direction = 'high_better'

                    records.append({'dataset': dataset_label, 'gene': gene, 'endpoint': endpoint, 'time_col': tcol,
                        'status_col': scol, 'cutoff': cutoff, 'pvalue': best_p, 'n_total': len(km_df),
                        'n_events': int(km_df[scol].sum()), 'n_high': n_high, 'n_low': n_low,
                        'median_surv_high': med_high, 'median_surv_low': med_low, 'direction': direction, })
                except Exception as e:
                    if verbose:
                        print(f"  [error] {dataset_label}/{gene}/{endpoint}: {e}")
                    continue

    result = pd.DataFrame.from_records(records)
    if not result.empty and add_fdr:
        result['fdr_bh'] = _bh_fdr(result['pvalue'].values)
        result = result.sort_values('pvalue').reset_index(drop=True)
    return result


# --------------------------------------------------------------------------- #
# 可选：约束版 best-cutoff（避免退化分组带来的虚假极小 p 值）
# 需要 from lifelines.statistics import logrank_test
# --------------------------------------------------------------------------- #
def km_best_cutoff_constrained(df, value, time='time', status='status', q_low=0.1, q_high=0.9, min_group=10):
    """只在 [q_low, q_high] 分位区间内搜索 cutoff，且要求两组各 >= min_group。

    返回 (best_cutoff, result_df)，result_df 含 cutoff / pvalue / n_high / n_low，
    与 screen_gene_survival_across_datasets 期望的 'pvalue' 列兼容。
    """
    from lifelines.statistics import logrank_test
    d = df[[value, time, status]].dropna().copy()
    lo, hi = d[value].quantile([q_low, q_high])
    candidates = sorted(d.loc[(d[value] >= lo) & (d[value] <= hi), value].unique())
    rows = []
    for cutoff in candidates:
        grp = d[value] > cutoff
        n_high, n_low = int(grp.sum()), int((~grp).sum())
        if min(n_high, n_low) < min_group:
            continue
        res = logrank_test(durations_A=d.loc[grp, time], durations_B=d.loc[~grp, time],
            event_observed_A=d.loc[grp, status], event_observed_B=d.loc[~grp, status], )
        rows.append({'cutoff': cutoff, 'pvalue': res.p_value, 'n_high': n_high, 'n_low': n_low})
    if not rows:
        return np.nan, pd.DataFrame(columns=['cutoff', 'pvalue', 'n_high', 'n_low'])
    out = pd.DataFrame(rows).sort_values('pvalue').reset_index(drop=True)
    return out.loc[0, 'cutoff'], out


# --------------------------------------------------------------------------- #
# 根据 screen 结果中的一行复现 KM 图
# --------------------------------------------------------------------------- #
def plot_km_for_row(row, cbioportal_data_dir, kaplan_meier, rna_filenames=DEFAULT_RNA_FILENAMES, tumor_only=True,
                    high_label='High', low_label='Low', savefig=None, **km_kwargs):
    """根据 screen 结果中的一行，重新读取数据、按 cutoff 分高/低组并画 KM 图。

    参数
    ----
    row : pandas.Series | dict
        screen_gene_survival_across_datasets 结果中的一行，需含
        dataset / gene / time_col / status_col / cutoff（endpoint 可选）。
    cbioportal_data_dir : str
        cBioPortal 数据根目录。
    kaplan_meier : callable
        你的 KM 绘图函数，签名
        ``kaplan_meier(grouped_df, groupby, time, status, groups=..., **kwargs)``，
        返回 ``(matrix, p_value, fig)``。
    high_label / low_label : str
        高/低表达组的显示名（分组顺序为 [low, high]，即以低表达为 cox 参照）。
    savefig : str | None
        若给定则把图保存到该路径。
    **km_kwargs
        透传给 kaplan_meier（如 figsize、color_dict、cox_analysis、cox_ref 等）。

    返回 (matrix, p_value, fig)。
    """
    dataset_label = row['dataset']
    gene = row['gene']
    tcol = row['time_col']
    scol = row['status_col']
    cutoff = row['cutoff']
    endpoint = row['endpoint'] if 'endpoint' in row else tcol.split('_')[0]

    pkg = _find_package_by_label(cbioportal_data_dir, dataset_label)
    if pkg is None:
        raise FileNotFoundError(f'未找到 dataset_label={dataset_label} 对应的 .tar.gz')

    with tarfile.open(pkg, 'r:gz') as handle:
        patient_info_df = pd.read_csv(handle.extractfile(f'{dataset_label}/data_clinical_patient.txt'), sep='\t',
            skiprows=4, index_col='PATIENT_ID')
        sur_df, _ = _build_survival_df(patient_info_df, endpoints=[endpoint])
        if sur_df is None:
            raise ValueError(f'{dataset_label} 中找不到终点 {endpoint}')
        rna_exp_file = _pick_rna_file(handle.getnames(), rna_filenames)
        if rna_exp_file is None:
            raise FileNotFoundError(f'{dataset_label} 中找不到 RNA 表达文件')
        sample2patient = _load_sample2patient(handle, dataset_label, tumor_only)
        expr_sample = _load_expression_sample_level(handle, rna_exp_file)

    if gene not in expr_sample.columns:
        raise KeyError(f'{gene} 不在 {dataset_label} 的表达矩阵中')

    # 病人级生存广播到样本级，再与表达 join（与 screen 完全一致）
    sample_surv = sur_df.reindex(sample2patient.values)
    sample_surv.index = sample2patient.index
    km_df = sample_surv[[tcol, scol]].copy()
    km_df[gene] = expr_sample[gene].reindex(km_df.index)
    km_df = km_df.dropna(how='any')

    # 按 cutoff 分组：value > cutoff -> 高组（与 km_best_cutoff 一致）
    grouped = km_df[[tcol, scol]].copy()
    grouped[gene] = np.where(km_df[gene] > cutoff, high_label, low_label)

    matrix, p_value, fig = kaplan_meier(grouped, groupby=gene, time=tcol, status=scol, groups=[low_label, high_label],
        **km_kwargs)

    if savefig:
        fig.savefig(savefig, dpi=300, bbox_inches='tight')
    return matrix, p_value, fig
