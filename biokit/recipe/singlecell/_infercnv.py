# coding: utf-8
"""InferCNV-based Cancer cell / Normal Epithelial classification pipeline.

Pipeline overview
-----------------
1. Read gene coordinates from a GTF file and map GTF ``gene_id`` values to
   ``adata.var[feature_id_col]``. In the current AE0058 data, the field named
   ``transcript`` contains GTF ``gene_id`` values and is therefore used as the
   default feature identifier.
2. For each ``sample``, use non-epithelial annotated cells as inferCNV
   references and run ``infercnvpy.tl.infercnv`` with each configured
   ``(window_size, step)`` pair.
3. For epithelial cellbins, define the cell-level CNV burden as the mean of
   the absolute values of the inferred CNV windows:

       cnv_score = mean(abs(X_cnv), axis=window)

4. For every sample and CNV-score setting, fit 1-, 2-, and 3-component GMMs
   to the raw epithelial ``cnv_score`` distribution. The two-component model
   (GMM-2) is the primary classification model. GMM-1 and GMM-3 are retained
   as diagnostics only: they report whether one or three components provide a
   better BIC fit, but never veto GMM-2 classification.
5. For GMM-2, order components by mean CNV burden. The weighted-density
   intersection between the low- and high-mean components is retained as the
   one-dimensional ``gmm2_threshold``. Each cellbin is classified directly by
   its GMM-2 posterior probability: posterior(high component) >= 0.5 is
   ``Cancer cell``; lower posterior is ``Normal Epithelial``.
6. Estimate the empirical score density with KDE and detect KDE modes only for
   visual audit. KDE valley location and peak topology are shown in all plots,
   but they do not determine eligibility, window selection, or cell labels.
7. Select the best CNV-score setting separately for every sample by GMM-2
   separation quality (component balance, Ashman's D, and posterior confidence)
   and assign ``Cancer cell`` / ``Normal Epithelial`` labels for every sample.

Method references
-----------------
1. Sturm G, et al. infercnvpy documentation. ``infercnvpy.tl.infercnv``
   calculates expression deviations from reference cells, smooths them along
   genomic position in windows, centres the signal by cell, and applies noise
   filtering. https://infercnvpy.readthedocs.io/
2. Liu T, et al. Single-cell profiling of primary and paired metastatic lymph
   node tumors in breast cancer patients. Nature Communications. 2022.
   doi:10.1038/s41467-022-34581-2. This study used a bimodal malignancy-score
   threshold estimated with ``scCancer::getBimodalThres``. The present pipeline
   follows the same bimodal-threshold principle but explicitly implements the
   KDE two-mode valley on each sample and CNV-score setting.
3. Schwarz G. Estimating the dimension of a model. Annals of Statistics.
   1978;6:461-464. BIC is used to compare the 1-, 2-, and 3-component GMMs.
4. Ashman KM, Bird CM, Zepf SE. Detecting bimodality in astronomical datasets.
   Astronomical Journal. 1994;108:2348-2351. Ashman's D is used as a scale-free
   measure of separation between the two fitted Gaussian components.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union
import os
import time

import anndata as ad
import gtfparse
import infercnvpy as cnv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.signal import find_peaks
from scipy.special import logsumexp
from scipy.stats import gaussian_kde, norm
from sklearn.mixture import GaussianMixture
from threadpoolctl import threadpool_limits

__version__ = '2026-07-07-gmm2-primary-v7'

OverwriteType = Union[bool, Mapping[str, bool]]

REFERENCES = {'infercnvpy': 'Sturm G, et al. infercnvpy documentation. https://infercnvpy.readthedocs.io/',
              'bimodal_threshold': 'Liu T, et al. Nat Commun. 2022. doi:10.1038/s41467-022-34581-2.',
              'bic': 'Schwarz G. Ann Stat. 1978;6:461-464.',
              'ashman_d': 'Ashman KM, Bird CM, Zepf SE. Astron J. 1994;108:2348-2351.'}


def _should_overwrite(overwrite: OverwriteType, step: str) -> bool:
    """Return whether one pipeline step must be recalculated.

    Parameters
    ----------
    overwrite
        ``True`` reruns every cached step. ``False`` reuses existing files.
        A mapping can target individual steps, for example
        ``{'cnv_score': True, 'cutoff': False}``.
    step
        One of ``'cnv_score'``, ``'cutoff'``, ``'selection'``, or ``'plot'``.

    Returns
    -------
    bool
        Whether the named step is explicitly requested for recalculation.
    """
    if isinstance(overwrite, bool):
        return overwrite
    return bool(overwrite.get(step, False))


def _safe_name(value: str) -> str:
    """Convert an arbitrary sample name into a filesystem-safe file stem.

    Parameters
    ----------
    value
        Sample identifier or any other string used in an output filename.

    Returns
    -------
    str
        Value with path separators and whitespace replaced by underscores.
    """
    return str(value).replace('/', '_').replace('\\', '_').replace(' ', '_')


def _savefig(fig: plt.Figure, outdir: Union[str, Path], name: str) -> None:
    """Save one matplotlib figure as both PNG and PDF, then close it.

    Parameters
    ----------
    fig
        Figure object to save. ``fig.tight_layout()`` should be called by the
        plotting function before this helper.
    outdir
        Directory receiving the image files. It is created when absent.
    name
        File stem without extension. The function writes ``{name}.png`` at
        300 dpi and ``{name}.pdf`` using the same bounding box.

    Returns
    -------
    None
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f'{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(outdir / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def _score_col(window_size: int) -> str:
    """Return the CNV-score column name used for one genomic window size.

    Parameters
    ----------
    window_size
        Number of genes included in each inferCNV smoothing window.

    Returns
    -------
    str
        Column name in the format ``cnv_score_w{window_size}``.

    Notes
    -----
    ``window_configs`` should contain one unique ``window_size`` per setting.
    The corresponding ``step`` remains recorded in the cutoff-result table.
    """
    return f'cnv_score_w{window_size}'


def add_gene_coordinates(adata: ad.AnnData, gtf_file: str, feature_id_col: str = 'transcript',
                         gtf_feature_attr='gene_id', chromosome_prefix: str = 'chr') -> np.ndarray:
    """Map GTF gene coordinates to ``adata.var`` and return mappable features.

    Parameters
    ----------
    adata
        AnnData object whose ``var[feature_id_col]`` values match GTF
        ``gene_id`` values. The default field is called ``transcript`` for
        historical reasons in AE0058, but its values are GTF ``gene_id`` values.
    gtf_file
        GTF annotation file containing at least ``feature``, ``seqname``,
        ``start``, ``end``, and ``gene_id`` columns.
    feature_id_col
        Column in ``adata.var`` used to join GTF ``gene_id`` coordinates.
    chromosome_prefix
        Prefix applied to GTF chromosome names. With the default ``'chr'``,
        GTF chromosome ``X`` becomes ``chrX`` and is compatible with
        ``exclude_chromosomes=('chrX', 'chrY')``.

    Returns
    -------
    numpy.ndarray
        Boolean vector aligned to ``adata.var``. ``True`` marks features with
        non-missing chromosome, start, and end coordinates after the join.

    Notes
    -----
    inferCNV smooths expression according to genomic order; therefore every
    feature used by inferCNV must have a valid chromosome and genomic position.
    """
    gtf_df = gtfparse.read_gtf(gtf_file, result_type='pandas')
    gene_pos_df = gtf_df.loc[gtf_df['feature'].eq('gene'), ['seqname', 'start', 'end', gtf_feature_attr]].dropna(
        subset=gtf_feature_attr).drop_duplicates(gtf_feature_attr).rename(columns={'seqname': 'chromosome'}).set_index(
        gtf_feature_attr)
    gene_pos_df['chromosome'] = gene_pos_df['chromosome'].astype(str)
    if chromosome_prefix:
        gene_pos_df['chromosome'] = np.where(gene_pos_df['chromosome'].str.startswith(chromosome_prefix),
                                             gene_pos_df['chromosome'], chromosome_prefix + gene_pos_df['chromosome'])
    adata.var = adata.var.drop(columns=['chromosome', 'start', 'end'], errors='ignore').join(
        gene_pos_df[['chromosome', 'start', 'end']], on=feature_id_col)
    return adata.var['chromosome'].notna().to_numpy()


def calculate_cnv_scores(adata: ad.AnnData, cnv_feature_mask: np.ndarray, window_configs: Sequence[Tuple[int, int]],
                         cell_type_col: str = 'cell_type', sample_col: str = 'sample', epi_name: str = 'Epithelial',
                         unassigned_name: str = 'Unassigned', chipid_col: Optional[str] = 'chipid',
                         layer: str = 'log1p', lfc_clip: float = 3, dynamic_threshold: float = 1.5,
                         exclude_chromosomes: Sequence[str] = ('chrX', 'chrY'), chunksize: int = 1000,
                         n_jobs: int = 8) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate epithelial cellbin CNV burden for every sample and window setting.

    Parameters
    ----------
    adata
        AnnData object with coordinates already added to ``adata.var`` and cell
        annotations in ``adata.obs[cell_type_col]``.
    cnv_feature_mask
        Boolean vector aligned to ``adata.var``. Only features with mapped GTF
        coordinates are used in inferCNV.
    window_configs
        Sequence of ``(window_size, step)`` pairs passed to
        ``infercnvpy.tl.infercnv``. ``window_size`` is the number of genes per
        genomic running window; ``step`` is the number of genes advanced before
        computing the next window. Each ``window_size`` should be unique because
        output score columns are named ``cnv_score_w{window_size}``.
    cell_type_col
        ``adata.obs`` column containing the main cell-type annotation.
    sample_col
        ``adata.obs`` column defining independent spatial samples. inferCNV is
        run separately for each sample.
    epi_name
        Cell-type label identifying the epithelial cellbins to classify.
    unassigned_name
        Cell-type label excluded before inferCNV. Set to ``None`` to retain all
        cells in each sample.
    chipid_col
        Optional chip identifier retained in output score tables.
    layer
        AnnData layer supplied to inferCNV. It should contain normalized,
        log-transformed expression values; the default is ``'log1p'``.
    lfc_clip
        Absolute clipping threshold for expression deviations from reference
        cells in infercnvpy.
    dynamic_threshold
        infercnvpy noise-filter multiplier. Values below
        ``dynamic_threshold × SD`` of smoothed expression are set to zero.
    exclude_chromosomes
        Chromosomes excluded from inferCNV, typically sex chromosomes.
    chunksize
        Number of cells processed per infercnvpy chunk.
    n_jobs
        Number of worker processes passed to infercnvpy.

    Returns
    -------
    cnv_score_df
        Per-epithelial-cellbin metadata plus one score column per window size.
        Each score is the mean absolute inferred CNV value across genomic
        windows for that epithelial cellbin.
    cnv_score_summary_df
        Sample-level mean, median, and standard deviation for every CNV score.

    Notes
    -----
    For each sample, all non-epithelial annotated cells are supplied as
    ``reference_cat``. infercnvpy subtracts reference expression, smooths the
    result along genomic position, centres each cell, and applies noise
    filtering before producing ``X_{key_added}``.

    References
    ----------
    See ``REFERENCES['infercnvpy']`` and the infercnvpy method documentation.
    """
    cnv_score_list = []
    score_columns = [_score_col(window_size) for window_size, _ in window_configs]
    samples = adata.obs[sample_col].dropna().astype(str).unique()
    for sample in samples:
        sample_mask = adata.obs[sample_col].astype(str).eq(sample)
        if unassigned_name is not None:
            sample_mask &= adata.obs[cell_type_col].astype(str).ne(unassigned_name)
        sample_adata = adata[sample_mask.to_numpy(), cnv_feature_mask].copy()
        reference_cat = sample_adata.obs.loc[
            sample_adata.obs[cell_type_col].astype(str).ne(epi_name), cell_type_col].astype(str).unique().tolist()
        epithelial_mask = sample_adata.obs[cell_type_col].astype(str).eq(epi_name).to_numpy()
        metadata_cols = [col for col in
                         [chipid_col, sample_col, cell_type_col, 'leiden', 'x', 'y', 'x_slide_mm', 'y_slide_mm'] if
                         col is not None and col in sample_adata.obs.columns]
        sample_score_df = sample_adata.obs.loc[epithelial_mask, metadata_cols].copy()
        for window_size, step in window_configs:
            key_added = f'cnv_w{window_size}'
            cnv.tl.infercnv(sample_adata, reference_key=cell_type_col, reference_cat=reference_cat, layer=layer,
                            lfc_clip=lfc_clip, window_size=window_size, step=step, dynamic_threshold=dynamic_threshold,
                            exclude_chromosomes=exclude_chromosomes, chunksize=chunksize, n_jobs=n_jobs,
                            key_added=key_added)
            cnv_matrix = sample_adata.obsm[f'X_{key_added}']
            sample_score_df[_score_col(window_size)] = np.asarray(
                np.abs(cnv_matrix[epithelial_mask]).mean(axis=1)).ravel()
            del sample_adata.obsm[f'X_{key_added}']
            del sample_adata.uns[key_added]
        cnv_score_list.append(sample_score_df)
        print(f'{sample}: {sample_score_df.shape[0]} epithelial cells completed')
    cnv_score_df = pd.concat(cnv_score_list, axis=0)
    summary_group_cols = [col for col in [chipid_col, sample_col] if col is not None and col in cnv_score_df.columns]
    cnv_score_summary_df = cnv_score_df.groupby(summary_group_cols)[score_columns].agg(['mean', 'median', 'std'])
    cnv_score_summary_df.columns = [f'{score}_{statistic}' for score, statistic in cnv_score_summary_df.columns]
    return cnv_score_df, cnv_score_summary_df


def _gmm_boundary(grid: np.ndarray, component_density: np.ndarray, gmm_means: np.ndarray) -> float:
    """Calculate the one-dimensional GMM two-component density intersection.

    Parameters
    ----------
    grid
        Ordered score values at which component densities were evaluated.
    component_density
        Two-row array containing weighted low- and high-component GMM density
        values over ``grid``.
    gmm_means
        Ordered two-element array of low- and high-component GMM means.

    Returns
    -------
    float
        Score where the weighted component densities cross between the two GMM
        means, or ``numpy.nan`` when no crossing is found.

    Notes
    -----
    This value is retained as a diagnostic only. It is not the final threshold,
    because it may differ from the empirical KDE valley when the two Gaussian
    components have unequal widths, unequal weights, or overlap strongly.
    """
    component_diff = component_density[0] - component_density[1]
    cross_idx = np.flatnonzero(np.diff(np.signbit(component_diff)))
    cross_idx = cross_idx[(grid[cross_idx] > gmm_means[0]) & (grid[cross_idx] < gmm_means[1])]
    if len(cross_idx) == 0:
        return np.nan
    i = cross_idx[0]
    x0, x1 = grid[i], grid[i + 1]
    y0, y1 = component_diff[i], component_diff[i + 1]
    return x0 - y0 * (x1 - x0) / (y1 - y0)


def fit_gmm_kde_cutoff(score_values: np.ndarray, sample: str, score_col: str, window_size: int, step: int,
                       n_grid: int = 5000, bw_scale: float = 1.0, peak_prominence_fraction: float = 0.03,
                       peak_distance_fraction: float = 0.05, gmm_n_init: int = 50, gmm_max_iter: int = 1000,
                       gmm_tol: float = 1e-6, gmm_reg_covar: float = 1e-8, random_state: int = 0,
                       timing: bool = False) -> Dict[str, object]:
    """Fit GMM diagnostics and obtain a KDE bimodal CNV cutoff for one sample.

    Parameters
    ----------
    score_values
        One-dimensional epithelial CNV-score vector from one sample and one
        ``(window_size, step)`` setting.
    sample
        Sample identifier written to the returned result dictionary.
    score_col
        Name of the CNV-score column, usually ``cnv_score_w{window_size}``.
    window_size
        Number of genes in the inferCNV smoothing window that generated this
        score distribution.
    step
        Number of genes advanced between adjacent inferCNV windows.
    n_grid
        Number of equally spaced points between score minimum and maximum used
        to evaluate KDE and GMM densities. This changes plot/valley resolution
        only; it does not alter the fitted GMM.
    bw_scale
        Multiplicative factor applied to Scott's KDE bandwidth. ``1.0`` uses
        Scott's default bandwidth. A larger value produces a smoother KDE and
        suppresses minor local peaks; the same value must be used for all
        samples and window settings in one analysis.
    peak_prominence_fraction
        Minimum peak prominence passed to ``scipy.signal.find_peaks``, expressed
        as a fraction of the maximum KDE density.
    peak_distance_fraction
        Minimum horizontal separation between KDE peaks, expressed as a fraction
        of ``n_grid``.
    gmm_n_init
        Number of GMM initializations. The best likelihood solution is retained;
        increasing this number reduces sensitivity to local EM optima.
    gmm_max_iter
        Maximum number of EM iterations per GMM initialization.
    gmm_tol
        EM convergence tolerance for the average lower-bound improvement.
    gmm_reg_covar
        Non-negative variance regularization added to GMM components. It avoids
        variance collapse; because CNV scores are small, this should remain much
        smaller than the empirical component variance.
    random_state
        Random seed for reproducible GMM initialization.
    timing
        If ``True``, print immediate start/end messages and elapsed seconds for
        GMM-1, GMM-2, GMM-3, BIC, posterior calculation, KDE construction, KDE
        evaluation, and peak detection. The same timing values are also retained
        in the returned dictionary as ``time_*_sec`` columns.

    Returns
    -------
    dict
        One row of diagnostics, including GMM BIC values, component means,
        standard deviations, weights, Ashman's D, posterior-confidence rates,
        KDE peak locations, KDE valley threshold, selection-relevant flags, and
        stage-level elapsed times.

    Method
    ------
    **Primary model: GMM-2**

    - Fit a two-component Gaussian mixture to the raw ``cnv_score`` vector.
    - Order its components by mean. The lower-mean component is CNV-low and the
      higher-mean component is CNV-high.
    - Calculate ``gmm2_threshold`` as the weighted-component density
      intersection between the ordered means. This is the primary one-dimensional
      cutoff displayed in figures.
    - Later cell-level classification uses GMM-2 posterior probability directly:
      ``P(CNV-high | score) >= 0.5`` is ``Cancer cell`` and lower probability is
      ``Normal Epithelial``.

    **Reference diagnostics**

    - GMM-1 and GMM-3 are fitted only to report BIC comparisons with GMM-2.
      They do not gate primary-window selection or labels.
    - KDE is calculated only to display the empirical density and optional KDE
      valley. KDE mode count and valley depth are retained for audit but do not
      determine the selected GMM-2 model or classification.

    References
    ----------
    - Liu T, et al. Nat Commun. 2022. doi:10.1038/s41467-022-34581-2.
    - Schwarz G. Ann Stat. 1978;6:461-464.
    - Ashman KM, Bird CM, Zepf SE. Astron J. 1994;108:2348-2351.
    """
    time_total_start = time.perf_counter()
    score_values = np.asarray(score_values, dtype=float)
    score_matrix = score_values[:, None]
    task_label = f'[{sample} | {score_col}]'

    gmm_1 = GaussianMixture(n_components=1, covariance_type='full', n_init=1, max_iter=gmm_max_iter, tol=gmm_tol,
                            reg_covar=gmm_reg_covar, random_state=random_state)
    gmm_2 = GaussianMixture(n_components=2, covariance_type='full', n_init=gmm_n_init, max_iter=gmm_max_iter,
                            tol=gmm_tol, reg_covar=gmm_reg_covar, random_state=random_state)
    gmm_3 = GaussianMixture(n_components=3, covariance_type='full', n_init=gmm_n_init, max_iter=gmm_max_iter,
                            tol=gmm_tol, reg_covar=gmm_reg_covar, random_state=random_state)

    if timing:
        print(f'{task_label} GMM1 start | n_cells={len(score_values):,}, n_init=1', flush=True)
    time_start = time.perf_counter()
    gmm_1.fit(score_matrix)
    time_gmm_1 = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} GMM1 done: {time_gmm_1:.2f}s', flush=True)

    if timing:
        print(f'{task_label} GMM2 start | n_cells={len(score_values):,}, n_init={gmm_n_init}', flush=True)
    time_start = time.perf_counter()
    gmm_2.fit(score_matrix)
    time_gmm_2 = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} GMM2 done: {time_gmm_2:.2f}s', flush=True)

    if timing:
        print(f'{task_label} GMM3 start | n_cells={len(score_values):,}, n_init={gmm_n_init}', flush=True)
    time_start = time.perf_counter()
    gmm_3.fit(score_matrix)
    time_gmm_3 = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} GMM3 done: {time_gmm_3:.2f}s', flush=True)

    if timing:
        print(f'{task_label} BIC start', flush=True)
    time_start = time.perf_counter()
    bic_1 = gmm_1.bic(score_matrix)
    bic_2 = gmm_2.bic(score_matrix)
    bic_3 = gmm_3.bic(score_matrix)
    time_bic = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} BIC done: {time_bic:.2f}s', flush=True)

    component_order = np.argsort(gmm_2.means_.ravel())
    gmm_means = gmm_2.means_.ravel()[component_order]
    gmm_stds = np.sqrt(gmm_2.covariances_.reshape(-1)[component_order])
    gmm_weights = gmm_2.weights_.ravel()[component_order]

    if timing:
        print(f'{task_label} posterior start', flush=True)
    time_start = time.perf_counter()
    posterior = gmm_2.predict_proba(score_matrix)[:, component_order]
    posterior_confidence = posterior.max(axis=1)
    time_posterior = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} posterior done: {time_posterior:.2f}s', flush=True)

    if timing:
        print(f'{task_label} grid/GMM density start | n_grid={n_grid:,}', flush=True)
    time_start = time.perf_counter()
    grid = np.linspace(score_values.min(), score_values.max(), n_grid)
    gmm_component_density = np.vstack(
        [weight * norm.pdf(grid, loc=mean, scale=std) for mean, std, weight in zip(gmm_means, gmm_stds, gmm_weights)])
    gmm_mixture_density = gmm_component_density.sum(axis=0)
    time_grid_gmm_density = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} grid/GMM density done: {time_grid_gmm_density:.2f}s', flush=True)

    if timing:
        print(f'{task_label} KDE build start', flush=True)
    time_start = time.perf_counter()
    kde = gaussian_kde(score_values, bw_method=lambda kde_obj: kde_obj.scotts_factor() * bw_scale)
    time_kde_build = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} KDE build done: {time_kde_build:.2f}s', flush=True)

    if timing:
        print(f'{task_label} KDE evaluate start | n_grid={n_grid:,}', flush=True)
    time_start = time.perf_counter()
    kde_density = kde(grid)
    time_kde_evaluate = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} KDE evaluate done: {time_kde_evaluate:.2f}s', flush=True)

    if timing:
        print(f'{task_label} peak detection start', flush=True)
    time_start = time.perf_counter()
    min_peak_distance = max(1, int(n_grid * peak_distance_fraction))
    peak_idx, _ = find_peaks(kde_density, prominence=kde_density.max() * peak_prominence_fraction,
                             distance=min_peak_distance)
    gmm_mode_idx, _ = find_peaks(gmm_mixture_density, prominence=gmm_mixture_density.max() * peak_prominence_fraction,
                                 distance=min_peak_distance)
    time_peak_detect = time.perf_counter() - time_start
    if timing:
        print(f'{task_label} peak detection done: {time_peak_detect:.2f}s', flush=True)

    gmm2_threshold = _gmm_boundary(grid, gmm_component_density, gmm_means)

    result = {'sample': sample, 'score_col': score_col, 'window_size': window_size, 'step': step,
              'n_cells': len(score_values), 'score_min': score_values.min(), 'score_max': score_values.max(),
              'bic_1component': bic_1, 'bic_2component': bic_2, 'bic_3component': bic_3,
              'bic_delta_1minus2': bic_1 - bic_2, 'bic_delta_2minus3': bic_2 - bic_3,
              'gmm_1_converged': gmm_1.converged_, 'gmm_2_converged': gmm_2.converged_,
              'gmm_3_converged': gmm_3.converged_, 'low_gmm_mean': gmm_means[0], 'high_gmm_mean': gmm_means[1],
              'low_gmm_std': gmm_stds[0], 'high_gmm_std': gmm_stds[1], 'low_gmm_weight': gmm_weights[0],
              'high_gmm_weight': gmm_weights[1],
              'ashman_d': np.sqrt(2) * abs(gmm_means[1] - gmm_means[0]) / np.sqrt(gmm_stds[0] ** 2 + gmm_stds[1] ** 2),
              'posterior_high_conf_fraction': (posterior_confidence >= 0.9).mean(),
              'posterior_uncertain_fraction': (posterior_confidence < 0.8).mean(), 'gmm_boundary': gmm2_threshold,
              'gmm2_threshold': gmm2_threshold, 'gmm2_primary_model': True, 'n_kde_peaks': len(peak_idx),
              'n_gmm_mixture_modes': len(gmm_mode_idx), 'time_gmm_1_sec': time_gmm_1, 'time_gmm_2_sec': time_gmm_2,
              'time_gmm_3_sec': time_gmm_3, 'time_gmm_total_sec': time_gmm_1 + time_gmm_2 + time_gmm_3,
              'time_bic_sec': time_bic, 'time_posterior_sec': time_posterior,
              'time_grid_gmm_density_sec': time_grid_gmm_density, 'time_kde_build_sec': time_kde_build,
              'time_kde_evaluate_sec': time_kde_evaluate, 'time_peak_detect_sec': time_peak_detect}

    if len(peak_idx) == 2:
        left_peak_idx, right_peak_idx = np.sort(peak_idx)
        valley_idx = left_peak_idx + 1 + np.argmin(kde_density[left_peak_idx + 1:right_peak_idx])
        low_kde_peak = grid[left_peak_idx]
        high_kde_peak = grid[right_peak_idx]
        kde_valley_threshold = grid[valley_idx]
        result['status'] = 'bimodal'
        result['low_kde_peak'] = low_kde_peak
        result['high_kde_peak'] = high_kde_peak
        result['kde_valley_threshold'] = kde_valley_threshold
        result['valley_depth'] = 1 - kde_density[valley_idx] / min(kde_density[left_peak_idx],
                                                                   kde_density[right_peak_idx])
        result['gmm_kde_boundary_offset_ratio'] = abs(result['gmm2_threshold'] - kde_valley_threshold) / (
                high_kde_peak - low_kde_peak) if np.isfinite(result['gmm_boundary']) else np.nan
    else:
        result['status'] = 'unimodal' if len(peak_idx) < 2 else 'multimodal'
        result['low_kde_peak'] = np.nan
        result['high_kde_peak'] = np.nan
        result['kde_valley_threshold'] = np.nan
        result['valley_depth'] = np.nan
        result['gmm_kde_boundary_offset_ratio'] = np.nan

    result['time_total_sec'] = time.perf_counter() - time_total_start
    if timing:
        print(f'{task_label} finished | status={result["status"]}, total={result["time_total_sec"]:.2f}s', flush=True)
    return result


def build_gmm_kde_tasks(cnv_score_df: pd.DataFrame, window_configs: Sequence[Tuple[int, int]],
                        sample_col: str = 'sample') -> List[Dict[str, object]]:
    """Create one independent GMM/KDE fitting task per ``sample × CNV score``.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin CNV-score table returned by
        :func:`calculate_cnv_scores`.
    window_configs
        Ordered sequence of ``(window_size, step)`` pairs used for CNV-score
        calculation.
    sample_col
        Column identifying samples in ``cnv_score_df``.

    Returns
    -------
    list of dict
        Each task contains one sample identifier, one CNV-score column, its
        ``window_size`` and ``step``, and only the one-dimensional score vector
        required for fitting. The full AnnData object and the full score table
        are deliberately not passed to worker processes.

    Notes
    -----
    The task is the natural parallel unit because GMM/KDE fitting for one
    sample and one CNV-score setting does not depend on any other sample or
    window. Passing only a float64 score vector keeps inter-process transfer
    small relative to the repeated GMM EM fits.
    """
    tasks = []
    task_index = 0
    for sample, sample_df in cnv_score_df.groupby(sample_col, sort=False):
        for window_size, step in window_configs:
            score_col = _score_col(window_size)
            score_values = np.ascontiguousarray(sample_df[score_col].dropna().to_numpy(dtype=np.float64))
            tasks.append(
                {'task_index': task_index, 'sample': str(sample), 'score_col': score_col, 'window_size': window_size,
                 'step': step, 'score_values': score_values})
            task_index += 1
    return tasks


def _fit_gmm_kde_task(task: Dict[str, object], n_grid: int, bw_scale: float, peak_prominence_fraction: float,
                      peak_distance_fraction: float, gmm_n_init: int, gmm_max_iter: int, gmm_tol: float,
                      gmm_reg_covar: float, random_state: int) -> Tuple[int, Dict[str, object]]:
    """Run one GMM/KDE task in a worker process with one BLAS/OpenMP thread.

    This helper is intentionally module-level so it is directly usable by
    ``ProcessPoolExecutor``. ``threadpool_limits(1)`` prevents every worker
    from independently expanding BLAS/OpenMP kernels to all CPU threads. With
    ``cutoff_n_jobs=8`` this keeps the cutoff stage at approximately eight CPU
    threads rather than creating an 8 × N-thread oversubscription.
    """
    with threadpool_limits(limits=1):
        result = fit_gmm_kde_cutoff(score_values=task['score_values'], sample=task['sample'],
                                    score_col=task['score_col'], window_size=task['window_size'], step=task['step'],
                                    n_grid=n_grid, bw_scale=bw_scale, peak_prominence_fraction=peak_prominence_fraction,
                                    peak_distance_fraction=peak_distance_fraction, gmm_n_init=gmm_n_init,
                                    gmm_max_iter=gmm_max_iter, gmm_tol=gmm_tol, gmm_reg_covar=gmm_reg_covar,
                                    random_state=random_state, timing=False)
    result['task_index'] = task['task_index']
    result['worker_pid'] = os.getpid()
    return task['task_index'], result


def calculate_gmm_kde_cutoffs(cnv_score_df: pd.DataFrame, window_configs: Sequence[Tuple[int, int]],
                              sample_col: str = 'sample', n_grid: int = 5000, bw_scale: float = 1.0,
                              peak_prominence_fraction: float = 0.03, peak_distance_fraction: float = 0.05,
                              gmm_n_init: int = 50, gmm_max_iter: int = 1000, gmm_tol: float = 1e-6,
                              gmm_reg_covar: float = 1e-8, random_state: int = 0, cutoff_n_jobs: int = 1,
                              timing: bool = False) -> pd.DataFrame:
    """Run sample-by-score GMM/KDE cutoff fitting sequentially or in processes.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin table returned by :func:`calculate_cnv_scores`.
    window_configs
        Ordered ``(window_size, step)`` settings used to create
        ``cnv_score_df``.
    sample_col
        Column identifying samples in ``cnv_score_df``.
    n_grid, bw_scale, peak_prominence_fraction, peak_distance_fraction
        KDE-grid and peak-detection parameters passed unchanged to
        :func:`fit_gmm_kde_cutoff`.
    gmm_n_init, gmm_max_iter, gmm_tol, gmm_reg_covar, random_state
        GMM fitting parameters passed unchanged to :func:`fit_gmm_kde_cutoff`.
    cutoff_n_jobs
        Number of operating-system processes used only for the GMM/KDE cutoff
        stage. ``1`` runs tasks sequentially. Values above ``1`` create one
        process pool and submit one independent task per ``sample × CNV score``.
        This parameter is intentionally separate from ``n_jobs`` in
        :func:`infercnv_cancer_normal_recipe`, which belongs to infercnvpy when
        the upstream CNV scores themselves are recomputed.
    timing
        If ``True``, print task completion order, worker PID, status, and
        stage-level elapsed times. Each task's timing values are retained in
        the returned table regardless of this switch.

    Returns
    -------
    pandas.DataFrame
        One diagnostic row per ``sample × window_config`` combination. The
        rows are returned in deterministic task order, even though parallel
        tasks complete in variable order.

    Parallelization strategy
    ------------------------
    Each task uses only one sample's one-dimensional ``cnv_score`` vector. Its
    1-, 2-, and 3-component GMM fits and KDE calculations are independent from
    all other tasks, so they are scheduled with ``ProcessPoolExecutor``.
    Processes, rather than Python threads, are used because the expensive
    operations are CPU-bound SciPy/scikit-learn calls and because separate
    processes avoid GIL contention and isolate BLAS/OpenMP thread limits.

    Worker processes call ``threadpool_limits(limits=1)``. This is essential:
    without it, e.g. 8 worker processes may each attempt to use all available
    BLAS threads, producing severe oversubscription and slower total runtime.

    Notes
    -----
    The implementation uses the Linux ``fork`` multiprocessing context. This
    matches the current server environment and lets the imported pipeline
    functions be used directly from a Jupyter kernel. The worker payload is
    only the relevant float64 score vector, not the multi-million-cell AnnData
    object.
    """
    tasks = build_gmm_kde_tasks(cnv_score_df, window_configs=window_configs, sample_col=sample_col)
    task_total = len(tasks)
    result_by_index = {}

    if cutoff_n_jobs == 1:
        for finished_count, task in enumerate(tasks, start=1):
            if timing:
                print(f'[{finished_count}/{task_total}] {task["sample"]} {task["score_col"]} cutoff fitting start',
                      flush=True)
            task_index, result = _fit_gmm_kde_task(task, n_grid=n_grid, bw_scale=bw_scale,
                                                   peak_prominence_fraction=peak_prominence_fraction,
                                                   peak_distance_fraction=peak_distance_fraction, gmm_n_init=gmm_n_init,
                                                   gmm_max_iter=gmm_max_iter, gmm_tol=gmm_tol,
                                                   gmm_reg_covar=gmm_reg_covar, random_state=random_state)
            result_by_index[task_index] = result
            if timing:
                print(f'[{finished_count}/{task_total}] {result["sample"]} {result["score_col"]} done | '
                      f'pid={result["worker_pid"]}, status={result["status"]}, '
                      f'GMM={result["time_gmm_total_sec"]:.2f}s, '
                      f'KDE={result["time_kde_evaluate_sec"]:.2f}s, total={result["time_total_sec"]:.2f}s', flush=True)
            else:
                print(f'{result["sample"]} {result["score_col"]}: {result["status"]}')
    else:
        if timing:
            print(f'GMM/KDE cutoff parallel mode: {task_total} sample × score tasks, '
                  f'{cutoff_n_jobs} processes, one BLAS/OpenMP thread per process.', flush=True)
        mp_context = get_context('fork')
        with ProcessPoolExecutor(max_workers=cutoff_n_jobs, mp_context=mp_context) as executor:
            future_to_task = {executor.submit(_fit_gmm_kde_task, task, n_grid, bw_scale, peak_prominence_fraction,
                                              peak_distance_fraction, gmm_n_init, gmm_max_iter, gmm_tol, gmm_reg_covar,
                                              random_state): task for task in tasks}
            for finished_count, future in enumerate(as_completed(future_to_task), start=1):
                task_index, result = future.result()
                result_by_index[task_index] = result
                if timing:
                    print(f'[{finished_count}/{task_total}] {result["sample"]} {result["score_col"]} done | '
                          f'pid={result["worker_pid"]}, status={result["status"]}, '
                          f'GMM1={result["time_gmm_1_sec"]:.2f}s, '
                          f'GMM2={result["time_gmm_2_sec"]:.2f}s, '
                          f'GMM3={result["time_gmm_3_sec"]:.2f}s, '
                          f'KDE={result["time_kde_evaluate_sec"]:.2f}s, total={result["time_total_sec"]:.2f}s',
                          flush=True)
                else:
                    print(f'{result["sample"]} {result["score_col"]}: {result["status"]}')

    result_rows = [result_by_index[task_index] for task_index in range(task_total)]
    return pd.DataFrame(result_rows).sort_values(['sample', 'window_size']).reset_index(drop=True)


def select_best_cnv_score(cutoff_df: pd.DataFrame, min_delta_bic_1vs2: float = 10, max_delta_bic_2vs3: float = 10,
                          min_ashman_d: float = 2, min_component_weight: float = 0.05,
                          min_valley_depth: float = 0.10) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Select one primary GMM-2 CNV-score model per sample.

    Parameters
    ----------
    cutoff_df
        One row per ``sample × window`` from :func:`calculate_gmm_kde_cutoffs`.
        It must contain GMM-2 component parameters, posterior-confidence
        metrics, and BIC diagnostics.
    min_delta_bic_1vs2
        Diagnostic threshold for flagging whether GMM-2 improves on GMM-1.
        This is recorded but is not a gate for classification or selection.
    max_delta_bic_2vs3
        Diagnostic threshold for flagging whether GMM-3 materially improves on
        GMM-2. This is recorded but is not a gate for classification or
        selection.
    min_ashman_d
        Diagnostic threshold for ``gmm2_quality_pass``. It does not veto final
        classification; every sample still receives its highest-quality GMM-2
        window.
    min_component_weight
        Diagnostic threshold for ``gmm2_quality_pass``. It does not veto final
        classification.
    min_valley_depth
        Diagnostic KDE-valley threshold retained for reporting only. KDE is not
        a criterion for GMM-2 selection.

    Returns
    -------
    selection_df
        ``cutoff_df`` plus GMM-2 quality metrics, reference-model diagnostics,
        ranks, and exactly one ``selected=True`` row per sample.
    best_df
        Exactly one GMM-2 primary model per sample. There is no
        ``no_valid_window`` state in this GMM-2-primary workflow.

    Selection method
    ----------------
    GMM-2 is the primary model for every tested window. Windows are ranked
    within each sample by a GMM-2-only composite quality score:

    ``component_balance × separation_score × posterior_high_conf_fraction``

    where ``component_balance = 2 × min(component weights)`` and
    ``separation_score = min(Ashman_D / 3, 1)``. Thus selection is driven by
    GMM-2 component balance, separation, and cell-level posterior certainty.
    GMM-1 and GMM-3 BIC comparisons, and KDE peak/valley metrics, are retained
    as reference diagnostics but do not alter the selected model or labels.
    """
    selection_df = cutoff_df.copy()
    if 'gmm2_threshold' not in selection_df.columns:
        selection_df['gmm2_threshold'] = selection_df['gmm_boundary']

    selection_df['component_balance'] = 2 * np.minimum(selection_df['low_gmm_weight'], selection_df['high_gmm_weight'])
    selection_df['separation_score'] = np.minimum(selection_df['ashman_d'] / 3, 1)
    selection_df['gmm2_quality_score'] = (
            selection_df['component_balance'] * selection_df['separation_score'] * selection_df[
        'posterior_high_conf_fraction'])
    selection_df['quality_score'] = selection_df['gmm2_quality_score']

    selection_df['gmm2_vs_gmm1_supported'] = selection_df['bic_delta_1minus2'] > min_delta_bic_1vs2
    selection_df['gmm3_vs_gmm2_supported'] = selection_df['bic_delta_2minus3'] > max_delta_bic_2vs3
    selection_df['two_component_better'] = selection_df['gmm2_vs_gmm1_supported']
    selection_df['two_components_sufficient'] = ~selection_df['gmm3_vs_gmm2_supported']
    selection_df['weight_pass'] = np.minimum(selection_df['low_gmm_weight'],
                                             selection_df['high_gmm_weight']) >= min_component_weight
    selection_df['separation_pass'] = selection_df['ashman_d'] >= min_ashman_d
    selection_df['valley_pass'] = selection_df['valley_depth'].fillna(-np.inf) >= min_valley_depth
    selection_df['gmm2_quality_pass'] = selection_df['weight_pass'] & selection_df['separation_pass']
    selection_df['eligible'] = selection_df['gmm2_quality_pass']

    selection_df['selection_rank'] = np.nan
    selection_df['selected'] = False
    selection_df['selection_method'] = 'GMM2 primary; GMM1/GMM3/KDE diagnostic only'

    best_rows = []
    for sample, sample_df in selection_df.groupby('sample', sort=False):
        ranking_df = sample_df.sort_values(
            ['gmm2_quality_score', 'ashman_d', 'posterior_high_conf_fraction', 'component_balance', 'window_size',
             'step'], ascending=[False, False, False, False, True, True])
        selection_df.loc[ranking_df.index, 'selection_rank'] = np.arange(1, len(ranking_df) + 1)
        selected_idx = ranking_df.index[0]
        selection_df.loc[selected_idx, 'selected'] = True
        best_row = selection_df.loc[selected_idx].to_dict()
        best_row['selection_status'] = 'selected_gmm2_primary'
        best_rows.append(best_row)

    best_df = pd.DataFrame(best_rows).sort_values('sample').reset_index(drop=True)
    return selection_df, best_df


def _gmm2_high_probability(score_values: np.ndarray, low_mean: float, low_std: float, low_weight: float,
                           high_mean: float, high_std: float, high_weight: float) -> np.ndarray:
    """Return posterior probability of the ordered high-CNV GMM-2 component.

    Parameters
    ----------
    score_values
        One-dimensional CNV-score vector.
    low_mean, low_std, low_weight
        Parameters of the lower-mean GMM-2 component.
    high_mean, high_std, high_weight
        Parameters of the higher-mean GMM-2 component.

    Returns
    -------
    numpy.ndarray
        Stable evaluation of ``P(high-CNV component | CNV score)`` for every
        score value.
    """
    log_low = np.log(low_weight) + norm.logpdf(score_values, loc=low_mean, scale=low_std)
    log_high = np.log(high_weight) + norm.logpdf(score_values, loc=high_mean, scale=high_std)
    return np.exp(log_high - logsumexp(np.vstack([log_low, log_high]), axis=0))


def assign_best_cnv_labels(cnv_score_df: pd.DataFrame, best_df: pd.DataFrame,
                           sample_col: str = 'sample') -> pd.DataFrame:
    """Assign Cancer/Normal labels directly from each sample's selected GMM-2.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin score table from :func:`calculate_cnv_scores`.
    best_df
        One selected GMM-2 model per sample from :func:`select_best_cnv_score`.
    sample_col
        Sample identifier column shared by the two input tables.

    Returns
    -------
    pandas.DataFrame
        Per-cellbin metadata plus ``best_score_col``, ``best_cnv_score``,
        ``best_cnv_threshold``, ``cnv_high_probability``,
        ``cnv_assignment_confidence``, and ``cnv_auto_state``.

    Method
    ------
    The higher-mean GMM-2 component is defined as the CNV-high component. For
    each cellbin, labels use the fitted posterior directly:

    ``P(CNV-high | score) >= 0.5 -> Cancer cell``
    ``P(CNV-high | score) < 0.5  -> Normal Epithelial``

    ``best_cnv_threshold`` stores the GMM-2 weighted-component intersection
    for one-dimensional display. It is not estimated from KDE.
    """
    metadata_cols = [col for col in ['chipid', sample_col, 'cell_type', 'leiden', 'x', 'y', 'x_slide_mm', 'y_slide_mm']
                     if col in cnv_score_df.columns]
    classification_df = cnv_score_df[metadata_cols].copy()
    classification_df['best_score_col'] = pd.NA
    classification_df['best_cnv_score'] = np.nan
    classification_df['best_cnv_threshold'] = np.nan
    classification_df['cnv_high_probability'] = np.nan
    classification_df['cnv_assignment_confidence'] = np.nan
    classification_df['cnv_auto_state'] = pd.NA

    for _, row in best_df.iterrows():
        sample_mask = cnv_score_df[sample_col].astype(str).eq(str(row['sample']))
        score_col = str(row['score_col'])
        score_values = cnv_score_df.loc[sample_mask, score_col].to_numpy(dtype=float)
        high_probability = _gmm2_high_probability(score_values, low_mean=float(row['low_gmm_mean']),
                                                  low_std=float(row['low_gmm_std']),
                                                  low_weight=float(row['low_gmm_weight']),
                                                  high_mean=float(row['high_gmm_mean']),
                                                  high_std=float(row['high_gmm_std']),
                                                  high_weight=float(row['high_gmm_weight']))
        threshold = row.get('gmm2_threshold', row.get('gmm_boundary', np.nan))
        classification_df.loc[sample_mask, 'best_score_col'] = score_col
        classification_df.loc[sample_mask, 'best_cnv_score'] = score_values
        classification_df.loc[sample_mask, 'best_cnv_threshold'] = threshold
        classification_df.loc[sample_mask, 'cnv_high_probability'] = high_probability
        classification_df.loc[sample_mask, 'cnv_assignment_confidence'] = np.maximum(high_probability,
                                                                                     1 - high_probability)
        classification_df.loc[sample_mask, 'cnv_auto_state'] = np.where(high_probability >= 0.5, 'Cancer cell',
                                                                        'Normal Epithelial')
    return classification_df


def _mark_selected(ax: plt.Axes, selection_df: pd.DataFrame, sample_order: Sequence[str],
                   score_order: Sequence[str]) -> None:
    """Mark the selected score setting with an asterisk on a heatmap.

    Parameters
    ----------
    ax
        Heatmap axis to annotate.
    selection_df
        Selection-result table containing Boolean column ``selected``.
    sample_order
        Ordered sample labels used as heatmap rows.
    score_order
        Ordered CNV-score labels used as heatmap columns.

    Returns
    -------
    None
    """
    selected_df = selection_df.loc[selection_df['selected']]
    for _, row in selected_df.iterrows():
        y = list(sample_order).index(row['sample']) + 0.5
        x = list(score_order).index(row['score_col']) + 0.5
        ax.text(x, y, '*', ha='center', va='center', fontsize=16)


def plot_window_selection_heatmaps(selection_df: pd.DataFrame, outdir: Union[str, Path]) -> None:
    """Plot GMM-2 quality and GMM-2 threshold heatmaps as PNG and PDF.

    Parameters
    ----------
    selection_df
        Output of :func:`select_best_cnv_score`.
    outdir
        Pipeline output directory. Figures are written to ``{outdir}/figures``.

    Notes
    -----
    An asterisk marks the selected primary GMM-2 model for each sample. KDE
    valley values are intentionally not used in these primary selection plots.
    """
    figure_dir = Path(outdir) / 'figures'
    sample_order = selection_df['sample'].drop_duplicates().tolist()
    score_order = selection_df.sort_values(['window_size', 'step'])['score_col'].drop_duplicates().tolist()

    quality_df = selection_df.pivot(index='sample', columns='score_col', values='gmm2_quality_score').reindex(
        index=sample_order, columns=score_order)
    fig, ax = plt.subplots(figsize=(max(7, len(score_order) * 2), max(5, len(sample_order) * 0.55)))
    sns.heatmap(quality_df, annot=True, fmt='.2f', ax=ax)
    _mark_selected(ax, selection_df, sample_order, score_order)
    ax.set_xlabel('CNV score')
    ax.set_ylabel('Sample')
    ax.set_title('GMM-2 primary model quality; * selected')
    fig.tight_layout()
    _savefig(fig, figure_dir, 'CNV_GMM2_primary_window_quality_heatmap')

    threshold_col = 'gmm2_threshold' if 'gmm2_threshold' in selection_df.columns else 'gmm_boundary'
    threshold_df = selection_df.pivot(index='sample', columns='score_col', values=threshold_col).reindex(
        index=sample_order, columns=score_order)
    fig, ax = plt.subplots(figsize=(max(7, len(score_order) * 2), max(5, len(sample_order) * 0.55)))
    sns.heatmap(threshold_df, annot=True, fmt='.5f', ax=ax)
    _mark_selected(ax, selection_df, sample_order, score_order)
    ax.set_xlabel('CNV score')
    ax.set_ylabel('Sample')
    ax.set_title('GMM-2 weighted-component thresholds; * selected')
    fig.tight_layout()
    _savefig(fig, figure_dir, 'CNV_GMM2_primary_threshold_heatmap')


# Legacy KDE-primary plotting wrapper removed in GMM2-primary v7.

# Legacy KDE-primary plotting wrapper removed in GMM2-primary v7.

def add_cnv_results_to_adata(adata: ad.AnnData, cnv_score_df: pd.DataFrame, classification_df: pd.DataFrame) -> None:
    """Write CNV scores and GMM-2-primary labels back to ``adata.obs``.

    Added columns include every tested score, selected score/threshold, GMM-2
    CNV-high posterior probability, posterior confidence, and Cancer/Normal
    state for epithelial cellbins.
    """
    score_cols = [col for col in cnv_score_df.columns if col.startswith('cnv_score_w')]
    for score_col in score_cols:
        adata.obs[score_col] = np.nan
        adata.obs.loc[cnv_score_df.index, score_col] = cnv_score_df[score_col]
    for col in ['best_score_col', 'best_cnv_score', 'best_cnv_threshold', 'cnv_high_probability',
                'cnv_assignment_confidence', 'cnv_auto_state']:
        adata.obs[col] = pd.NA if col in ['best_score_col', 'cnv_auto_state'] else np.nan
        adata.obs.loc[classification_df.index, col] = classification_df[col]


def select_best_plot_window(selection_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Use the selected primary GMM-2 window as the diagnostic plot window.

    Parameters
    ----------
    selection_df
        Output of :func:`select_best_cnv_score`.

    Returns
    -------
    selection_df
        Input table augmented with ``plot_selected`` and
        ``plot_selection_reason``.
    best_plot_df
        Exactly one selected GMM-2 primary model per sample.

    Notes
    -----
    In the GMM-2-primary workflow every sample has a selected model, so no KDE
    or fallback selection is required. KDE remains visible inside the plot as a
    diagnostic curve and optional valley line.
    """
    selection_df = selection_df.copy()
    selection_df['plot_selected'] = selection_df['selected'].astype(bool)
    selection_df['plot_selection_reason'] = np.where(selection_df['selected'], 'selected_gmm2_primary', pd.NA)
    best_plot_df = selection_df.loc[selection_df['plot_selected']].copy().sort_values('sample').reset_index(drop=True)
    return selection_df, best_plot_df


def _make_density_curves(score_values: np.ndarray, result_row: pd.Series, n_grid: int = 500, bw_scale: float = 1.0,
                         peak_prominence_fraction: float = 0.03, peak_distance_fraction: float = 0.05) -> Dict[
    str, np.ndarray]:
    """Recalculate plot-ready KDE and GMM density curves for one fitted result.

    Parameters
    ----------
    score_values
        One-dimensional CNV-score vector for one sample and one window.
    result_row
        One row from ``selection_df`` or ``cutoff_df`` containing the two-GMM
        component parameters saved during cutoff fitting.
    n_grid
        Number of score positions used for density curves and figure rendering.
    bw_scale
        Same Scott-bandwidth multiplier used in the cutoff stage.
    peak_prominence_fraction
        Relative KDE peak-prominence threshold used only to display detected
        modes consistently with cutoff fitting.
    peak_distance_fraction
        Relative KDE peak-distance threshold used only to display detected
        modes consistently with cutoff fitting.

    Returns
    -------
    dict
        Arrays for the x-grid, KDE density, two weighted GMM component curves,
        GMM mixture curve, and KDE peak indices.
    """
    score_values = np.asarray(score_values, dtype=float)
    grid = np.linspace(score_values.min(), score_values.max(), n_grid)
    kde = gaussian_kde(score_values, bw_method=lambda kde_obj: kde_obj.scotts_factor() * bw_scale)
    kde_density = kde(grid)
    gmm_means = np.array([result_row['low_gmm_mean'], result_row['high_gmm_mean']], dtype=float)
    gmm_stds = np.array([result_row['low_gmm_std'], result_row['high_gmm_std']], dtype=float)
    gmm_weights = np.array([result_row['low_gmm_weight'], result_row['high_gmm_weight']], dtype=float)
    component_order = np.argsort(gmm_means)
    gmm_means = gmm_means[component_order]
    gmm_stds = gmm_stds[component_order]
    gmm_weights = gmm_weights[component_order]
    gmm_component_density = np.vstack(
        [weight * norm.pdf(grid, loc=mean, scale=std) for mean, std, weight in zip(gmm_means, gmm_stds, gmm_weights)])
    gmm_mixture_density = gmm_component_density.sum(axis=0)
    min_peak_distance = max(1, int(n_grid * peak_distance_fraction))
    peak_idx, _ = find_peaks(kde_density, prominence=kde_density.max() * peak_prominence_fraction,
                             distance=min_peak_distance)
    return {'grid': grid, 'kde_density': kde_density, 'gmm_component_density': gmm_component_density,
            'gmm_mixture_density': gmm_mixture_density, 'peak_idx': peak_idx}


def _plot_gmm_kde_axis(ax: plt.Axes, score_values: np.ndarray, result_row: pd.Series, n_grid: int = 500,
                       bw_scale: float = 1.0, peak_prominence_fraction: float = 0.03,
                       peak_distance_fraction: float = 0.05, bins: int = 100, show_legend: bool = True) -> None:
    """Draw histogram, KDE audit curve, GMM-2 components, and GMM-2 threshold.

    The dashed vertical line is the primary GMM-2 threshold. When available,
    the dotted vertical line is the KDE valley diagnostic and does not drive
    selection or labels.
    """
    curve_dict = _make_density_curves(score_values, result_row=result_row, n_grid=n_grid, bw_scale=bw_scale,
                                      peak_prominence_fraction=peak_prominence_fraction,
                                      peak_distance_fraction=peak_distance_fraction)
    grid = curve_dict['grid']
    kde_density = curve_dict['kde_density']
    gmm_component_density = curve_dict['gmm_component_density']
    gmm_mixture_density = curve_dict['gmm_mixture_density']
    peak_idx = curve_dict['peak_idx']
    sns.histplot(score_values, bins=bins, stat='density', element='step', fill=False, ax=ax, label='Histogram')
    ax.plot(grid, kde_density, label='KDE')
    ax.plot(grid, gmm_component_density[0], label='GMM2 low component')
    ax.plot(grid, gmm_component_density[1], label='GMM2 high component')
    ax.plot(grid, gmm_mixture_density, label='GMM2 mixture')
    if len(peak_idx) > 0:
        ax.scatter(grid[peak_idx], kde_density[peak_idx], s=22, label='KDE modes')

    gmm2_threshold = result_row.get('gmm2_threshold', result_row.get('gmm_boundary', np.nan))
    if pd.notna(gmm2_threshold):
        ax.axvline(float(gmm2_threshold), linestyle='--', label=f'GMM2 cutoff = {float(gmm2_threshold):.6f}')
    kde_valley_threshold = result_row.get('kde_valley_threshold', np.nan)
    if pd.notna(kde_valley_threshold):
        ax.axvline(float(kde_valley_threshold), linestyle=':', label=f'KDE valley = {float(kde_valley_threshold):.6f}')
    ax.set_xlabel(str(result_row['score_col']))
    ax.set_ylabel('Density')
    if show_legend:
        ax.legend(fontsize=8)


def _plot_title(result_row: pd.Series, include_sample: bool = True) -> str:
    """Build a compact GMM-2-primary diagnostic title for one density plot."""
    gmm2_threshold = result_row.get('gmm2_threshold', result_row.get('gmm_boundary', np.nan))
    gmm2_text = f'{float(gmm2_threshold):.6f}' if pd.notna(gmm2_threshold) else 'NA'
    kde_threshold = result_row.get('kde_valley_threshold', np.nan)
    kde_text = f'{float(kde_threshold):.6f}' if pd.notna(kde_threshold) else 'NA'
    selected = 'yes' if bool(result_row.get('selected', False)) else 'no'
    quality_pass = 'yes' if bool(result_row.get('gmm2_quality_pass', False)) else 'no'
    kde_status = result_row.get('status', 'NA')
    title_lines = [
        f'{result_row["sample"]} | {result_row["score_col"]}' if include_sample else str(result_row['score_col']),
        f'GMM2 primary; selected={selected}; quality_pass={quality_pass}; KDE={kde_status}',
        f'ΔBIC(1-2)={float(result_row["bic_delta_1minus2"]):.1f}; '
        f'ΔBIC(2-3)={float(result_row["bic_delta_2minus3"]):.1f}; '
        f'Ashman D={float(result_row["ashman_d"]):.2f}', f'GMM2 cutoff={gmm2_text}; KDE valley={kde_text}']
    return '\n'.join(title_lines)


def plot_all_windows_for_each_sample(cnv_score_df: pd.DataFrame, selection_df: pd.DataFrame, sample_col: str = 'sample',
                                     outdir: Union[str, Path] = 'analysis_result/infercnv', n_grid: int = 500,
                                     bw_scale: float = 1.0, peak_prominence_fraction: float = 0.03,
                                     peak_distance_fraction: float = 0.05, bins: int = 100) -> None:
    """Save one KDE/GMM diagnostic figure for every ``sample × CNV score`` pair.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin CNV-score table.
    selection_df
        Selection table containing every score setting, fitted GMM parameters,
        KDE diagnostics, strict eligibility, and ``plot_selected`` flags.
    sample_col
        Sample identifier column in ``cnv_score_df``.
    outdir
        Pipeline output directory. Files are written under
        ``{outdir}/figures/all_windows`` as both PNG and PDF.
    n_grid, bw_scale, peak_prominence_fraction, peak_distance_fraction
        Plotting settings matched to the cutoff stage.
    bins
        Number of histogram bins.

    Returns
    -------
    None

    Notes
    -----
    This function draws every tested window regardless of strict eligibility,
    KDE status, or whether a finite valley cutoff exists. It is therefore the
    complete visual audit set for window selection.
    """
    figure_dir = Path(outdir) / 'figures' / 'all_windows'
    ordered_df = selection_df.sort_values(['sample', 'window_size', 'step'])
    total = len(ordered_df)
    for plot_index, (_, result_row) in enumerate(ordered_df.iterrows(), start=1):
        sample = str(result_row['sample'])
        score_col = str(result_row['score_col'])
        score_values = cnv_score_df.loc[cnv_score_df[sample_col].astype(str).eq(sample), score_col].dropna().to_numpy(
            dtype=float)
        fig, ax = plt.subplots(figsize=(7, 5))
        _plot_gmm_kde_axis(ax, score_values=score_values, result_row=result_row, n_grid=n_grid, bw_scale=bw_scale,
                           peak_prominence_fraction=peak_prominence_fraction,
                           peak_distance_fraction=peak_distance_fraction, bins=bins, show_legend=True)
        ax.set_title(_plot_title(result_row, include_sample=True), fontsize=10)
        fig.tight_layout()
        _savefig(fig, figure_dir, f'{_safe_name(sample)}_{score_col}_KDE_GMM')
        print(f'[{plot_index}/{total}] plotted {sample} {score_col}', flush=True)


def plot_panel_for_each_sample(cnv_score_df: pd.DataFrame, selection_df: pd.DataFrame, sample_col: str = 'sample',
                               outdir: Union[str, Path] = 'analysis_result/infercnv', n_grid: int = 500,
                               bw_scale: float = 1.0, peak_prominence_fraction: float = 0.03,
                               peak_distance_fraction: float = 0.05, bins: int = 100, ncols: int = 2) -> None:
    """Save a multi-window density panel for every sample.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin CNV-score table.
    selection_df
        Selection table containing every score setting and fitted diagnostics.
    sample_col
        Sample identifier column in ``cnv_score_df``.
    outdir
        Pipeline output directory. Panel figures are saved under
        ``{outdir}/figures/panel_by_sample`` as PNG and PDF.
    n_grid, bw_scale, peak_prominence_fraction, peak_distance_fraction
        Plotting settings matched to the cutoff stage.
    bins
        Number of histogram bins per panel.
    ncols
        Number of panel columns. With the default four tested windows, ``2``
        creates a 2 × 2 layout.

    Returns
    -------
    None

    Notes
    -----
    Unlike the all-window figures, this view makes within-sample changes in
    score distribution, KDE peak topology, and fitted GMM components directly
    comparable across window sizes.
    """
    figure_dir = Path(outdir) / 'figures' / 'panel_by_sample'
    sample_groups = list(selection_df.groupby('sample', sort=False))
    for sample_index, (sample, sample_result_df) in enumerate(sample_groups, start=1):
        sample_result_df = sample_result_df.sort_values(['window_size', 'step'])
        n_panels = len(sample_result_df)
        nrows = int(np.ceil(n_panels / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows))
        axes = np.atleast_1d(axes).ravel()
        sample_score_df = cnv_score_df.loc[cnv_score_df[sample_col].astype(str).eq(str(sample))]
        for ax, (_, result_row) in zip(axes, sample_result_df.iterrows()):
            score_col = str(result_row['score_col'])
            score_values = sample_score_df[score_col].dropna().to_numpy(dtype=float)
            _plot_gmm_kde_axis(ax, score_values=score_values, result_row=result_row, n_grid=n_grid, bw_scale=bw_scale,
                               peak_prominence_fraction=peak_prominence_fraction,
                               peak_distance_fraction=peak_distance_fraction, bins=bins, show_legend=False)
            ax.set_title(_plot_title(result_row, include_sample=False), fontsize=9)
        for ax in axes[n_panels:]:
            ax.axis('off')
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=3, fontsize=9, bbox_to_anchor=(0.5, 1.01))
        fig.suptitle(str(sample), y=1.08, fontsize=15)
        fig.tight_layout()
        _savefig(fig, figure_dir, f'{_safe_name(sample)}_all_windows_panel')
        print(f'[{sample_index}/{len(sample_groups)}] plotted panel {sample}', flush=True)


def plot_best_window_for_each_sample(cnv_score_df: pd.DataFrame, best_plot_df: pd.DataFrame, sample_col: str = 'sample',
                                     outdir: Union[str, Path] = 'analysis_result/infercnv', n_grid: int = 500,
                                     bw_scale: float = 1.0, peak_prominence_fraction: float = 0.03,
                                     peak_distance_fraction: float = 0.05, bins: int = 100) -> None:
    """Save one fallback-aware best-window density figure for every sample.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin CNV-score table.
    best_plot_df
        One-row-per-sample output of :func:`select_best_plot_window`.
    sample_col
        Sample identifier column in ``cnv_score_df``.
    outdir
        Pipeline output directory. Figures are saved under
        ``{outdir}/figures/best_windows`` as PNG and PDF.
    n_grid, bw_scale, peak_prominence_fraction, peak_distance_fraction
        Plotting settings matched to the cutoff stage.
    bins
        Number of histogram bins.

    Returns
    -------
    None

    Notes
    -----
    This function always produces one representative plot per sample. In a
    strictly eligible sample it matches the classification window; otherwise it
    visualizes the highest-ranked fallback window without assigning cell labels.
    """
    figure_dir = Path(outdir) / 'figures' / 'best_windows'
    ordered_df = best_plot_df.sort_values('sample')
    for plot_index, (_, result_row) in enumerate(ordered_df.iterrows(), start=1):
        sample = str(result_row['sample'])
        score_col = str(result_row['score_col'])
        score_values = cnv_score_df.loc[cnv_score_df[sample_col].astype(str).eq(sample), score_col].dropna().to_numpy(
            dtype=float)
        fig, ax = plt.subplots(figsize=(7, 5))
        _plot_gmm_kde_axis(ax, score_values=score_values, result_row=result_row, n_grid=n_grid, bw_scale=bw_scale,
                           peak_prominence_fraction=peak_prominence_fraction,
                           peak_distance_fraction=peak_distance_fraction, bins=bins, show_legend=True)
        ax.set_title(f'{_plot_title(result_row, include_sample=True)}\n'
                     f'plot selection={result_row["plot_selection_reason"]}', fontsize=10)
        fig.tight_layout()
        _savefig(fig, figure_dir, f'{_safe_name(sample)}_{score_col}_best_window_KDE_GMM')
        print(f'[{plot_index}/{len(ordered_df)}] plotted best window {sample} {score_col}', flush=True)


def infercnv_cancer_normal_recipe(adata: ad.AnnData, gtf_file: str,
                                  window_configs: Sequence[Tuple[int, int]] = ((100, 10), (300, 30), (500, 50),
                                                                               (1000, 100)),
                                  cell_type_col: str = 'cell_type', sample_col: str = 'sample',
                                  epi_name: str = 'Epithelial', outdir: Union[str, Path] = 'analysis_result/infercnv',
                                  feature_id_col: str = 'transcript', gtf_feature_attr='gene_name',
                                  chipid_col: Optional[str] = 'chipid', unassigned_name: str = 'Unassigned',
                                  layer: str = 'log1p', lfc_clip: float = 3, dynamic_threshold: float = 1.5,
                                  exclude_chromosomes: Sequence[str] = ('chrX', 'chrY'), chunksize: int = 1000,
                                  n_jobs: int = 8, cutoff_n_jobs: int = 1, n_grid: int = 5000, bw_scale: float = 1.0,
                                  peak_prominence_fraction: float = 0.03, peak_distance_fraction: float = 0.05,
                                  gmm_n_init: int = 50, gmm_max_iter: int = 1000, gmm_tol: float = 1e-6,
                                  gmm_reg_covar: float = 1e-8, random_state: int = 0, min_delta_bic_1vs2: float = 10,
                                  max_delta_bic_2vs3: float = 10, min_ashman_d: float = 2,
                                  min_component_weight: float = 0.05, min_valley_depth: float = 0.10,
                                  plot_all_windows: bool = True, plot_panels: bool = True,
                                  plot_best_windows: bool = True, plot_bins: int = 100,
                                  overwrite: OverwriteType = False, timing: bool = False) -> Dict[str, object]:
    """Run inferCNV score calculation, GMM-2-primary selection/classification, and diagnostic plotting.

    Parameters
    ----------
    adata
        Cellbin-level AnnData object. Required fields are
        ``adata.var[feature_id_col]``, ``adata.obs[cell_type_col]``, and
        ``adata.obs[sample_col]``. The object is modified in place with genomic
        coordinates in ``var`` and final GMM-2-primary CNV labels in ``obs``.
    gtf_file
        GTF file used to map genomic coordinates by GTF ``gene_id``.
    window_configs
        Sequence of unique ``(window_size, step)`` settings for inferCNV.
    cell_type_col, sample_col, epi_name, feature_id_col, chipid_col,
    unassigned_name, layer, lfc_clip, dynamic_threshold, exclude_chromosomes,
    chunksize, n_jobs
        Upstream inferCNV score-calculation settings. ``n_jobs`` belongs only
        to infercnvpy and is used only when ``cnv_score`` is recalculated.
    cutoff_n_jobs
        Number of processes for downstream independent
        ``sample × CNV-score`` GMM/KDE fitting tasks. Each worker is limited to
        one BLAS/OpenMP thread by :func:`_fit_gmm_kde_task`.
    n_grid, bw_scale, peak_prominence_fraction, peak_distance_fraction
        Shared KDE and peak-detection settings used consistently for cutoff
        fitting and all diagnostic figures.
    gmm_n_init, gmm_max_iter, gmm_tol, gmm_reg_covar, random_state
        GMM fitting controls. GMM-1 always uses one initialization because a
        one-component fit has no multi-start local-optimum issue; GMM-2 and
        GMM-3 use ``gmm_n_init`` starts.
    min_delta_bic_1vs2, max_delta_bic_2vs3, min_ashman_d,
    min_component_weight, min_valley_depth
        Diagnostic thresholds recorded by :func:`select_best_cnv_score`; GMM-2 selection and labels remain available for every sample.
    plot_all_windows
        Save every ``sample × score`` KDE/GMM figure under
        ``figures/all_windows``. This is the complete visual audit set.
    plot_panels
        Save one multi-window panel per sample under ``figures/panel_by_sample``.
    plot_best_windows
        Save the selected primary GMM-2 window per sample under
        ``figures/best_windows``. Every sample has one selected window.
    plot_bins
        Number of histogram bins in all density figures.
    overwrite
        Cache policy. ``False`` reuses saved results. ``True`` recalculates all
        analysis and plotting stages. A mapping supports individual steps, e.g.
        ``{'cnv_score': False, 'cutoff': True, 'selection': True, 'plot': True}``.
        Recalculation propagates downstream dependencies automatically.
    timing
        Print GMM/KDE task completion timing. Timing columns are always saved
        in the cutoff result table when cutoff fitting is run.

    Returns
    -------
    dict
        Contains modified ``adata``, ``cnv_score_df``, ``cnv_score_summary_df``,
        ``cutoff_df``, ``selection_df``, strict ``best_df``, GMM-2-primary
        ``best_plot_df``, and ``classification_df``.

    Output files
    ------------
    ``outdir`` contains score and cutoff CSV files, GMM-2 primary window
    selections, GMM-2 epithelial labels, timing CSV, heatmaps, every
    ``sample × score`` density figure, one panel per sample, and one best-window
    diagnostic figure per sample. All figures are saved as PNG and PDF through
    :func:`_savefig`.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    coordinate_path = outdir / 'feature_gene_coordinates.csv'
    cnv_score_path = outdir / 'epithelial_cell_CNV_scores_all_window_sizes.csv'
    cnv_summary_path = outdir / 'epithelial_CNV_score_summary_by_sample.csv'
    cutoff_path = outdir / 'CNV_GMM_KDE_thresholds_all_samples_windows.csv'
    timing_path = outdir / 'CNV_GMM_KDE_timing_all_samples_windows.csv'
    selection_path = outdir / 'CNV_score_window_selection_all_samples.csv'
    best_path = outdir / 'CNV_best_score_threshold_by_sample.csv'
    best_plot_path = outdir / 'CNV_best_plot_window_by_sample.csv'
    classification_path = outdir / 'epithelial_best_CNV_classification.csv'
    plot_manifest_path = outdir / 'CNV_plot_manifest.csv'
    plot_complete_path = outdir / 'figures' / '.all_plots_complete'
    selection_schema_ready = False
    if selection_path.exists():
        selection_schema_ready = 'gmm2_quality_score' in pd.read_csv(selection_path, nrows=1).columns

    cnv_feature_mask = add_gene_coordinates(adata, gtf_file=gtf_file, feature_id_col=feature_id_col,
                                            gtf_feature_attr=gtf_feature_attr)
    coordinate_cols = [col for col in [feature_id_col, 'real_gene_name', 'chromosome', 'start', 'end'] if
                       col in adata.var.columns]
    if _should_overwrite(overwrite, 'cnv_score') or not coordinate_path.exists():
        adata.var[coordinate_cols].to_csv(coordinate_path)
    print(f'总 feature 数：{adata.n_vars}')
    print(f'成功匹配 gene_id 坐标的 feature 数：{cnv_feature_mask.sum()}')
    print(f'未匹配 feature 数：{(~cnv_feature_mask).sum()}')

    plot_config_text = (f'version={__version__}\nall_windows={plot_all_windows}\npanels={plot_panels}\n'
                        f'best_windows={plot_best_windows}\nplot_bins={plot_bins}\nn_grid={n_grid}\n'
                        f'bw_scale={bw_scale}\npeak_prominence_fraction={peak_prominence_fraction}\n'
                        f'peak_distance_fraction={peak_distance_fraction}\nselection_mode=GMM2_primary\n')
    current_plot_config = plot_complete_path.read_text(encoding='utf-8') if plot_complete_path.exists() else ''
    run_cnv_score = _should_overwrite(overwrite, 'cnv_score') or not cnv_score_path.exists()
    run_cutoff = _should_overwrite(overwrite, 'cutoff') or run_cnv_score or not cutoff_path.exists()
    run_selection = (_should_overwrite(overwrite,
                                       'selection') or run_cutoff or not selection_path.exists() or not best_path.exists() or not best_plot_path.exists() or not selection_schema_ready)
    run_plot = _should_overwrite(overwrite, 'plot') or run_selection or current_plot_config != plot_config_text

    if run_cnv_score:
        cnv_score_df, cnv_score_summary_df = calculate_cnv_scores(adata, cnv_feature_mask=cnv_feature_mask,
                                                                  window_configs=window_configs,
                                                                  cell_type_col=cell_type_col, sample_col=sample_col,
                                                                  epi_name=epi_name, unassigned_name=unassigned_name,
                                                                  chipid_col=chipid_col, layer=layer, lfc_clip=lfc_clip,
                                                                  dynamic_threshold=dynamic_threshold,
                                                                  exclude_chromosomes=exclude_chromosomes,
                                                                  chunksize=chunksize, n_jobs=n_jobs)
        cnv_score_df.to_csv(cnv_score_path, index_label='obs_name')
        cnv_score_summary_df.to_csv(cnv_summary_path)
    else:
        cnv_score_df = pd.read_csv(cnv_score_path, index_col=0)
        summary_index_col = [0, 1] if chipid_col is not None and chipid_col in cnv_score_df.columns else 0
        cnv_score_summary_df = pd.read_csv(cnv_summary_path, index_col=summary_index_col)

    if run_cutoff:
        cutoff_df = calculate_gmm_kde_cutoffs(cnv_score_df, window_configs=window_configs, sample_col=sample_col,
                                              n_grid=n_grid, bw_scale=bw_scale,
                                              peak_prominence_fraction=peak_prominence_fraction,
                                              peak_distance_fraction=peak_distance_fraction, gmm_n_init=gmm_n_init,
                                              gmm_max_iter=gmm_max_iter, gmm_tol=gmm_tol, gmm_reg_covar=gmm_reg_covar,
                                              random_state=random_state, cutoff_n_jobs=cutoff_n_jobs, timing=timing)
        cutoff_df.to_csv(cutoff_path, index=False)
        timing_cols = ['sample', 'score_col', 'window_size', 'step', 'n_cells'] + [col for col in cutoff_df.columns if
                                                                                   col.startswith('time_')]
        cutoff_df[timing_cols].to_csv(timing_path, index=False)
    else:
        cutoff_df = pd.read_csv(cutoff_path)
        if not timing_path.exists() and any(col.startswith('time_') for col in cutoff_df.columns):
            timing_cols = ['sample', 'score_col', 'window_size', 'step', 'n_cells'] + [col for col in cutoff_df.columns
                                                                                       if col.startswith('time_')]
            cutoff_df[timing_cols].to_csv(timing_path, index=False)

    if 'gmm2_threshold' not in cutoff_df.columns:
        cutoff_df['gmm2_threshold'] = cutoff_df['gmm_boundary']
        cutoff_df['gmm2_primary_model'] = True
        cutoff_df.to_csv(cutoff_path, index=False)

    if run_selection:
        selection_df, best_df = select_best_cnv_score(cutoff_df, min_delta_bic_1vs2=min_delta_bic_1vs2,
                                                      max_delta_bic_2vs3=max_delta_bic_2vs3, min_ashman_d=min_ashman_d,
                                                      min_component_weight=min_component_weight,
                                                      min_valley_depth=min_valley_depth)
        selection_df, best_plot_df = select_best_plot_window(selection_df)
        selection_df.to_csv(selection_path, index=False)
        best_df.to_csv(best_path, index=False)
        best_plot_df.to_csv(best_plot_path, index=False)
        classification_df = assign_best_cnv_labels(cnv_score_df, best_df, sample_col=sample_col)
        classification_df.to_csv(classification_path, index_label='obs_name')
    else:
        selection_df = pd.read_csv(selection_path)
        best_df = pd.read_csv(best_path)
        best_plot_df = pd.read_csv(best_plot_path)
        classification_df = pd.read_csv(classification_path, index_col=0)

    if run_plot:
        plot_window_selection_heatmaps(selection_df, outdir=outdir)
        if plot_all_windows:
            plot_all_windows_for_each_sample(cnv_score_df, selection_df, sample_col=sample_col, outdir=outdir,
                                             n_grid=n_grid, bw_scale=bw_scale,
                                             peak_prominence_fraction=peak_prominence_fraction,
                                             peak_distance_fraction=peak_distance_fraction, bins=plot_bins)
        if plot_panels:
            plot_panel_for_each_sample(cnv_score_df, selection_df, sample_col=sample_col, outdir=outdir, n_grid=n_grid,
                                       bw_scale=bw_scale, peak_prominence_fraction=peak_prominence_fraction,
                                       peak_distance_fraction=peak_distance_fraction, bins=plot_bins)
        if plot_best_windows:
            plot_best_window_for_each_sample(cnv_score_df, best_plot_df, sample_col=sample_col, outdir=outdir,
                                             n_grid=n_grid, bw_scale=bw_scale,
                                             peak_prominence_fraction=peak_prominence_fraction,
                                             peak_distance_fraction=peak_distance_fraction, bins=plot_bins)
        plot_manifest_cols = [col for col in
                              ['sample', 'score_col', 'window_size', 'step', 'status', 'selected', 'plot_selected',
                               'plot_selection_reason', 'gmm2_threshold', 'gmm2_quality_score', 'gmm2_quality_pass',
                               'gmm2_vs_gmm1_supported', 'gmm3_vs_gmm2_supported', 'kde_valley_threshold', 'ashman_d',
                               'component_balance', 'posterior_high_conf_fraction'] if col in selection_df.columns]
        selection_df[plot_manifest_cols].to_csv(plot_manifest_path, index=False)
        plot_complete_path.parent.mkdir(parents=True, exist_ok=True)
        plot_complete_path.write_text(plot_config_text, encoding='utf-8')

    add_cnv_results_to_adata(adata, cnv_score_df, classification_df)
    return {'adata': adata, 'cnv_score_df': cnv_score_df, 'cnv_score_summary_df': cnv_score_summary_df,
            'cutoff_df': cutoff_df, 'selection_df': selection_df, 'best_df': best_df, 'best_plot_df': best_plot_df,
            'classification_df': classification_df}
