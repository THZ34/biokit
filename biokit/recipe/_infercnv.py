# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

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
   to the raw epithelial ``cnv_score`` distribution. GMM is used to quantify
   evidence for a two-component distribution and the quality of its separation;
   it is not used as the final cutoff.
5. Estimate the empirical score density with KDE, detect two dominant KDE
   modes, and use the global density minimum between these modes as the final
   CNV cutoff. This avoids treating the GMM component intersection or the
   minimum of the overall GMM mixture density as the threshold. Either can be
   displaced from the visibly observed valley when the fitted Gaussian
   components are asymmetric or overlap substantially.
6. Select the best CNV-score setting separately for every sample, assign
   ``Cancer cell`` / ``Normal Epithelial`` labels, and retain
   ``Epithelial uncertain`` when no setting passes the specified quality rules.

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

from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

import anndata as ad
import gtfparse
import infercnvpy as cnv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, norm
from sklearn.mixture import GaussianMixture

__version__ = '2026-07-07-doc-savefig-v2'

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
                         chromosome_prefix: str = 'chr') -> np.ndarray:
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
    gene_pos_df = gtf_df.loc[gtf_df['feature'].eq('gene'), ['seqname', 'start', 'end', 'gene_id']].dropna(
        subset=['gene_id']).drop_duplicates('gene_id').rename(columns={'seqname': 'chromosome'}).set_index('gene_id')
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
                       gmm_tol: float = 1e-6, gmm_reg_covar: float = 1e-8, random_state: int = 0) -> Dict[str, object]:
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

    Returns
    -------
    dict
        One row of diagnostics, including GMM BIC values, component means,
        standard deviations, weights, Ashman's D, posterior-confidence rates,
        KDE peak locations, KDE valley threshold, and selection-relevant flags.

    Method
    ------
    **GMM diagnostics**

    - Fit 1-, 2-, and 3-component Gaussian mixtures to the raw ``cnv_score``
      vector. ``BIC(1) - BIC(2) > 0`` indicates that two components explain the
      distribution better than one; ``BIC(2) - BIC(3)`` checks whether a third
      component materially improves fit.
    - Order the two GMM components by mean, then report their means, standard
      deviations, weights, Ashman's D, posterior confidence, and the weighted
      component-density intersection.
    - These GMM values are diagnostics for evidence and separation quality, not
      the final cutoff.

    **Final cutoff**

    - Estimate the empirical density using Gaussian KDE on the full score range.
    - Detect KDE modes using prominence and minimum-distance rules.
    - When exactly two KDE modes are detected, find the global KDE minimum in
      the open interval between the low and high modes.
    - Store this minimum as ``kde_valley_threshold``. It is the actual score
      cutoff used later for ``Cancer cell`` / ``Normal Epithelial`` assignment.

    The KDE-valley rule follows the bimodal malignancy-score threshold concept
    used by Liu et al., but does not reproduce ``scCancer::getBimodalThres``
    internals. The deliberate use of the empirical valley rather than the GMM
    intersection addresses cases where the fitted GMM mixture remains formally
    unimodal although the observed KDE has two visible modes.

    References
    ----------
    - Liu T, et al. Nat Commun. 2022. doi:10.1038/s41467-022-34581-2.
    - Schwarz G. Ann Stat. 1978;6:461-464.
    - Ashman KM, Bird CM, Zepf SE. Astron J. 1994;108:2348-2351.
    """
    score_values = np.asarray(score_values, dtype=float)
    score_matrix = score_values[:, None]
    gmm_1 = GaussianMixture(n_components=1, covariance_type='full', n_init=gmm_n_init, max_iter=gmm_max_iter,
                            tol=gmm_tol, reg_covar=gmm_reg_covar, random_state=random_state)
    gmm_2 = GaussianMixture(n_components=2, covariance_type='full', n_init=gmm_n_init, max_iter=gmm_max_iter,
                            tol=gmm_tol, reg_covar=gmm_reg_covar, random_state=random_state)
    gmm_3 = GaussianMixture(n_components=3, covariance_type='full', n_init=gmm_n_init, max_iter=gmm_max_iter,
                            tol=gmm_tol, reg_covar=gmm_reg_covar, random_state=random_state)
    gmm_1.fit(score_matrix)
    gmm_2.fit(score_matrix)
    gmm_3.fit(score_matrix)
    component_order = np.argsort(gmm_2.means_.ravel())
    gmm_means = gmm_2.means_.ravel()[component_order]
    gmm_stds = np.sqrt(gmm_2.covariances_.reshape(-1)[component_order])
    gmm_weights = gmm_2.weights_.ravel()[component_order]
    posterior = gmm_2.predict_proba(score_matrix)[:, component_order]
    posterior_confidence = posterior.max(axis=1)
    grid = np.linspace(score_values.min(), score_values.max(), n_grid)
    gmm_component_density = np.vstack(
        [weight * norm.pdf(grid, loc=mean, scale=std) for mean, std, weight in zip(gmm_means, gmm_stds, gmm_weights)])
    gmm_mixture_density = gmm_component_density.sum(axis=0)
    kde = gaussian_kde(score_values, bw_method=lambda kde_obj: kde_obj.scotts_factor() * bw_scale)
    kde_density = kde(grid)
    min_peak_distance = max(1, int(n_grid * peak_distance_fraction))
    peak_idx, _ = find_peaks(kde_density, prominence=kde_density.max() * peak_prominence_fraction,
                             distance=min_peak_distance)
    gmm_mode_idx, _ = find_peaks(gmm_mixture_density, prominence=gmm_mixture_density.max() * peak_prominence_fraction,
                                 distance=min_peak_distance)
    result = {'sample': sample, 'score_col': score_col, 'window_size': window_size, 'step': step,
        'n_cells': len(score_values), 'score_min': score_values.min(), 'score_max': score_values.max(),
        'bic_1component': gmm_1.bic(score_matrix), 'bic_2component': gmm_2.bic(score_matrix),
        'bic_3component': gmm_3.bic(score_matrix),
        'bic_delta_1minus2': gmm_1.bic(score_matrix) - gmm_2.bic(score_matrix),
        'bic_delta_2minus3': gmm_2.bic(score_matrix) - gmm_3.bic(score_matrix), 'gmm_1_converged': gmm_1.converged_,
        'gmm_2_converged': gmm_2.converged_, 'gmm_3_converged': gmm_3.converged_, 'low_gmm_mean': gmm_means[0],
        'high_gmm_mean': gmm_means[1], 'low_gmm_std': gmm_stds[0], 'high_gmm_std': gmm_stds[1],
        'low_gmm_weight': gmm_weights[0], 'high_gmm_weight': gmm_weights[1],
        'ashman_d': np.sqrt(2) * abs(gmm_means[1] - gmm_means[0]) / np.sqrt(gmm_stds[0] ** 2 + gmm_stds[1] ** 2),
        'posterior_high_conf_fraction': (posterior_confidence >= 0.9).mean(),
        'posterior_uncertain_fraction': (posterior_confidence < 0.8).mean(),
        'gmm_boundary': _gmm_boundary(grid, gmm_component_density, gmm_means), 'n_kde_peaks': len(peak_idx),
        'n_gmm_mixture_modes': len(gmm_mode_idx)}
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
        result['gmm_kde_boundary_offset_ratio'] = abs(result['gmm_boundary'] - kde_valley_threshold) / (
                high_kde_peak - low_kde_peak) if np.isfinite(result['gmm_boundary']) else np.nan
    else:
        result['status'] = 'unimodal' if len(peak_idx) < 2 else 'multimodal'
        result['low_kde_peak'] = np.nan
        result['high_kde_peak'] = np.nan
        result['kde_valley_threshold'] = np.nan
        result['valley_depth'] = np.nan
        result['gmm_kde_boundary_offset_ratio'] = np.nan
    return result


def calculate_gmm_kde_cutoffs(cnv_score_df: pd.DataFrame, window_configs: Sequence[Tuple[int, int]],
                              sample_col: str = 'sample', n_grid: int = 5000, bw_scale: float = 1.0,
                              peak_prominence_fraction: float = 0.03, peak_distance_fraction: float = 0.05,
                              gmm_n_init: int = 50, gmm_max_iter: int = 1000, gmm_tol: float = 1e-6,
                              gmm_reg_covar: float = 1e-8, random_state: int = 0) -> pd.DataFrame:
    """Run the GMM/KDE cutoff procedure for every sample and score setting.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin table returned by :func:`calculate_cnv_scores`.
    window_configs
        Same ordered ``(window_size, step)`` settings used to create
        ``cnv_score_df``.
    sample_col
        Column identifying samples in ``cnv_score_df``.
    n_grid, bw_scale, peak_prominence_fraction, peak_distance_fraction
        KDE-grid and peak-detection parameters passed unchanged to
        :func:`fit_gmm_kde_cutoff`.
    gmm_n_init, gmm_max_iter, gmm_tol, gmm_reg_covar, random_state
        GMM fitting parameters passed unchanged to :func:`fit_gmm_kde_cutoff`.

    Returns
    -------
    pandas.DataFrame
        One diagnostic row for every ``sample × window_config`` combination.

    Notes
    -----
    The returned ``kde_valley_threshold`` is only non-missing when exactly two
    KDE modes are detected under the common density and peak-detection settings.
    """
    result_rows = []
    for window_size, step in window_configs:
        score_col = _score_col(window_size)
        for sample, sample_df in cnv_score_df.groupby(sample_col, sort=False):
            result = fit_gmm_kde_cutoff(sample_df[score_col].dropna().to_numpy(), sample=str(sample),
                                        score_col=score_col, window_size=window_size, step=step, n_grid=n_grid,
                                        bw_scale=bw_scale, peak_prominence_fraction=peak_prominence_fraction,
                                        peak_distance_fraction=peak_distance_fraction, gmm_n_init=gmm_n_init,
                                        gmm_max_iter=gmm_max_iter, gmm_tol=gmm_tol, gmm_reg_covar=gmm_reg_covar,
                                        random_state=random_state)
            result_rows.append(result)
            print(f'{sample} {score_col}: {result["status"]}')
    return pd.DataFrame(result_rows).sort_values(['sample', 'window_size']).reset_index(drop=True)


def select_best_cnv_score(cutoff_df: pd.DataFrame, min_delta_bic_1vs2: float = 10, max_delta_bic_2vs3: float = 10,
                          min_ashman_d: float = 2, min_component_weight: float = 0.05,
                          min_valley_depth: float = 0.10) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Select one CNV-score window and one KDE-valley threshold per sample.

    Parameters
    ----------
    cutoff_df
        Output of :func:`calculate_gmm_kde_cutoffs` containing one diagnostic
        row per sample and score setting.
    min_delta_bic_1vs2
        Minimum required ``BIC(1) - BIC(2)``. Values above zero favour two GMM
        components; ``10`` is a conservative evidence threshold.
    max_delta_bic_2vs3
        Maximum allowed ``BIC(2) - BIC(3)``. A value above this limit indicates
        that three components improve BIC substantially and a simple binary
        partition may be insufficient.
    min_ashman_d
        Minimum Ashman's D for the two fitted GMM components. ``D >= 2`` is
        conventionally interpreted as clearly separated components.
    min_component_weight
        Minimum required GMM weight for both low and high components. This
        prevents a tiny tail component from being selected as a Cancer group.
    min_valley_depth
        Minimum required relative KDE valley depth. The depth is defined as
        ``1 - density_valley / min(density_low_peak, density_high_peak)``.

    Returns
    -------
    selection_df
        ``cutoff_df`` augmented with quality components, Boolean criteria,
        eligibility, rank, and the selected-window flag.
    best_df
        One selected row per sample. Samples with no eligible score setting are
        retained with ``selection_status='no_valid_window'``.

    Selection logic
    ---------------
    A score setting is eligible only when:

    1. KDE detects exactly two modes;
    2. the two-component GMM is supported over one component;
    3. both fitted components meet the minimum weight;
    4. Ashman's D meets the requested separation threshold; and
    5. the KDE valley is sufficiently deep.

    Eligible settings are ranked by whether two components are sufficient versus
    three, then by a composite quality score, valley depth, Ashman's D, and
    posterior confidence. The final classification threshold remains the KDE
    valley of the highest-ranked eligible setting.
    """
    selection_df = cutoff_df.copy()
    selection_df['component_balance'] = 2 * np.minimum(selection_df['low_gmm_weight'], selection_df['high_gmm_weight'])
    selection_df['separation_score'] = np.minimum(selection_df['ashman_d'] / 3, 1)
    selection_df['quality_score'] = selection_df['valley_depth'] * selection_df['separation_score'] * selection_df[
        'posterior_high_conf_fraction'] * selection_df['component_balance']
    selection_df['two_component_better'] = selection_df['bic_delta_1minus2'] > min_delta_bic_1vs2
    selection_df['two_components_sufficient'] = selection_df['bic_delta_2minus3'] <= max_delta_bic_2vs3
    selection_df['weight_pass'] = np.minimum(selection_df['low_gmm_weight'],
                                             selection_df['high_gmm_weight']) >= min_component_weight
    selection_df['separation_pass'] = selection_df['ashman_d'] >= min_ashman_d
    selection_df['valley_pass'] = selection_df['valley_depth'] >= min_valley_depth
    selection_df['eligible'] = selection_df['status'].eq('bimodal') & selection_df['two_component_better'] & \
                               selection_df['weight_pass'] & selection_df['separation_pass'] & selection_df[
                                   'valley_pass']
    selection_df['selection_rank'] = np.nan
    selection_df['selected'] = False
    best_rows = []
    for sample, sample_df in selection_df.groupby('sample', sort=False):
        eligible_df = sample_df.loc[sample_df['eligible']].sort_values(
            ['two_components_sufficient', 'quality_score', 'valley_depth', 'ashman_d', 'posterior_high_conf_fraction'],
            ascending=False)
        if eligible_df.empty:
            best_rows.append({'sample': sample, 'selection_status': 'no_valid_window'})
            continue
        selection_df.loc[eligible_df.index, 'selection_rank'] = np.arange(1, len(eligible_df) + 1)
        selected_idx = eligible_df.index[0]
        selection_df.loc[selected_idx, 'selected'] = True
        best_row = selection_df.loc[selected_idx].to_dict()
        best_row['selection_status'] = 'selected'
        best_rows.append(best_row)
    best_df = pd.DataFrame(best_rows).sort_values('sample').reset_index(drop=True)
    return selection_df, best_df


def assign_best_cnv_labels(cnv_score_df: pd.DataFrame, best_df: pd.DataFrame,
                           sample_col: str = 'sample') -> pd.DataFrame:
    """Assign epithelial Cancer/Normal labels using each sample's selected cutoff.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin score table from :func:`calculate_cnv_scores`.
    best_df
        Selected-window table from :func:`select_best_cnv_score`.
    sample_col
        Sample identifier column shared by the two input tables.

    Returns
    -------
    pandas.DataFrame
        Per-epithelial-cellbin metadata and the columns ``best_score_col``,
        ``best_cnv_score``, ``best_cnv_threshold``, and ``cnv_auto_state``.

    Notes
    -----
    For selected samples, cells with ``best_cnv_score >= best_cnv_threshold``
    are labeled ``Cancer cell``; lower cells are labeled ``Normal Epithelial``.
    If no score setting passes quality selection for a sample, its epithelial
    cellbins remain ``Epithelial uncertain``.
    """
    metadata_cols = [col for col in ['chipid', sample_col, 'cell_type', 'leiden', 'x', 'y', 'x_slide_mm', 'y_slide_mm']
                     if col in cnv_score_df.columns]
    classification_df = cnv_score_df[metadata_cols].copy()
    classification_df['best_score_col'] = np.nan
    classification_df['best_cnv_score'] = np.nan
    classification_df['best_cnv_threshold'] = np.nan
    classification_df['cnv_auto_state'] = 'Epithelial uncertain'
    selected_df = best_df.loc[best_df['selection_status'].eq('selected')]
    for _, row in selected_df.iterrows():
        sample_mask = cnv_score_df[sample_col].astype(str).eq(str(row['sample']))
        score_col = row['score_col']
        classification_df.loc[sample_mask, 'best_score_col'] = score_col
        classification_df.loc[sample_mask, 'best_cnv_score'] = cnv_score_df.loc[sample_mask, score_col].to_numpy()
        classification_df.loc[sample_mask, 'best_cnv_threshold'] = row['kde_valley_threshold']
        classification_df.loc[sample_mask, 'cnv_auto_state'] = np.where(
            classification_df.loc[sample_mask, 'best_cnv_score'].to_numpy() >= row['kde_valley_threshold'],
            'Cancer cell', 'Normal Epithelial')
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
    """Plot window-selection quality and KDE-threshold heatmaps as PNG and PDF.

    Parameters
    ----------
    selection_df
        Output of :func:`select_best_cnv_score`.
    outdir
        Pipeline output directory. Heatmaps are written to ``{outdir}/figures``
        through :func:`_savefig` in both PNG and PDF formats.

    Returns
    -------
    None

    Notes
    -----
    An asterisk marks the selected score setting for each sample. The quality
    heatmap shows only eligible settings; the threshold heatmap shows all KDE-
    bimodal settings, whether or not they passed final eligibility filters.
    """
    figure_dir = Path(outdir) / 'figures'
    sample_order = selection_df['sample'].drop_duplicates().tolist()
    score_order = selection_df.sort_values('window_size')['score_col'].drop_duplicates().tolist()
    quality_df = selection_df.loc[selection_df['eligible']].pivot(index='sample', columns='score_col',
                                                                  values='quality_score').reindex(index=sample_order,
        columns=score_order)
    fig, ax = plt.subplots(figsize=(max(7, len(score_order) * 2), max(5, len(sample_order) * 0.55)))
    sns.heatmap(quality_df, annot=True, fmt='.2f', ax=ax)
    _mark_selected(ax, selection_df, sample_order, score_order)
    ax.set_xlabel('CNV score')
    ax.set_ylabel('Sample')
    ax.set_title('Eligible CNV-score window quality; * selected')
    fig.tight_layout()
    _savefig(fig, figure_dir, 'CNV_window_selection_quality_heatmap')
    threshold_df = selection_df.loc[selection_df['status'].eq('bimodal')].pivot(index='sample', columns='score_col',
        values='kde_valley_threshold').reindex(index=sample_order, columns=score_order)
    fig, ax = plt.subplots(figsize=(max(7, len(score_order) * 2), max(5, len(sample_order) * 0.55)))
    sns.heatmap(threshold_df, annot=True, fmt='.5f', ax=ax)
    _mark_selected(ax, selection_df, sample_order, score_order)
    ax.set_xlabel('CNV score')
    ax.set_ylabel('Sample')
    ax.set_title('KDE valley thresholds; * selected')
    fig.tight_layout()
    _savefig(fig, figure_dir, 'CNV_window_KDE_threshold_heatmap')


def plot_kde_bimodal_threshold(score_values: np.ndarray, result_row: pd.Series, outdir: Union[str, Path],
                               n_grid: int = 5000, bw_scale: float = 1.0) -> None:
    """Plot one sample's histogram, KDE, GMM diagnostics, and selected cutoff.

    Parameters
    ----------
    score_values
        Epithelial CNV-score values from the selected sample and window setting.
    result_row
        One selected row from ``best_df``. It must contain GMM parameters,
        KDE-mode positions, and ``kde_valley_threshold``.
    outdir
        Pipeline output directory. Figures are saved below
        ``{outdir}/figures/KDE_bimodal_threshold`` in PNG and PDF formats.
    n_grid
        Density-evaluation grid length. Must match the general analysis setting
        for visual comparability, although it does not change the stored cutoff.
    bw_scale
        KDE bandwidth multiplier used for the original cutoff calculation.

    Returns
    -------
    None

    Notes
    -----
    The dashed vertical line is the KDE valley threshold used for classification.
    The GMM component curves and mixture are displayed as diagnostics only.
    """
    figure_dir = Path(outdir) / 'figures' / 'KDE_bimodal_threshold'
    score_values = np.asarray(score_values, dtype=float)
    grid = np.linspace(score_values.min(), score_values.max(), n_grid)
    kde = gaussian_kde(score_values, bw_method=lambda kde_obj: kde_obj.scotts_factor() * bw_scale)
    kde_density = kde(grid)
    gmm_means = np.array([result_row['low_gmm_mean'], result_row['high_gmm_mean']])
    gmm_stds = np.array([result_row['low_gmm_std'], result_row['high_gmm_std']])
    gmm_weights = np.array([result_row['low_gmm_weight'], result_row['high_gmm_weight']])
    gmm_component_density = np.vstack(
        [weight * norm.pdf(grid, loc=mean, scale=std) for mean, std, weight in zip(gmm_means, gmm_stds, gmm_weights)])
    gmm_mixture_density = gmm_component_density.sum(axis=0)
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.histplot(score_values, bins=100, stat='density', element='step', fill=False, ax=ax, label='Histogram')
    ax.plot(grid, kde_density, label='KDE')
    ax.plot(grid, gmm_component_density[0], label='GMM low component')
    ax.plot(grid, gmm_component_density[1], label='GMM high component')
    ax.plot(grid, gmm_mixture_density, label='GMM mixture')
    ax.scatter([result_row['low_kde_peak'], result_row['high_kde_peak']],
               kde([result_row['low_kde_peak'], result_row['high_kde_peak']]), s=35, label='KDE modes')
    ax.axvline(result_row['kde_valley_threshold'], linestyle='--',
               label=f'KDE valley = {result_row["kde_valley_threshold"]:.6f}')
    ax.set_xlabel(result_row['score_col'])
    ax.set_ylabel('Density')
    ax.set_title(f'{result_row["sample"]} | BIC(1)={result_row["bic_1component"]:.1f}, '
                 f'BIC(2)={result_row["bic_2component"]:.1f}, Ashman D={result_row["ashman_d"]:.2f}')
    ax.legend()
    fig.tight_layout()
    _savefig(fig, figure_dir, f'{_safe_name(result_row["sample"])}_{result_row["score_col"]}_KDE_bimodal_threshold')


def plot_best_kde_bimodal_thresholds(cnv_score_df: pd.DataFrame, best_df: pd.DataFrame, outdir: Union[str, Path],
                                     sample_col: str = 'sample', n_grid: int = 5000, bw_scale: float = 1.0) -> None:
    """Create one KDE/GMM threshold figure for every sample with a selected score.

    Parameters
    ----------
    cnv_score_df
        Per-epithelial-cellbin score table.
    best_df
        Output of :func:`select_best_cnv_score`.
    outdir
        Pipeline output directory.
    sample_col
        Sample identifier column in ``cnv_score_df``.
    n_grid
        Density-evaluation grid length passed to
        :func:`plot_kde_bimodal_threshold`.
    bw_scale
        KDE bandwidth multiplier passed to :func:`plot_kde_bimodal_threshold`.

    Returns
    -------
    None
    """
    selected_df = best_df.loc[best_df['selection_status'].eq('selected')]
    for _, row in selected_df.iterrows():
        score_values = cnv_score_df.loc[
            cnv_score_df[sample_col].astype(str).eq(str(row['sample'])), row['score_col']].dropna().to_numpy()
        plot_kde_bimodal_threshold(score_values, row, outdir=outdir, n_grid=n_grid, bw_scale=bw_scale)


def add_cnv_results_to_adata(adata: ad.AnnData, cnv_score_df: pd.DataFrame, classification_df: pd.DataFrame) -> None:
    """Write CNV scores and automatic epithelial labels back to ``adata.obs``.

    Parameters
    ----------
    adata
        Original AnnData object. It is modified in place.
    cnv_score_df
        Per-epithelial-cellbin CNV-score table.
    classification_df
        Output of :func:`assign_best_cnv_labels`.

    Returns
    -------
    None

    Added columns
    -------------
    - ``cnv_score_w{window_size}`` for every tested score setting;
    - ``best_score_col``;
    - ``best_cnv_score``;
    - ``best_cnv_threshold``;
    - ``cnv_auto_state``.
    """
    score_cols = [col for col in cnv_score_df.columns if col.startswith('cnv_score_w')]
    for score_col in score_cols:
        adata.obs[score_col] = np.nan
        adata.obs.loc[cnv_score_df.index, score_col] = cnv_score_df[score_col]
    for col in ['best_score_col', 'best_cnv_score', 'best_cnv_threshold', 'cnv_auto_state']:
        adata.obs[col] = pd.NA if col in ['best_score_col', 'cnv_auto_state'] else np.nan
        adata.obs.loc[classification_df.index, col] = classification_df[col]


def infercnv_cancer_normal_recipe(adata: ad.AnnData, gtf_file: str,
                                  window_configs: Sequence[Tuple[int, int]] = ((100, 10), (300, 30), (500, 50),
                                                                               (1000, 100)),
                                  cell_type_col: str = 'cell_type', sample_col: str = 'sample',
                                  epi_name: str = 'Epithelial', outdir: Union[str, Path] = 'analysis_result/infercnv',
                                  feature_id_col: str = 'transcript', chipid_col: Optional[str] = 'chipid',
                                  unassigned_name: str = 'Unassigned', layer: str = 'log1p', lfc_clip: float = 3,
                                  dynamic_threshold: float = 1.5, exclude_chromosomes: Sequence[str] = ('chrX', 'chrY'),
                                  chunksize: int = 1000, n_jobs: int = 8, n_grid: int = 5000, bw_scale: float = 1.0,
                                  peak_prominence_fraction: float = 0.03, peak_distance_fraction: float = 0.05,
                                  gmm_n_init: int = 50, gmm_max_iter: int = 1000, gmm_tol: float = 1e-6,
                                  gmm_reg_covar: float = 1e-8, random_state: int = 0, min_delta_bic_1vs2: float = 10,
                                  max_delta_bic_2vs3: float = 10, min_ashman_d: float = 2,
                                  min_component_weight: float = 0.05, min_valley_depth: float = 0.10,
                                  overwrite: OverwriteType = False) -> Dict[str, object]:
    """Run the complete per-sample inferCNV Cancer/Normal Epithelial pipeline.

    Parameters
    ----------
    adata
        Cellbin-level AnnData object. Required fields are:

        - ``adata.var[feature_id_col]``: GTF ``gene_id`` values;
        - ``adata.obs[cell_type_col]``: main cell-type labels;
        - ``adata.obs[sample_col]``: sample labels.

        The object is modified in place: genomic coordinates and final CNV
        fields are added to ``adata.var`` and ``adata.obs`` respectively.
    gtf_file
        GTF annotation file used to map genomic coordinates.
    window_configs
        Unique ``(window_size, step)`` pairs used for inferCNV smoothing. The
        default evaluates 100/10, 300/30, 500/50, and 1000/100 gene windows.
    cell_type_col
        Main cell-type annotation column in ``adata.obs``.
    sample_col
        Sample identifier column in ``adata.obs``. All CNV calculation,
        GMM/KDE fitting, window selection, and thresholds are sample-specific.
    epi_name
        Label identifying epithelial cellbins to divide into Cancer/Normal.
    outdir
        Output directory for CSV intermediates and figures.
    feature_id_col
        ``adata.var`` column containing GTF ``gene_id`` values. In the AE0058
        object this column is named ``transcript`` despite containing gene IDs.
    chipid_col
        Optional chip identifier copied into output tables.
    unassigned_name
        Cell-type label removed before inferCNV. Set to ``None`` to retain it.
    layer
        Normalized log-expression layer supplied to infercnvpy.
    lfc_clip
        infercnvpy log-fold-change clipping threshold.
    dynamic_threshold
        infercnvpy dynamic noise-filter multiplier.
    exclude_chromosomes
        Chromosomes excluded from CNV inference. This must use the same naming
        scheme generated by ``chromosome_prefix='chr'`` in
        :func:`add_gene_coordinates`; default is ``('chrX', 'chrY')``.
    chunksize
        infercnvpy cell chunk size.
    n_jobs
        Number of infercnvpy worker processes.
    n_grid
        Shared density-evaluation grid length used for KDE, GMM diagnostics,
        cutoff detection, and threshold plots.
    bw_scale
        Shared KDE bandwidth multiplier. Do not tune this separately for each
        sample, because sample-specific smoothing changes peak detectability and
        makes selected thresholds non-comparable.
    peak_prominence_fraction
        Minimum KDE peak prominence relative to maximum density.
    peak_distance_fraction
        Minimum KDE peak separation relative to ``n_grid``.
    gmm_n_init
        Number of GMM initializations per component count.
    gmm_max_iter
        Maximum EM iterations per GMM fit.
    gmm_tol
        GMM convergence tolerance.
    gmm_reg_covar
        GMM variance regularization.
    random_state
        Random seed used in all GMM fits.
    min_delta_bic_1vs2
        Minimum BIC improvement required to favour two components over one.
    max_delta_bic_2vs3
        Maximum allowed BIC improvement of three components over two.
    min_ashman_d
        Minimum two-component separation for a selectable setting.
    min_component_weight
        Minimum weight required for both fitted components.
    min_valley_depth
        Minimum relative KDE valley depth for a selectable setting.
    overwrite
        Cache control. Supported values:

        - ``False``: use existing files whenever possible;
        - ``True``: rerun every analysis and plotting step;
        - mapping, for example ``{'cnv_score': True, 'cutoff': False}``.

        Valid mapping keys are ``'cnv_score'``, ``'cutoff'``, ``'selection'``,
        and ``'plot'``. Dependency propagation is automatic: rerunning
        ``cnv_score`` also reruns cutoff, selection, labels, and figures;
        rerunning cutoff also reruns selection, labels, and figures.

    Returns
    -------
    dict
        Dictionary containing modified ``adata`` and the following tables:

        - ``cnv_score_df``;
        - ``cnv_score_summary_df``;
        - ``cutoff_df``;
        - ``selection_df``;
        - ``best_df``;
        - ``classification_df``.

    Output files
    ------------
    ``outdir`` receives:

    - ``feature_gene_coordinates.csv``;
    - ``epithelial_cell_CNV_scores_all_window_sizes.csv``;
    - ``epithelial_CNV_score_summary_by_sample.csv``;
    - ``CNV_GMM_KDE_thresholds_all_samples_windows.csv``;
    - ``CNV_score_window_selection_all_samples.csv``;
    - ``CNV_best_score_threshold_by_sample.csv``;
    - ``epithelial_best_CNV_classification.csv``;
    - ``figures/CNV_window_selection_quality_heatmap.png/.pdf``;
    - ``figures/CNV_window_KDE_threshold_heatmap.png/.pdf``;
    - ``figures/KDE_bimodal_threshold/{sample}_{score}_KDE_bimodal_threshold.png/.pdf``.

    References
    ----------
    See module-level ``REFERENCES`` and the detailed method description in
    :func:`fit_gmm_kde_cutoff`.


    Example
    ---------
    infercnv_outdir = 'analysis_result/infercnv'
    result = infercnv_cancer_normal_recipe(adata=adata, gtf_file='/data/pipeline/reference/human/genes/genes.gtf',
    window_configs=[(100, 10), (300, 30), (500, 50), (1000, 100)], cell_type_col='cell_type', sample_col='sample',
    epi_name='Epithelial', feature_id_col='transcript', chipid_col='chipid', unassigned_name='Unassigned',
    layer='log1p', lfc_clip=3, dynamic_threshold=1.5, exclude_chromosomes=('chrX', 'chrY'), chunksize=1000, n_jobs=8,
    n_grid=5000, bw_scale=1.0, peak_prominence_fraction=0.03, peak_distance_fraction=0.05, gmm_n_init=50,
    gmm_max_iter=1000, gmm_tol=1e-6, gmm_reg_covar=1e-8, random_state=0, min_delta_bic_1vs2=10, max_delta_bic_2vs3=10,
    min_ashman_d=2, min_component_weight=0.05, min_valley_depth=0.10, outdir=infercnv_outdir,
    overwrite={'cnv_score': False, 'cutoff': True, 'selection': True, 'plot': True})

    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    coordinate_path = outdir / 'feature_gene_coordinates.csv'
    cnv_score_path = outdir / 'epithelial_cell_CNV_scores_all_window_sizes.csv'
    cnv_summary_path = outdir / 'epithelial_CNV_score_summary_by_sample.csv'
    cutoff_path = outdir / 'CNV_GMM_KDE_thresholds_all_samples_windows.csv'
    selection_path = outdir / 'CNV_score_window_selection_all_samples.csv'
    best_path = outdir / 'CNV_best_score_threshold_by_sample.csv'
    classification_path = outdir / 'epithelial_best_CNV_classification.csv'
    figure_path = outdir / 'figures' / 'CNV_window_selection_quality_heatmap.png'
    cnv_feature_mask = add_gene_coordinates(adata, gtf_file=gtf_file, feature_id_col=feature_id_col)
    coordinate_cols = [col for col in [feature_id_col, 'real_gene_name', 'chromosome', 'start', 'end'] if
                       col in adata.var.columns]
    if _should_overwrite(overwrite, 'cnv_score') or not coordinate_path.exists():
        adata.var[coordinate_cols].to_csv(coordinate_path)
    print(f'总 feature 数：{adata.n_vars}')
    print(f'成功匹配 gene_id 坐标的 feature 数：{cnv_feature_mask.sum()}')
    print(f'未匹配 feature 数：{(~cnv_feature_mask).sum()}')
    run_cnv_score = _should_overwrite(overwrite, 'cnv_score') or not cnv_score_path.exists()
    run_cutoff = _should_overwrite(overwrite, 'cutoff') or run_cnv_score or not cutoff_path.exists()
    run_selection = _should_overwrite(overwrite,
                                      'selection') or run_cutoff or not selection_path.exists() or not best_path.exists()
    run_plot = _should_overwrite(overwrite, 'plot') or run_selection or not figure_path.exists()
    if run_cnv_score:
        cnv_score_df, cnv_score_summary_df = calculate_cnv_scores(adata, cnv_feature_mask=cnv_feature_mask,
            window_configs=window_configs, cell_type_col=cell_type_col, sample_col=sample_col, epi_name=epi_name,
            unassigned_name=unassigned_name, chipid_col=chipid_col, layer=layer, lfc_clip=lfc_clip,
            dynamic_threshold=dynamic_threshold, exclude_chromosomes=exclude_chromosomes, chunksize=chunksize,
            n_jobs=n_jobs)
        cnv_score_df.to_csv(cnv_score_path, index_label='obs_name')
        cnv_score_summary_df.to_csv(cnv_summary_path)
    else:
        cnv_score_df = pd.read_csv(cnv_score_path, index_col=0)
        summary_index_col = [0, 1] if chipid_col is not None and chipid_col in cnv_score_df.columns else 0
        cnv_score_summary_df = pd.read_csv(cnv_summary_path, index_col=summary_index_col)
    if run_cutoff:
        cutoff_df = calculate_gmm_kde_cutoffs(cnv_score_df, window_configs=window_configs, sample_col=sample_col,
            n_grid=n_grid, bw_scale=bw_scale, peak_prominence_fraction=peak_prominence_fraction,
            peak_distance_fraction=peak_distance_fraction, gmm_n_init=gmm_n_init, gmm_max_iter=gmm_max_iter,
            gmm_tol=gmm_tol, gmm_reg_covar=gmm_reg_covar, random_state=random_state)
        cutoff_df.to_csv(cutoff_path, index=False)
    else:
        cutoff_df = pd.read_csv(cutoff_path)
    if run_selection:
        selection_df, best_df = select_best_cnv_score(cutoff_df, min_delta_bic_1vs2=min_delta_bic_1vs2,
            max_delta_bic_2vs3=max_delta_bic_2vs3, min_ashman_d=min_ashman_d, min_component_weight=min_component_weight,
            min_valley_depth=min_valley_depth)
        selection_df.to_csv(selection_path, index=False)
        best_df.to_csv(best_path, index=False)
        classification_df = assign_best_cnv_labels(cnv_score_df, best_df, sample_col=sample_col)
        classification_df.to_csv(classification_path, index_label='obs_name')
    else:
        selection_df = pd.read_csv(selection_path)
        best_df = pd.read_csv(best_path)
        classification_df = pd.read_csv(classification_path, index_col=0)
    if run_plot:
        plot_window_selection_heatmaps(selection_df, outdir=outdir)
        plot_best_kde_bimodal_thresholds(cnv_score_df, best_df, outdir=outdir, sample_col=sample_col, n_grid=n_grid,
                                         bw_scale=bw_scale)
    add_cnv_results_to_adata(adata, cnv_score_df, classification_df)
    return {'adata': adata, 'cnv_score_df': cnv_score_df, 'cnv_score_summary_df': cnv_score_summary_df,
            'cutoff_df': cutoff_df, 'selection_df': selection_df, 'best_df': best_df,
            'classification_df': classification_df}
