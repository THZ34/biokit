# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import logging
import numpy as np
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import os
from multiprocessing import cpu_count
import glob
import logging
import os
import glob
import traceback
from pathlib import Path
from multiprocessing import Pool, cpu_count

import anndata as ad
import pandas as pd
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

from multiprocessing import cpu_count, get_context
from queue import Empty
import traceback


def _run_cpdb_single_sample(task, result_queue):
    sample, outdir, meta_file_path, counts_file_path, cpdb_file_path, threads = task
    try:
        cpdb_statistical_analysis_method.call(cpdb_file_path=cpdb_file_path, meta_file_path=meta_file_path,
                                              counts_file_path=counts_file_path, counts_data='hgnc_symbol',
                                              output_path=outdir, threads=threads)
        result = {'sample': sample, 'status': 'completed', 'outdir': outdir, 'error': None}
    except Exception:
        result = {'sample': sample, 'status': 'failed', 'outdir': outdir, 'error': traceback.format_exc()}
    result_queue.put(result)


def _run_cpdb_tasks_non_daemon(tasks, processes):
    if len(tasks) == 0:
        return {}

    context = get_context('fork')
    result_queue = context.Queue()
    task_iter = iter(tasks)
    active_processes = {}
    run_status = {}

    def start_next_task():
        try:
            task = next(task_iter)
        except StopIteration:
            return False

        sample = task[0]
        process = context.Process(target=_run_cpdb_single_sample, args=(task, result_queue), name=f'CPDB_{sample}')
        process.daemon = False
        process.start()
        active_processes[sample] = process
        logging.info(f'CellPhoneDB started: {sample}, PID={process.pid}')
        return True

    for _ in range(min(processes, len(tasks))):
        start_next_task()

    while active_processes:
        try:
            result = result_queue.get(timeout=5)
            sample = result['sample']
            process = active_processes.pop(sample, None)

            if process is not None:
                process.join()

            run_status[sample] = result

            if result['status'] == 'completed':
                logging.info(f'CellPhoneDB completed: {sample}')
            else:
                logging.error(f'CellPhoneDB failed: {sample}\n{result["error"]}')

            start_next_task()

        except Empty:
            crashed_samples = []

            for sample, process in active_processes.items():
                if not process.is_alive() and process.exitcode not in (0, None):
                    crashed_samples.append((sample, process))

            for sample, process in crashed_samples:
                process.join()
                active_processes.pop(sample)

                run_status[sample] = {'sample': sample, 'status': 'failed', 'outdir': None,
                                      'error': f'CellPhoneDB child process exited unexpectedly, exitcode={process.exitcode}.'}

                logging.error(f'CellPhoneDB child process failed: {sample}, exitcode={process.exitcode}')
                start_next_task()

    result_queue.close()
    result_queue.join_thread()
    return run_status


def _write_cpdb_input_files(adata, cell_mask, cell_type_col, sample, outdir):
    sample_adata = adata[cell_mask].copy()
    meta_file_path = outdir / f'{sample}_meta.txt'
    counts_file_path = outdir / f'{sample}_expression.h5ad'

    meta_df = sample_adata.obs[[cell_type_col]].copy()
    meta_df.index.name = 'Cell'
    meta_df.to_csv(meta_file_path, sep='\t')

    counts_adata = ad.AnnData(X=sample_adata.X, obs=pd.DataFrame(index=sample_adata.obs_names.copy()),
                              var=pd.DataFrame(index=sample_adata.var_names.copy()), )
    counts_adata.obs_names.name = 'Cell'
    counts_adata.var_names.name = 'Gene'
    counts_adata.write_h5ad(counts_file_path)

    del sample_adata, counts_adata
    return str(meta_file_path), str(counts_file_path)


def ccc_recipe(adata, cell_type_col, sample_col=None, samples=None, outdir_root=None, cpdb_file_path=None,
               processes=None, threads=1, overwrite=False):
    if cell_type_col not in adata.obs.columns:
        raise KeyError(f'cell_type_col "{cell_type_col}" is not present in adata.obs.')

    if sample_col is None:
        sample_col = '_ccc_sample'
        sample_series = pd.Series('Sample1', index=adata.obs_names, name=sample_col)
    else:
        if sample_col not in adata.obs.columns:
            raise KeyError(f'sample_col "{sample_col}" is not present in adata.obs.')
        sample_series = adata.obs[sample_col]

    if samples is None:
        samples = sample_series.dropna().unique().tolist()
    else:
        samples = list(samples)

    if len(samples) == 0:
        raise ValueError('No samples were found for CellPhoneDB analysis.')

    if outdir_root is None:
        outdir_root = Path.cwd() / 'cellphonedb'
    else:
        outdir_root = Path(outdir_root).resolve()

    if cpdb_file_path is None:
        cpdb_file_path = '/data/project/tanghongzhen/data/database/cellphonedb/cellphonedb_v5.zip'

    cpdb_file_path = str(Path(cpdb_file_path).resolve())
    if not os.path.isfile(cpdb_file_path):
        raise FileNotFoundError(f'CellPhoneDB database was not found: {cpdb_file_path}')

    if processes is None:
        processes = min(len(samples), max(cpu_count() // max(threads, 1), 1))

    if not isinstance(processes, int) or processes < 1:
        raise ValueError('processes must be a positive integer.')

    if not isinstance(threads, int) or threads < 1:
        raise ValueError('threads must be a positive integer.')

    total_cpu_requested = processes * threads
    if total_cpu_requested > cpu_count():
        logging.warning(f'processes × threads = {total_cpu_requested}, larger than available CPUs ({cpu_count()}). '
                        'This may cause CPU oversubscription.')

    outdir_root.mkdir(parents=True, exist_ok=True)
    tasks = []
    run_status = {}

    for sample in samples:
        sample = str(sample)
        outdir = outdir_root / sample
        mean_files = glob.glob(str(outdir / '*statistical_analysis_means*'))

        if mean_files and not overwrite:
            logging.warning(f'CellPhoneDB results for sample {sample} already exist. Skipping.')
            run_status[sample] = {'sample': sample, 'status': 'skipped', 'outdir': str(outdir), 'error': None, }
            continue

        cell_mask = sample_series.astype(str).eq(sample)
        n_cells = int(cell_mask.sum())

        if n_cells == 0:
            logging.warning(f'No cells were found for sample {sample}. Skipping.')
            run_status[sample] = {'sample': sample, 'status': 'skipped', 'outdir': str(outdir),
                                  'error': 'No cells found for this sample.', }
            continue

        logging.info(f'Preparing CellPhoneDB input for {sample}: {n_cells:,} cells.')
        outdir.mkdir(parents=True, exist_ok=True)

        meta_file_path, counts_file_path = _write_cpdb_input_files(adata=adata, cell_mask=cell_mask.to_numpy(),
                                                                   cell_type_col=cell_type_col, sample=sample,
                                                                   outdir=outdir, )

        tasks.append((sample, str(outdir), meta_file_path, counts_file_path, cpdb_file_path, threads,))

    if len(tasks) == 0:
        logging.warning('No CellPhoneDB tasks need to be run.')
        return run_status

    logging.info(f'Running {len(tasks)} CellPhoneDB tasks with processes={processes}, '
                 f'threads_per_process={threads}.')
    run_status.update(_run_cpdb_tasks_non_daemon(tasks=tasks, processes=min(processes, len(tasks)), ))
    return run_status


def read_ccc_recipe(outdir_root, samples=None):
    if not outdir_root:
        raise ValueError("Please provide a valid outdir_root path.")
    if not samples:
        samples = os.listdir(outdir_root)

    mean_df_dict = {}
    pvalue_df_dict = {}


# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import logging

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import os
from multiprocessing import cpu_count
import glob
import logging
import os
import glob
import traceback
from pathlib import Path
from multiprocessing import Pool, cpu_count

import anndata as ad
import pandas as pd
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

from multiprocessing import cpu_count, get_context
from queue import Empty
import traceback


def _run_cpdb_single_sample(task, result_queue):
    sample, outdir, meta_file_path, counts_file_path, cpdb_file_path, threads = task
    try:
        cpdb_statistical_analysis_method.call(cpdb_file_path=cpdb_file_path, meta_file_path=meta_file_path,
                                              counts_file_path=counts_file_path, counts_data='hgnc_symbol',
                                              output_path=outdir, threads=threads)
        result = {'sample': sample, 'status': 'completed', 'outdir': outdir, 'error': None}
    except Exception:
        result = {'sample': sample, 'status': 'failed', 'outdir': outdir, 'error': traceback.format_exc()}
    result_queue.put(result)


def _run_cpdb_tasks_non_daemon(tasks, processes):
    if len(tasks) == 0:
        return {}

    context = get_context('fork')
    result_queue = context.Queue()
    task_iter = iter(tasks)
    active_processes = {}
    run_status = {}

    def start_next_task():
        try:
            task = next(task_iter)
        except StopIteration:
            return False

        sample = task[0]
        process = context.Process(target=_run_cpdb_single_sample, args=(task, result_queue), name=f'CPDB_{sample}')
        process.daemon = False
        process.start()
        active_processes[sample] = process
        logging.info(f'CellPhoneDB started: {sample}, PID={process.pid}')
        return True

    for _ in range(min(processes, len(tasks))):
        start_next_task()

    while active_processes:
        try:
            result = result_queue.get(timeout=5)
            sample = result['sample']
            process = active_processes.pop(sample, None)

            if process is not None:
                process.join()

            run_status[sample] = result

            if result['status'] == 'completed':
                logging.info(f'CellPhoneDB completed: {sample}')
            else:
                logging.error(f'CellPhoneDB failed: {sample}\n{result["error"]}')

            start_next_task()

        except Empty:
            crashed_samples = []

            for sample, process in active_processes.items():
                if not process.is_alive() and process.exitcode not in (0, None):
                    crashed_samples.append((sample, process))

            for sample, process in crashed_samples:
                process.join()
                active_processes.pop(sample)

                run_status[sample] = {'sample': sample, 'status': 'failed', 'outdir': None,
                                      'error': f'CellPhoneDB child process exited unexpectedly, exitcode={process.exitcode}.'}

                logging.error(f'CellPhoneDB child process failed: {sample}, exitcode={process.exitcode}')
                start_next_task()

    result_queue.close()
    result_queue.join_thread()
    return run_status


def _write_cpdb_input_files(adata, cell_mask, cell_type_col, sample, outdir):
    sample_adata = adata[cell_mask].copy()
    meta_file_path = outdir / f'{sample}_meta.txt'
    counts_file_path = outdir / f'{sample}_expression.h5ad'

    meta_df = sample_adata.obs[[cell_type_col]].copy()
    meta_df.index.name = 'Cell'
    meta_df.to_csv(meta_file_path, sep='\t')

    counts_adata = ad.AnnData(X=sample_adata.X, obs=pd.DataFrame(index=sample_adata.obs_names.copy()),
                              var=pd.DataFrame(index=sample_adata.var_names.copy()), )
    counts_adata.obs_names.name = 'Cell'
    counts_adata.var_names.name = 'Gene'
    counts_adata.write_h5ad(counts_file_path)

    del sample_adata, counts_adata
    return str(meta_file_path), str(counts_file_path)


def ccc_recipe(adata, cell_type_col, sample_col=None, samples=None, outdir_root=None, cpdb_file_path=None,
               processes=None, threads=1, overwrite=False):
    if cell_type_col not in adata.obs.columns:
        raise KeyError(f'cell_type_col "{cell_type_col}" is not present in adata.obs.')

    if sample_col is None:
        sample_col = '_ccc_sample'
        sample_series = pd.Series('Sample1', index=adata.obs_names, name=sample_col)
    else:
        if sample_col not in adata.obs.columns:
            raise KeyError(f'sample_col "{sample_col}" is not present in adata.obs.')
        sample_series = adata.obs[sample_col]

    if samples is None:
        samples = sample_series.dropna().unique().tolist()
    else:
        samples = list(samples)

    if len(samples) == 0:
        raise ValueError('No samples were found for CellPhoneDB analysis.')

    if outdir_root is None:
        outdir_root = Path.cwd() / 'cellphonedb'
    else:
        outdir_root = Path(outdir_root).resolve()

    if cpdb_file_path is None:
        cpdb_file_path = '/data/project/tanghongzhen/data/database/cellphonedb/cellphonedb_v5.zip'

    cpdb_file_path = str(Path(cpdb_file_path).resolve())
    if not os.path.isfile(cpdb_file_path):
        raise FileNotFoundError(f'CellPhoneDB database was not found: {cpdb_file_path}')

    if processes is None:
        processes = min(len(samples), max(cpu_count() // max(threads, 1), 1))

    if not isinstance(processes, int) or processes < 1:
        raise ValueError('processes must be a positive integer.')

    if not isinstance(threads, int) or threads < 1:
        raise ValueError('threads must be a positive integer.')

    total_cpu_requested = processes * threads
    if total_cpu_requested > cpu_count():
        logging.warning(f'processes × threads = {total_cpu_requested}, larger than available CPUs ({cpu_count()}). '
                        'This may cause CPU oversubscription.')

    outdir_root.mkdir(parents=True, exist_ok=True)
    tasks = []
    run_status = {}

    for sample in samples:
        sample = str(sample)
        outdir = outdir_root / sample
        mean_files = glob.glob(str(outdir / '*statistical_analysis_means*'))

        if mean_files and not overwrite:
            logging.warning(f'CellPhoneDB results for sample {sample} already exist. Skipping.')
            run_status[sample] = {'sample': sample, 'status': 'skipped', 'outdir': str(outdir), 'error': None, }
            continue

        cell_mask = sample_series.astype(str).eq(sample)
        n_cells = int(cell_mask.sum())

        if n_cells == 0:
            logging.warning(f'No cells were found for sample {sample}. Skipping.')
            run_status[sample] = {'sample': sample, 'status': 'skipped', 'outdir': str(outdir),
                                  'error': 'No cells found for this sample.', }
            continue

        logging.info(f'Preparing CellPhoneDB input for {sample}: {n_cells:,} cells.')
        outdir.mkdir(parents=True, exist_ok=True)

        meta_file_path, counts_file_path = _write_cpdb_input_files(adata=adata, cell_mask=cell_mask.to_numpy(),
                                                                   cell_type_col=cell_type_col, sample=sample,
                                                                   outdir=outdir, )

        tasks.append((sample, str(outdir), meta_file_path, counts_file_path, cpdb_file_path, threads,))

    if len(tasks) == 0:
        logging.warning('No CellPhoneDB tasks need to be run.')
        return run_status

    logging.info(f'Running {len(tasks)} CellPhoneDB tasks with processes={processes}, '
                 f'threads_per_process={threads}.')
    run_status.update(_run_cpdb_tasks_non_daemon(tasks=tasks, processes=min(processes, len(tasks)), ))
    return run_status


def read_ccc_recipe(outdir_root, samples=None, alpha=0.05, save=True, save_dir=None):
    if not outdir_root:
        raise ValueError("Please provide a valid outdir_root path.")
    outdir_root = Path(outdir_root)
    if not outdir_root.exists():
        raise FileNotFoundError(f'outdir_root does not exist: {outdir_root}')
    if not samples:
        samples = sorted([sample for sample in os.listdir(outdir_root) if (outdir_root / sample).is_dir()])
    else:
        samples = [str(sample) for sample in samples]

    def get_latest_file(sample_dir, pattern):
        files = glob.glob(str(sample_dir / pattern))
        if len(files) == 0:
            return None
        return max(files, key=os.path.getmtime)

    def get_cellpair_cols(mean_df, pvalue_df):
        cellpair_cols = [col for col in mean_df.columns if col in pvalue_df.columns and '|' in col]
        if len(cellpair_cols) == 0:
            raise ValueError('No cell-pair columns were found. Please check CellPhoneDB output columns.')
        return cellpair_cols

    def add_cellpair_cols(df, cell_pair_col='cell_pair'):
        df = df.copy()
        if df.shape[0] == 0:
            df['source_cell_type'] = pd.Series(dtype=object)
            df['target_cell_type'] = pd.Series(dtype=object)
            return df
        split_values = df[cell_pair_col].astype(str).str.split('|', n=1, regex=False)
        df['source_cell_type'] = split_values.str[0]
        df['target_cell_type'] = split_values.str[1]
        return df

    mean_df_dict = {}
    pvalue_df_dict = {}
    significant_mean_df_dict = {}
    summary_df_list = []
    detail_df_list = []
    audit_rows = []

    for sample in samples:
        sample_dir = outdir_root / sample
        if not sample_dir.is_dir():
            audit_rows.append({'sample': sample, 'status': 'missing_sample_dir', 'mean_file': None, 'pvalue_file': None,
                               'significant_mean_file': None, 'error': f'Missing sample directory: {sample_dir}'})
            continue
        mean_file = get_latest_file(sample_dir, 'statistical_analysis_means*.txt')
        pvalue_file = get_latest_file(sample_dir, 'statistical_analysis_pvalues*.txt')
        significant_mean_file = get_latest_file(sample_dir, 'statistical_analysis_significant_means*.txt')
        if mean_file is None or pvalue_file is None:
            audit_rows.append({'sample': sample, 'status': 'missing_required_file', 'mean_file': mean_file,
                               'pvalue_file': pvalue_file, 'significant_mean_file': significant_mean_file,
                               'error': 'Missing means or pvalues file.'})
            continue
        mean_df = pd.read_csv(mean_file, sep='\t', dtype=object)
        pvalue_df = pd.read_csv(pvalue_file, sep='\t', dtype=object)
        if mean_df.shape[0] != pvalue_df.shape[0]:
            audit_rows.append({'sample': sample, 'status': 'failed', 'mean_file': mean_file, 'pvalue_file': pvalue_file,
                               'significant_mean_file': significant_mean_file,
                               'error': f'means and pvalues row count mismatch: {mean_df.shape[0]} vs {pvalue_df.shape[0]}'})
            continue
        if 'interacting_pair' not in mean_df.columns:
            audit_rows.append({'sample': sample, 'status': 'failed', 'mean_file': mean_file, 'pvalue_file': pvalue_file,
                               'significant_mean_file': significant_mean_file,
                               'error': 'Missing interacting_pair column in means file.'})
            continue
        mean_df_dict[sample] = mean_df
        pvalue_df_dict[sample] = pvalue_df
        if significant_mean_file is not None:
            significant_mean_df_dict[sample] = pd.read_csv(significant_mean_file, sep='\t', dtype=object)

        cellpair_cols = get_cellpair_cols(mean_df, pvalue_df)
        mean_num_df = mean_df.loc[:, cellpair_cols].apply(pd.to_numeric, errors='coerce')
        pvalue_num_df = pvalue_df.loc[:, cellpair_cols].apply(pd.to_numeric, errors='coerce')
        significant_mask_df = pvalue_num_df.lt(alpha) & mean_num_df.gt(0)
        strength_sum = mean_num_df.where(significant_mask_df).sum(axis=0, skipna=True)
        lr_count = significant_mask_df.sum(axis=0)
        strength_mean = strength_sum / lr_count.replace(0, np.nan)

        sample_summary_df = pd.DataFrame({'sample': sample, 'cell_pair': cellpair_cols,
                                          'significant_strength_sum': strength_sum.reindex(cellpair_cols).to_numpy(),
                                          'significant_strength_mean': strength_mean.reindex(cellpair_cols).to_numpy(),
                                          'significant_lr_count': lr_count.reindex(cellpair_cols).to_numpy()})
        sample_summary_df = add_cellpair_cols(sample_summary_df)
        sample_summary_df = sample_summary_df.loc[
            :, ['sample', 'cell_pair', 'source_cell_type', 'target_cell_type', 'significant_strength_sum',
                'significant_strength_mean', 'significant_lr_count']]
        summary_df_list.append(sample_summary_df)

        id_cols = [col for col in mean_df.columns if col not in cellpair_cols]
        mean_long_df = mean_num_df.where(significant_mask_df).copy()
        mean_long_df.insert(0, 'row_id', np.arange(mean_long_df.shape[0]))
        mean_long_df = mean_long_df.melt(id_vars='row_id', var_name='cell_pair',
                                         value_name='communication_strength').dropna(subset=['communication_strength'])
        pvalue_long_df = pvalue_num_df.copy()
        pvalue_long_df.insert(0, 'row_id', np.arange(pvalue_long_df.shape[0]))
        pvalue_long_df = pvalue_long_df.melt(id_vars='row_id', var_name='cell_pair', value_name='pvalue')
        id_df = mean_df.loc[:, id_cols].copy()
        id_df.insert(0, 'row_id', np.arange(id_df.shape[0]))
        sample_detail_df = mean_long_df.merge(pvalue_long_df, on=['row_id', 'cell_pair'], how='left').merge(id_df,
                                                                                                            on='row_id',
                                                                                                            how='left')
        sample_detail_df.insert(0, 'sample', sample)
        sample_detail_df = add_cellpair_cols(sample_detail_df)
        front_cols = ['sample', 'cell_pair', 'source_cell_type', 'target_cell_type', 'communication_strength', 'pvalue']
        sample_detail_df = sample_detail_df.loc[
            :, front_cols + [col for col in sample_detail_df.columns if col not in front_cols + ['row_id']]]
        detail_df_list.append(sample_detail_df)

        audit_rows.append({'sample': sample, 'status': 'completed', 'mean_file': mean_file, 'pvalue_file': pvalue_file,
                           'significant_mean_file': significant_mean_file, 'n_interactions': mean_df.shape[0],
                           'n_cell_pairs': len(cellpair_cols), 'n_significant_records': sample_detail_df.shape[0],
                           'error': None})

    audit_df = pd.DataFrame(audit_rows)
    failed_df = audit_df.loc[audit_df['status'] != 'completed'] if 'status' in audit_df.columns else pd.DataFrame()
    if len(summary_df_list) == 0:
        raise RuntimeError(f'No valid CellPhoneDB results were read. Audit:\n{audit_df.to_string(index=False)}')

    summary_df = pd.concat(summary_df_list, axis=0, ignore_index=True)
    if len(detail_df_list) == 0:
        detail_df = pd.DataFrame()
    else:
        detail_df = pd.concat(detail_df_list, axis=0, ignore_index=True)
    sample_order = sorted(summary_df['sample'].unique())
    cellpair_order = sorted(summary_df['cell_pair'].unique())
    strength_sum_df = summary_df.pivot(index='sample', columns='cell_pair', values='significant_strength_sum').reindex(
        index=sample_order, columns=cellpair_order).fillna(0)
    strength_mean_df = summary_df.pivot(index='sample', columns='cell_pair',
                                        values='significant_strength_mean').reindex(index=sample_order,
                                                                                    columns=cellpair_order)
    lr_count_df = summary_df.pivot(index='sample', columns='cell_pair', values='significant_lr_count').reindex(
        index=sample_order, columns=cellpair_order).fillna(0).astype(int)

    result = {'mean_df_dict': mean_df_dict, 'pvalue_df_dict': pvalue_df_dict,
              'significant_mean_df_dict': significant_mean_df_dict, 'summary_df': summary_df,
              'sample_cellpair_significant_strength_sum_df': strength_sum_df,
              'sample_cellpair_significant_strength_mean_df': strength_mean_df,
              'sample_cellpair_significant_lr_count_df': lr_count_df,
              'sample_cellpair_lr_strength_detail_df': detail_df, 'audit_df': audit_df}

    if save:
        if save_dir is None:
            save_dir = outdir_root.parent / f'{outdir_root.name}_integrated'
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        strength_sum_df.to_csv(save_dir / 'sample_cellpair_significant_strength_sum.csv', encoding='utf-8-sig')
        strength_mean_df.to_csv(save_dir / 'sample_cellpair_significant_strength_mean.csv', encoding='utf-8-sig')
        lr_count_df.to_csv(save_dir / 'sample_cellpair_significant_lr_count.csv', encoding='utf-8-sig')
        detail_df.to_csv(save_dir / 'sample_cellpair_lr_strength_detail.csv.gz', index=False, encoding='utf-8-sig')
        summary_df.to_csv(save_dir / 'sample_cellpair_significant_summary_long.csv', index=False, encoding='utf-8-sig')
        audit_df.to_csv(save_dir / 'read_ccc_recipe_audit.csv', index=False, encoding='utf-8-sig')
        with pd.ExcelWriter(save_dir / 'CellPhoneDB_cell_type_integrated.xlsx') as writer:
            strength_sum_df.to_excel(writer, sheet_name='01_strength_sum')
            strength_mean_df.to_excel(writer, sheet_name='02_strength_mean')
            lr_count_df.to_excel(writer, sheet_name='03_lr_count')
            summary_df.to_excel(writer, sheet_name='04_summary_long', index=False)
            detail_df.to_excel(writer, sheet_name='05_lr_strength_detail', index=False)
            audit_df.to_excel(writer, sheet_name='06_audit', index=False)
        result['save_dir'] = str(save_dir)

    if len(failed_df) > 0:
        print('以下样本未成功读取，请检查 audit_df：')
        print(failed_df.to_string(index=False))
    print(
        f'读取完成：成功样本数={len(mean_df_dict)}，细胞对数={len(cellpair_order)}，显著样本-细胞对-LR记录数={detail_df.shape[0]}')
    return result
