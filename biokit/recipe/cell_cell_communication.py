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


def read_ccc_recipe(outdir_root, samples=None):
    if not outdir_root:
        raise ValueError("Please provide a valid outdir_root path.")
    if not samples:
        samples = os.listdir(outdir_root)

    mean_df_dict = {}
    pvalue_df_dict = {}
