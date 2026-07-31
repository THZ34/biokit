import os
import re
import json
import tarfile
import hashlib
import shutil
import subprocess
import argparse
from pathlib import Path

import pandas as pd
import nibabel as nib
from tqdm import tqdm


def package_mri_batches(batch_csv, batch_start, batch_end, output_dir, dx_upload_dir, bulk_root='/mnt/project',
                        dcm2niix='/usr/bin/dcm2niix', mri_index_csv='mri_file_df.csv'):
    WORKDIR = Path(output_dir) / 'workdir'
    PACKAGE_DIR = Path(output_dir) / 'packages'
    WORKDIR.mkdir(parents=True, exist_ok=True)
    PACKAGE_DIR.mkdir(parents=True, exist_ok=True)

    MRI_ROOT = f'{bulk_root}/Bulk/Whole Body MRI/Dixon'

    def get_existing_batches_dx(dx_dir):
        out = subprocess.check_output(['dx', 'ls', dx_dir], text=True)
        batches = set()
        for line in out.splitlines():
            m = re.match(r'batch_(\d+)\.tar\.gz$', line)
            if m: batches.add(int(m.group(1)))
        return batches

    def calc_md5(fp):
        h = hashlib.md5()
        with open(fp, 'rb') as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b''):
                h.update(chunk)
        return h.hexdigest()

    def infer_signal_type(meta):
        it = " ".join(meta.get("ImageType", [])).upper()
        sd = meta.get("SeriesDescription", "").upper()
        if "WATER" in it or sd.endswith("_W"): return "W"
        if "FAT" in it or sd.endswith("_F"): return "F"
        if "INPHASE" in it or "IN_PHASE" in it or sd.endswith("_IN"): return "in"
        if "OUTPHASE" in it or "OUT_OF_PHASE" in it or "OPPOSED" in it or sd.endswith("_OPP"): return "opp"
        return None

    def extract_spatial_info(nii_path):
        img = nib.load(nii_path)
        shape = img.shape[:3]
        voxel = img.header.get_zooms()[:3]
        return dict(shape_x=shape[0], shape_y=shape[1], shape_z=shape[2], voxel_size_x=voxel[0], voxel_size_y=voxel[1],
                    voxel_size_z=voxel[2], fov_x_mm=shape[0] * voxel[0], fov_y_mm=shape[1] * voxel[1],
                    fov_z_mm=shape[2] * voxel[2])

    sample_df = pd.read_csv(batch_csv)
    assert {'sample_id', 'download_batch'}.issubset(sample_df.columns)

    if not os.path.exists(mri_index_csv):
        import glob
        mri_paths = glob.glob(f'{MRI_ROOT}/*/*zip')
        mri_file_df = pd.DataFrame(mri_paths, columns=['full_path'])
        mri_file_df['file_name'] = mri_file_df['full_path'].apply(os.path.basename)
        mri_file_df['sample_id'] = mri_file_df['file_name'].str.split('_').str[0]
        mri_file_df['field'] = mri_file_df['file_name'].str.split('_').str[1]
        mri_file_df['instance'] = mri_file_df['file_name'].str.split('_').str[2]
        mri_file_df['i'] = mri_file_df['file_name'].str.split('.zip').str[0].str.split('_').str[3]
        mri_file_df.index = mri_file_df['sample_id'] + '_' + mri_file_df['instance']
        mri_file_df.to_csv(mri_index_csv)

    mri_df = pd.read_csv(mri_index_csv, index_col=0)
    mri_df = mri_df[mri_df['instance'] == 2]

    EXISTING_BATCHES = get_existing_batches_dx(dx_upload_dir)
    print('[INFO] Existing DX batches:', sorted(EXISTING_BATCHES))

    for batch_id, batch_df in sample_df.groupby('download_batch'):
        batch_id = int(batch_id)
        if batch_id < batch_start or batch_id > batch_end: continue
        if batch_id in EXISTING_BATCHES:
            print(f'[SKIP] batch_{batch_id:02d} exists on DX')
            continue

        batch_name = f'batch_{batch_id:02d}'
        print(f'\n[RUN] Processing {batch_name} | n={len(batch_df)}')
        batch_root = WORKDIR / batch_name
        batch_root.mkdir(exist_ok=True)
        included = []

        for sid in tqdm(batch_df['sample_id'], desc=batch_name):
            sid = str(sid)
            key = f'{sid}_2'
            if key not in mri_df.index: continue

            zip_path = mri_df.loc[key, 'full_path']
            sample_dir = batch_root / key
            dicom_dir = sample_dir / 'dicom'
            nii_dir = sample_dir / 'nii'

            dicom_dir.mkdir(parents=True, exist_ok=True)
            subprocess.run(['unzip', '-q', zip_path, '-d', dicom_dir], check=True)

            nii_dir.mkdir(exist_ok=True)
            subprocess.run([dcm2niix, '-z', 'y', '-b', 'y', '-ba', 'n', '-f', '%p_%s', '-o', nii_dir, dicom_dir],
                           check=True)

            records = []
            for nii in nii_dir.glob('*.nii.gz'):
                js = nii.with_suffix('').with_suffix('.json')
                meta = {}
                if js.exists():
                    with open(js) as f: meta = json.load(f)
                row = dict(sample_id=sid, instance=2, nii_name=nii.name, signal_type=infer_signal_type(meta))
                row.update(extract_spatial_info(nii))
                records.append(row)

            pd.DataFrame(records).to_csv(nii_dir / 'nii_info.csv', index=False)
            shutil.rmtree(dicom_dir)
            included.append(key)

        if not included:
            shutil.rmtree(batch_root)
            continue

        tar_path = PACKAGE_DIR / f'{batch_name}.tar.gz'
        with tarfile.open(tar_path, 'w:gz') as tar:
            for k in included:
                tar.add(batch_root / k, arcname=f'{batch_name}/{k}')

        md5 = calc_md5(tar_path)
        md5_path = tar_path.with_suffix('.tar.gz.md5')
        md5_path.write_text(f'{md5}  {tar_path.name}\n')

        subprocess.run(['dx', 'upload', tar_path, '--destination', dx_upload_dir], check=True)
        subprocess.run(['dx', 'upload', md5_path, '--destination', dx_upload_dir], check=True)

        shutil.rmtree(batch_root)
        tar_path.unlink()
        md5_path.unlink()

    print('\n[OK] ALL DONE')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch-csv', required=True)
    parser.add_argument('--batch-start', type=int, required=True)
    parser.add_argument('--batch-end', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--dx-upload-dir', required=True)
    parser.add_argument('--bulk-root', default='/mnt/project')
    parser.add_argument('--dcm2niix', default='/usr/bin/dcm2niix')
    parser.add_argument('--mri-index-csv', default='mri_file_df.csv')

    args = parser.parse_args()
    dx_upload_dir = args.dx_upload_dir
    if not dx_upload_dir.endswith('/'):
        dx_upload_dir += '/'

    package_mri_batches(batch_csv=args.batch_csv, batch_start=args.batch_start, batch_end=args.batch_end,
                        output_dir=args.output_dir, dx_upload_dir=dx_upload_dir, bulk_root=args.bulk_root,
                        dcm2niix=args.dcm2niix, mri_index_csv=args.mri_index_csv)
