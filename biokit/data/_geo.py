# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Thread
import shutil
import subprocess
import tempfile
import re
import time
import gzip
import pandas as pd
import os


def _soft_url(gseid):
    prefix = f"{gseid[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gseid}/soft/{gseid}_family.soft.gz"


def _norm_key(s):
    s = s.strip().lower().replace(' ', '_')
    s = re.sub(r'[^0-9a-z_]+', '_', s)
    return s.strip('_')


def _parse_characteristics(lines):
    out = {}
    for s in lines:
        if s is None: continue
        s = str(s).strip()
        if not s: continue
        if ':' in s:
            k, v = s.split(':', 1)
            out[_norm_key(k)] = v.strip()
        else:
            out.setdefault('characteristics_ch1_raw', []).append(s)
    if 'characteristics_ch1_raw' in out:
        out['characteristics_ch1_raw'] = '; '.join(out['characteristics_ch1_raw'])
    return out


def get_gse_sampleinfo(gseid, dest='dataset', save_tsv=True):
    import requests
    os.makedirs(f"{dest}/{gseid}", exist_ok=True)
    url = _soft_url(gseid)
    gz_path = f"{dest}/{gseid}/{gseid}_family.soft.gz"
    if not os.path.exists(gz_path):
        r = requests.get(url, timeout=600)
        r.raise_for_status()
        open(gz_path, 'wb').write(r.content)
    samples = {}
    with gzip.open(gz_path, 'rt', encoding='utf-8', errors='ignore') as fh:
        cur = None
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE = '):
                cur = line.split('=', 1)[1].strip()
                samples[cur] = {}
                samples[cur]['_char_lines'] = []
                continue
            if cur is None: continue
            if line.startswith('!Sample_') and '=' in line:
                k, v = line.split('=', 1)
                k = k.strip()[1:]  # drop leading '!'
                v = v.strip()
                if k.lower().startswith('sample_characteristics_ch1'):
                    samples[cur]['_char_lines'].append(v)
                else:
                    key = _norm_key(k.replace('Sample_', ''))
                    samples[cur].setdefault(key, []).append(v)
    rows = []
    for gsm, d in samples.items():
        row = {}
        for k, vals in d.items():
            if k == '_char_lines': continue
            if isinstance(vals, list): row[k] = '; '.join([str(x) for x in vals if x is not None])
            else: row[k] = str(vals)
        ch = _parse_characteristics(d.get('_char_lines', []))
        row.update(ch)
        row['gsm'] = gsm
        rows.append(row)
    df = pd.DataFrame(rows).set_index('gsm', drop=True)
    # 常用字段别名补充
    if 'title' not in df.columns and 'title_ch1' in df.columns: df['title'] = df['title_ch1']
    if 'platform_id' not in df.columns and 'gpl' in df.columns: df['platform_id'] = df['gpl']
    # 关系字段中提取 SRA/BioProject/BioSample
    def pick_relation(x, pat):
        if not isinstance(x, str): return ''
        m = re.findall(pat, x, flags=re.I)
        return m[0] if m else ''
    rel = df['relation'] if 'relation' in df.columns else pd.Series('', index=df.index)
    df['sra_study'] = rel.apply(lambda s: pick_relation(s, r'(SRP[0-9]+)'))
    df['bioproject'] = rel.apply(lambda s: pick_relation(s, r'(PRJ[EN][A-Z0-9]+)'))
    df['biosample'] = rel.apply(lambda s: pick_relation(s, r'(SAMN?[A-Z0-9]+)'))
    if save_tsv:
        df.to_csv(f"{dest}/{gseid}/clinical_from_soft.tsv", sep='\t')
    return df

# def get_gse_sampleinfo(gseid):
#     url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gseid}&targ=self&view=brief&form=text'
#     series_info = requests.get(url).text
#     pattern = re.compile(r'!Series_sample_id = (GSM\d+)')
#     gsmids = pattern.findall(series_info)
#     t_list = []
#     sample_info = {}
#     for gsmid in gsmids:
#         t = Thread(target=get_gsm_sampleinfo, args=(gsmid, sample_info))
#         t.start()
#         time.sleep(0.5)
#         t_list.append(t)
#     n_completed = 0
#     n_samples = len(t_list)
#     for t in t_list:
#         t.join()
#         n_completed += 1
#         print(f'已爬取: {n_completed}/{n_samples}\r', end='')
#
#     # 补漏
#     for gsmid in gsmids:
#         if gsmid not in sample_info:
#             try:
#                 get_gsm_sampleinfo(gsmid, sample_info)
#             except:
#                 continue
#
#     sample_info_df = pd.DataFrame(sample_info).T
#     return sample_info_df


def get_gsm_sampleinfo(gsmid, sampleinfo_dict):
    import requests
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(requests.get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsmid}', timeout=5).text,
                         'html.parser')
    tr_tags = soup.find_all('tr', valign='top')
    gsm_info = {}
    for tr_tag in tr_tags:
        td_list = tr_tag.find_all('td')
        if len(td_list) == 2:
            if td_list[0].text == 'Characteristics':
                content = td_list[1].get_text(separator="\n")
                for line in content.split('\n'):
                    if line:
                        key = line.split(': ')[0]
                        value = line.split(': ')[1] if len(line.split(': ')) == 2 else line.split(': ')[1:]
                        gsm_info[key] = value
                break
            else:
                gsm_info[td_list[0].text] = td_list[1].text
    gsm_info['gsmid'] = gsmid
    sampleinfo_dict[gsmid] = gsm_info


def download_gse(gseid, output=None, max_workers=4):
    """下载GSE数据集的补充文件

    :param gseid: GSE编号，如GSE132465
    :param output: 输出目录，如None则不下载，只返回文件链接列表
    :param max_workers: 并发下载任务数
    :return:
    """
    from bs4 import BeautifulSoup
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gseid[:-3]}nnn/{gseid}/suppl'
    downloader = shutil.which('aria2c') or shutil.which('wget')
    if not downloader:
        raise EnvironmentError('download_gse requires aria2c or wget to be installed')

    if os.path.basename(downloader).lower().startswith('aria2c') and shutil.which('wget'):
        list_cmd = [shutil.which('wget'), '-q', '-O', '-', url]
        list_proc = subprocess.Popen(list_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        list_stdout, list_stderr = list_proc.communicate()
        if list_proc.returncode != 0:
            raise RuntimeError(
                'download_gse failed to fetch listing for {} with {}:\n{}'.format(
                    url,
                    list_cmd[0],
                    list_stderr.decode('utf-8', errors='ignore') if list_stderr else '',
                )
            )
    elif os.path.basename(downloader).lower().startswith('aria2c'):
        tmp_dir = tempfile.mkdtemp(prefix='biokit_gse_')
        try:
            list_name = '{}_suppl_listing.html'.format(gseid)
            list_path = os.path.join(tmp_dir, list_name)
            list_cmd = [
                downloader,
                '--allow-overwrite=true',
                '--file-allocation=none',
                '--dir',
                tmp_dir,
                '--out',
                list_name,
                url,
            ]
            list_proc = subprocess.Popen(list_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            _, list_stderr = list_proc.communicate()
            if list_proc.returncode != 0:
                raise RuntimeError(
                    'download_gse failed to fetch listing for {} with {}:\n{}'.format(
                        url,
                        downloader,
                        list_stderr.decode('utf-8', errors='ignore') if list_stderr else '',
                    )
                )
            if not os.path.exists(list_path):
                raise RuntimeError('download_gse failed to fetch listing for {}'.format(url))
            with open(list_path, 'rb') as f:
                list_stdout = f.read()
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        list_cmd = [downloader, '-q', '-O', '-', url]
        list_proc = subprocess.Popen(list_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        list_stdout, list_stderr = list_proc.communicate()
        if list_proc.returncode != 0:
            raise RuntimeError(
                'download_gse failed to fetch listing for {} with {}:\n{}'.format(
                    url,
                    downloader,
                    list_stderr.decode('utf-8', errors='ignore') if list_stderr else '',
                )
            )

    soup = BeautifulSoup(list_stdout.decode('utf-8', errors='ignore'), 'html.parser')
    a_tags = soup.find_all('a')
    urls = []
    files = []
    for a_tag in a_tags:
        if not a_tag.text == 'Parent Directory' and not a_tag.text == 'HHS Vulnerability Disclosure':
            filename = a_tag.get('href')
            if not filename:
                continue
            file_url = f'{url}/{filename}'
            files.append((filename, file_url))
            urls.append(file_url)

    if output and files:
        os.makedirs(f'{output}/{gseid}', exist_ok=True)

        def _download_one(item):
            filename, file_url = item
            dst = os.path.join(output, gseid, filename)
            if os.path.exists(dst):
                return file_url

            if os.path.basename(downloader).lower().startswith('aria2c'):
                cmd = [
                    downloader,
                    '--allow-overwrite=true',
                    '--continue=true',
                    '--file-allocation=none',
                    '--split=16',
                    '--max-connection-per-server=16',
                    '--min-split-size=1M',
                    '--dir',
                    os.path.join(output, gseid),
                    '--out',
                    filename,
                    file_url,
                ]
            else:
                cmd = [
                    downloader,
                    '-c',
                    '-O',
                    dst,
                    file_url,
                ]

            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                raise RuntimeError(
                    'download_gse failed for {} with {}:\n{}'.format(
                        file_url, downloader, stderr.decode('utf-8', errors='ignore') if stderr else ''
                    )
                )
            return file_url

        max_workers = max(1, min(max_workers, len(files)))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_download_one, item) for item in files]
            for future in as_completed(futures):
                future.result()
    return urls


def get_gse_info(gseid, info_dict):
    import requests
    from bs4 import BeautifulSoup
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gseid}'
    soup = BeautifulSoup(requests.get(url).text, 'html.parser')
    temp_info_dict = {}
    status_tag = soup.find_all('td', string='Status')[0]
    table_tag = status_tag.parent.parent
    for tr_line in table_tag.find_all('tr', attrs={'valign': 'top'}):
        try:
            key = tr_line.find_all('td')[0].text
            value = tr_line.find_all('td')[1].text
            if f'Series {gseid}' in key:
                continue

            temp_info_dict[key] = value
        except IndexError:
            continue
    info_dict[gseid] = temp_info_dict

def geoparse_sample_info(gseid):
    import GEOparse
    gse = GEOparse.get_GEO(geo=gseid, destdir="./GEO")

    sample_info = {}
    for gsm_name, gsm in gse.gsms.items():
        sample_info[gsm_name] = {}
        for key_value in gsm.metadata.get('characteristics_ch1', 'No characteristics'):
            key, value = key_value.split(': ')[0], ': '.join(key_value.split(': ')[1:])
            sample_info[gsm_name][key] = value
        sample_info[gsm_name]['Title'] = gsm.metadata.get('title', ['No title'])[0]
    sample_info = pd.DataFrame(sample_info).T
    return sample_info