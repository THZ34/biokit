# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

from threading import Thread
import requests
import re
import time
import pandas as pd
from bs4 import BeautifulSoup
import os


def get_gse_sampleinfo(gseid):
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gseid}&targ=self&view=brief&form=text'
    series_info = requests.get(url).text
    pattern = re.compile(r'!Series_sample_id = (GSM\d+)')
    gsmids = pattern.findall(series_info)
    t_list = []
    sample_info = {}
    for gsmid in gsmids:
        t = Thread(target=get_gsm_sampleinfo, args=(gsmid, sample_info))
        t.start()
        time.sleep(0.5)
        t_list.append(t)
    n_completed = 0
    n_samples = len(t_list)
    for t in t_list:
        t.join()
        n_completed += 1
        print(f'已爬取: {n_completed}/{n_samples}\r', end='')

    # 补漏
    for gsmid in gsmids:
        if gsmid not in sample_info:
            get_gsm_sampleinfo(gsmid, sample_info)

    sample_info_df = pd.DataFrame(sample_info).T
    return sample_info_df


def get_gsm_sampleinfo(gsmid, sampleinfo_dict):
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
                        key, value = line.split(': ')
                        gsm_info[key] = value
                break
            else:
                gsm_info[td_list[0].text] = td_list[1].text
    gsm_info['gsmid'] = gsmid
    sampleinfo_dict[gsmid] = gsm_info


def download_gse(gseid, output):
    os.makedirs(f'{output}/{gseid}', exist_ok=True)
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gseid[:-3]}nnn/{gseid}/suppl'
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    a_tags = soup.find_all('a')
    for a_tag in a_tags:
        if not a_tag.text == 'Parent Directory' and not a_tag.text == 'HHS Vulnerability Disclosure':
            filename = a_tag.get('href')
            file_url = f'{url}/{filename}'
            file_response = requests.get(f'{file_url}')
            with open(f'datasets/{gseid}/{filename}', 'wb') as f:
                f.write(file_response.content)


def get_gse_info(gseid, info_dict):
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gseid}'
    soup = BeautifulSoup(requests.get(url).text, 'html.parser')
    temp_info_dict = {}
    status_tag = soup.find_all('td', string='Status')[0]
    table_tag = status_tag.parent.parent
    for tr_line in table_tag.find_all('tr', attrs={'valign': 'top'}):
        try:
            key = tr_line.find_all('td')[0].text
            value = tr_line.find_all('td')[1].text
            temp_info_dict[key] = value
        except IndexError:
            continue
    info_dict[gseid] = temp_info_dict
