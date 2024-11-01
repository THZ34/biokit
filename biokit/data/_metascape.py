# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import pandas as pd


def load_metascape(file):
    """

    :param file: xlsx�ļ���zip�ļ�·��
    :return:
    """
    if isinstance(file, str):
        if file.endswith('.xlsx'):
            df = load_metascape_xlsx(file)
        elif file.endswith('.zip'):
            df = load_metascape_zip(file)
    df = df.loc[df['GroupID'].str.contains('Member')]
    df[['Input number', 'Background number']] = df['InTerm_InList'].str.split('/', expand=True).astype(int)
    df['Ratio'] = df['Input number'] / df['Background number']
    df.rename({'Category': 'Database', 'Term': 'PathwayID', 'Description': 'PathwayName', 'LogP': 'log10(pvalue)',
               'Log(q-value)': 'log10(padj)'}, axis=1, inplace=True)
    df['-log10(padj)'] = -df['log10(padj)']
    df['-log10(pvalue)'] = -df['log10(pvalue)']
    return df


def load_metascape_zip(file):
    if isinstance(file, str):
        if file.endswith('.zip'):
            import zipfile
            return load_metascape_xlsx(zipfile.ZipFile(file).extract('metascape_result.xlsx'))


def load_metascape_xlsx(file):
    return pd.read_excel(file, sheet_name='Enrichment')
