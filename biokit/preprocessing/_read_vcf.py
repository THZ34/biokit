import pandas as pd


def read_vcf(filename):
    header = -2
    with open(filename) as file:
        while True:
            header += 1
            if not file.readline().startswith('#'):
                break
    return pd.read_csv(filename, header=header, sep='\t', index_col=None)
