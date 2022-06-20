import pandas as pd
from scipy.sparse import coo_matrix


def vcf_to_mutation(vcfs, filter=None):
    if filter is None:
        filter = {
            'ExonicFunc.refGene': ['nonsynonymous SNV', 'frameshift deletion', 'frameshift insertion', 'stopgain',
                                   'stoploss']}
    merge_df = pd.DataFrame()
    mut_type_counts = pd.DataFrame()
    for vcf in vcfs:
        sample = vcf.split('.')[0]
        df = pd.read_csv(vcf)
        # 过滤
        for col in filter:
            df = df[[(i in filter[col]) for i in df[col]]]
        # 突变类型
        df['muts'] = df['Ref'] + '>' + df['Alt']
        # 统计突变数量
        mut_type_count = pd.value_counts(df['muts'])
        mut_type_count = pd.DataFrame([mut_type_count.to_list(), [sample.split('/')[-1]] * len(mut_type_count)],
                                      columns=mut_type_count.index).T
        mut_type_count['MutType'] = mut_type_count.index
        mut_type_count.index = mut_type_count[1]
        mut_type_count.index.name = 'Sample'
        mut_type_count.drop(1, axis=1, inplace=True)
        mut_type_count.columns = ['Num', 'MutType']
        mut_type_count = mut_type_count[['MutType', 'Num']]
        mut_type_counts = pd.concat([mut_type_counts, mut_type_count], axis=0)
        #
        df = df[['Gene.refGene', 'ExonicFunc.refGene']]
        df.columns = ['gene', 'type']
        df.drop_duplicates(inplace=True)
        df['type'] = df['type'] + ';'
        df['type'] = df['type'].str.replace(' ', '_')
        df.index = df['gene']
        df.index.name = None
        gene_count = pd.value_counts(df.index)
        dup_genes = gene_count[gene_count > 1].index
        nodup_genes = gene_count[gene_count == 1].index
        df = pd.concat([df.loc[nodup_genes], df.loc[dup_genes].groupby('gene').sum()], axis=0, sort=True)
        df['gene'] = df.index
        df.index = range(df.shape[0])
        df['sample'] = sample.split('/')[-1]
        merge_df = pd.concat([merge_df, df], axis=0)

    merge_df['type'] = list(map(lambda x: ';'.join(sorted(list(set(x.split(';'))))).lstrip(';'), merge_df['type']))

    samples = sorted(list(set(merge_df['sample'])))
    genes = sorted(list(set(merge_df['gene'])))
    values = sorted(list(set(merge_df['type'])))
    sample_dict = dict(zip(samples, range(len(samples))))
    gene_dict = dict(zip(genes, range(len(genes))))
    value_dict = dict(zip(values, range(1, len(values) + 1)))
    for map_dict in [sample_dict, gene_dict, value_dict]:
        merge_df.replace(map_dict, inplace=True)

    sparse_df = coo_matrix((merge_df['type'], (merge_df['gene'], merge_df['sample'])))
    sparse_csr_df = sparse_df.tocsr()
    mutation_df = pd.DataFrame(sparse_df.todense(), index=genes, columns=samples)
    reverse_value_dict = dict(zip(value_dict.values(), [i.strip(';') for i in value_dict.keys()]))
    mutation_df.replace(reverse_value_dict, inplace=True)
    mutation_df.replace({0: None}, inplace=True)

    tmb_df = pd.read_csv('sstmb_correlation.csv', index_col=0)
    sample_info = pd.read_excel('YConeSelectSamplesInfo.xlsx')
    sample_info.index = sample_info['SampleID']
    sample_info = sample_info.loc[tmb_df.index]
    sample_info = pd.concat(
        [tmb_df, sample_info.drop(['Normal', 'TMB', 'Workdir', 'Time', 'QCbed', 'SampleID'], axis=1)],
        axis=1)
    sample_info.index.name = 'Sample'
