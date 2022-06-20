# %%
import pandas as pd

# %%
from scipy.sparse import coo_matrix


def sparse_mutation(snv_df):
    """稀疏矩阵转换密集矩阵"""
    mutation_df = snv_df.copy()
    genes = sorted(list(set(mutation_df['gene'])))
    gene_dict = dict(zip(genes, range(len(genes))))
    samples = sorted(list(set(mutation_df['sample'])))
    sample_dict = dict(zip(samples, range(len(samples))))
    mutation_df['gene'].replace(gene_dict, inplace=True)
    mutation_df['sample'].replace(sample_dict, inplace=True)
    mutation_sparse = coo_matrix((mutation_df['effect'], (mutation_df['sample'], mutation_df['gene'])))
    mutation_df = pd.DataFrame(mutation_sparse.todense(), index=samples, columns=genes)
    return mutation_df


def read_aachange(patients, files):
    """读取Vardict的AAchange.xls文件，生成突变矩阵"""

    # 读取数据
    snv = pd.DataFrame()
    no_mutate_patients = []
    if isinstance(files, pd.DataFrame):
        snv = files
    else:
        for patient, file in zip(patients, files):
            try:
                aa_df = pd.read_csv(file, sep='\t')
            except pd.errors.EmptyDataError:
                no_mutate_patients.append(patient)
                continue

            temp_snv = aa_df[['Gene', 'Effect']]
            temp_snv['sample'] = patient
            snv = pd.concat([snv, temp_snv], axis=0)
        snv.columns = ['gene', 'effect', 'sample']

    # 识别多个位点突变的基因标记为Multi_Hits
    multi_marked_snv = pd.DataFrame()
    for patient in patients:
        temp_snv = snv[snv['sample'] == patient]
        temp_snv.index = temp_snv['gene']
        counts = pd.value_counts(temp_snv['gene'])
        multi_hit_genes = counts[counts > 1].index
        single_hit_genes = counts[counts == 1].index
        multi_hit_snvs = temp_snv.loc[multi_hit_genes]
        multi_hit_snvs['effect'] = 'multi_hits'
        multi_hit_snvs.drop_duplicates(inplace=True, keep='first')
        multi_marked_snv = pd.concat([multi_marked_snv, multi_hit_snvs, temp_snv.loc[single_hit_genes]], axis=0)

    # effect包含 ',' 替换为 multi_hits
    multi_marked_snv.index = range(multi_marked_snv.shape[0])
    indexer = multi_marked_snv.index[multi_marked_snv['effect'].str.contains(',')]
    multi_marked_snv.loc[indexer, 'effect'] = 'multi_hits'
    variants = multi_marked_snv['effect'].unique()
    variants.sort()

    # 稀疏矩阵转换为密集矩阵
    variant_dict = dict(zip(variants, range(1, len(variants) + 1)))
    variant_dict_reverse = dict(zip(range(1, len(variants) + 1), variants))
    variant_dict_reverse[0] = 'no mutate'

    multi_marked_snv['effect'].replace(variant_dict, inplace=True)
    mutations = sparse_mutation(multi_marked_snv)
    mutations = mutations.T
    # mutations.replace({0: 'no_mutate'}, inplace=True)
    mutation_stat = mutations.astype(bool).astype(int)
    mutation_stat['sum'] = mutation_stat.sum(1)
    mutation_stat.sort_values(by='sum', inplace=True, ascending=False)
    mutation_stat = mutation_stat.drop('sum', axis=1).sort_values(by=list(mutation_stat.index), axis=1, ascending=False)
    mutations = mutations.loc[mutation_stat.index][mutation_stat.columns]
    mutations.replace(variant_dict_reverse, inplace=True)
    mutations[no_mutate_patients] = 'no mutate'
    return mutations
