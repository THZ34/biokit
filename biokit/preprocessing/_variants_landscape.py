# %%
import pandas as pd
from scipy.sparse import coo_matrix

def sort_mutation(mut_df):
    mutation_stat = (mut_df != 'no mutate').astype(bool).astype(int)
    mutation_stat['sum'] = mutation_stat.sum(1)
    mutation_stat.sort_values(by='sum', inplace=True, ascending=False)
    mutation_stat = mutation_stat.drop('sum', axis=1).sort_values(by=list(mutation_stat.index), axis=1, ascending=False)
    mut_df = mut_df.loc[mutation_stat.index][mutation_stat.columns]
    return mut_df


def read_aachange(patients, files, allow_multi_hits=False):
    """读取Vardict的AAchange.xls文件，生成突变矩阵"""

    # 读取数据
    snv = pd.DataFrame()
    if isinstance(files, pd.DataFrame):
        snv = files
    else:
        for sample, file in zip(patients, files):
            try:
                aa_df = pd.read_csv(file, sep='\t')
            except pd.errors.EmptyDataError:
                continue
            temp_snv = aa_df[['Gene', 'Effect']]
            temp_snv['sample'] = sample
            snv = pd.concat([snv, temp_snv], axis=0)
        snv.columns = ['gene', 'effect', 'sample']

    snv = snv[['gene', 'effect', 'sample']]
    no_mutate_patients = sorted(list(set(patients) - set(snv['sample'])))

    snv.drop_duplicates(inplace=True, ignore_index=True)

    # 合并多个位点突变的基因标记
    multi_marked_snv = pd.DataFrame()
    for sample in patients:
        temp_snv = snv[snv['sample'] == sample]
        temp_snv.index = temp_snv['gene']
        counts = pd.value_counts(temp_snv['gene'])
        multi_hit_genes = counts[counts > 1].index
        single_hit_genes = counts[counts == 1].index
        multi_hit_snvs = temp_snv.loc[multi_hit_genes]
        for gene in multi_hit_genes:
            effects = ','.join(multi_hit_snvs[multi_hit_snvs['gene'] == gene]['effect'])
            effects = ','.join(sorted(list(set(effects.split(',')))))
            row_indexer = multi_hit_snvs.index[multi_hit_snvs['gene'] == gene]
            multi_hit_snvs.loc[row_indexer]['effect'] = effects
            multi_hit_snvs.drop(gene, axis=0, inplace=True)
            multi_hit_snvs.loc[gene] = [gene, effects, sample]
        multi_hit_snvs.drop_duplicates(inplace=True, keep='first')
        multi_marked_snv = pd.concat([multi_marked_snv, multi_hit_snvs, temp_snv.loc[single_hit_genes]], axis=0)

    # effect包含 ',' 替换为 multi_hits
    if not allow_multi_hits:
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
    mutations = multi_marked_snv.pivot(index='gene', columns='sample', values='effect')
    mutation_stat = mutations.astype(bool).astype(int)
    mutation_stat['sum'] = mutation_stat.sum(1)
    mutation_stat.sort_values(by='sum', inplace=True, ascending=False)
    mutation_stat = mutation_stat.drop('sum', axis=1).sort_values(by=list(mutation_stat.index), axis=1, ascending=False)
    mutations = mutations.loc[mutation_stat.index][mutation_stat.columns]
    mutations.replace(variant_dict_reverse, inplace=True)
    mutations[no_mutate_patients] = 'no mutate'
    return mutations
