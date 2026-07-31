import tqdm
from cellphonedb.utils import db_utils


def extract_ligand_receptor_genes(cpi_ids, cpdb_file_path):
    interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name, receptor2tfs = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)
    lr_dict = {}
    for cpi_id in tqdm.tqdm(cpi_ids):
        interaction_row = interactions[interactions['id_cp_interaction'] == cpi_id].iloc[0]
        prot1 = interaction_row['id_multidata_1']
        prot2 = interaction_row['id_multidata_2']
        iscomplex1 = interaction_row['is_complex_1']
        iscomplex2 = interaction_row['is_complex_2']

        if iscomplex1:
            complex_df = complex_composition.set_index('complex_multidata_id').loc[[prot1]]
            protids = complex_df['protein_multidata_id'].values.tolist()
            genenames1 = genes[genes['id_multidata'].isin(protids)]['gene_name'].values.tolist()
        else:
            genenames1 = [genes[genes['id_multidata'] == prot1]['gene_name'].values[0]]
        if iscomplex2:
            complex_df = complex_composition.set_index('complex_multidata_id').loc[[prot2]]
            protids = complex_df['protein_multidata_id'].values.tolist()
            genenames2 = genes[genes['id_multidata'].isin(protids)]['gene_name'].values.tolist()
        else:
            genenames2 = [genes[genes['id_multidata'] == prot2]['gene_name'].values[0]]

        lr_list = []
        for gene1 in genenames1:
            for gene2 in genenames2:
                lr_list.append([gene1, gene2])
        lr_dict[cpi_id] = lr_list
    return lr_dict
