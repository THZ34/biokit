# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%

def extract_cell_types_recursive(tree_dict, level, current_level=1, prefix=''):
    cell_types = []
    for cell_type, sub_dict in tree_dict.items():
        full_cell_type = f"{prefix}{cell_type}"
        if isinstance(sub_dict, dict):
            if current_level < level:
                sub_cell_types = extract_cell_types_recursive(sub_dict, level, current_level + 1,
                                                              prefix=f"{full_cell_type}.")
                cell_types.extend(sub_cell_types)
            elif current_level == level:
                cell_types.append(full_cell_type)
        else:
            cell_types.append(full_cell_type)
    return cell_types


def get_max_level(cell_tree):
    max_level = 1

    def traverse(node, level):
        nonlocal max_level
        if level > max_level:
            max_level = level
        if isinstance(node, dict):
            for child in node.values():
                traverse(child, level + 1)

    traverse(cell_tree, 1)
    return max_level - 1


def anno_level_cells(adata, cutoff_dict, cell_tree, marker_dict, all_cluster_marker_gene_mean, cluster_key='leiden',
                     level=1):
    cell_chains = extract_cell_types_recursive(cell_tree, level=level)
    cell_chains = [i.split('.') for i in cell_chains]
    cell_chains = [i for i in cell_chains if len(i) == level]
    colname = f'level{level} cellanno'
    for cells in cell_chains:
        if level > 1:
            # 根据父类筛选本级需要注释的细胞
            parent_celltype = cells[level - 2]
            clusters = adata.obs.loc[adata.obs[f'level{level - 1} cellanno'] == parent_celltype, cluster_key].unique()
        else:
            clusters = adata.obs[cluster_key].unique()
        cell_type = cells[-1]
        markers = marker_dict[cell_type]
        for cluster in clusters:
            if all_cluster_marker_gene_mean.loc[cluster, markers].values[0] > cutoff_dict[level]:
                adata.obs.loc[adata.obs[cluster_key] == cluster, colname] = cell_type


def anno_cells(adata, cutoff_dict, cell_tree, marker_dict, all_cluster_marker_gene_ratio, cluster_key='leiden'):
    max_level = get_max_level(cell_tree)
    for level in range(1, max_level + 1):
        adata.obs[f'level{level} cellanno'] = None
        anno_level_cells(adata, cutoff_dict, cell_tree, marker_dict, all_cluster_marker_gene_ratio,
                         cluster_key, level=level)
        # 任意一级没有注释到子类则沿用上一级
        if level == 1:
            adata.obs[f'level{level} cellanno'].fillna('Non-Immune', inplace=True)
        if level > 1:
            adata.obs[f'level{level} cellanno'] = np.where(adata.obs[f'level{level} cellanno'].isna(),
                                                           adata.obs[f'level{level - 1} cellanno'],
                                                           adata.obs[f'level{level} cellanno'])


def get_all_subclasses(celltree, parent):
    subclasses = []
    if parent in celltree:
        children = celltree[parent]
        if isinstance(children, dict):
            for child in children:
                subclasses.append(child)
                subclasses.extend(get_all_subclasses(children, child))
    for key in celltree:
        if isinstance(celltree[key], dict):
            subclasses.extend(get_all_subclasses(celltree[key], parent))
    return subclasses

