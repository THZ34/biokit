# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%

def load_singlecell_marker():
    marker_dict = {'Immune': ['PTPRC'],
                   'B-cell': ['CD79A', 'MS4A1'],
                   'Germinal center B': ['RGS13'],
                   'Naive B': ["IGHD", "IGHM", "FCER2", "CD38", "CD27"],
                   'Memory B': ["CD38", "CD27", "IGHD"],
                   'Activate B': ["GPR183", "TNFRSF1B"],
                   'Plasma cell': ['MS4A1', "CD38", "MZB1", "JCHAIN", "IGHA1", "IGHG1"],
                   'T-cell': ['CD3E', 'CD3D'],
                   'CD8 T': ["CD8A", "CD8B", "GZMH"],
                   'Memory CD8T': ["SELL", "CCR7"],
                   'effector memory CD8T': ['GZMK'],
                   'effector CD8T': ["GZMK", "FGFBP2", "FCGR3A", "SPON2", "GNLY"],
                   'exhausted CD8T': ["HAVCR2", "LAG3", "PDCD1", "ENTPD1", "GNLY", "CXCL13"],
                   'CD4 T': ['CD4'],
                   "Naive CD4T/CD8T": ["IL7R", "SELL", "CCR7"],
                   "Memory CD4T": ["SELLlow", "CCR7low"],
                   "CD4 Treg": ["FOXP3", "IL2RA", "CTLA4", "TIGIT"],
                   "CD4 Th1 & Th2": ["BTLA", "CD200", "CXCL13", "PTPN133", "CXCR5"],
                   "CD4 Th17": ["CCL20", "CTSH", "KLRB1"],
                   'NK': ['GNLY'],
                   'NK1': ["NCAM1", "FCGR3A", "KLRC1", "TRDC"],
                   'NK2': ["NCAM1", "FCGR3A", "KLRC1", "TRDC", "GZMA", "GZMB", "GZMH", "NKG7"],
                   'Monocyte': ["CD14", "FCGR3A", "S100A8", "VCAN", "FCN1"],
                   'Macrophage': ["CD68", "CD163", "C1QA", "APOC1"],
                   'cDC1': ['CD1C', 'CLEC10A', 'FCER1A'],
                   'cDC2': ['LAMP3', 'CCR7', 'CCL19', 'CCL22'],
                   'cDC3': ['CLEC9A', 'IDO1', 'CADM1', 'CAMK2D'],
                   'pDC': ['LILRA4', 'IRF7', 'IL3RA', 'SELL'],
                   'Mast cell': ['KIT', 'MS4A2', 'TPSAB1', 'TPSB2'],
                   'Fibroblast': ['DCN', 'THY1', 'COL6A1', 'COL1A1'],
                   'smooth muscle cells': ['TAGLN', 'ACTA2', 'ACTG2', 'MYH11', 'GNN1'],
                   'Pericyte': ['TAGLN', 'ACTA2', 'RGS5', 'PDGFRB'],
                   'Epithelial': ['EPCAM', 'KRT8', 'KRT18'],
                   'Endothelial': ['PECAM1', 'CD34', 'VWF']
                   }
    return marker_dict


def load_celltree():
    celltree = cell_tree_dict = {
        'Immune': {
            'B cells': {
                'Germinal Center B': 0,
                'Memory B': 0,
                'Naive B': 0,
                'Activate B': 0,
                'Plasma B': 0},
            'T cells': {
                'CD8 T': {
                    'Memory CD8 T': 0,
                    'effector memory CD8 T': 0,
                    'effector CD8 T': 0,
                    'exhausted CD8 T': 0,
                    "Naive CD8 T": 0},
                'CD4 T': {
                    "Naive CD4 T": 0,
                    "Memory CD4 T": 0,
                    "CD4 Treg": 0,
                    "CD4 Th1 & Th2": 0,
                    "CD4 Th17": 0}, },
            'NK cells': {
                'NK1': 0,
                'NK2': 0, },
            'Monocytes': 0,
            'Macrophages': 0,
            'DCs': {
                'cDC1': 0,
                'CDC2': 0,
                'cDC3': 0,
                'pDC': 0},
        },
        'Non-Immune': {
            'Fibroblast': 0,
            'smooth muscle cells': 0,
            'Pericyte': 0,
            'Epithelial': 0,
            'Endothelial': 0}
    }
