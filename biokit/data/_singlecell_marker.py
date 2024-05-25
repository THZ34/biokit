# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%

def load_singlecell_marker():
    marker_dict = {'Immune': ['PTPRC'],
                   'ILC-Innate Lymphoid Cell': ['GNLY', 'IL7R'],
                   'NK': ['GNLY', 'KLRB1C'],
                   'NK1': ['NCAM1', 'FCGR3A', 'KLRC1', 'TRDC'],
                   'NK2': ['NCAM1', 'FCGR3A', 'KLRC1', 'TRDC', 'GZMA', 'GZMB', 'GZMH', 'NKG7'],
                   'ILC1': ['NCR2', 'ITAGE', 'CD160', 'KLRD1'],
                   'Intraepithelial ILC1': ['NCR2', 'ITGAE', 'IL7R-low', 'CD160', 'KLRD1'],
                   'Lamnia propria ILC1': ['NCAM1-low', 'KIT-low', 'NCR2-low', 'IL7R', 'KLRB1'],
                   'ILC2': ['KLRB1', 'IL1RL1', 'PTGDR2', 'IL7R', 'KLRB1C', 'NCR2-low', 'IL2RA', 'KIT', 'THY1'],
                   'LTi': ['CCR6', 'THY1', 'KIT', 'IL7R', 'CD4-low', 'NCAM1-low'],
                   'NCR- ILC3': ['CCR6', 'NCR2-low', 'IL2RA', 'KIT', 'IL7R', 'KLRB1'],
                   'NCR+ ILC3': ['CCR6', 'NCR2', 'NCR3', 'KIT', 'IL7R', 'KLRB1'],
                   'ILC regs': ['IL2RA', 'THY1', 'IL10', 'CD4-low', 'FOXP3-low'],
                   'NCAM1+ ILC': ['NCAM1', 'FCGR3A', 'NCR1', 'XCL1', 'KLRC1', 'GZMK'],
                   'B cells': ['CD79A', 'MS4A1'],
                   'Germinal Center B': ['RGS13'],
                   'Naive B': ['IGHD', 'IGHM', 'FCER2', 'CD38', 'CD27'],
                   'Memory B': ['CD38', 'CD27', 'IGHD'],
                   'Activate B': ['GPR183', 'TNFRSF1B'],
                   'Plasma B': ['MS4A1', 'CD38', 'MZB1', 'JCHAIN', 'IGHA1', 'IGHG1'],
                   'T cells': ['CD3E', 'CD3D'],
                   'CD8 T': ['CD8A', 'CD8B', 'GZMH'],
                   'Memory CD8 T': ['SELL', 'CCR7'],
                   'Effector memory CD8 T': ['GZMK'],
                   'Effector CD8 T': ['GZMK', 'FGFBP2', 'FCGR3A', 'SPON2', 'GNLY'],
                   'Exhausted CD8 T': ['HAVCR2', 'LAG3', 'PDCD1', 'ENTPD1', 'GNLY', 'CXCL13'],
                   'Naive CD8 T': ['IL7R', 'SELL', 'CCR7'],
                   'CD4 T': ['CD4'],
                   'Naive CD4 T': ['IL7R', 'SELL', 'CCR7'],
                   'Memory CD4 T': ['SELLlow', 'CCR7low'],
                   'CD4 Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TIGIT'],
                   'CD4 Th1 & Th2': ['BTLA', 'CD200', 'CXCL13', 'PTPN133', 'CXCR5'],
                   'CD4 Th17': ['CCL20', 'CTSH', 'KLRB1'],
                   'Monocytes': ['CD14', 'FCGR3A', 'S100A8', 'VCAN', 'FCN1', 'CCR2', 'TYROBP', 'LYZ', 'S100A9'],
                   'Macrophages': ['CD68', 'CD163', 'C1QA', 'APOC1', 'CD14', 'CD80'],
                   'M1 Macrophages': ['HLA-DRA', 'STAT1', 'ITGAX', 'CD86', 'IL10', 'MRC1'],
                   'M2 Macrophages': ['CD163', 'IL10', 'MRC1'],
                   'DCs': ['CLEC9A'],
                   'cDC1': ['CD1C', 'CLEC10A', 'FCER1A'],
                   'cDC2': ['LAMP3', 'CCR7', 'CCL19', 'CCL22'],
                   'cDC3': ['CLEC9A', 'IDO1', 'CADM1', 'CAMK2D'],
                   'pDC': ['LILRA4', 'IRF7', 'IL3RA', 'SELL'],
                   'Non-Immune': ['PTPRC-low'],
                   'Mast cell': ['KIT', 'MS4A2', 'TPSAB1', 'TPSB2'],
                   'Fibroblast': ['DCN', 'THY1', 'COL6A1', 'COL1A1', 'MMP2', 'PDGFRA'],
                   'CAF-Cancer Associated Fibroblast': ['FAP', 'TWIST1', 'WNT2', 'MMP11', 'MMP1', 'PDGFRL', 'GREM1'],
                   'Smooth muscle cells': ['ACTA2', 'TAGLN', 'MYL9', 'MYH11', 'CSPG4', 'ABCC9', 'KCNJ8', 'CNN1', 'DES',
                                           'SYNPO2'],
                   'Enteric glial cells': ['CXCL8', 'CXCL14', 'SOX10', 'S100B', 'SOX2', 'HAND2', 'PLP1'],
                   'Pericyte': ['RGS5', 'ACTA2', 'TAGLN', 'MYL9', 'CSPG4', 'ABCC9', 'KCNJ8', 'CNN1', 'DES', 'SYNPO2', ],
                   'Epithelial': ['EPCAM', 'KRT8', 'KRT18'],
                   'Epithelial ciliated': ["C1orf194", "FAM183A", "OMG", "FOXJ1"],
                   'Epithelial AT2': ['SFTPC', 'LAMP3'],
                   'Epithelial AT1': ['AGER', 'PDPN'],
                   'Epithelial Club': ['SCGB1A1'],
                   'Endothelial': ['PECAM1', 'CD34', 'VWF', 'CD34', 'ACTA2', 'MCAM', 'NID1', 'NID2',
                                   'HECW2', 'GRB10', 'CA2'],
                   'Endo arterial': ['LTC4S', 'PCSK5', 'BMX', 'HEY1'],
                   'Endo capillary': ['ACE', 'CD320', 'ENPP2', 'SLC14A1', 'TMEM88', 'CD36', 'CA4', 'PASK', 'HIGD1B',
                                      'COX4I2', 'NOTCH4', 'PDGFRB', 'RGS5', 'MYL9'],
                   'Endo venous': ['CPE', 'ADGRG6', 'MADCAM1', 'ACKR1', 'IL1R1', 'CCL14', 'CCL23'],
                   'Endo lymphatic': ['CCL21', 'LYVE1', 'MMRN1'],
                   'Endo capillary-like': ['RBP7', 'PPARG', 'FABP4', 'BTNL9', 'CD300LG'],
                   'Endo arterial-like': ['GJA4', 'GJA5', 'VEGFC'],
                   'Endo tip cells': ['ANGPT2', 'ESM1', 'APLN', 'PDGFB', 'PGF', 'CXCR4'],
                   'Endo proif': ['STMN1', 'HMGN2', 'HMGB2', 'MKI67', 'CSF3', 'HIF1A', 'S100A3', 'SELE'],
                   'Endo venous-like': ['SELP', 'IL33', 'VCAM1']
                   }
    return marker_dict


def load_celltree():
    celltree = {
        'Immune': {
            'NK cells': {
                'NK1': 0,
                'NK2': 0,
            },
            'ILC-Innate Lymphoid Cell': {
                'ILC2': 0,
                'ILC regs': 0,
                'ILC1': {
                    'Intraepithelial ILC1': 0,
                    'Lamnia propria ILC1': 0,

                },
                'ILC3': {
                    'LTi': 0,
                    'NCR- ILC3': 0,
                    'NCR+ ILC3': 0,
                }
            },
            'B cells': {
                'Germinal Center B': 0,
                'Memory B': 0,
                'Naive B': 0,
                'Activate B': 0,
                'Plasma B': 0
            },
            'T cells': {
                'CD8 T': {
                    'Memory CD8 T': 0,
                    'Effector memory CD8 T': 0,
                    'Effector CD8 T': 0,
                    'Exhausted CD8 T': 0,
                    'Naive CD8 T': 0
                },
                'CD4 T': {
                    'Naive CD4 T': 0,
                    'Memory CD4 T': 0,
                    'CD4 Treg': 0,
                    'CD4 Th1 & Th2': 0,
                    'CD4 Th17': 0
                },
            },
            'Monocytes': 0,
            'Macrophages': {
                'M1 Macrophages': 0,
                'M2 Macrophages': 0,
            },
            'DCs': {
                'cDC1': 0,
                'cDC2': 0,
                'cDC3': 0,
                'pDC': 0
            },
        },
        'Non-Immune': {
            'Fibroblast': {'CAF-Cancer Associated Fibroblast': 0},
            'Smooth muscle cells': 0,
            'Pericyte': 0,
            'Epithelial': {
                'Epithelial subtype1': 0,
                'Epithelial subtype2': 0,
                'Epithelial subtype3': 0,
                'Epithelial subtype4': 0,
                'Epithelial subtype5': 0,
                'Epithelial subtype6': 0,
                'Epithelial subtype7': 0,
                'Epithelial subtype8': 0,
                'Epithelial subtype9': 0,
                'Epithelial ciliated': 0,
                'Epithelial AT2': 0,
                'Epithelial AT1': 0,
                'Epithelial Club': 0,
            },
            'Endothelial': {
                'Endo arterial': 0,
                'Endo capillary': 0,
                'Endo venous': 0,
                'Endo lymphatic': 0,
                'Endo capillary-like': 0,
                'Endo arterial-like': 0,
                'Endo tip cells': 0,
                'Endo proif': 0,
                'Endo venous-like': 0
            },
        }
    }
    return celltree

# 记录数据来源
# ILC: 优宁维ILC (Innate Lymphoid Cell)固有免疫细胞解决方案 https://univ-shop.oss-cn-shanghai.aliyuncs.com/pc/images/upload/Image/2019060410.jpg


# %%  下载并提取cellxgene数据库
