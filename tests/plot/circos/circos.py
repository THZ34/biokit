import pandas as pd
from biokit.plot import Circos
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from biokit.data import load_hg19_ref

# %%
ref_df = load_hg19_ref()
fusion_df = pd.read_excel('tests/plot/circos/27基因检测范围选择依据.xlsx')
fusion_df = fusion_df[fusion_df['融合形式'] != 'Fusion']
fusion_df = fusion_df[fusion_df['融合形式'].str.contains('-')]

bezier_df = []
for fusion in fusion_df['融合形式']:
    gene1, gene2 = fusion.strip().split('-')
    gene1 = gene1.strip()
    gene2 = gene2.strip()
    chr1, start1, end1 = ref_df[ref_df['name2'] == gene1][['chrom', 'txStart', 'txEnd']].to_numpy()[0]
    chr2, start2, end2 = ref_df[ref_df['name2'] == gene2][['chrom', 'txStart', 'txEnd']].to_numpy()[0]
    bezier_df.append([chr1, start1, end1, chr2, start2, end2, gene1, gene2])

bezier_df = pd.DataFrame(bezier_df, columns=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'gene1', 'gene2'])
bezier_df['type'] = fusion_df['证据等级'].to_list()
bezier_df.drop(14, axis=0, inplace=True)

# %%
circos = Circos(bottom=10)

ring = circos.add_ring()
bezier_layer = circos.add_layer(data=bezier_df, value_col='type', kind='bezier', ring=ring, linewidth=8,
                                linepoints=100)

ring = circos.add_ring()
text_df1 = bezier_df[['chr1', 'start1', 'end1', 'gene1']].drop_duplicates(ignore_index=True)
text_df1.columns = ['chr', 'start', 'end', 'gene']
gene1_layer = circos.add_layer(data=text_df1, value_col='gene', kind='text', ring=ring, fontsize=5)

ring = circos.add_ring()
text_df2 = bezier_df[['chr2', 'start2', 'end2', 'gene2']].drop_duplicates(ignore_index=True)
text_df2.columns = ['chr', 'start', 'end', 'gene']
gene2_layer = circos.add_layer(data=text_df2, value_col='gene', kind='text', ring=ring, fontsize=5)

circos.show()

pdf = PdfPages('circos.pdf')
pdf.savefig()
plt.close()
pdf.close()
