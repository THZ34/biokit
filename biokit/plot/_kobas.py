# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def kobas_barplot(file):
    pathway_df = pd.read_excel(file)
    if 'select' in pathway_df.columns:
        pathway_df = pathway_df[pathway_df['select'] == 1]
    # replace cluster
    clusters = pathway_df['Cluster'].dropna().unique()
    clusters.sort()
    cluster_replace_dict = dict(zip(clusters, range(1, len(clusters) + 1)))
    pathway_df['Cluster'].replace(cluster_replace_dict, inplace=True)
    pathway_df['Cluster'].fillna('Other', inplace=True)
    color_dict = dict(zip(range(1, len(clusters) + 1), sns.husl_palette(n_colors=len(clusters), l=0.75, s=1)))
    color_dict.update({'Other': 'silver'})
    clusters = range(1, len(clusters) + 1)

    yticklabels = []
    fig, ax = plt.subplots(figsize=(8, 4))
    i = pathway_df.shape[0]
    for cluster in clusters:
        temp_pathway_df = pathway_df[pathway_df['Cluster'] == cluster].copy()
        temp_pathway_df.sort_values('Enrich_ratio', ascending=False, inplace=True)
        yticklabels.extend(temp_pathway_df['Term'].to_list())
        ax.barh(y=range(i, i - temp_pathway_df.shape[0], -1), width=temp_pathway_df['Enrich_ratio'],
                color=color_dict[cluster], height=0.3, left=0, label=cluster)
        i -= temp_pathway_df.shape[0]

    cluster = 'Other'
    temp_pathway_df = pathway_df[pathway_df['Cluster'] == cluster].copy()
    temp_pathway_df.sort_values('Enrich_ratio', ascending=False, inplace=True)
    yticklabels.extend(temp_pathway_df['Term'].to_list())
    ax.barh(y=range(i, i - temp_pathway_df.shape[0], -1), width=temp_pathway_df['Enrich_ratio'],
            color=color_dict[cluster], height=0.3, left=0, label=cluster)
    i -= temp_pathway_df.shape[0]

    yticklabels.reverse()
    ax.set_yticks(range(1, pathway_df.shape[0] + 1))
    ax.set_yticklabels(yticklabels)
    ax.xaxis.grid(True)
    ax.set_xlabel('Enrich ratio')
    ax.legend(bbox_to_anchor=(1, 0.5), loc='center left', title='Cluster')
    ax.set_title(f'MS {file.split("/")[-1].split(".")[0]} pathway')
    plt.subplots_adjust(0.5, 0.11, 0.85, 0.9)
    return pathway_df, ax
