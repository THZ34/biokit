# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%


# %%
def pathway_enrichment(genes, prefix, outputdir):
    with open(f'{prefix}.txt', 'w') as f:
        f.write('\n'.join(genes))
    command = f'docker run -v "$(pwd)":/workdir -w /workdir metadocker8/msbio2 python /msbio/mylib/ms/msbio2.py /workdir/{prefix}.txt' \
              f' -o /workdir/{outputdir}/prefix -t Symbol -s -u --license /workdir/license'
    return command
