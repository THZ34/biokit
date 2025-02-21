# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import os


# %%
def pathway_enrichment(genes, prefix, outputdir, license=None):
    os.makedirs(outputdir, exist_ok=True)
    with open(f'{outputdir}/{prefix}.txt', 'w') as f:
        f.write('\n'.join(genes))
    if license is None:
        license = "$(pwd)/license"
    command = (f'docker run -u "$(id -u)" '
               f'-v "$(pwd)":/workdir ',
               f'-v "{license}":/workdir/license ',
               f'-w /workdir metadocker8/msbio2 python /msbio/mylib/ms/msbio2.py "/workdir/{outputdir}/{prefix}.txt" ',
               f'-o "/workdir/{outputdir}/{prefix}" -t Symbol -s -u --license /workdir/license')
    return command
