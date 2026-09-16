# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import os
import json
import subprocess
from urllib.request import urlopen
from biokit.data import load_metascape


# %%
def detect_cytoscape(port=1234, timeout=2):
    try:
        with urlopen(f'http://127.0.0.1:{port}/v1/', timeout=timeout) as response:
            status = json.loads(response.read().decode())
    except Exception:
        return False
    return response.status == 200 and status.get('apiVersion') == 'v1' and 'memoryStatus' in status


# %%
def pathway_enrichment(genes, prefix, outputdir, license=None, run=False, cytoscape_port=1234):
    os.makedirs(outputdir, exist_ok=True)
    input_file = f'{outputdir}/{prefix}.txt'
    result_dir = f'{outputdir}/{prefix}'
    with open(input_file, 'w') as f:
        f.write('\n'.join(genes))
    if license is None:
        license = '$(pwd)/license'
    # docker_options = ''
    # msbio_options = ''
    # if detect_cytoscape(cytoscape_port):
    docker_options = '--network host '
    msbio_options = f'--cytoport {cytoscape_port} '
    command = (f'docker run --rm {docker_options}'
               f'-u "$(id -u)" '
               f'-v "$(pwd)":/workdir '
               f'-v "{license}":/workdir/license '
               f'-w /workdir '
               f'metadocker8/msbio2 '
               f'python /msbio/mylib/ms/msbio2.py '
               f'"/workdir/{input_file}" '
               f'-o "/workdir/{result_dir}" '
               f'-t Symbol -s -u {msbio_options}'
               f'--license /workdir/license')
    if run:
        subprocess.run(command, shell=True, check=True)
        return load_metascape(f'{result_dir}/metascape_result.xlsx')
    return command
