# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %% 
import getpass
from copy import deepcopy

import json5
from matplotlib import rcParams

config = {"font.family": 'Microsoft YaHei', 'pdf.fonttype': 42}
rcParams.update(config)


# %%
def refactor_json(seq_df, template_json, outdir, platform='MGI', project='CD14000', queue='prj.q', owner=None):
    """

    :param seq_df: 下机表
    :param template_json:
    :param tumor_type:
    :param owner:
    :param platform:
    :param outdir:
    :param sample_info:
    :param oss_path:
    :return:
    """
    if not owner:
        owner = getpass.getuser()
    template_json = json5.load(open(template_json))
    product_head = list(template_json.keys())[0].split('.')[0]

    # 检查数据列
    required_cols = ['run_path', 'i7', 'gender', 'tumorid', 'tumor_type', 'tumor_read1', 'tumor_read2']
    lost_cols = sorted(list(set(required_cols) - set(seq_df.columns)))
    if len(lost_cols) > 0:
        raise ValueError(f'seq_df中缺少 {",".join(lost_cols)}')

    # DNA
    if 'tumorid' in seq_df.columns:
        # 单样本DNA
        if not 'normalid' in seq_df.columns:
            for sample in seq_df.index:
                target_json = deepcopy(template_json)
                run_path, gender, tumorid, tumor_type, tumor_read1, tumor_read2, product_id = seq_df.loc[sample][
                    ['run_path', 'gender', 'tumorid', 'tumor_type', 'tumor_read1', 'tumor_read2', 'product']]

                target_json[f'{product_head}.prefix'] = f'{tumorid}'
                target_json[f'{product_head}.queue'] = queue
                target_json[f'{product_head}.project'] = project
                target_json[f'{product_head}.patient_id'] = f'{tumorid[2:9]}'
                target_json[f'{product_head}.project_id'] = f'{tumorid}'
                target_json[f'{product_head}.product_id'] = product_id
                target_json[f'{product_head}.tumor_name'] = tumor_type
                target_json[f'{product_head}.owner'] = owner
                target_json[f'{product_head}.platform'] = platform

                target_json[f'{product_head}.gender'] = gender
                target_json[f'{product_head}.tumor_id'] = tumorid
                target_json[f'{product_head}.tumor_read1'] = [tumor_read1]
                target_json[f'{product_head}.tumor_read2'] = [tumor_read2]

                with open(f'{outdir}/{tumorid}.json5', 'w') as f:
                    f.write(json5.dumps(target_json, indent=4, separators=(',', ':')))
        # 配对样本DNA
        else:
            for sample in seq_df.index:
                target_json = deepcopy(template_json)
                run_path, gender, tumorid, tumor_type, tumor_read1, tumor_read2, product_id, normalid, normal_read1, normal_read2 = \
                    seq_df.loc[sample][
                        ['run_path', 'gender', 'tumorid', 'tumor_type', 'tumor_read1', 'tumor_read2', 'product',
                         'normalid', 'normal_read1', 'normal_read2']]
                target_json[f'{product_head}.prefix'] = f'{tumorid}'
                target_json[f'{product_head}.queue'] = queue
                target_json[f'{product_head}.project'] = project
                target_json[f'{product_head}.patient_id'] = f'{tumorid[2:9]}'
                target_json[f'{product_head}.project_id'] = f'{tumorid}-VS-{normalid}'
                target_json[f'{product_head}.product_id'] = product_id
                target_json[f'{product_head}.tumor_name'] = tumor_type
                target_json[f'{product_head}.owner'] = owner
                target_json[f'{product_head}.platform'] = platform

                target_json[f'{product_head}.gender'] = gender
                target_json[f'{product_head}.tumor_id'] = tumorid
                target_json[f'{product_head}.tumor_read1'] = [tumor_read1]
                target_json[f'{product_head}.tumor_read2'] = [tumor_read2]
                target_json[f'{product_head}.normal_id'] = normalid
                target_json[f'{product_head}.normal_read1'] = [normal_read1]
                target_json[f'{product_head}.normal_read2'] = [normal_read2]

                with open(f'{outdir}/{tumorid}-VS-{normalid}.json5', 'w') as f:
                    f.write(json5.dumps(target_json, indent=4, separators=(',', ':')))

    elif 'rnaid' in seq_df.columns:
        pass
