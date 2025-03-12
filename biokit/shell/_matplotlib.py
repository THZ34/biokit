# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


import subprocess
import os
import shutil
from pathlib import Path


def get_matplotlibrc_path():
    result = subprocess.run(
        'conda info --envs',
        capture_output=True,
        text=True,
        shell=True  # Windows 需要设置 shell=True
    )
    envs = [i.split()[0] for i in result.stdout.split('\n') if 'envs' in i]

    conda_path = os.path.dirname(
        os.path.dirname(subprocess.run('which conda', capture_output=True, text=True, shell=True).stdout.strip()))
    envs = [os.path.join(conda_path, 'envs', i) for i in envs]

    pythons = []
    for env in envs:
        python_dir = [i for i in os.listdir(f'{env}/lib') if i.startswith('python')][0]
        pythons.append(python_dir)

    matplotlibrc_paths = [f'{env}/lib/{python}/site-packages/matplotlib/mpl-data/matplotlibrc' for env, python in
                          zip(envs, pythons)]

    return matplotlibrc_paths


def modify_matplotlibrc(config_path, font_name="Microsoft YaHei", font_path=None):
    """
    修改 matplotlib 配置文件以设置字体

    参数：
    - config_path: matplotlibrc 文件路径
    - font_name: 要使用的主字体名称（如 "SimHei"）
    - font_path: 自定义字体目录路径（可选）
    """
    # 1. 备份原始文件
    backup_path = config_path + ".bak"
    if not os.path.exists(backup_path):
        shutil.copyfile(config_path, backup_path)
        print(f"已创建备份文件: {backup_path}")

    # 2. 读取配置文件内容
    with open(config_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 3. 定义要修改的配置项
    target_configs = {
        "font.family": "sans-serif",  # 设置为 sans-serif/serif 等通用族
        "font.sans-serif": f"{font_name}, DejaVu Sans, Arial",  # 添加中文字体到首位
        "pdf.fonttype": 42,  # 优化 PDF 字体显示效果
    }

    # 可选：添加自定义字体路径
    if font_path:
        target_configs["font.dir"] = f"{font_path}:..."

    # 4. 执行修改
    modified = False
    new_lines = []
    for line in lines:
        line_strip = line.strip()

        # 匹配目标配置项
        for key, value in target_configs.items():
            if line_strip.startswith(f"{key} :") or line_strip.startswith(f"# {key} :"):
                new_line = f"{key} : {value}\n"
                if new_line != line:
                    new_lines.append(new_line)
                    modified = True
                break
        else:
            new_lines.append(line)

    # 5. 写入修改后的文件
    if modified:
        with open(config_path, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        print("配置文件修改成功！")
    else:
        print("未检测到需要修改的配置项")


def set_fontpath(fontpath='/usr/share/fonts/truetype/msttcorefonts'):
    config_paths = get_matplotlibrc_path()
    for config_path in config_paths:
        modify_matplotlibrc(config_path, font_path=fontpath)
