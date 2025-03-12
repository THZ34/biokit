# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

# %%
import os
import shutil


# %%

def set_index_url():
    os.makedirs('~/.pip', exist_ok=True)
    if os.path.exists('~/.pip/pip.conf'):
        shutil.copy('~/.pip/pip.conf', '~/.pip/pip.conf.bak')
    with open('~/.pip/pip.conf', 'w') as f:
        content = """[global] 
index-url = https://pypi.tuna.tsinghua.edu.cn/simple/"""
        f.write(content)
