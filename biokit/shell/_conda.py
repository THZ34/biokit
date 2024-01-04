# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%

def init_conda():
    command = """
conda init --
# 配置国内镜像源
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/main/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/msys2/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/menpo/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/

conda config --add channels https://mirrors.aliyun.com/pypi/simple/
conda config --add channels http://pypi.douban.com/simple/ 
conda config --add channels http://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/win-64/

# 在安装时总是显示包源
conda config --set show_channel_urls yes
conda config --set always_yes True
    """