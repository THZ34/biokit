# biokit/cli.py
# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


import argparse
import fire
from shell import init_conda
from shell import set_fontpath
from shell import set_index_url


class BioKitCLI:
    """BioKit 命令行工具"""

    def __init__(self):
        self.conda = CondaCommands()
        self.matplotlib = MatplotlibCommands()


class CondaCommands:
    def __init__(self):
        self.init = init_conda


class MatplotlibCommands:
    def __init__(self):
        self.set_fontpath = set_fontpath


class PythonCommands:
    def __init__(self):
        self.set_index_url = set_index_url


if __name__ == "__main__":
    fire.Fire(BioKitCLI)
