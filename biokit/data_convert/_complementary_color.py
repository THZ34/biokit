# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import matplotlib.colors as mcolors


def complementary_color(color: str | tuple) -> str:
    """
    ��ȡ����ɫ
    :param color: 16������ɫ��RGB��ɫ
    :return: color�Ļ���ɫ
    """
    if isinstance(color, str):
        rgb = mcolors.hex2color(color)
    cmy = [(1 - x) for x in rgb]
    return mcolors.rgb2hex(cmy)
