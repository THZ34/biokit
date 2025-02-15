# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
from typing import Union

import matplotlib.colors as mcolors


def complementary_color(color: Union[str, tuple]) -> str:
    """
    获取互补色
    :param color: 16进制颜色或RGB颜色
    :return: color的互补色
    """
    if isinstance(color, str):
        rgb = mcolors.hex2color(color)
    cmy = [(1 - x) for x in rgb]
    return mcolors.rgb2hex(cmy)
