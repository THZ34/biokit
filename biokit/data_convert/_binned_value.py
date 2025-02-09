# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import numpy as np


def binned_value(x, bins=50):
    """
    将x中的值按照bins中的值进行分箱
    """
    vmin = np.min(x)
    vmax = np.max(x)
    interval_width = (vmax - vmin) / bins
    var_exp_binned = np.digitize(x, np.arange(vmin, vmax, interval_width))
    var_exp_binned = (var_exp_binned / bins) * (vmax - vmin) + vmin
    return var_exp_binned
