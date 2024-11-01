# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


from scipy.special import comb
import numpy as np


def bezier_curve_S(points, n_points=1000):
    """
    生成S形贝塞尔曲线
    :param points: 起点, 控制点1, 控制点2, 终点
    :param n_points: S形曲线的点数
    :return: S形曲线的点
    """

    def bernstein_poly(i, n, t):
        return comb(n, i) * (t ** (n - i)) * (1 - t) ** i

    n = len(points) - 1
    curve = np.zeros((n_points, 2))
    for i in range(n_points):
        t = i / (n_points - 1)
        curve[i] = np.sum([point * bernstein_poly(j, n, t) for j, point in enumerate(points)], axis=0)
    return curve
