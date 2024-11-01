# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


from scipy.special import comb
import numpy as np


def bezier_curve_S(points, n_points=1000):
    """
    ����S�α���������
    :param points: ���, ���Ƶ�1, ���Ƶ�2, �յ�
    :param n_points: S�����ߵĵ���
    :return: S�����ߵĵ�
    """

    def bernstein_poly(i, n, t):
        return comb(n, i) * (t ** (n - i)) * (1 - t) ** i

    n = len(points) - 1
    curve = np.zeros((n_points, 2))
    for i in range(n_points):
        t = i / (n_points - 1)
        curve[i] = np.sum([point * bernstein_poly(j, n, t) for j, point in enumerate(points)], axis=0)
    return curve
