# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import pandas as pd
import numpy as np


def grid_average(df, n_grid_x=5, n_grid_y=5):
    """
    将二维平面的点 (x, y, value) 划分为 n x n 网格，并计算每个网格中 value 的平均值。

    参数:
    df: pd.DataFrame, 包含 3 列, 分别为 [x, y, value]
    n: int, 网格划分的行列数，例如 n=5 则划分为 5x5=25 个网格子区域

    返回:
    pd.DataFrame, 行索引为 y_bin, 列索引为 x_bin, 值为对应子网格的平均 value
    """
    # 1. 计算 x, y 的最小值和最大值，用于确定网格边界
    min_x, max_x = df['x'].min(), df['x'].max()
    min_y, max_y = df['y'].min(), df['y'].max()

    # 2. 生成网格边界
    x_bins = np.linspace(min_x, max_x, n_grid_x + 1)
    y_bins = np.linspace(min_y, max_y, n_grid_y + 1)

    # 3. 利用 pd.cut 将 x, y 分到对应的网格编号（0 ~ n-1）
    df['x_bin'] = pd.cut(df['x'], bins=x_bins, labels=False, include_lowest=True)
    df['y_bin'] = pd.cut(df['y'], bins=y_bins, labels=False, include_lowest=True)

    # 4. 按照 (x_bin, y_bin) 分组，然后计算 value 的平均值
    grouped = df.groupby(['x_bin', 'y_bin'])['value'].mean().reset_index()

    # 5. 将结果透视为一个 n x n 的二维表
    #    - 行索引: y_bin
    #    - 列索引: x_bin
    #    - 值: mean(value)
    # result = grouped.pivot(index='y_bin', columns='x_bin', values='value')

    # 如果有网格没有任何点落在其中，会显示为 NaN，可根据需求进行填充
    # 例如使用 fillna(0) 或者保留 NaN
    # result = result.fillna(0)

    return grouped
