# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
import pandas as pd
import numpy as np


def grid_average(df, n_grid_x=5, n_grid_y=5):
    """
    ����άƽ��ĵ� (x, y, value) ����Ϊ n x n ���񣬲�����ÿ�������� value ��ƽ��ֵ��

    ����:
    df: pd.DataFrame, ���� 3 ��, �ֱ�Ϊ [x, y, value]
    n: int, ���񻮷ֵ������������� n=5 �򻮷�Ϊ 5x5=25 ������������

    ����:
    pd.DataFrame, ������Ϊ y_bin, ������Ϊ x_bin, ֵΪ��Ӧ�������ƽ�� value
    """
    # 1. ���� x, y ����Сֵ�����ֵ������ȷ������߽�
    min_x, max_x = df['x'].min(), df['x'].max()
    min_y, max_y = df['y'].min(), df['y'].max()

    # 2. ��������߽�
    x_bins = np.linspace(min_x, max_x, n_grid_x + 1)
    y_bins = np.linspace(min_y, max_y, n_grid_y + 1)

    # 3. ���� pd.cut �� x, y �ֵ���Ӧ�������ţ�0 ~ n-1��
    df['x_bin'] = pd.cut(df['x'], bins=x_bins, labels=False, include_lowest=True)
    df['y_bin'] = pd.cut(df['y'], bins=y_bins, labels=False, include_lowest=True)

    # 4. ���� (x_bin, y_bin) ���飬Ȼ����� value ��ƽ��ֵ
    grouped = df.groupby(['x_bin', 'y_bin'])['value'].mean().reset_index()

    # 5. �����͸��Ϊһ�� n x n �Ķ�ά��
    #    - ������: y_bin
    #    - ������: x_bin
    #    - ֵ: mean(value)
    # result = grouped.pivot(index='y_bin', columns='x_bin', values='value')

    # ���������û���κε��������У�����ʾΪ NaN���ɸ�������������
    # ����ʹ�� fillna(0) ���߱��� NaN
    # result = result.fillna(0)

    return grouped
