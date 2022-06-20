import pandas as pd


def group4(df, var1, var2, border1, border2) -> pd.DataFrame:
    """按两个变量分为4组

    :param df: pd.DataFrame
    :param var1: 变量1
    :param var2: 变量2
    :param border1: 阈值1
    :param border2: 阈值2
    :return: pd.DataFrame
    """
    matrix = []
    columns = [f'{var1}>{str(border1)} and {var2}>{str(border2)}', f'{var1}≤{str(border1)} and {var2}>{str(border2)}',
               f'{var1}>{str(border1)} and {var2}≤{str(border2)}', f'{var1}≤{str(border1)} and {var2}≤{str(border2)}']
    for i, j in df[[var1, var2]].to_numpy():
        if i > border1 and j > border2:
            matrix.append([1, 0, 0, 0])
        elif i <= border1 and j > border2:
            matrix.append([0, 1, 0, 0])
        elif i > border1 and j <= border2:
            matrix.append([0, 0, 1, 0])
        elif i <= border1 and j <= border2:
            matrix.append([0, 0, 0, 1])
        else:
            print(i, j)
    group4_df = pd.DataFrame(matrix, columns=columns, index=df.index)
    return group4_df
