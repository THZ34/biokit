import warnings
from copy import copy

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.outliers_influence import variance_inflation_factor


# %%
def cox(df, time='time', status='status', mod='single', drop_by_vif=True,
        ref_dict=None) -> pd.DataFrame:
    """单因素和多因素 Cox-regression
    使用多因素分析时候应对高相关性变量进行筛除，否则会导致函数不收敛

    :param ref_dict: 离散变量的参考值
    :param drop_by_vif: VIF==np.inf时自动去除以避免cox不收敛
    :param df: pd.DataFrame
    :param time: time
    :param status: status
    :param mod: single or multiple
    :return: pd.DataFrame
    """
    if df.isna().sum().sum() > 0:
        raise ValueError('数据中包含缺省值')

    sur_df = df[[time, status]]
    raw_variables = df.columns.drop(time).drop(status)
    continuous_index = df.dtypes[(df.dtypes == int) | (df.dtypes == 'int64') | (df.dtypes == float)].index.intersection(
        raw_variables)
    discrete_index = df.dtypes[df.dtypes == object].index.intersection(raw_variables)
    continuous_df = df[continuous_index]
    continuous_df.columns = [continuous_index, continuous_index]

    # 拆分离散变量
    discrete_df = pd.DataFrame()
    if not ref_dict:
        ref_dict = {}
    multi_ref_variables = []
    if len(discrete_index) > 0:
        for groupby in discrete_index:
            group_df = pd.DataFrame()
            groups = sorted(list(df[groupby].unique()))
            groupdf_columns = [[], []]
            for group in groups:
                group_df[groupby, group] = (df[groupby] == group).replace({True: 1, False: 0})
                groupdf_columns[0].append(groupby)
                groupdf_columns[1].append(group)
            group_df.columns = groupdf_columns
            discrete_df = pd.concat([discrete_df, group_df], axis=1)

            multi_ref_variables.append((groupby, ref_dict.get(groupby, groups[0])))

    cox_input = pd.concat([sur_df, continuous_df, discrete_df], axis=1)
    # 找出重复列
    dup_cols = cox_input.T[cox_input.T.duplicated()].index
    dup_col_pair = dict([(i, cox_input.corr()[cox_input.corr()[i] == 1].drop(i).index[0]) for i in dup_cols])
    dup_refs = dup_cols.intersection(multi_ref_variables)

    # cox分析
    cph = CoxPHFitter()
    if mod == 'single':
        single_cox_result = pd.DataFrame()
        variables = cox_input.columns.drop(multi_ref_variables).drop([time, status])
        for var in variables:
            tmp_df = cox_input[[time, status, var]]
            cph.fit(tmp_df, duration_col=time, event_col=status)
            single_cox_result = pd.concat([single_cox_result, cph.summary], axis=0)
        cox_result = single_cox_result[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
        cox_result.index = variables

    elif mod == 'multiple':
        # 去掉重复列和离散变量的参考值
        for dup_col in dup_cols:
            if dup_col in multi_ref_variables:
                multi_ref_variables.remove(dup_col)

        cox_input.drop(multi_ref_variables, axis=1, inplace=True)
        cox_input.drop(dup_cols, axis=1, inplace=True)

        # 检查变量VIF(方差膨胀系数)
        if cox_input.shape[1] > 3:
            vif_result = {}
            vif_input = cox_input.iloc[:, 2:].to_numpy()
            vif_drop = []
            print_vif = False
            for i, variable in zip(range(cox_input.shape[0] - 2), cox_input.columns[2:]):
                vif_result[variable] = variance_inflation_factor(vif_input, i)
                if np.isinf(vif_result[variable]):
                    print_vif = True
                    vif_drop.append(variable)
            if print_vif:
                warnings.warn(f'以下变量VIF过高，{vif_drop}')

            if drop_by_vif:
                # 删除VIF过高的变量
                cox_input.drop(vif_drop, axis=1, inplace=True)

                vif_drop_discrete = sorted(list(set([i[0] for i in vif_drop]) & set(discrete_index)))
                # 删除值完全相同的重复变量
                for col in copy(dup_cols):
                    if dup_col_pair[col][0] in vif_drop_discrete:
                        dup_cols = dup_cols.drop(col)
                for col in copy(dup_refs):
                    if col[0] in vif_drop_discrete:
                        dup_refs.remove(col)
                # 如果是离散变量，删除同一变量的其他子类
                for col in copy(multi_ref_variables):
                    if col[0] in vif_drop_discrete:
                        multi_ref_variables.remove(col)

        # 拟合cox模型
        cph.fit(cox_input, duration_col=time, event_col=status)
        cox_result = cph.summary
        cox_result = cox_result[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]

    else:
        raise ValueError('mod parameter support "single" and "multiple"')

    # 填充参考值和重复值
    ref_df = pd.DataFrame([[1, 1, 1, 1]] * len(multi_ref_variables), index=multi_ref_variables,
                          columns=['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p'])
    cox_result = pd.concat([cox_result, ref_df], axis=0)
    if mod == 'multiple':
        dup_df = cox_result.loc[[dup_col_pair[dup_col] for dup_col in dup_cols]]
        dup_df.index = dup_cols
        multi_ref_variables.extend(dup_refs)
    else:
        dup_df = pd.DataFrame()
    cox_result = pd.concat([cox_result, dup_df], axis=0)

    cox_result.sort_index(inplace=True)
    cox_result.index = pd.MultiIndex.from_tuples(cox_result.index)
    cox_result.columns = ['HR', 'HR(95CI-Low)', 'HR(95CI-High)', 'p-value']
    cox_result.time_col = time
    cox_result.status_col = status
    cox_result.df = df
    cox_result.multi_ref_variables = multi_ref_variables
    return cox_result
