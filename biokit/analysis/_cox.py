import copy
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceError
from statsmodels.stats.outliers_influence import variance_inflation_factor


# %%
def cox(df, time='time', status='status', variables=None, mod='single', drop_by_vif=True,
        ref_dict=None) -> pd.DataFrame:
    """单因素和多因素 Cox-regression
    使用多因素分析时候应对高相关性变量进行筛除，否则会导致函数不收敛

    :param variables:
    :param ref_dict: 离散变量的参考值
    :param drop_by_vif: VIF==np.inf时自动去除以避免cox不收敛
    :param df: pd.DataFrame
    :param time: time
    :param status: status
    :param mod: single or multiple
    :return: pd.DataFrame
    """
    df = df.copy()
    # 排除缺省值
    if df.isna().sum().sum() > 0:
        raise ValueError('数据中包含缺省值')

    sur_df = df[[time, status]]

    if not variables:
        raw_variables = df.columns.drop(time).drop(status)
    else:
        raw_variables = variables
    continuous_index = df.dtypes[(df.dtypes == int) | (df.dtypes == 'int64') | (df.dtypes == float)].index.intersection(
        raw_variables)
    discrete_index = df.dtypes[(df.dtypes == object) | (df.dtypes == bool)].index.intersection(raw_variables)
    continuous_df = df[continuous_index]
    continuous_df.columns = [continuous_index, continuous_index]

    # 离散变量哑编码
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
    corr_df = cox_input.corr()

    # 找出重复列
    factor_df = cox_input.drop([time, status], axis=1).copy()
    dup_cols = factor_df.T[factor_df.T.duplicated()].index
    dup_col_pair = dict([(i, corr_df[corr_df[i] == 1].drop(i).index[0]) for i in dup_cols])
    dup_refs = dup_cols.intersection(multi_ref_variables)

    # 过滤0方差特征
    # variance_zero_vars = df.columns[df.std() == 0]
    # df.drop(variance_zero_vars, axis=1, inplace=True)
    # print(f'变量 {variance_zero_vars} 方差为0')

    # cox分析
    cph = CoxPHFitter()
    if mod == 'single':
        variables = cox_input.columns.drop(multi_ref_variables).drop([time, status])
        single_cox_result = []
        for var in variables:
            try:
                tmp_df = cox_input[[time, status, var]]
                cph.fit(tmp_df, duration_col=time, event_col=status)
                single_cox_result.append(cph.summary)
            except ConvergenceError:
                single_cox_result.append(pd.DataFrame(index=[var]))

        single_cox_result = pd.concat(single_cox_result, axis=0)
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
                for col in copy.copy(dup_cols):
                    if dup_col_pair[col][0] in vif_drop_discrete:
                        dup_cols = dup_cols.drop(col)
                for col in copy.copy(dup_refs):
                    if col[0] in vif_drop_discrete:
                        dup_refs.remove(col)
                # 如果是离散变量，删除同一变量的其他子类
                for col in copy.copy(multi_ref_variables):
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

    n_samples = []
    index_df = cox_result.index.to_frame()
    for groupby in index_df[0].unique():
        groups = index_df.loc[groupby].index.to_list()
        if len(groups) > 1:
            for group in groups:
                n_samples.append(df[df[groupby] == group].shape[0])
        else:
            n_samples.append(df.shape[0])
    cox_result['n_sample'] = n_samples

    return cox_result
