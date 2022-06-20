import numpy as np
import pandas as pd
from lifelines import CoxPHFitter


# %%
def cox(df, time='time', status='status', mod='single') -> pd.DataFrame:
    """单因素和多因素 Cox-regression
    使用多因素分析时候应对高相关性变量进行筛除，否则会导致函数不收敛

    :param df: pd.DataFrame
    :param time: time
    :param status: status
    :param mod: single or multiple
    :return: pd.DataFrame
    """
    raw_variables = df.columns.drop(time).drop(status)
    continuous_index = df.dtypes[(df.dtypes == int) | (df.dtypes == 'int64') | (df.dtypes == float)].index.intersection(
        raw_variables)
    group_colname = df.dtypes[df.dtypes == object].index.intersection(raw_variables)
    single_variables = list(map(lambda x: list(x), zip(continuous_index, continuous_index)))
    multi_variables = list(map(lambda x: list(x), zip(continuous_index, continuous_index)))

    group_df = pd.DataFrame()
    if len(group_colname) > 0:
        for groupby in group_colname:
            groups = df[groupby].unique()
            group_df[groups[0]] = (df[groupby] == groups[0]).replace({True: 1, False: 0})
            single_variables.append([groupby, groups[0]])
            for group in groups[1:]:
                group_df[group] = (df[groupby] == group).replace({True: 1, False: 0})
                single_variables.append([groupby, group])
                multi_variables.append([groupby, group])
        df = pd.concat([df.drop(group_colname, axis=1), group_df], axis=1)
    single_variables = np.array(single_variables).T.tolist()
    multi_variables = np.array(multi_variables).T.tolist()
    print(single_variables, multi_variables)
    cph = CoxPHFitter()
    if mod == 'single':
        variables = single_variables[1]
        univariate_df = pd.DataFrame()
        for var in variables:
            tmp_df = df[[time, status, var]]
            cph.fit(tmp_df, duration_col=time, event_col=status)
            univariate_df = pd.concat([univariate_df, cph.summary], axis=0)
        univariate_df = univariate_df[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
        cox_df = univariate_df

    elif mod == 'multiple':
        variables = multi_variables[1]
        cph.fit(pd.concat([df[[time, status]], df[variables]], axis=1), duration_col=time, event_col=status)
        multivariate_df = cph.summary
        multivariate_df = multivariate_df[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
        unknown_var = [i for i in variables if i not in variables]
        unknown_df = pd.DataFrame([[1, 1, 1, 1]] * len(unknown_var), index=unknown_var,
                                  columns=['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p'])
        cox_df = pd.concat([multivariate_df, unknown_df], axis=0)

    empty_df = pd.DataFrame(index=single_variables[1])
    cox_df = pd.concat([cox_df, empty_df], axis=1).loc[single_variables[1]].fillna(1)
    cox_df = pd.DataFrame(cox_df.to_numpy(), columns=['HR', 'HR(95CI-Low)', 'HR(95CI-High)', 'p-value'],
                          index=single_variables)
    cox_df.time_col = time
    cox_df.status_col = status
    return cox_df
