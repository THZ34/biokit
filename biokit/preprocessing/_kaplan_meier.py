# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt


# %%
def km_best_cutoff(df, value, time='time', status='status', min_samples=None):
    """找到p值最低的cutoff

    :param min_samples:
    :param df:
    :param value:
    :param time:
    :param status:
    :return:
    """
    df = df.copy()[[value, time, status]]
    pvalues = []
    n_samples = []
    for cutoff in df[value]:
        df['group'] = df[value] > cutoff
        duration_A = df[df['group'] == True][time]
        duration_B = df[df['group'] == False][time]
        event_A = df[df['group'] == True][status]
        event_B = df[df['group'] == False][status]
        pvalue = logrank_test(durations_A=duration_A, durations_B=duration_B, event_observed_A=event_A,
                              event_observed_B=event_B).p_value
        pvalues.append(pvalue)
        n_sample = df['group'].value_counts().to_dict().get(True, 0)
        n_samples.append(n_sample)

    df['pvalue'] = pvalues
    df['high_samples'] = n_samples
    df['low_samples'] = df.shape[0] - df['high_samples']
    df.sort_values(by=[value], inplace=True)
    print(df)
    if min_samples is not None:
        cutoff = df[(df['high_samples'] >= min_samples) & (df['low_samples'] >= min_samples)].sort_values(by='pvalue')[value].to_numpy()[0]
    else:
        cutoff = df.sort_values(by='pvalue')[value].to_numpy()[0]
    return cutoff, df


def km_base_cutoff_plot(df, value):
    fig, ax = plt.subplots()
    ax.plot(range(df.shape[0]), df.sort_values(by=value)['pvalue'])
    ax.plot([0, df.shape[0]], [0.05, 0.05])
    return ax
