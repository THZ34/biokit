# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
from lifelines import CoxPHFitter
from sklearn.linear_model import LassoCV
from biokit.analysis import cox


def lassocox(df, time, status, variants=None, alpha=None, cv=5, random_state=0):
    if not variants:
        variants = [i for i in df.columns if i not in [time, status]]
    X = df[variants]
    y = df[time]
    lasso = LassoCV(cv=cv, alphas=alpha, random_state=0)
    lasso.fit(X, y[:, 0])  # ֻʹ������ʱ����Ϊ�����
    selected_features = X[:, lasso.coef_ != 0]  # ѡ��ϵ����Ϊ���Ԥ������

    cph = CoxPHFitter()
    cph.fit(selected_features, y)
    cph.print_summary()  # ��ӡģ��ժҪ��Ϣ

