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
    lasso.fit(X, y[:, 0])  # 只使用生存时间作为因变量
    selected_features = X[:, lasso.coef_ != 0]  # 选择系数不为零的预测因素

    cph = CoxPHFitter()
    cph.fit(selected_features, y)
    cph.print_summary()  # 打印模型摘要信息

