import random

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import ensemble
from sklearn import linear_model
from sklearn import naive_bayes
from sklearn import neighbors
from sklearn import neural_network
from sklearn import svm
from sklearn import tree
from sklearn.feature_selection import RFE
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
from sklearn.metrics import r2_score
from sklearn.metrics import recall_score
from sklearn.model_selection import KFold
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import train_test_split


def performance_evaluation(X, y, model):
    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    accuracy = accuracy_score(y, y_pred)
    precision = precision_score(y, y_pred)
    recall = recall_score(y, y_pred)
    f1 = f1_score(y, y_pred)
    return [r2, accuracy, precision, recall, f1]


def fit_models(X, y, models: dict, n_splits=10, n_repeats=5, random_state=0) -> tuple:
    """模型训练及测试

    :param X: 特征 m行 n列
    :param y: 目标值 m个元素
    :param models: 字典，键：模型名称 值：实例化的模型
    :param n_splits: n折
    :param n_repeats: 随即重复次数
    :param random_state: 随机数种
    :return: 训练集分数，测试集分数
    """
    if n_splits > 1:
        rkf = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
        indices = rkf.split(X=X, y=y)
        iterator = [(X.iloc[train_indices, :], X.iloc[test_indices, :], y[train_indices], y[test_indices]) for
                    train_indices, test_indices in indices]
    elif 0 < n_splits < 1:
        random.seed(random_state)
        random_state_list = random.sample(range(10000), n_repeats)
        train_size = round(X.shape[0] * n_splits)
        test_size = X.shape[0] - train_size
        iterator = [train_test_split(X, y, train_size=train_size, test_size=test_size, random_state=random_state) for
                    random_state in random_state_list]
    else:
        raise ValueError('unexpected n_splits value')

    train_scores = {}
    test_scores = {}
    for model_name in models:
        train_score = []
        test_score = []
        for train_X, test_X, train_y, test_y in iterator:
            try:
                models[model_name].fit(X=train_X, y=train_y)
                train_score.append(performance_evaluation(train_X, train_y, models[model_name]))
                test_score.append(performance_evaluation(test_X, test_y, models[model_name]))
            except:
                break
        if len(train_score) != 0 and len(test_score) != 0:
            train_scores[model_name] = pd.DataFrame(train_score,
                                                    columns=['R2', 'accuracy', 'precision', 'recall', 'f1'])
            test_scores[model_name] = pd.DataFrame(test_score, columns=['R2', 'accuracy', 'precision', 'recall', 'f1'])
    return train_scores, test_scores


def sklearn_models():
    models_1 = {'ARDRegression': linear_model.ARDRegression,
                'BayesianRidge': linear_model.BayesianRidge,
                'ElasticNet': linear_model.ElasticNet,
                'ElasticNetCV': linear_model.ElasticNetCV,
                'HuberRegressor': linear_model.HuberRegressor,
                'Lars': linear_model.Lars,
                'LarsCV': linear_model.LarsCV,
                'Lasso': linear_model.Lasso,
                'LassoCV': linear_model.LassoCV,
                'LassoLars': linear_model.LassoLars,
                'LassoLarsCV': linear_model.LassoLarsCV,
                'LassoLarsIC': linear_model.LassoLarsIC,
                'LinearRegression': linear_model.LinearRegression,
                'OrthogonalMatchingPursuit': linear_model.OrthogonalMatchingPursuit,
                'OrthogonalMatchingPursuitCV': linear_model.OrthogonalMatchingPursuitCV,
                'PassiveAggressiveRegressor': linear_model.PassiveAggressiveRegressor,
                'PoissonRegressor': linear_model.PoissonRegressor,
                'Ridge': linear_model.Ridge,
                'RidgeCV': linear_model.RidgeCV,
                'SGDRegressor': linear_model.SGDRegressor,
                'SGDClassifier': linear_model.SGDClassifier,
                'TheilSenRegressor': linear_model.TheilSenRegressor,
                'TweedieRegressor': linear_model.TweedieRegressor}

    models_2 = {'LinearSVC': svm.LinearSVC,
                'LinearSVR': svm.LinearSVR,
                'NuSVC': svm.NuSVC,
                'NuSVR': svm.NuSVR,
                'SVC': svm.SVC,
                'SVR': svm.SVR}

    models_3 = {'KNeighborsClassifier': neighbors.KNeighborsClassifier,
                'KNeighborsRegressor': neighbors.KNeighborsRegressor,
                'KNeighborsTransformer': neighbors.KNeighborsTransformer, 'KernelDensity': neighbors.KernelDensity,
                'LocalOutlierFactor': neighbors.LocalOutlierFactor, 'NearestCentroid': neighbors.NearestCentroid,
                'NearestNeighbors': neighbors.NearestNeighbors,
                'NeighborhoodComponentsAnalysis': neighbors.NeighborhoodComponentsAnalysis,
                'RadiusNeighborsClassifier': neighbors.RadiusNeighborsClassifier,
                'RadiusNeighborsRegressor': neighbors.RadiusNeighborsRegressor,
                'RadiusNeighborsTransformer': neighbors.RadiusNeighborsTransformer}

    models_4 = {'BernoulliNB': naive_bayes.BernoulliNB,
                'CategoricalNB': naive_bayes.CategoricalNB,
                'ComplementNB': naive_bayes.ComplementNB,
                'GaussianNB': naive_bayes.GaussianNB,
                'MultinomialNB': naive_bayes.MultinomialNB}

    models_5 = {'DecisionTreeClassifier': tree.DecisionTreeClassifier,
                'DecisionTreeRegressor': tree.DecisionTreeRegressor,
                'ExtraTreeClassifier': tree.ExtraTreeClassifier,
                'ExtraTreeRegressor': tree.ExtraTreeRegressor}

    models_6 = {'AdaBoostClassifier': ensemble.AdaBoostClassifier,
                'AdaBoostRegressor': ensemble.AdaBoostRegressor,
                'BaggingClassifier': ensemble.BaggingClassifier,
                'BaggingRegressor': ensemble.BaggingRegressor,
                'ExtraTreesClassifier': ensemble.ExtraTreesClassifier,
                'ExtraTreesRegressor': ensemble.ExtraTreesRegressor,
                'GradientBoostingClassifier': ensemble.GradientBoostingClassifier,
                'GradientBoostingRegressor': ensemble.GradientBoostingRegressor,
                'IsolationForest': ensemble.IsolationForest,
                'RandomForestClassifier': ensemble.RandomForestClassifier,
                'RandomForestRegressor': ensemble.RandomForestRegressor,
                'RandomTreesEmbedding': ensemble.RandomTreesEmbedding}

    models_7 = {'MLPClassifier': neural_network.MLPClassifier(hidden_layer_sizes=(220, 110, 50), batch_size=40)}

    models = dict()
    for i in [models_1, models_2, models_3, models_4, models_5, models_6, models_7]:
        models.update(i)
    return models


def leave_one_out(model, X, y):
    """留一法"""
    n_sample = X.shape[0]
    kf = KFold(n_splits=n_sample)
    predict_y = []
    for train_idx, test_idx in kf.split(X, y):
        train_X = X.iloc[train_idx, :]
        train_y = y[train_idx]
        test_X = X.iloc[test_idx, :]
        model.fit(train_X, train_y)
        predict_y.append(model.predict(test_X)[0])
    return predict_y


def feature_selection(fold_result, model_names=None, pdf=None):
    if not model_names:
        model_names = sorted(list(fold_result.keys()))
    if not pdf:
        pdf = PdfPages('feature_selection.pdf')
    for model_name in model_names:
        train_scores = []
        test_scores = []
        scores_list = []
        x = []
        for n_feature in sorted(list(fold_result[model_name].keys())):
            train_score = np.mean(fold_result[model_name][n_feature]['train'])
            test_score = np.mean(fold_result[model_name][n_feature]['test'])
            train_scores.append(train_score)
            test_scores.append(test_score)
            scores_list.append(train_score + test_score)
            x.append(n_feature)

        fig, ax = plt.subplots(figsize=(12, 4))
        ax.plot(x, train_scores, label='Fitting score', zorder=0)
        ax.plot(x, test_scores, label='Generalization score', zorder=0)
        ax.plot(x, scores_list, label='The sum of scores', zorder=0)
        ax.scatter(x=scores_list.index(max(scores_list)) + 1, y=max(scores_list), s=20, color='red', zorder=2)
        ax.text(x=scores_list.index(max(scores_list)) + 1, y=max(scores_list), s='number of feature with the max score',
                ha='left', zorder=2)

        ax.set_ylim(0, 2)
        ax.set_xlabel('Max depth' if 'Tree' in model_name or 'Forest' else 'number of features')
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        ax.set_title(f'{model_name} feature selection')

        plt.subplots_adjust(left=0.05, right=0.8, bottom=0.15, top=0.9)
        pdf.savefig()
    pdf.close()
    plt.close()


def feature_ranking(X, y, model_names, models, random_state=0):
    ranking_dict = {}
    for model_name in model_names:
        try:
            model = models[model_name](random_state=random_state)
        except TypeError:
            model = models[model_name]()
        ranking_dict[model_name] = rfe_features(X, y, model)

    rank_df = pd.DataFrame()
    for model_name in ranking_dict:
        rank_df[model_name] = ranking_dict[model_name]['rank']

    rank_df.to_csv('rank_df.csv')
    rank_df = rank_df.astype(int)
    # sns.clustermap(rank_df, cmap='viridis_r', method='ward')
    # plt.savefig('rank_cluster.png', dpi=1200)
    # plt.close()
    return rank_df


def rfe_features(X, y, model):
    """为了节约时间和计算资源，我们在应用递归特征消除的时候分3次进行
    每次消除100个特征至10000特征
    每次消除20个特征至2000特征
    每次消除1个特征至1个特征"""
    # 第一轮
    n_features_to_select = 10000
    step = 100
    rfe = RFE(estimator=model, n_features_to_select=n_features_to_select, step=step)
    rfe.fit(X, y)
    ranking_1 = pd.DataFrame([X.columns, rfe.ranking_]).T.sort_values(by=1)
    ranking_1[1] = ranking_1[1][ranking_1[1] > 1] * step + n_features_to_select
    X = X.T[rfe.ranking_ == 1].T

    # 第二轮
    n_features_to_select = 2000
    step = 20
    rfe = RFE(estimator=model, n_features_to_select=n_features_to_select, step=step)
    rfe.fit(X, y)
    ranking_2 = pd.DataFrame([X.columns, rfe.ranking_]).T.sort_values(by=1)
    ranking_2[1] = ranking_2[1][ranking_2[1] > 1] * step + n_features_to_select
    X = X.T[rfe.ranking_ == 1].T

    # 第三轮
    n_features_to_select = 1
    step = 1
    rfe = RFE(estimator=model, n_features_to_select=n_features_to_select, step=step)
    rfe.fit(X, y)
    ranking_3 = pd.DataFrame([X.columns, rfe.ranking_]).T.sort_values(by=1)

    ranking_1.index = ranking_1[0]
    ranking_2.index = ranking_2[0]
    ranking_3.index = ranking_3[0]
    ranking = pd.concat([ranking_1, ranking_2, ranking_3]).dropna().sort_values(by=1)
    ranking.columns = ['gene', 'rank']
    return ranking
