import random

import pandas as pd
from sklearn import linear_model, svm, neighbors, naive_bayes, tree, ensemble, neural_network
from sklearn.metrics import f1_score, recall_score, precision_score, accuracy_score, r2_score
from sklearn.model_selection import RepeatedKFold, train_test_split, KFold


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
