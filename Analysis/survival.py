import os
import sys
from sklearn import set_config
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.preprocessing import OneHotEncoder
from sksurv.nonparametric import kaplan_meier_estimator
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

def getKaplanSurvivalProbsOverTime(y, event='event', time='time'):
    return kaplan_meier_estimator(
        y[event], y[time]
    )

def plotKaplanSurvivalProbsOverTime(y, event='event', time='time'):
    fig, ax = plt.subplots(figsize=(8, 5))
    time, survival_prob = getKaplanSurvivalProbsOverTime(y, event, time)
    ax.step(time, survival_prob, where="post")
    ax.set_ylim(0, 1)
    ax.set_ylabel("est. probability of survival $\hat{S}(t)$")
    ax.set_xlabel("time $t$")
    ax.set_title(event)

def getOneHot(X):
    for colname in X.columns:
        try:
            X[colname] = X[colname].astype(float)
        except:
            X[colname] = X[colname].astype('category')
    return OneHotEncoder().fit_transform(X)

def getTrainedEstimator(X_train, y_train):
    set_config(display="text")
    #estimator = CoxPHSurvivalAnalysis()
    estimator = CoxnetSurvivalAnalysis()
    estimator.fit(X_train, y_train)
    return estimator

def fit_and_score_features(X_train, y_train, X_test, y_test, feature_list):
    n_features = len(feature_list)
    scores = np.empty(n_features)
    #m = CoxPHSurvivalAnalysis()
    m = CoxnetSurvivalAnalysis()
    for j in range(n_features):
        feature = feature_list[j]
        Xj_train = X_train.loc[:, X_train.columns.str.startswith(feature)].values
        Xj_test = X_test.loc[:, X_test.columns.str.startswith(feature)].values
        #try:
        m.fit(Xj_train, y_train)
        #except:
        #    return None
        scores[j] = m.score(Xj_test, y_test)
    return pd.Series(scores, index=feature_list)
    #.sort_values(ascending=False)

def DFtoStrucArray(y, event='event', time='time'):
    return np.array(y[[event, time]].apply(tuple, axis=1), dtype=[(event, 'bool'), (time, 'float')])

#scaler = StandardScaler()
#survival_x_scaled = pd.DataFrame(data=scaler.fit_transform(survival_x), index=survival_x.index, columns=survival_x.columns)