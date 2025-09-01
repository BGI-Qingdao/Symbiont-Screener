#!/usr/bin/python
# -*- coding: UTF-8 -*-
from sklearn.decomposition import PCA
from sklearn import mixture

def doPCA(X):
    pca=PCA(n_components=X.shape[1],whiten=True)
    X2=pca.fit_transform(X)
    return X2 , pca.explained_variance_ratio_

def doBGM(X2 , n_components=5):
    dpgmm = mixture.BayesianGaussianMixture(n_components=n_components,covariance_type='full').fit(X2)
    predict_Y=dpgmm.predict(X2)
    return predict_Y , dpgmm.covariances_


