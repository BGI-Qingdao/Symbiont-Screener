#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import math
import numpy as np

def Counter(Y):
    d = {}
    for y in Y:
        if y in d :
            d[y] += 1;
        else :
            d[y] = 1;
    return d

def takeSecond(elem):
    return elem[1]

def ReverseSortCounters(D):
    l = []
    for k  in D.keys():
        l.append([k,D[k]])
    l.sort(key=takeSecond,reverse = True)
    return l

class BestHitFinder:

    def BestHit(self,predict_Y , formula_predict , covariances):
        self.covariances = covariances
        self.all_counter = Counter(predict_Y)
        total_predict=np.sum(formula_predict==1)
        #print(self.all_counter)
        PY=predict_Y[formula_predict==1]
        self.pf_counter = Counter(PY)
        #print(self.pf_counter)
        sorted_counters = ReverseSortCounters(self.pf_counter)
        best_var = 0.0 ;
        best_hit = -1
        for i in range(3):
            if i >= len(sorted_counters):
                break
            cluster_id , counter = sorted_counters[i]
            if counter <= total_predict / 10 :
                break 
            if counter <= 0 :
                break
            cov = covariances[cluster_id][0][0]
            if best_var == 0.0 or math.fabs(cov) < best_var :
                best_var = math.fabs(cov)
                best_hit = cluster_id
        self.best_hit = best_hit
        self.best_var = best_var 
        self.valid = bool( best_hit != -1 )
        if best_hit == -1 :
            return False , -1
        else :
            return True , best_hit

