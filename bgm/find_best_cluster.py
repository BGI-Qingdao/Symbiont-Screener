#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import math

def Counter(Y):
    d = {}
    for y in Y:
        if y in d :
            d[y] += 1;
        else :
            d[y] = 1;
    return d

def ReverseSortCounters(D):
    l = []
    for k , v in D:
        l.append([k,v])
    l.sort(key=takeSecond,reverse = True)
    return l

class BestHitFinder:

    def BestHit(self,predict_Y , formula_predict , covariances):
        self.covariances = covariances
        self.all_counter = Counter(predict_Y)

        PY=predict_Y[formula_predict==1]
        self.pf_counter = Counter(PY)
        sorted_counters = ReverseSortCounters(self.pf_counter)
        best_var = 0.0 ;
        best_hit = -1
        for i in range(3):
            cluster_id , counter = sorted_counters[i]
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

