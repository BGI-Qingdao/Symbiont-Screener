#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np

def save_gmm_results(results, best, second_best):
    data=np.vstack(results)
    np.savetxt("gmm.result.txt",data.T,"%d")
    data=np.vstack((best,second_best))
    np.savetxt("gmm.besthit.txt",data.T,"%d")


def print_result(predict_Y,consensus_Y,Y_score,Ys_score):
    print("formula_predict\tconsensus_result\tbest_count\tsecond_best_count")
    for i in range(len(predict_Y)):
        print("%d\t%d\t%d\t%d"%(predict_Y[i],consensus_Y[i],Y_score[i],Ys_score[i]))

