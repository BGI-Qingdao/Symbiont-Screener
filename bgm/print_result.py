#!/usr/bin/python3
# -*- coding: UTF-8 -*-

def print_result(predict_Y,consensus_Y,Y_score):
    print("formula_predict\tconsensus_result\thit_count")
    for i in range(len(predict_Y)):
        print("%d\t%d\t%d"%(predict_Y[i],consensus_Y[i],Y_score[i]))
