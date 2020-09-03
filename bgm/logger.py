#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
import numpy as np
from sklearn.metrics import classification_report

import parse_opts
import load_matrix
import find_best_cluster

def log_parameters(opts):
        print("Info : load trio from %s "%self.trio_file,file=sys.stderr);
        print("Info : load 2mer from %s "%self.mer2_file,file=sys.stderr);
        print("Info : loop num is %d"%self.loop_num,file=sys.stderr);

def log_matrix(X):
    print("Info : the shape of matrix is %s " % str(X.shape) ,file=sys.stderr)

def log_pca_ration(ratio):
    print("Info : the pca ration is %s "%str(ratio) , file=sys.stderr)

def log_bgm_result(i,finder):
    print("################################  Round %d"%i,file=sys.stderr)
    print("Info : covariances is \n%s\n"% finder.covariances,file=sys.stderr)
    print("Info : cluster size is\n%s\n"% str(finder.all_counter),file=sys.stderr)
    print("Info : formula_predict in cluster is\n%s\n"% str(finder.pf_counter),file=sys.stderr)
    print("Info : best-hit is %d with smallest covariances %f "%(finder.best_hit ,finder.best_var),file=sys.stderr)

def log_meterics(data,new_Y):
    if not data.debug :
        print("Info : no Y in non-debug mode.",file=sys.stderr)
    else :
        print("Info : in debug mode now : ",file=sys.stderr)
        print(classification_report(data.Y,new_Y),file=sys.stderr)
