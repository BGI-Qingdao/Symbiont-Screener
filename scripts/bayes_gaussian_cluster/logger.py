#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
import numpy as np
from sklearn.metrics import classification_report

import parse_opts
import load_matrix
import find_best_cluster

def log_parameters(opts):
        print("Info : load trio from %s "%self.trio_file,file=sys.stderr ,flush = True);
        print("Info : load 2mer from %s "%self.mer2_file,file=sys.stderr ,flush = True);
        print("Info : loop num is %d"%self.loop_num,file=sys.stderr,flush = True);

def log_matrix(X):
    print("Info : the shape of matrix is %s " % str(X.shape) ,file=sys.stderr,flush = True)

def log_pca_ration(ratio):
    print("Info : the pca ration is %s "%str(ratio) , file=sys.stderr,flush = True)

def log_bgm_result(i,finder):
    print("################################  Round %d"%i,file=sys.stderr,flush = True)
    print("Info : covariances is \n%s\n"% finder.covariances,file=sys.stderr,flush = True)
    print("Info : cluster size is\n%s\n"% str(finder.all_counter),file=sys.stderr,flush = True)
    print("Info : formula_predict in cluster is\n%s\n"% str(finder.pf_counter),file=sys.stderr,flush = True)
    print("Info : best-hit is %d with smallest covariances %f "%(finder.best_hit ,finder.best_var),file=sys.stderr,flush = True)
    print("Info : second-best-hit is %d "%(finder.second_best_hit),file=sys.stderr,flush = True)

def log_cluster(data,predict_Y):
    if not data.debug :
        print("Info : no Y in non-debug mode.",file=sys.stderr,flush = True)
    else :
        print("Info : in debug mode now : ",file=sys.stderr,flush = True)
        for j in range(1,10):
            print("Y==%d:"%j ,file=sys.stderr,flush = True)
            t1=predict_Y[data.Y==j]
            ll={}
            for i in range(5):
                tt=np.sum(t1==i)
                ll[i]=tt
            print(ll,file=sys.stderr,flush = True)

def log_meterics(data,new_Y):
    if not data.debug :
        print("Info : no Y in non-debug mode.",file=sys.stderr,flush = True)
    else :
        print("Info : in debug mode now : ",file=sys.stderr,flush = True)
        print(classification_report(data.Y,new_Y),file=sys.stderr,flush = True)


def log_hit(i,hits):
    ll={}
    for j in range(i):
        tt=np.sum(hits==j)
        ll[j]=tt
    print("Info : hit frequency %s"% str(ll),file=sys.stderr,flush = True)

