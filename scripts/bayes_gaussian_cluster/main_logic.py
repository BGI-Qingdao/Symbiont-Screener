#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import numpy as np
import parse_opts
import load_matrix
import consensus_result
import machine_learning
import find_best_cluster
import print_result
import logger
import mkl

if __name__ == '__main__':
    opts = parse_opts.OPTs()
    opts.Parse(sys.argv[1:])
    mkl.set_num_threads(opts.threads)
    np.random.seed(opts.rseed)
    data = load_matrix.MatrixLoader(trio_file=opts.trio_file
            ,mer2_file=opts.mer2_file,debug=opts.debug)
    data.Load()
    logger.log_matrix(data.X)

    X, pca_ratio = machine_learning.doPCA(data.X)
    logger.log_pca_ration(pca_ratio);

    results = consensus_result.ConsensusResult(X.shape[0],opts.loop_num)
    for i in range(opts.loop_num):
        predict_Y , covariances = machine_learning.doBGM(X,opts.clusters_number)
        logger.log_cluster(data,predict_Y)
        best_hit_finder = find_best_cluster.BestHitFinder()
        valid , best_hit, second_best_hit = best_hit_finder.BestHit(predict_Y,data.formula_predict,covariances)
        if not valid :
            continue
        logger.log_bgm_result(i,best_hit_finder);
        results.PushOneResult(predict_Y,best_hit,second_best_hit)

    results.GenFinalResult()
    logger.log_meterics(data,results.hosts)
    logger.log_hit(opts.loop_num,results.scores)
    print_result.save_gmm_results(results.results,results.best_hits,results.second_best_hits)
    print_result.print_result(data.formula_predict,results.hosts,results.scores,results.second_scores);

print("ALL DONE ... ", file=sys.stderr)


