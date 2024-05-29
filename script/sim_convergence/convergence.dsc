#!/usr/bin/env dsc
# This file produces results for optimization behaviours for TED, ED and FA.  


simulate: sim_init.R
    K: 10
    par: c(p = 5, n_train = 1e4),
         c(p = 50, n_train = 1e3)
    $dat: dat
    $f0: f0
    $p: p
   
fit: fit.R
    # module input and variables
    dat: $dat
    p: $p
    f0: $f0
    prior_cov_update: "ted", "ed", "fa", "ed_fa", "fa_ted", "ted.iw", "ed.iw", "ed_ted", "ed.iw_ted.iw"
    # module output
    $loglik: loglik
    $loglik_pen: loglik_pen
    $progress: progress
    $progress_pen: progress_pen   

DSC:
    run: simulate * fit
    replicate: 100
    R_libs: udr,mashr
    exec_path: bin
    output: dsc_result
