#!/usr/bin/env dsc

simulate: simdata.R
    K: 10
    n_test: 1e4
    U_pattern: "varied", "rank1"
    par: c(p = 5, n_train = 1e4),
         c(p = 50, n_train = 1e3)
    $dat_train: dat_train
    $dat_test: dat_test
    $param: param
    $U_true_pattern: U_pattern
    $p: p
   

fit: fit.R
    # module input and variables
    dat_train: $dat_train
    param: $param
    p: $p
    prior_cov_update: "oracle", "ted", "ed", "fa", "ted.iw", "ted.nu", "ed.iw", "ted.rank1", "fa.rank1" 
    # module output
    $fit: fit
    
posterior_train: posterior.R
    # module input and variables
    dat: $dat_train
    param: $param
    fit: $fit
    # module output
    $lfsr: lfsr
    $theta: theta
    $posterior_mean: posterior_mean

    
loglik_train: loglik.R
    # module input and variables
    dat: $dat_train
    param: $param
    fit: $fit
    # module output
    $loglik: loglik

loglik_test: loglik.R
    # module input and variables
    dat: $dat_test
    param: $param
    fit: $fit
    # module output
    $loglik: loglik



DSC:
    define:
      evaluate: loglik_train, loglik_test, posterior_train
    run: simulate * fit * evaluate
    replicate: 20
    R_libs: udr,mashr
    exec_path: bin
    output: dsc_result

