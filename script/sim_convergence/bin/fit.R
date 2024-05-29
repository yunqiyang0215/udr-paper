##### Helper functions #####

# Function to create the fit object "g", which is the same as in mashr
# @param w: the vector for mixture component weights
# @param U: a list of prior covariance matrices
create_g <- function(w, U){
  g = list(w, U, c(1), FALSE)
  names(g) = c("pi", "Ulist", "grid", "usepointmass")
  return(g)
}


##### Fit: TED/ED/FA/TED.iw/TED.nu/ED.nu/TED.rank1 #####
# Load packages
library(udr)
library(mashr)
library(mvtnorm)
library(LaplacesDemon)


# Fitting parameters
lambda <- p
maxiter = 2e3

# Warm start. Run ED or ED.iw or FA for 20 iterations, and then run other algorithms.
f1 <- ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                             tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 20), verbose=TRUE)

f1.reg <- ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                    tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 20), verbose=TRUE)
  
f.fa <- ud_fit(f0, control = list(unconstrained.update = "fa", resid.update = 'none',
                             tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 20), verbose=TRUE)

  
# Algorithms without regularization
if (prior_cov_update == "ed"){
  fit = ud_fit(f1, control = list(unconstrained.update = "ed", resid.update = 'none',
                                     tol = 0, tol.lik = 0, lambda = 0, maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
	
}

if (prior_cov_update == "ted"){
  fit = ud_fit(f1, control = list(unconstrained.update = "ted", resid.update = 'none',
                                      tol = 0, tol.lik = 0, lambda = 0, maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}

if (prior_cov_update == "fa"){
  fit = ud_fit(f.fa, control = list(unconstrained.update = "fa", resid.update = 'none',
                                     tol = 0, tol.lik = 0, lambda = 0, maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}

if (prior_cov_update == "fa_ted"){
  fit = ud_fit(f.fa, control = list(unconstrained.update = "ted", resid.update = 'none',
                                      tol = 0, tol.lik = 0, lambda = 0, maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}


if (prior_cov_update == "ed_fa"){
  fit = ud_fit(f1, control = list(unconstrained.update = "fa", resid.update = 'none',
                                    tol = 0, tol.lik = 0, lambda = 0, maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}



# Algorithms with regularization 
if (prior_cov_update == "ted.iw"){
  fit = ud_fit(f1.reg, control = list(unconstrained.update = "ted", resid.update = 'none',
                                         tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}

if (prior_cov_update == "ed.iw"){
  fit = ud_fit(f1.reg, control = list(unconstrained.update = "ed", resid.update = 'none',
                                  tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}

# Run TED from ED solution
if (prior_cov_update == "ed_ted"){
  fit = ud_fit(f1, control = list(unconstrained.update = "ed", resid.update = 'none',
                                  tol = 0, tol.lik = 0, lambda = 0, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
  fit = ud_fit(fit, control = list(unconstrained.update = "ted", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}

if (prior_cov_update == "ed.iw_ted.iw"){
  fit = ud_fit(f1.reg, control = list(unconstrained.update = "ed", resid.update = 'none',
                                  tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
  fit = ud_fit(fit, control = list(unconstrained.update = "ted", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
  res <- tail(fit$progress, 1)
  loglik = res$loglik
  loglik_pen = res$loglik.pen
  progress = fit$progress$loglik
  progress_pen = fit$progress$loglik.pen
}

