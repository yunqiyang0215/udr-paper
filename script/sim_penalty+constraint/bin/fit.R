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

# Initialization
set.seed(888)
U.init = c()
K = length(param$w)
V = param$V
maxiter = 2e3

for (k in 1:K){
  U.init[[k]] <- udr:::sim_unconstrained(p)
}
f0 = ud_init(dat_train$X, V = V, U_scaled = NULL, U_unconstrained = U.init, n_rank1 = 0)


if (prior_cov_update == "oracle"){
  fit = NULL
}

# Algorithms without regularization
if (prior_cov_update == "ed"){
  fit = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                     tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = maxiter), verbose=TRUE)
}

if (prior_cov_update == "ted"){
  fit = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                      tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = maxiter), verbose=TRUE)
}

if (prior_cov_update == "fa"){
  fit = ud_fit(f0, control = list(unconstrained.update = "fa", resid.update = 'none',
                                  tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = maxiter), verbose=TRUE)
}


# Algorithms with regularization 
if (prior_cov_update == "ted.iw"){
  fit = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                         tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
}

if (prior_cov_update == "ted.nu"){
  fit = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                  tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "nu", maxiter = maxiter), verbose=TRUE)
}

if (prior_cov_update == "ed.iw"){
  fit = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                  tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = maxiter), verbose=TRUE)
}

if (prior_cov_update == "ted.rank1"){
  # Initialize rank1 Us
  for (k in 1:K){
    U.init[[k]] <- udr:::sim_rank1(p)
  }
  f0 = ud_init(dat_train$X, V = V, U_scaled = NULL, n_unconstrained = 0, U_rank1 = U.init)
  fit = ud_fit(f0, control = list(unconstrained.update = "ted", rank1.update = "ted", resid.update = 'none',
                                  tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = maxiter), verbose=TRUE)
}

if (prior_cov_update == "fa.rank1"){
  # Initialize rank1 Us
  for (k in 1:K){
    U.init[[k]] <- udr:::sim_rank1(p)
  }
  f0 = ud_init(dat_train$X, V = V, U_scaled = NULL, n_unconstrained = 0, U_rank1 = U.init)
  fit = ud_fit(f0, control = list(unconstrained.update = "ted", rank1.update = "fa", resid.update = 'none',
                                  tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = maxiter), verbose=TRUE)
}

