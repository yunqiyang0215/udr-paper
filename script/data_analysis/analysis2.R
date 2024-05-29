## Description: Compare old workflow with new workflow
setwd('/project2/mstephens/yunqiyang/udr-data-application')

#### Load packages

library(Matrix)
library(mvtnorm)
library(ashr)
library(udr)
library(mashr)
library(flashr)
library(LaplacesDemon)


smart_initialization <- function(mash_data){
  X.center = apply(mash_data$Bhat, 2, function(x) x - mean(x))
  XtX <- t(X.center) %*% X.center / nrow(X.center)
  # From mashr vignette
  U.f = cov_flash(mash_data, factors="nonneg", tag="non_neg", var_type="constant")
  U.pca = cov_pca(mash_data, 5)
  U.init = c(U.f, U.pca, list(XtX))
  return(U.init)
}


# set directory
dat = readRDS('/project2/mstephens/yunqiyang/udr-data-application/mvSuSiE_output/GTEx_V8_strong_z.rds')

# Fixed parameters
tol = 0.01
tol.lik = 0.01
R = nrow(dat$strong.z)


# split data
set.seed(999)
n = nrow(dat$strong.z)
output = list()

nfold <- 5
split_points <- round(seq(from = 0, to = n, length.out = nfold + 1))


for (i in 1:nfold) {
  # split data
  indx = c((split_points[i]+1):split_points[i+1])
  dat.test = dat$strong.z[indx, ]
  dat.train = dat$strong.z[-indx, ]
  
  # set mashr data object
  mash_data = mashr::mash_set_data(dat.train, V = dat$null.cor)
  
  # initialization
  V = dat$null.cor
  t1 = proc.time()
  U.smart <- smart_initialization(mash_data)
  t2 = proc.time()
  K = length(U.smart)
  R <- nrow(U.smart[[1]])
  lambda <- R
  
  U.random <- list()
  for (k in 1:K){
    U.random[[k]] <- udr:::sim_unconstrained(R)
  }
  
  
  # Smart initialization 
  ## Run ED & ED.iw
  f0 = ud_init(dat.train, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.smart, n_rank1 = 0)
  ed = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                 tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
  ed.iw = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none', penalty.type = "iw",
                                    tol = 1e-02, tol.lik = 1e-2, lambda = lambda, maxiter = 5e3), verbose=TRUE)
  
  
  ## Run TED & TED.iw
  f0 = ud_init(dat.train, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.smart, n_rank1 = 0)
  ted = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none', penalty.type = "iw",
                                  tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
  ted.iw = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                     tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 5e3), verbose=TRUE)
  
  res1 = list(ed, ed.iw, ted, ted.iw)
  
  # Random initialization
  ## Run ED & ED.iw
  f0 = ud_init(dat.train, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.random, n_rank1 = 0)
  ed = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                 tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
  ed.iw = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none', penalty.type = "iw",
                                    tol = 1e-02, tol.lik = 1e-2, lambda = lambda, maxiter = 5e3), verbose=TRUE)
  
  ## Run TED & TED.iw
  f0 = ud_init(dat.train, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.random, n_rank1 = 0)
  ted = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none', penalty.type = "iw",
                                  tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
  ted.iw = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                     tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 5e3), verbose=TRUE)
  
  time.init = t2 - t1
  res2 = list(ed, ed.iw, ted, ted.iw)
  result = list(res1, res2, dat.test, time.init)
  names(result) = c("smart", "random", "test", "time")
  output[[i]] = result
}

saveRDS(output, "/project2/mstephens/yunqiyang/udr-paper/result202306/data_analysis/result2.rds")





