
# Description: GTEx data application
# 1. Compare train log-likelihood between ED + specialized initialization vs. TED + random initialization
# 2. Compare regularized vs. unregularized by init with K=100 and see how many components are selected out..

smart_initialization <- function(mash_data){
  X.center = apply(mash_data$Bhat, 2, function(x) x - mean(x))
  XtX <- t(X.center) %*% X.center / nrow(X.center)
  # From mashr vignette
  U.f = cov_flash(mash_data, factors="nonneg", tag="non_neg", var_type="constant")
  U.pca = cov_pca(mash_data, 5)
  U.init = c(U.f, U.pca, list(XtX))
  return(U.init)
}

# load package
library(Matrix)
library(mvtnorm)
library(ashr)
library(udr)
library(mashr)
library(flashr)
library(LaplacesDemon)


# set directory
dat = readRDS('/project2/mstephens/yunqiyang/udr-data-application/mvSuSiE_output/GTEx_V8_strong_z.rds')


# set mashr data object
bhat = dat$strong.z
sbhat = bhat
sbhat[!is.na(sbhat)] = 1
mash_data = mashr::mash_set_data(bhat,sbhat)

U.smart <- smart_initialization(mash_data)
K = length(U.smart)
R <- nrow(U.smart[[1]])
# penalty strength
lambda <- R
print(K)

set.seed(215)
U.random <- list()
for (k in 1:K){
  U.random[[k]] <- udr:::sim_unconstrained(R)
}


# Smart initialization 
## Run ED & ED.iw
f0 = ud_init(dat$strong.z, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.smart, n_rank1 = 0)
ed = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                               tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
ed.iw = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none', penalty.type = "iw",
                                  tol = 1e-02, tol.lik = 1e-2, lambda = lambda, maxiter = 5e3), verbose=TRUE)


## Run TED & TED.iw
f0 = ud_init(dat$strong.z, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.smart, n_rank1 = 0)
ted = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none', penalty.type = "iw",
                                tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
ted.iw = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                   tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 5e3), verbose=TRUE)

res1 = list(ed, ed.iw, ted, ted.iw)

# Random initialization
## Run ED & ED.iw
f0 = ud_init(dat$strong.z, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.random, n_rank1 = 0)
ed = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                               tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
ed.iw = ud_fit(f0, control = list(unconstrained.update = "ed", resid.update = 'none', penalty.type = "iw",
                                  tol = 1e-02, tol.lik = 1e-2, lambda = lambda, maxiter = 5e3), verbose=TRUE)


## Run TED & TED.iw
f0 = ud_init(dat$strong.z, V = dat$null.cor, U_scaled = NULL, U_unconstrained = U.random, n_rank1 = 0)
ted = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none', penalty.type = "iw",
                                tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 5e3), verbose=TRUE)
ted.iw = ud_fit(f0, control = list(unconstrained.update = "ted", resid.update = 'none',
                                   tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 5e3), verbose=TRUE)

res2 = list(ed, ed.iw, ted, ted.iw)
result = list(res1, res2)
names(result) = c("smart", "random")
saveRDS(result, '/project2/mstephens/yunqiyang/udr-paper/result202306/data_analysis/analysis1.rds')
