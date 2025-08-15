# This script performs a simple cross-validation on the GTEx data
# to choose the penalty strength for the IW and NN penalties.
#
# This is an updated version of the script from here:
# /project2/mstephens/yunqiyang/udr-paper/result202306/data_analysis/analysis2.R
#
# sinteractive -p mstephens --account=pi-mstephens -c 4 --mem=8G \
#   --time=48:00:00
# module load R/4.1.0
# R --no-save
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_10"
#

# Load the packages needed to perform the analyses.
library(tools)
library(udr) # 0.3.158

# Load the z-scores stored as a 15,636 x 49 matrix (genes x tissues).
dat <- readRDS("../../data/data_analysis/GTEx_V8_strong_z.rds")
n <- nrow(dat$strong.z)
R <- ncol(dat$strong.z)

# Tuning parameters.
LAMBDA  <- 10^seq(-1,2.7,length.out = 20)
K       <- 30
tol     <- 0.01
tol_lik <- 0.01
maxiter <- 200

# Split the data into a "training set" (80%) and a "test set" (20%).
set.seed(1)
i <- sort(sample(n,12509))
j <- sort(setdiff(1:n,i))
dat_train <- dat$strong.z[i,]
dat_test  <- dat$strong.z[j,]

# Randomly initialize the U matrices.
U_rand <- vector("list",K)
names(U_rand) <- paste0("U",1:K)
for (k in 1:K) {
  U_rand[[k]] <- udr:::sim_unconstrained(R)
}

# Initialize the model.
fit0 <- ud_init(dat_train,V = dat$null.cor,U_scaled = NULL,
                U_unconstrained = U_rand,n_rank1 = 0)

# Repeat for each setting of the penalty strength.
# This step is expected to take roughly 5 h.
n <- length(LAMBDA)
cv_iw_fits <- vector("list",n)
cv_nn_fits <- vector("list",n)
cv_iw_loglik_test <- rep(0,n)
cv_nn_loglik_test <- rep(0,n)
for (i in 1:n) {
  lambda <- LAMBDA[i]

  # Run the TED updates with the IW penalty.
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  maxiter = maxiter,lambda = lambda,penalty.type = "iw"))
  fit["X"] <- NULL
  fit["P"] <- NULL
  cv_iw_fits[[i]] <- fit

  # Compute the test-set log-likelihood.
  cv_iw_loglik_test[i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                          version = "Rcpp")
  
  # Run the TED updates with the nuclear norm penalty.
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  maxiter = maxiter,lambda = lambda,penalty.type = "nu"))
  fit["X"] <- NULL
  fit["P"] <- NULL
  cv_nn_fits[[i]] <- fit

  # Compute the test-set log-likelihood.
  cv_nn_loglik_test[i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                          version = "Rcpp")
}

# Save the results to an .RData file.
save(list = c("K","LAMBDA","tol","tol_lik","maxiter",
              "cv_iw_fits","cv_nn_fits","cv_iw_loglik_test",
              "cv_nn_loglik_test"),
     file = "gtex_cv_out.RData")
resaveRdaFiles("gtex_cv_out.RData")
