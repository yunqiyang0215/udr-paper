# This script performs the main analysis of the GTEx data.
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
library(udr)   # 0.3.158
library(mashr) # 0.2.81
#
# library(Matrix)        # 1.3.3
# library(mvtnorm)       # 1.3.3
# library(ashr)          # 2.2.66
# library(flashier)      # 1.0.56
# library(LaplacesDemon) # 16.1.6
#

# Load the z-scores stored as a 15,636 x 49 matrix (genes x tissues).
dat <- readRDS("../../data/data_analysis/GTEx_V8_strong_z.rds")
n <- nrow(dat$strong.z)
R <- ncol(dat$strong.z)

# Tuning parameters.
lambda_iw <- 83.4
lambda_nn <- 34.0
nfold     <- 5
tol       <- 0.01
tol_lik   <- 0.01
maxiter   <- 50 # 1000

smart_initialization <- function (mash_data) {
  X.center <- apply(mash_data$Bhat,2,function(x) x - mean(x))
  XtX      <- t(X.center) %*% X.center / nrow(X.center)
  U.f      <- cov_flash(mash_data,factors = "nonneg",tag = "non_neg")
  U.pca    <- cov_pca(mash_data,5)
  U.init   <- c(U.f,U.pca,list(XtX = XtX))
  return(U.init)
}

# Repeat for each "fold".
split_points <- round(seq(from = 0,to = n,length.out = nfold + 1))
variants     <- c("smart_ed","smart_ed_iw","smart_ted","smart_ted_iw",
                   "smart_ted_nn")
nvar         <- length(variants)
fits         <- vector("list",nfold)
loglik_test  <- matrix(0,nvar,nfold)
timings      <- matrix(0,nvar,nfold)
names(fits) <- paste0("k",1:nfold)
rownames(loglik_test) <- variants
colnames(loglik_test) <- paste0("k",1:nfold)
rownames(timings) <- variants
colnames(timings) <- paste0("k",1:nfold)
for (i in 1:nfold) {
  f <- vector("list",nvar)
  names(f) <- variants
  
  # Split the data into a "training set" and a "test set".
  indx      <- c(seq(split_points[i] + 1,split_points[i + 1]))
  dat_train <- dat$strong.z[-indx,]
  dat_test  <- dat$strong.z[indx,]
  mash_data <- mashr::mash_set_data(dat_train,V = dat$null.cor)

  # Use the "smart" initialization.
  U_smart <- smart_initialization(mash_data)
  K <- length(U_smart)
  fit0 <- ud_init(dat_train,V = dat$null.cor,U_scaled = NULL,
                  U_unconstrained = U_smart,n_rank1 = 0)

  # Run the ED updates.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ed",
                  resid.update = "none",lambda = 0,tol = tol,
                  tol.lik = tol_lik,maxiter = maxiter))
  t1 <- proc.time()
  timings["smart_ed",i] <- (t1 - t0)["elapsed"]
  loglik_test["smart_ed",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                               version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$smart_ed <- fit

  # Run the ED updates using the IW penalty.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ed",
                  resid.update = "none",penalty.type = "iw",
                  tol = tol,tol.lik = tol_lik,lambda = lambda_iw,
                  maxiter = maxiter))
  t1 <- proc.time()
  timings["smart_ed_iw",i] <- (t1 - t0)["elapsed"]
  loglik_test["smart_ed_iw",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                                  version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$smart_ed_iw <- fit

  # Run the TED updates.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  lambda = 0,maxiter = maxiter))
  t1 <- proc.time()
  timings["smart_ted",i] <- (t1 - t0)["elapsed"]
  loglik_test["smart_ted",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                                version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$smart_ted <- fit
  
  # Run the TED updates using the IW penalty.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  lambda = lambda_iw,penalty.type = "iw",
                  maxiter = maxiter))
  t1 <- proc.time()
  timings["smart_ted_iw",i] <- (t1 - t0)["elapsed"]
  loglik_test["smart_ted_iw",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                                   version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$smart_ted_iw <- fit

  # Run the TED updates with the nuclear norm penalty.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  lambda = lambda_nn,penalty.type = "nu",
                  maxiter = maxiter))
  t1 <- proc.time() 
  timings["smart_ted_nn",i] <- (t1 - t0)["elapsed"]
  loglik_test["smart_ted_nn",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                                   version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$smart_ted_nn <- fit
  
  # Randomly initialize the U matrices.
  U_rand <- vector("list",K)
  names(U_rand) <- paste0("U",1:K)
  for (k in 1:K) {
    U_rand[[k]] <- udr:::sim_unconstrained(R)
  }

  # Store the results for this "fold".
  fits[[i]] <- f

  stop()
  
  # Random initialization
  ## Run ED & ED.iw
  f0 <- ud_init(dat.train,V = dat$null.cor,U_scaled = NULL,
                U_unconstrained = U.random,n_rank1 = 0)
  ed <- ud_fit(f0,verbose = TRUE,
               control = list(unconstrained.update = "ed",
                 resid.update = "none",tol = 1e-02,tol.lik = 1e-2,
                 lambda = 0,maxiter = 5e3))
  ed.iw <- ud_fit(f0,verbose = TRUE,
                  control = list(unconstrained.update = "ed",
                    resid.update = "none",penalty.type = "iw",
                    tol = 1e-02,tol.lik = 1e-2,lambda = lambda,
                    maxiter = 5e3))
  
  ## Run TED & TED.iw
  f0 <- ud_init(dat.train,V = dat$null.cor,U_scaled = NULL,
                U_unconstrained = U.random,n_rank1 = 0)
  ted <- ud_fit(f0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",penalty.type = "iw",
                  tol = 1e-02,tol.lik = 1e-2,lambda = 0,maxiter = 5e3))
  ted.iw <- ud_fit(f0,verbose = TRUE,
                   control = list(unconstrained.update = "ted",
                     resid.update = "none",tol = 1e-02,tol.lik = 1e-2,
                     lambda = lambda,penalty.type = "iw",maxiter = 5e3))
  
  time.init <- t2 - t1
  res2      <- list(ed,ed.iw,ted,ted.iw)
  result    <- list(res1,res2,dat.test,time.init)
  names(result) <- c("smart","random","test","time")
  output[[i]]   <- result
}

# Save the results to an .RData file.
# TO DO.
