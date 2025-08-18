# This script performs the main analysis of the GTEx data.
#
# This is an updated version of the script from here:
# /project2/mstephens/yunqiyang/udr-paper/result202306/data_analysis/analysis2.R
#
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
maxiter   <- 1000

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
for (i in 1:nfold) {
  
  # Split the data into a "training set" and a "test set".
  indx      <- c(seq(split_points[i] + 1,split_points[i + 1]))
  dat_test  <- dat$strong.z[indx,]
  dat_train <- dat$strong.z[-indx,]
  mash_data <- mashr::mash_set_data(dat_train,V = dat$null.cor)

  # "Smart" initialization.
  U_smart <- smart_initialization(mash_data)
  K <- length(U_smart)

  stop()
  
  U.random <- list()
  for (k in 1:K) {
    U.random[[k]] <- udr:::sim_unconstrained(R)
  }

  # Smart initialization 
  ## Run ED & ED.iw
  f0 <- ud_init(dat.train,V = dat$null.cor,U_scaled = NULL,
                U_unconstrained = U.smart,n_rank1 = 0)
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
                U_unconstrained = U.smart,n_rank1 = 0)
  ted <- ud_fit(f0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",penalty.type = "iw",
                  tol = 1e-02,tol.lik = 1e-2,lambda = 0,maxiter = 5e3))
  ted.iw <- ud_fit(f0,verbose = TRUE,
                   control = list(unconstrained.update = "ted",
                     resid.update = "none",tol = 1e-02,tol.lik = 1e-2,
                     lambda = lambda,penalty.type = "iw",maxiter = 5e3))
  
  res1 <- list(ed,ed.iw,ted,ted.iw)
  
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

saveRDS(output,"result2.rds")
