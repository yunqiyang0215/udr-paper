# This script performs the main analysis of the GTEx data.
#
# This is an updated version of the script from here:
# /project2/mstephens/yunqiyang/udr-paper/result202306/data_analysis/analysis2.R
#
# sinteractive -p mstephens --account=pi-mstephens -c 4 --mem=8G \
#   --time=36:00:00
# module load R/4.1.0
# R --no-save
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_10"
#

# Load the packages needed to perform the analyses.
library(tools)
library(ashr)     # 2.2.66
library(udr)      # 0.3.158
library(mashr)    # 0.2.81
library(flashier) # 1.0.56

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
variants <- c("smart_ed","smart_ed_iw","smart_ted","smart_ted_iw",
              "smart_ted_nn","rand_ed","rand_ed_iw","rand_ted",
              "rand_ted_iw","rand_ted_nn")
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

  # "SMART" INITIALIZATION
  # ----------------------
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

  # RANDOM INITIALIZATION
  # ---------------------
  # Randomly initialize the U matrices.
  U_rand <- vector("list",K)
  names(U_rand) <- paste0("U",1:K)
  for (k in 1:K) {
    U_rand[[k]] <- udr:::sim_unconstrained(R)
  }
  fit0 <- ud_init(dat_train,V = dat$null.cor,U_scaled = NULL,
                  U_unconstrained = U_rand,n_rank1 = 0)

  # Run the ED updates.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ed",
                  resid.update = "none",lambda = 0,tol = tol,
                  tol.lik = tol_lik,maxiter = maxiter))
  t1 <- proc.time()
  timings["rand_ed",i] <- (t1 - t0)["elapsed"]
  loglik_test["rand_ed",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                              version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$rand_ed <- fit

  # Run the ED updates using the IW penalty.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ed",
                  resid.update = "none",penalty.type = "iw",
                  tol = tol,tol.lik = tol_lik,lambda = lambda_iw,
                  maxiter = maxiter))
  t1 <- proc.time()
  timings["rand_ed_iw",i] <- (t1 - t0)["elapsed"]
  loglik_test["rand_ed_iw",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                                 version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$rand_ed_iw <- fit

  # Run the TED updates.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  lambda = 0,maxiter = maxiter))
  t1 <- proc.time()
  timings["rand_ted",i] <- (t1 - t0)["elapsed"]
  loglik_test["rand_ted",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                               version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$rand_ted <- fit
  
  # Run the TED updates using the IW penalty.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  lambda = lambda_iw,penalty.type = "iw",
                  maxiter = maxiter))
  t1 <- proc.time()
  timings["rand_ted_iw",i] <- (t1 - t0)["elapsed"]
  loglik_test["rand_ted_iw",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                                  version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$rand_ted_iw <- fit

  # Run the TED updates with the nuclear norm penalty.
  t0 <- proc.time()
  fit <- ud_fit(fit0,verbose = TRUE,
                control = list(unconstrained.update = "ted",
                  resid.update = "none",tol = tol,tol.lik = tol_lik,
                  lambda = lambda_nn,penalty.type = "nu",
                  maxiter = maxiter))
  t1 <- proc.time() 
  timings["rand_ted_nn",i] <- (t1 - t0)["elapsed"]
  loglik_test["rand_ted_nn",i] <- udr:::loglik_ud(dat_test,fit$w,fit$U,fit$V,
                                                  version = "Rcpp")
  fit["X"] <- NULL
  fit["P"] <- NULL
  f$rand_ted_nn <- fit
  
  # Store the results for this "fold".
  fits[[i]] <- f
}

# Save the results to an .RData file.
save(list = c("lambda_iw","lambda_nn","tol","tol_lik","maxiter",
              "fits","loglik_test","timings"),
     file = "gtex_main_out.RData")
resaveRdaFiles("gtex_main_out.RData")

