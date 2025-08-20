# This script performs final analysis of the GTEx data
# informed by the results of the "gtex_main" analysis.
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
lambda  <- 83.4
tol     <- 0.0001
tol_lik <- 0.0001
maxiter <- 50 # 1000

# "Smart" initialization
# ----------------------
smart_initialization <- function (mash_data) {
  X.center <- apply(mash_data$Bhat,2,function(x) x - mean(x))
  XtX      <- t(X.center) %*% X.center / nrow(X.center)
  U.f      <- cov_flash(mash_data,factors = "nonneg",tag = "non_neg")
  U.pca    <- cov_pca(mash_data,5)
  U.init   <- c(U.f,U.pca,list(XtX = XtX))
  return(U.init)
}
mash_data <- mashr::mash_set_data(dat$strong.z,V = dat$null.cor)
U_smart <- smart_initialization(mash_data)
K <- length(U_smart)
fit0 <- ud_init(dat$strong.z,V = dat$null.cor,U_scaled = NULL,
                U_unconstrained = U_smart,n_rank1 = 0)

# Run the ED updates with no penalty.
timings <- rep(0,2)
names(timings) <- c("ed","ted_iw")
t0 <- proc.time()
fit_ed <- ud_fit(fit0,verbose = TRUE,
                 control = list(unconstrained.update = "ed",
                   resid.update = "none",lambda = 0,tol = tol,
                   tol.lik = tol_lik,maxiter = maxiter))
t1 <- proc.time()
timings["ed"] <- (t1 - t0)["elapsed"]
fit_ed["X"] <- NULL
fit_ed["P"] <- NULL

# Run the TED updates using the IW penalty.
t0 <- proc.time()
fit_ted_iw <- ud_fit(fit0,verbose = TRUE,
                     control = list(unconstrained.update = "ted",
                       resid.update = "none",tol = tol,tol.lik = tol_lik,
                       lambda = lambda,penalty.type = "iw",
                       maxiter = maxiter))
t1 <- proc.time()
timings["ted_iw"] <- (t1 - t0)["elapsed"]
fit_ted_iw["X"] <- NULL
fit_ted_iw["P"] <- NULL

stop()

# Save the results to an .RData file.
# TO DO.
