#' Function to simulate data from multivariate EBNM models
#' with true mean (\theta) stored.
#' @param n: number of observations
#' @param w: component weight for each mixture
#' @param U: a list of prior covariance matrices
#' @param V: residual covariance matrix
simulate_mixture_ebnm <- function (n, w, U, V) {

  # Get the dimension of the data points (m) and the number of
  # mixture components (k).
  if (is.matrix(U))
    U <- list(U = U)
  m <- nrow(U[[1]])
  k <- length(U)
  if (n < 2)
    stop("n should be 2 or greater")
  if (m < 2)
    stop("The univariate case (m = 1) is not implemented")

  # Check the residual covariance matrix.
  if (missing(V))
    V <- diag(nrow(U[[1]]))
  if (!udr:::issemidef(V))
    stop("Input argument \"V\" should be a positive semi-definite matrix")

  # Check the prior covariance matrices.
  for (i in 1:k)
    if (!udr:::issemidef(U[[i]]))
      stop("All \"U\" matrices should be positive semi-definite")

  # Check the mixture weights.
  if (missing(w))
    w <- rep(1/k,k)
  if (!(is.numeric(w) & length(w) == k & all(w >= 0)))
    stop("Input argument \"w\" should be a vector of length \"k\" ",
         "containing non-negative weights")
  w <- w/sum(w)

  # Draw mixture components according to the mixture weights.
  z <- sample(k,n,replace = TRUE,prob = w)

  # Draw the data points.
  X <- matrix(0,n,m)
  theta <- matrix(0, n, m)
  for (j in 1:k) {
    i <- which(z == j)
    if (length(i) > 0) {
      theta[i, ] = rmvnorm(length(i),sigma = U[[j]])
      for (s in i){
        X[s, ] = rmvnorm(1, theta[s, ], sigma = V)
      }
    }
  }
  rownames(X) <- paste0("s",1:n)
  colnames(X) <- rownames(V)
  rownames(theta) <- paste0("s",1:n)
  colnames(theta) <- rownames(V)
  return(list(X=X, theta = theta, z = z))
}

#' Function to simulate data from EBNM model for one component.
#' It first simulate the true mean and then observed data.
#' @param n Number of data points to simulate.
#' @param U Prior covariance matrix
#' @param V Residual covariance matrix
simulate_ebnm_data <- function(n, U, V){
  m = ncol(U)
  X = matrix(NA, nrow = n, ncol = m)
  theta <- rmvnorm(n, sigma = U)
  for (i in 1:n)
    X[i, ] <- rmvnorm(1, theta[i, ], sigma = V)
  return(list(X = X, theta = theta))
}


#' Function to create the fit object "g", which is the same as in mashr
#' @param w: the vector for mixture component weights
#' @param U: a list of prior covariance matrices
create_g <- function(w, U){
  g = list(w, U, c(1), FALSE)
  names(g) = c("pi", "Ulist", "grid", "usepointmass")
  return(g)
}


#' @param X: observed data matrix of n by R
#' @param theta: data matrix of true means of n by R.
#' @param W: the posterior weight matrix
#' @param U: a list of prior covariance matrices estimated from udr
#' @param V: the residual covariance
compute_posterior_loglik2 <- function(X, theta, W, U, V){
  K = ncol(W)
  m = ncol(X)
  n = nrow(X)
  posterior_mean <- array(0, dim = c(n, m, K))
  posterior_cov <- array(0, dim = c(m, m, K))
  log_post_mat = matrix(NA, nrow = n, ncol = K)

  for (k in 1:K){
    # For each mixture component, compute posterior mean for each observation,
    # and compute posterior covariance.
    posterior_mean[,,k] <- t(apply(X, 1, function(x) udr:::compute_posterior_mvtnorm(x, U[[k]], V)$mu1))
    posterior_cov[,,k] <- solve(U[[k]] %*% solve(V) + diag(m)) %*% U[[k]]
    for (i in 1:n){
      log_post_mat[i,k] <- log(W[i, k]) + dmvnorm(theta[i, ], mean = posterior_mean[i,,k],
                                                  sigma = posterior_cov[,,k], log = TRUE)
    }
  }
  log_post <- sum(apply(log_post_mat, 1, udr:::log_sum_exp))
  return(log_post)
}


#' @param X: observed data matrix of n by R
#' @param theta: data matrix of true means of n by R.
#' @param W: the posterior weight matrix
#' @param U: a list of prior covariance matrices estimated from udr
#' @param V: the residual covariance
compute_posterior_loglik <- function(X, theta, W, U, V){
  K = ncol(W)
  m = ncol(X)
  n = nrow(X)
  posterior_mean <- array(0, dim = c(n, m, K))
  posterior_cov <- array(0, dim = c(m, m, K))
  log_post_mat = matrix(NA, nrow = n, ncol = K)

  for (k in 1:K){
    # For each mixture component, compute posterior mean for each observation,
    # and compute posterior covariance.
    posterior_mean[,,k] <- t(apply(X, 1, function(x) udr:::compute_posterior_mvtnorm(x, U[[k]], V)$mu1))
    posterior_cov[,,k] <- solve(U[[k]] %*% solve(V) + diag(m)) %*% U[[k]]

    log_post_mat[ ,k] <- sapply(1:n, function(i) log(W[i, k]) + dmvnorm(theta[i, ], mean = posterior_mean[i,,k],
                                                                        sigma = posterior_cov[,,k], log = TRUE))
  }
  log_post <- sum(apply(log_post_mat, 1, udr:::log_sum_exp))
  return(log_post)
}


#' Function for computing evaluation metrics
#' @param fit: fit object from udr
#' @param dat: data object contains observed data X and true mean theta
#' @param param: a list contains true parameter values for simulating object dat.
compute_metrics <- function(fit, dat, param){
  # compute data log-likelihood difference between truth
  U <- lapply(fit$U,function (e) "[["(e,"mat"))
  loglik <- udr:::loglik_ud_iid_helper(dat$X, fit$w, simplify2array(U), param$V) -
    udr:::loglik_ud_iid_helper(dat$X, param$w, simplify2array(param$U), param$V)
  # Set mash data
  g.fitted <- create_g(fit$w, U)
  g.true <- create_g(param$w, param$U)
  mc <- mash_set_data(dat$X, Shat = 1)
  mc.fit <- mash(mc, g = g.fitted, fixg = TRUE)
  mc.true <- mash(mc, g = g.true, fixg = TRUE)

  # Compare lfsr
  res.lfsr <- compare_lfsr_fitted_vs_truth(dat, g.fitted, mc, threshold = 0.95)

  # Compute posterior fit
  logpost <- compute_posterior_loglik(dat$X, dat$theta, mc.fit$posterior_weights, U, param$V) -
    compute_posterior_loglik(dat$X, dat$theta, mc.true$posterior_weights, param$U, param$V)

  res = c(loglik, logpost, res.lfsr[1], res.lfsr[2])
  names(res) = c("loglik", "logpost", names(res.lfsr)[1], names(res.lfsr)[2])
  return(res)
}


#' Function to compare lfsr obtained under true parameters vs fitted parameters.
#' @param g.true: an object created for mash input under true parameters
#' @param g.fitted: an object created for mash input under fitted parameters
#' @param mc: a mash object
#' @returns ppv: positive predictive value; tpr: true positive rate. This is different from the traditional definition
compare_lfsr_fitted_vs_true_param <- function(g.true, g.fitted, mc){

  res.true <- mashr::mash_compute_posterior_matrices(g.true, mc)
  res.fitted <- mashr::mash_compute_posterior_matrices(g.fitted, mc)

  # Calculate probability that an effect is positive
  res.true$PositiveProb = 1 - res.true$NegativeProb - res.true$lfdr
  res.fitted$PositiveProb = 1 - res.fitted$NegativeProb - res.fitted$lfdr

  # number of true positives, true discoveries with correct signs
  tp = sum(sum((res.true$PositiveProb >= 0.95) & (res.fitted$PositiveProb >= 0.95)) +
             sum((res.true$NegativeProb >= 0.95) & (res.fitted$NegativeProb >= 0.95)))
  # ns.true: non-signals under truth
  # ns.fitted: non-signals under fitted param
  ns.true = (res.true$NegativeProb < 0.95) & (res.true$PositiveProb < 0.95)
  ns.fitted = (res.fitted$NegativeProb < 0.95) & (res.fitted$PositiveProb < 0.95)

  # fp(false positives): real non-signals classified as signals under fitted param + signals with wrong signs compared to truth
  fp = sum(ns.true) - sum(ns.true & ns.fitted) +
    sum((res.true$PositiveProb >= 0.95) & (res.fitted$NegativeProb >= 0.95)) +
    sum((res.true$NegativeProb >= 0.95) & (res.fitted$PositiveProb >= 0.95))

  tpr =  tp/ sum(sum(res.true$PositiveProb >= 0.95) + sum(res.true$NegativeProb >= 0.95))
  ppv = tp/(tp + fp)
  res <- c(tpr, ppv)
  names(res) = c("tpr", "ppv")
  return(res)
}


#' Function to compare lfsr under fitted parameters vs. true means, theta.
#' @param dat: a data object from simulate_mixture_ebnm() containing X, theta, Z.
#' @param g.fitted: an object created for mash input under fitted parameters
#' @param mc: a mash object
#' @param threshold: a threshold for lfsr. Usually it's 0.95.
#' @returns ppv: positive predictive value; tpr: true positive rate. This is different from the traditional definition
compare_lfsr_fitted_vs_truth <- function(dat, g.fitted, mc, threshold){


  res.fitted <- mashr::mash_compute_posterior_matrices(g.fitted, mc)

  # Calculate probability that an effect is positive
  res.fitted$PositiveProb = 1 - res.fitted$NegativeProb - res.fitted$lfdr

  # number of true positives, true discoveries with correct signs
  tp = sum((dat$theta > 0) & (res.fitted$PositiveProb >= threshold)) +
             sum((dat$theta < 0) & (res.fitted$NegativeProb >= threshold))


  # ns.true: non-signals
  # ns.fitted: non-signals under fitted param
  ns.true = dat$theta == 0
  ns.fitted = (res.fitted$NegativeProb < threshold) & (res.fitted$PositiveProb < threshold)

  # fp(false positives): real non-signals classified as signals under fitted param + signals with wrong signs compared to truth
  fp = sum(ns.true) - sum(ns.true & ns.fitted) +
    sum((dat$theta > 0) & (res.fitted$NegativeProb >= threshold)) +
    sum((dat$theta < 0) & (res.fitted$PositiveProb >= threshold))

  tpr =  tp/(sum(dat$theta > 0) + sum(dat$theta<0))
  ppv = tp/(tp + fp)
  res <- c(tpr, ppv)
  names(res) = c("tpr", "ppv")
  return(res)
}



# Simulate single condition covariance matrix.
# @param R: covariance dimension
# @param s: scale of the matrix
# @param num: number of matrices to simulate
cov_singletons = function(R, s, num){
  Ulist = list()
  if (num == 0) {
    return(Ulist)
  }
  if (num >= R) {num = R}
  for (i in 1:num){
    U = matrix(0, ncol = R, nrow = R)
    U[i, i] = s
    Ulist[[i]] = U
  }
  return(Ulist)
}



# Create partial sharings where sharing only occurs in half of the conditions.
# @param R: number of conditions.
# @param s: the scaling factor
# @param corr: strength of sharing
cov_structured_sharing = function(R, s, corr=c(0.25,0.5,0.75)){
  U <- list()
  r = round(R/2)
  for(i in 1:length(corr)){
    start <- sample(1:(R-r), 1)
    mat <- matrix(0, ncol = R, nrow =R)
    mat[start:(start+r),  start:(start+r)] <- corr[i]
    diag(mat[start:(start+r),  start:(start+r)]) <-1
    U[[i]] <- mat*s
  }
  return(U)
}


# Function to simulate various Us.
# @param R: data dimension
# @param s: scale of the matrix
sim_U_true <- function(R, s, null.mat = FALSE, identity = FALSE, partial_sharing = FALSE,
                       equal_effect = FALSE, num_singleton = 0, num_unconstrained = 0){

  U.unconstrained <- list()
  U.c <- list()
  U.structured <- list()
  U.singletons <- cov_singletons(R, s, num_singleton)

  if (num_unconstrained != 0){
    for (i in 1:num_unconstrained){
      U.unconstrained[[i]] <- sim_invwishart(R, nu = R + 2, s = s)
    }
  }
  if (null.mat == TRUE){
    U.null <- matrix(0, ncol = R, nrow = R)
    U.c = c(U.c, list(U.null))
  }
  if (identity == TRUE){
    U.identity = s*diag(R)
    U.c = c(U.c, list(U.identity))
  }
  if (equal_effect == TRUE){
    U1 <- s*matrix(1, ncol = R, nrow = R)
    U.c = c(U.c, list(U1))
  }
  if (partial_sharing == TRUE){
    U.structured <- cov_structured_sharing(R, s)
  }

  Ulist = c(U.unconstrained, U.singletons, U.c, U.structured)
  return(Ulist)
}





# Function to simulate from inverse Wishart distribution
# @param nu: degree of freedom
# @param R: matrix dimension
# @param s: scale of the matrix
# @return U
sim_invwishart <- function(R, nu = R + 2, s = 5){
  # With the standard parameterization, U <- rinvwishart(nu, S), E(U) = S/(nu-R - 1).
  # If we want E(U) <- s*I, then
  S <- s*(nu- R - 1)*diag(R)
  U <- rinvwishart(nu, S)
  return(U)
}

## Function to compute empirical lfsr. That is, around a threshold t for lfsr,
## the proportion of units with wrong signs.
# @param theta: a matrix of underlying true means
# @param posterior_mean: a matrix of posterior means
# @param lfsr: a matrix contains the lfsr for each unit.
# @param range: a vector of length 2 containing the local range around a threshold t.
compute_empirical_lfsr <- function(theta, posterior_mean, lfsr, range){
  lfsr.empirical <- NA
  indx <- (lfsr >= range[1] & lfsr <= range[2])
  n <- sum(indx)   # total number of units in the range.
  if (n != 0){
    theta.lc <- theta[indx]
    posterior_mean.lc <- posterior_mean[indx]
    # Units with wrong signs.
    # If the estimated sign is wrong, then theta*posterior_mean <= 0.
    # (except when theta == posterior_mean == 0)
    m <- sum( theta.lc * posterior_mean.lc <= 0) - sum(theta.lc == 0 & posterior_mean.lc == 0)
    lfsr.empirical <- m/n
  }
  return(lfsr.empirical)
}


## Function to compute empirical fsr. That is, for all lfsr <= t, the proportion of
## lfsr with wrong signs.
# @param theta: a matrix of underlying true means
# @param posterior_mean: a matrix of posterior means
# @param lfsr: a matrix contains the lfsr for each unit.
# @param t: a threshold for lfsr
compute_empirical_fsr <- function(theta, posterior_mean, lfsr, t){
  fsr.empirical <- NA
  indx <- (lfsr <= t)
  n <- sum(indx)   # total number of units in the range.
  if (n != 0){
    theta.tail <- theta[indx]
    posterior_mean.tail <- posterior_mean[indx]
    # Units with wrong signs.
    # If the estimated sign is wrong, then theta*posterior_mean <= 0.
    # (except when theta == posterior_mean == 0)
    m <- sum( theta.tail * posterior_mean.tail <= 0) - sum(theta.tail == 0 & posterior_mean.tail == 0)
    fsr.empirical <- m/n
  }
  return(fsr.empirical)
}


## Function to compute nominal fsr. That is, mean(lfsr <= t)
# @param lfsr: a matrix contains the lfsr for each unit.
# @param t: a threshold for lfsr
compute_nominal_fsr <- function(lfsr, t){
  fsr.nominal <- NA
  if (sum(lfsr <= t) != 0){
    fsr.nominal <- mean(lfsr[lfsr <= t])
  }
  return(fsr.nominal)
}


## Function to compute true positive rate & false discovery rate
## at a given threshold t.
# @param theta: a matrix of underlying true means
# @param posterior_mean: a matrix of posterior means
# @param lfsr: a matrix contains the lfsr for each unit.
# @param t: a threshold for lfsr
create_tpr_vs_fdr_curve <- function(theta, posterior_mean, lfsr, t){
  res <- c(0, 0)
  names(res) = c("tpr", "fdr")
  indx <- (lfsr <= t)
  n <- sum(indx) # total number of units in the range.
  if (n != 0){
    theta.subset <- theta[indx]
    posterior_mean.subset <- posterior_mean[indx]
    ## Number of true positives
    tp <- sum(theta.subset * posterior_mean.subset > 0)
    tpr <-  tp/(sum(theta > 0) + sum(theta<0))
    # Units with wrong signs within all significant units at t.
    fp <- sum(theta.subset * posterior_mean.subset <= 0) - sum(theta.subset == 0 & posterior_mean.subset == 0)
    fdr <- fp / (fp + tp)
    res[1] <- tpr
    res[2] <- fdr
  }
  return(res)
}



