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
#' @param threshold: a threshold for lfsr. Usually it's 0.95, but can be more stringent.
#' @returns ppv: positive predictive value; tpr: true positive rate. This is different from the traditional definition
compare_lfsr_fitted_vs_truth <- function(dat, g.fitted, mc, threshold){


  res.fitted <- mashr::mash_compute_posterior_matrices(g.fitted, mc)

  # Calculate probability that an effect is positive
  res.fitted$PositiveProb = 1 - res.fitted$NegativeProb - res.fitted$lfdr

  # number of true positives, true discoveries with correct signs
  tp = sum(sum((dat$theta > 0) & (res.fitted$PositiveProb >= threshold)) +
             sum((dat$theta < 0) & (res.fitted$NegativeProb >= threshold)))


  # ns.true: non-signals
  # ns.fitted: non-signals under fitted param
  ns.true = dat$theta == 0
  ns.fitted = (res.fitted$NegativeProb < threshold) & (res.fitted$PositiveProb < threshold)

  # fp(false positives): real non-signals classified as signals under fitted param + signals with wrong signs compared to truth
  fp = sum(ns.true) - sum(ns.true & ns.fitted) +
    sum((dat$theta > 0) & (res.fitted$NegativeProb >= threshold)) +
    sum((dat$theta < 0) & (res.fitted$PositiveProb >= threshold))

  tpr =  tp/ sum(sum(dat$theta > 0) + sum(dat$theta<0))
  ppv = tp/(tp + fp)
  res <- c(tpr, ppv)
  names(res) = c("tpr", "ppv")
  return(res)
}



# Simulate single condition covariance matrix.
# @param R: covariance dimension
# @param num: number of matrices to simulate
cov_singletons = function(R, num){
  Ulist = list()
  if (num == 0) {
    return(Ulist)
  }
  if (num >= R) {num = R}
  for (i in 1:num){
    U = matrix(0, ncol = R, nrow = R)
    U[i, i] = 1
    Ulist[[i]] = U
  }
  return(Ulist)
}



# Function to simulate various Us.
# @param R: data dimension
sim_U_true <- function(R, null.mat = TRUE, identity = TRUE, num_singleton, num_unconstrained){
  U_unconstrained = list()
  U_singletons <- cov_singletons(R, num_singleton)

  for (i in 1:num_unconstrained){
    U <- udr:::sim_unconstrained(R)
    U <- U/max(U)
    U_unconstrained[[i]] <- U
  }

  U.c = list()
  if (null.mat == TRUE){
    U.null <- matrix(0, ncol = R, nrow = R)
    U.c = c(U.c, list(U.null))
  }
  if (identity == TRUE){
    U.identity = diag(R)
    U.c = c(U.c, list(U.identity))
  }

  Ulist = c(U_unconstrained, U_singletons, U.c)
  return(Ulist)
}






