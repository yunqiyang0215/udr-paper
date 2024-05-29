##### Load helper functions #####

# Function to simulate data from multivariate EBNM models
# with true mean (\theta) stored.
# @param n: number of observations
# @param w: component weight for each mixture
# @param U: a list of prior covariance matrices
# @param V: residual covariance matrix
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


##### Simulate data #####
# Load packages
library(udr)
library(mvtnorm)
library(LaplacesDemon)

# Simulation parameters 
p <- par["p"] # data dimension
n_train <- par["n_train"]
s = 5 # scalar for true Us
w = rep(1/K, K)
param = list(w = rep(1/K, K), U = list(), V = diag(p))

if (U_pattern == "varied"){
  # True Us: 1 singleton + all.one + identity + 7 rinvWishart
  param$U <- sim_U_true(p, s, identity = TRUE, partial_sharing = FALSE,
                        equal_effect = TRUE, num_singleton = 1, num_unconstrained = 7)
  dat_train = simulate_mixture_ebnm(n_train, param$w, param$U, param$V)
  dat_test = simulate_mixture_ebnm(n_test, param$w, param$U, param$V)
}

if (U_pattern == "rank1"){
  # True Us: 5 singletons and 5 randomly simulated rank1s
  U.singleton <- cov_singletons(p, s, 5)
  U.rank1 <- list()
  for (i in 1:5){
    U.rank1[[i]] <- udr:::sim_rank1(p)
  }
  param$U <- c(U.singleton, U.rank1)
  dat_train <- simulate_mixture_ebnm(n_train, param$w,param$U,param$V)
  dat_test <- simulate_mixture_ebnm(n_test, param$w, param$U, param$V)
}


