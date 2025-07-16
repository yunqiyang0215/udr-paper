library(mvtnorm)

# The first objective is to minimize
#
#   ||xx^T - S||^2
#
# for each data point x, where S = U + I.
f <- function (params) {
  y <- 0
  R <- matrix(c(params[1],0,params[2:3]),2,2)
  U <- crossprod(R)
  S <- U + I
  for (i in 1:n) {
    x  <- X[i,]
    xx <- tcrossprod(x)
    y  <- y + sum((xx - S)^2)
  }
  return(y)
}

# The second objective is to minimize
#
#   ||xx^T - S||^2
#
# for each data point x, where S = w1*U1 + w2*U2 + I
g <- function (params) {
  y <- 0
  R1 <- matrix(c(params[1],0,params[2:3]),2,2)
  R2 <- matrix(c(params[4],0,params[5:6]),2,2)
  U1 <- crossprod(R1)
  U2 <- crossprod(R2)
  S  <- w[1]*U1 + w[2]*U2 + I
  for (i in 1:n) {
    x  <- X[i,]
    xx <- tcrossprod(x)
    y  <- y + sum((xx - S)^2)
  }
  return(y)
}

# Simulate a data set.
set.seed(1)
n <- 1000
I <- diag(2)
S <- matrix(c(4,-1,-1,2),2,2)
X <- rmvnorm(n,rep(0,2),S)

# Fit U using the first objective.
params0 <- rep(1,3)
res <- optim(params0,f,method = "Nelder-Mead",
             control = list(maxit = 1000,trace = TRUE))
R <- matrix(c(res$par[1],0,res$par[2:3]),2,2)
U <- crossprod(R)
S1 <- U + I
print(S1)

# Fit U using the second objective.
k <- 2
w <- c(0.8,0.2)
params0 <- rep(1,6)
res <- optim(params0,g,method = "Nelder-Mead",
             control = list(maxit = 1000,trace = TRUE))
R1 <- matrix(c(res$par[1],0,res$par[2:3]),2,2)
R2 <- matrix(c(res$par[4],0,res$par[5:6]),2,2)
U1 <- crossprod(R1)
U2 <- crossprod(R2)
S2 <- w[1]*U1 + w[2]*U2 + I
print(S2)
