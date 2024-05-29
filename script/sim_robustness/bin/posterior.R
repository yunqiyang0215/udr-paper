##### Helper functions #####

# Function to create the fit object "g", which is the same as in mashr
# @param w: the vector for mixture component weights
# @param U: a list of prior covariance matrices
create_g <- function(w, U){
  g = list(w, U, c(1), FALSE)
  names(g) = c("pi", "Ulist", "grid", "usepointmass")
  return(g)
}



##### Save posterior quantities for comparison: lfsr(mat), theta(mat) #####
library(udr)
library(mashr)
library(mvtnorm)

if (is.null(fit)){
  U.est = param$U
  w.est = param$w
}else{
  U.est <- lapply(fit$U,function (e) "[["(e,"mat"))
  w.est <- fit$w
}


mc <- mash_set_data(dat$X, Shat = 1)
g.fitted <- create_g(w.est, U.est)
res <- mash(mc, g = g.fitted, fixg = TRUE)
theta <- dat$theta
lfsr <- res$result$lfsr
posterior_mean <- res$result$PosteriorMean

