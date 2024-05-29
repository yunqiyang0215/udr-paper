library(udr)
library(mashr)
library(mvtnorm)

if (is.null(fit)){
  loglik <- udr:::loglik_ud_iid_helper(dat$X, param$w, simplify2array(param$U), param$V)
}else{
  U <- lapply(fit$U,function (e) "[["(e,"mat"))
  loglik <- udr:::loglik_ud_iid_helper(dat$X, fit$w, simplify2array(U), param$V) 
}



