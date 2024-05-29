# Script to get data example for Fig1
library(udr)
library(mashr)
library(mvtnorm)


setwd("/project2/mstephens/yunqiyang/udr-paper/result202306/Fig1.202309/convergence_dat_example")
dat1 = readRDS("/project2/mstephens/yunqiyang/udr-paper/result202306/Fig1.202309/dsc_result/simulate/simulate_101.rds")
dat2 = readRDS("/project2/mstephens/yunqiyang/udr-paper/result202306/Fig1.202309/dsc_result/simulate/simulate_2.rds")
dat3 = readRDS("/project2/mstephens/yunqiyang/udr-paper/result202306/Fig1.202309/dsc_result/simulate/simulate_103.rds")
dat4 = readRDS("/project2/mstephens/yunqiyang/udr-paper/result202306/Fig1.202309/dsc_result/simulate/simulate_104.rds")


### Data example 1
# Fitting parameters
lambda <- dat1$p

# Warm start. Run ED or ED.iw or FA for 20 iterations, and then run other algorithms.
f1 <- ud_fit(dat1$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                     tol = 0, tol.lik = 0, lambda = 0, maxiter = 20), verbose=TRUE)

f1.reg <- ud_fit(dat1$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 20), verbose=TRUE)

fit.ed = ud_fit(f1, control = list(unconstrained.update = "ed", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)
fit.ted = ud_fit(f1, control = list(unconstrained.update = "ted", resid.update = 'none',
                                    tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.fa = ud_fit(f1, control = list(unconstrained.update = "fa", resid.update = 'none',
                                    tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.ediw = ud_fit(f1.reg, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)
fit.tediw = ud_fit(f1.reg, control = list(unconstrained.update = "ted", resid.update = 'none',
                                          tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)


example1 = list(fit.ed, fit.ted, fit.fa, fit.ediw, fit.tediw)
names(example1) = c("ed", "ted", "fa", "ed.iw", "ted.iw")
saveRDS(example1, "example1.rds")



### Data example 2
# Fitting parameters
lambda <- dat2$p

# Warm start. Run ED or ED.iw or FA for 20 iterations, and then run other algorithms.
f1 <- ud_fit(dat2$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                     tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 20), verbose=TRUE)

f1.reg <- ud_fit(dat2$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 20), verbose=TRUE)

fit.ed = ud_fit(f1, control = list(unconstrained.update = "ed", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)
fit.ted = ud_fit(f1, control = list(unconstrained.update = "ted", resid.update = 'none',
                                    tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.fa = ud_fit(f1, control = list(unconstrained.update = "fa", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.ediw = ud_fit(f1.reg, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)
fit.tediw = ud_fit(f1.reg, control = list(unconstrained.update = "ted", resid.update = 'none',
                                          tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)


example2 = list(fit.ed, fit.ted, fit.fa, fit.ediw, fit.tediw)
names(example2) = c("ed", "ted", "fa", "ed.iw", "ted.iw")
saveRDS(example2, "example2.rds")


### Data example 3
# Fitting parameters
lambda <- dat3$p

# Warm start. Run ED or ED.iw or FA for 20 iterations, and then run other algorithms.
f1 <- ud_fit(dat3$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                     tol = 1e-02, tol.lik = 1e-2, lambda = 0, maxiter = 20), verbose=TRUE)

f1.reg <- ud_fit(dat3$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 1e-02, tol.lik = 1e-2, lambda = lambda, penalty.type = "iw", maxiter = 20), verbose=TRUE)

fit.ed = ud_fit(f1, control = list(unconstrained.update = "ed", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)
fit.ted = ud_fit(f1, control = list(unconstrained.update = "ted", resid.update = 'none',
                                    tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.fa = ud_fit(f1, control = list(unconstrained.update = "fa", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.ediw = ud_fit(f1.reg, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)
fit.tediw = ud_fit(f1.reg, control = list(unconstrained.update = "ted", resid.update = 'none',
                                          tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)


example3 = list(fit.ed, fit.ted, fit.fa, fit.ediw, fit.tediw)
names(example3) = c("ed", "ted", 'fa', "ed.iw", "ted.iw")
saveRDS(example3, "example3.rds")


### Data example 4
# Fitting parameters
lambda <- dat4$p

# Warm start. Run ED or ED.iw or FA for 20 iterations, and then run other algorithms.
f1 <- ud_fit(dat4$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                     tol = 0, tol.lik = 0, lambda = 0, maxiter = 20), verbose=TRUE)

f1.reg <- ud_fit(dat4$f0, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 20), verbose=TRUE)

fit.ed = ud_fit(f1, control = list(unconstrained.update = "ed", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)
fit.ted = ud_fit(f1, control = list(unconstrained.update = "ted", resid.update = 'none',
                                    tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.fa = ud_fit(f1, control = list(unconstrained.update = "fa", resid.update = 'none',
                                   tol = 0, tol.lik = 0, lambda = 0, maxiter = 2e3), verbose=TRUE)

fit.ediw = ud_fit(f1.reg, control = list(unconstrained.update = "ed", resid.update = 'none',
                                         tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)
fit.tediw = ud_fit(f1.reg, control = list(unconstrained.update = "ted", resid.update = 'none',
                                          tol = 0, tol.lik = 0, lambda = lambda, penalty.type = "iw", maxiter = 2e3), verbose=TRUE)


example4 = list(fit.ed, fit.ted, fit.fa, fit.ediw, fit.tediw)
names(example4) = c("ed", "ted", 'fa', "ed.iw", "ted.iw")
saveRDS(example4, "example4.rds")

