library(dscrutils)
loglik <- dscquery("dsc_result",
                   c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","evaluate","evaluate.loglik"))

dat.lfsr <- dscquery("dsc_result",
                     c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","evaluate", "evaluate.lfsr"))

dat.theta <- dscquery("dsc_result",
                      c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","evaluate", "evaluate.theta"))

dat.post_mean <- dscquery("dsc_result",
                          c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","evaluate", "evaluate.posterior_mean"))



# Scenario 1: varied + p = 5
# Scenario 2: rank1 + p = 5
# Scenario 3: varied + p = 50
# Scenario 4: rank1 + p = 50

U_pattern = c("varied", "rank1")
p = c(5, 50)
file_indx = 0 

for (i in 1:length(p)){
  for (j in 1:length(U_pattern)){
    file_indx = file_indx + 1
    dsc <- dat.lfsr$DSC[dat.lfsr$evaluate == "posterior_train" & dat.lfsr$simulate.p == p[i] & dat.lfsr$simulate.U_true_pattern == U_pattern[j]]
    lfsr <- dat.lfsr$evaluate.lfsr[dat.lfsr$evaluate == "posterior_train" & dat.lfsr$simulate.p == p[i] & dat.lfsr$simulate.U_true_pattern == U_pattern[j]]
    theta <- dat.theta$evaluate.theta[dat.theta$evaluate == "posterior_train" & dat.theta$simulate.p == p[i] & dat.theta$simulate.U_true_pattern == U_pattern[j]]
    post_mean <- dat.post_mean$evaluate.posterior_mean[dat.post_mean$evaluate == "posterior_train" & dat.post_mean$simulate.p == p[i] & dat.post_mean$simulate.U_true_pattern == U_pattern[j]]
    post <- list(dsc, lfsr, theta, post_mean)
    names(post) = c("dsc", "lfsr", "theta", "posterior_mean")
    saveRDS(post, paste0("post",file_indx, ".rds"))
  }
}

saveRDS(loglik, "loglik.rds")

