

library(dscrutils)
#loglik1 <- dscquery("dsc_result1",
#                    c("simulate.p", "simulate.U_pattern", "fit.prior_cov_update", "fit.K_init", "evaluate","evaluate.loglik"))

#loglik2 <- dscquery("dsc_result2",
#                    c("simulate.p", "simulate.U_pattern", "fit.prior_cov_update", "fit.K_init", "evaluate","evaluate.loglik"))

#loglik3 <- dscquery("dsc_result3",
#                    c("simulate.p", "simulate.U_pattern", "fit.prior_cov_update", "fit.K_init", "evaluate","evaluate.loglik"))

#loglik <- rbind(loglik1, loglik2, loglik3)
#saveRDS(loglik, "loglik.rds")



lfsr1 <- dscquery("dsc_result1",
                     c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","fit.K_init", "evaluate", "evaluate.lfsr"))
lfsr2 <- dscquery("dsc_result2",
                  c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","fit.K_init", "evaluate", "evaluate.lfsr"))
lfsr3 <- dscquery("dsc_result3",
                  c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","fit.K_init", "evaluate", "evaluate.lfsr"))


theta1<- dscquery("dsc_result1",
                      c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update", "fit.K_init", "evaluate", "evaluate.theta"))
theta2<- dscquery("dsc_result2",
                  c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update", "fit.K_init", "evaluate", "evaluate.theta"))
theta3<- dscquery("dsc_result3",
                  c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update", "fit.K_init", "evaluate", "evaluate.theta"))


post_mean1 <- dscquery("dsc_result1",
                          c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update","fit.K_init", "evaluate", "evaluate.posterior_mean"))
post_mean2 <- dscquery("dsc_result2",
                       c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update", "fit.K_init", "evaluate", "evaluate.posterior_mean"))
post_mean3 <- dscquery("dsc_result3",
                       c("simulate.p", "simulate.U_true_pattern",  "fit.prior_cov_update", "fit.K_init", "evaluate", "evaluate.posterior_mean"))



# Scenario 1: varied + p = 5
# Scenario 2: rank1 + p = 5
# Scenario 3: varied + p = 50
# Scenario 4: rank1 + p = 50

U_pattern = c("varied", "rank1")
p = c(5, 50)
K.init = c(2, 4, 8)

for (i in 1:length(p)){
  for (j in 1:length(U_pattern)){
    for (k in 1:length(K.init)){
      dsc <- lfsr1$DSC[lfsr1$evaluate == "posterior_train" & lfsr1$simulate.p == p[i] & lfsr1$simulate.U_true_pattern == U_pattern[j] & lfsr1$fit.K_init == K.init[k]]
      lfsr <- lfsr1$evaluate.lfsr[lfsr1$evaluate == "posterior_train" & lfsr1$simulate.p == p[i] & lfsr1$simulate.U_true_pattern == U_pattern[j] & lfsr1$fit.K_init == K.init[k]]
      theta <- theta1$evaluate.theta[theta1$evaluate == "posterior_train" & theta1$simulate.p == p[i] & theta1$simulate.U_true_pattern == U_pattern[j] & theta1$fit.K_init == K.init[k]]
      post_mean <- post_mean1$evaluate.posterior_mean[post_mean1$evaluate == "posterior_train" & post_mean1$simulate.p == p[i] & post_mean1$simulate.U_true_pattern == U_pattern[j] & post_mean1$fit.K_init == K.init[k]]
      post <- list(dsc, lfsr, theta, post_mean)
      names(post) = c("dsc", "lfsr", "theta", "posterior_mean")
      file_name = paste0(U_pattern[j],"_R", p[i], "_K", K.init[k])
      saveRDS(post, paste0("./post/",file_name, ".rds"))
    }
  }
}


K.init = c(16, 32)

for (i in 1:length(p)){
  for (j in 1:length(U_pattern)){
    for (k in 1:length(K.init)){
      dsc <- lfsr2$DSC[lfsr2$evaluate == "posterior_train" & lfsr2$simulate.p == p[i] & lfsr2$simulate.U_true_pattern == U_pattern[j] & lfsr2$fit.K_init == K.init[k]]
      lfsr <- lfsr2$evaluate.lfsr[lfsr2$evaluate == "posterior_train" & lfsr2$simulate.p == p[i] & lfsr2$simulate.U_true_pattern == U_pattern[j] & lfsr2$fit.K_init == K.init[k]]
      theta <- theta2$evaluate.theta[theta2$evaluate == "posterior_train" & theta2$simulate.p == p[i] & theta2$simulate.U_true_pattern == U_pattern[j] & theta2$fit.K_init == K.init[k]]
      post_mean <- post_mean2$evaluate.posterior_mean[post_mean2$evaluate == "posterior_train" & post_mean2$simulate.p == p[i] & post_mean2$simulate.U_true_pattern == U_pattern[j] & post_mean2$fit.K_init == K.init[k]]
      post <- list(dsc, lfsr, theta, post_mean)
      names(post) = c("dsc", "lfsr", "theta", "posterior_mean")
      file_name = paste0(U_pattern[j],"_R", p[i], "_K", K.init[k])
      saveRDS(post, paste0("./post/",file_name, ".rds"))
    }
  }
}

U_pattern = c("varied", "rank1")
p = c(5, 50)
K.init = c(64, 128)

for (i in 1:length(p)){
  for (j in 1:length(U_pattern)){
    for (k in 1:length(K.init)){
      dsc <- lfsr3$DSC[lfsr3$evaluate == "posterior_train" & lfsr3$simulate.p == p[i] & lfsr3$simulate.U_true_pattern == U_pattern[j] & lfsr3$fit.K_init == K.init[k]]
      lfsr <- lfsr3$evaluate.lfsr[lfsr3$evaluate == "posterior_train" & lfsr3$simulate.p == p[i] & lfsr3$simulate.U_true_pattern == U_pattern[j] & lfsr3$fit.K_init == K.init[k]]
      theta <- theta3$evaluate.theta[theta3$evaluate == "posterior_train" & theta3$simulate.p == p[i] & theta3$simulate.U_true_pattern == U_pattern[j] & theta3$fit.K_init == K.init[k]]
      post_mean <- post_mean3$evaluate.posterior_mean[post_mean3$evaluate == "posterior_train" & post_mean3$simulate.p == p[i] & post_mean3$simulate.U_true_pattern == U_pattern[j] & post_mean3$fit.K_init == K.init[k]]
      post <- list(dsc, lfsr, theta, post_mean)
      names(post) = c("dsc", "lfsr", "theta", "posterior_mean")
      file_name = paste0(U_pattern[j],"_R", p[i], "_K", K.init[k])
      saveRDS(post, paste0("./post/",file_name, ".rds"))
    }
  }
}


