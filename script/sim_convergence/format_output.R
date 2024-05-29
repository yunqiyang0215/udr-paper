# Description: This file read and post-process simulation results from DSC
# Before running it, you have to activate DSC environment, for example, conda activate dsc
# R 4.2.0 is used when running DSC 


library(dscrutils)
res <- dscquery("dsc_result",
                c("simulate.p", "fit.prior_cov_update", "fit.loglik","fit.loglik_pen", "fit.progress", "fit.progress_pen"))
saveRDS(res, "res202309.rds")

