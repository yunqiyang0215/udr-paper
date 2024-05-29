library(udr)

# Summarize cross-validation result
compute_loglik_per_test_dat <- function(model.fit, test.data, V){
  n = nrow(test.data)
  U <- lapply(fit$U,function (e) "[["(e,"mat"))
  loglik <- udr:::loglik_ud_iid_helper(test.data, fit$w, simplify2array(U), V)
  return(loglik/n)
}

res = readRDS("/project2/mstephens/yunqiyang/udr-paper/result202306/data_analysis/result2.rds")

test.loglik.smart <- matrix(NA, ncol = 4, nrow = length(res))
colnames(test.loglik.smart) = c("ed", "ed.iw", "ted", "ted.iw")
test.loglik.random <- matrix(NA, ncol = 4, nrow = length(res))
colnames(test.loglik.random) = c("ed", "ed.iw", "ted", "ted.iw")

iteration.random <- matrix(NA, ncol = 4, nrow = length(res))
iteration.smart <- matrix(NA, ncol = 4, nrow = length(res))

time.random <- matrix(NA, ncol = 4, nrow = length(res))
time.smart <- matrix(NA, ncol = 4, nrow = length(res))



for (i in 1:length(res)){
  res.fold <- res[[i]]
  smart = res.fold[[1]]
  random = res.fold[[2]]
  test.data = res.fold[[3]]
  for (j in 1:length(smart)){
    fit = smart[[j]]
    iteration.smart[i, j] = nrow(fit$progress)
    time.smart[i, j] = mean(fit$progress$timing)
    test.loglik.smart[i, j] = compute_loglik_per_test_dat(fit, test.data, fit$V)
  }
  
  for (j in 1:length(random)){
    fit = random[[j]]
    iteration.random[i, j] = nrow(fit$progress)
    time.random[i, j] = mean(fit$progress$timing)
    test.loglik.random[i, j] = compute_loglik_per_test_dat(fit, test.data, fit$V)
  }
}

average_test_loglik <- c(apply(test.loglik.smart, 2, mean), apply(test.loglik.random, 2, mean))
average_iteration <- c(apply(iteration.smart, 2, mean), apply(iteration.random, 2, mean))
average_time <- c(apply(time.smart, 2, mean), apply(time.random, 2, mean))


