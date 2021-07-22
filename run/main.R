rm(list = ls())
graphics.off()
library(magrittr)
library(pbapply)
source('reduced_form.R')

start_time = Sys.time()

# \beginn{config} ---------------------------------------------------------
M     = 2
Core  = 64
NGML  = TRUE
Jaki  = FALSE
# \end{config} ------------------------------------------------------------

if (M == 1) {
  Rcpp::sourceCpp('../fun/Kernel_ICA.cpp')
  source('../fun/ng_ica.R')
} else Rcpp::sourceCpp('../fun/MS_ICA.cpp')

if (Jaki) require(kdensity)

Nrot  = choose(3,2)*M
K     = 3

# variable(oil production, real activities, oil price)
Sig   = crossprod(u)/(var3$obs - 3 * var3$p - 1)
C     = Sig %>% chol %>% t
u_st  = t(solve(C) %*% t(u))

if (M == 2){
  grid_prob  = seq(0.5, 1, by = 0.005) 
  start_prob = expand.grid(grid_prob, grid_prob) %>% as.matrix()
  startpars = list()
  set.seed(1234)
  if (NGML){
    for (i in 1:nrow(start_prob)) {
      startpars[[i]] = c(runif(Nrot, 0, pi/2), 
                         sort(rchisq(K*2,3), decreasing = T),
                         runif(K*2, 2, 20),
                         start_prob[i,])
    }
    
    # run
    erg_list = pblapply(startpars, nlm, f = loglike_ngMS_ICA, u_st = u_st, init = c(.5, .5),
                        gradtol = 1e-10, iterlim = 10000, cl = Core)
    saveRDS(erg_list, '../out/a_M2_ngml.rds')
  } else {
    for (i in 1:nrow(start_prob)) {
      startpars[[i]] = c(runif(Nrot, 0, pi/2), start_prob[i,])
    }
    # run
    erg_list = pblapply(startpars, optim, fn = loglike_MS_ICA, r = u, C = C, init = c(.5, .5), cl = Core)
    saveRDS(erg_list, '../out/a_grid_fein.rds')
  }
} else if (M == 3){
  grid_prob  = seq(0, 0.6, by = 0.1) 
  start_prob = expand.grid(grid_prob, grid_prob, grid_prob, grid_prob, grid_prob, grid_prob) %>% as.matrix()
  startpars = list()
  set.seed(1234)
  for (i in 1:nrow(start_prob)) {
    startpars[[i]] = c(runif(Nrot, 0, pi/2), start_prob[i,])
  }
  # run
  erg_list = pblapply(startpars, optim, fn = loglike_MS_ICA_M3, r = u, C = C, init = c(1/3, 1/3, 1/3), cl = Core)
  saveRDS(erg_list, '../out/a_M3.rds')
} else if (M == 1){
  set.seed(1234)
  Start_N = 100
  startpars = list()
  for (i in 1:Start_N) {
    startpars[[i]] = runif(Nrot, 0, pi/2)
  }
  if (Jaki) {
    erg_list = pblapply(startpars, optim, fn = Jaki_kernel_ICA, u = u, C = C,
                        control = list(maxit = 25000, tmax = 100, abstol = 1e-30), 
                        cl = 7)
    saveRDS(erg_list, '../out/a_M1_Jaki.rds')
  } else {
    erg_list = ng_ica(u)
    saveRDS(erg_list, '../out/a_M1_ng.rds')
  }
}


#loglike_MS_ICA_M3(startpars[[sample(1:nrow(start_prob), 1)]], r = u, C = C, init = c(1/3, 1/3, 1/3))

end_time = Sys.time()
elapsed_time = difftime(end_time, start_time)
Journal = data.frame(start_time, end_time, elapsed_time)
write.csv(Journal, paste0("../out/a_M", M, "_", format(Sys.Date(), "%m_%d"), "_journal.csv"))

