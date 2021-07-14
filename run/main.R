rm(list = ls())
library(magrittr)
library(pbapply)
source('reduced_form.R')

start_time = Sys.time()

# \beginn{config} ---------------------------------------------------------
M     = 1
Core  = 108
Jaki  = TRUE
# \end{config} ------------------------------------------------------------

if (M == 1) Rcpp::sourceCpp('../fun/Kernel_ICA.cpp') else Rcpp::sourceCpp('../fun/MS_ICA.cpp')
if (Jaki) require(kdensity)

Nrot  = choose(3,2)*M

# variable(oil production, real activities, oil price)
Sig = crossprod(u)/(var3$obs - 3 * var3$p - 1)
C   = Sig %>% chol %>% t

if (M == 2){
  grid_prob  = seq(0, 1, by = 0.01) 
  start_prob = expand.grid(grid_prob, grid_prob) %>% as.matrix()
  startpars = list()
  set.seed(1234)
  for (i in 1:nrow(start_prob)) {
    startpars[[i]] = c(runif(Nrot, 0, pi/2), start_prob[i,])
  }
  # run
  erg_list = pblapply(startpars, optim, fn = loglike_MS_ICA, r = u, C = C, init = c(.5, .5), cl = Core)
  saveRDS(erg_list, '../out/a_grid_fein.rds')
  
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
    erg_list = pblapply(startpars, optim, fn = fast_kernel_ICA, u = u, C = C,
                      control = list(maxit = 25000, tmax = 100, abstol = 1e-30), 
                      cl = 7)
    saveRDS(erg_list, '../out/a_M1.rds')
  }
}


#loglike_MS_ICA_M3(startpars[[sample(1:nrow(start_prob), 1)]], r = u, C = C, init = c(1/3, 1/3, 1/3))

end_time = Sys.time()
elapsed_time = difftime(end_time, start_time)
Journal = data.frame(start_time, end_time, elapsed_time)
write.csv(Journal, paste0("../out/a_M", M, "_", format(Sys.Date(), "%m_%d"), "_journal.csv"))

