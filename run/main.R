rm(list = ls())
library(magrittr)
library(kdensity)
library(pbapply)
source('reduced_form.R')

start_time = Sys.time()

# \beginn{config} ---------------------------------------------------------
M     = 1
Core  = 108
# \end{config} ------------------------------------------------------------

if (M == 1) Rcpp::sourceCpp('../fun/Kernel_ICA.cpp') else Rcpp::sourceCpp('../fun/MS_ICA.cpp')

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
  startpars = runif(Nrot, 0, pi/2)
  Rcpp::sourceCpp('../fun/Kernel_ICA.cpp')
  fast_kernel_ICA(startpars, u, C)
  erg = optim(startpars, fn = fast_kernel_ICA, u = u, C = C)
  saveRDS(erg, '../out/a_M1.rds')
}


#loglike_MS_ICA_M3(startpars[[sample(1:nrow(start_prob), 1)]], r = u, C = C, init = c(1/3, 1/3, 1/3))

end_time = Sys.time()
elapsed_time = difftime(end_time, start_time)
Journal = data.frame(start_time, end_time, elapsed_time)
write.csv(Journal, paste0("../out/a_M", M, "_", format(Sys.Date(), "%m_%d"), "_journal.csv"))

