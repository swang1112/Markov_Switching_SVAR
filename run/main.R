rm(list = ls())
library(magrittr)
library(kdensity)
library(pbapply)
source('reduced_form.R')
Rcpp::sourceCpp('../fun/MS_ICA.cpp')
start_time = Sys.time()

Core = 108

# variable(oil production, real activities, oil price)
Sig = crossprod(u)/(var3$obs - 3 * var3$p - 1)
C   = Sig %>% chol %>% t

grid_prob  = seq(0, 1, by = 0.1) 
start_prob = expand.grid(grid_prob, grid_prob) %>% as.matrix()
startpars = list()
for (i in 1:nrow(start_prob)) {
  startpars[[i]] = c(runif(6, 0, pi/2), start_prob[i,])
}

#startpars = c(runif(6, 0, pi/2), .9, .9)
#loglike_MS_ICA(startpars[[109]], r = u, C = C, init = c(.5, .5))

erg_list = pblapply(startpars, optim, fn = loglike_MS_ICA, r = u, C = C, init = c(.5, .5), cl = Core)
saveRDS(erg_list, '../out/a_grid.rds')

end_time = Sys.time()
elapsed_time = difftime(end_time, start_time)
Journal = data.frame(start_time, end_time, elapsed_time)
write.csv(Journal, "../out/erg_list.csv")

# erg_list = optim(startpars,fn = loglike_MS_ICA, r = u, C = C, init = c(.5, .5))
# erg_list$value
# erg_list$par
# par_optim = erg_list$par

# states distribution -----------------------------------------------------
# TProb = filter_MS_ICA(par_optim, u, C, init = c(.5, .5))
# par(mfrow = c(2,1))
# par(mar = c(2,2,1,1))
# plot(data_time[-1], TProb[-1,1], type = 'l', main = 'state1')
# plot(data_time[-1], TProb[-1,2], type = 'l', main = 'state2')
# par(mfrow = c(1,1))


# structural multiplier ---------------------------------------------------
# theta_optim1 = par_optim[1:3]
# theta_optim2 = par_optim[4:6]
# (B1 = getB(theta_optim1, C))
# (B2 = getB(theta_optim2, C))
# 


