rm(list = ls())
library(magrittr)
library(kdensity)
library(pbapply)
source('reduced_form.R')
Rcpp::sourceCpp('../fun/MS_ICA.cpp')


# variable(oil production, real activities, oil price)
Sig = crossprod(u)/(var3$obs - 3 * var3$p - 1)
C   = Sig %>% chol %>% t

erg_list  = readRDS( '../out/a_grid.rds')
erg_vals  = erg_list %>% lapply('[[', 2) %>% unlist
erg_optim = erg_list[[which.min(erg_vals)]]
par_optim = erg_optim$par


# states distribution -----------------------------------------------------
TProb = filter_MS_ICA(par_optim, u, C, init = c(.5, .5))
par(mfrow = c(2,1))
par(mar = c(2,2,1,1))
plot(data_time[-1], TProb[-1,1], type = 'l', main = 'state1')
plot(data_time[-1], TProb[-1,2], type = 'l', main = 'state2')
par(mfrow = c(1,1))


# structural multiplier ---------------------------------------------------
# theta_optim1 = par_optim[1:3]
# theta_optim2 = par_optim[4:6]
# (B1 = getB(theta_optim1, C))
# (B2 = getB(theta_optim2, C))
# 

