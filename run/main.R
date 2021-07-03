rm(list = ls())
library(magrittr)
library(kdensity)

source('reduced_form.R')
Rcpp::sourceCpp('../fun/MS_ICA.cpp')
#Rcpp::sourceCpp('../fun/kl_fast.cpp')

# variable(oil production, real activities, oil price)
Sig = crossprod(u)/(var3$obs - 3 * var3$p - 1)
C   = Sig %>% chol %>% t

startpars = c(runif(6, 0, pi/2), .9, .9)
loglike_MS_ICA(startpars, r = u, C = C, init = c(.5, .5))

erg_list = optim(startpars,fn = loglike_MS_ICA, r = u, C = C, init = c(.5, .5))
erg_list$value
erg_list$par
par_optim = erg_list$par

theta_optim1 = par_optim[1:3]
theta_optim2 = par_optim[4:6]
(B1 = getB(theta_optim1, C))
(B2 = getB(theta_optim2, C))

TProb = filter_MS_ICA(par_optim, u, C, init = c(.5, .5))
par(mfrow = c(2,1))
par(mar = c(2,2,1,1))
plot(data_time[-1], TProb[-1,1], type = 'l', main = 'state1')
plot(data_time[-1], TProb[-1,2], type = 'l', main = 'state2')
par(mfrow = c(1,1))

