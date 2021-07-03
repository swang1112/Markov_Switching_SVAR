rm(list = ls())
library(magrittr)
library(kdensity)
library(optimParallel)
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


par_optim = c( 0.09787152,  0.55440439,  0.09156761,  1.47105007, -0.11356540,  1.06635796,
               0.39450519,0.85060587)

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


# rest --------------------------------------------------------------------
cores = 7

cl = makeCluster(cores)  
setDefaultCluster(cl=cl)
startpars = c(runif(6, 0, pi/2), .5, .5)
erg_list = optimParallel(startpars,fn = loglike_MS_ICA, r = u, 
                         C = C,
                         lower = c(rep(0,6), rep(0,2)), 
                         upper = c(rep(pi/2, 6), rep(1,2)), 
                         parallel = cl)
stopCluster(cl)