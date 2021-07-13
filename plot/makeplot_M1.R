rm(list = ls())
library(magrittr)
setwd('../run/')
source('reduced_form.R')
source('../fun/trans_irf.R')
Rcpp::sourceCpp('../fun/IRF.cpp')
Rcpp::sourceCpp('../fun/Kernel_ICA.cpp')

# variable(oil production, real activities, oil price)
Sig = crossprod(u)/(var3$obs - 3 * var3$p - 1)
C   = Sig %>% chol %>% t

## M = 1
erg_optim = readRDS( '../out/a_M1.rds')
par_optim = erg_optim$par

# shock(oil supply, aggregate demand, oil stock demand)
(B = getB(par_optim, C))
B[,c(2,3)] = B[,c(2,3)]*-1
B[,c(1,3)] = B[,c(3,1)]
B

irf = IRF_fast(Bcoef, B, 36) %>% trans_irf 
ip  = irf %>% plot_irf() 
ggsave('../plot/ip_M1.pdf', ip, width = 12, height = 8)
