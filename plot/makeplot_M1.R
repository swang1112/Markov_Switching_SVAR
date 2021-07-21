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
# erg_list = readRDS( '../out/a_M1.rds')
# erg_list = readRDS( '../out/a_M1_Jaki.rds')
# erg_vals  = erg_list %>% lapply('[[', 2) %>% unlist
# erg_optim = erg_list[[which.min(erg_vals)]]
# par_optim = erg_optim$par

# shock(oil supply, aggregate demand, oil stock demand)
# (B = getB(par_optim, C))

erg_list = readRDS( '../out/a_M1_ng.rds')
D = erg_list$D
(B = erg_list$B)

B[,1] = B[,1]*-1

# check
B %*% D %*% t(B)
Sig

irf = IRF_fast(Bcoef, B, 36) %>% trans_irf 
ip  = irf %>% plot_irf() 
ip
# ggsave('../plot/ip_M1.pdf', ip, width = 12, height = 8)
# ggsave('../plot/ip_M1_Jaki.pdf', ip, width = 12, height = 8)
ggsave('../plot/ip_M1_ng.pdf', ip, width = 12, height = 8)

# dcov --------------------------------------------------------------------
dcov = var3 %>% svars::id.dc()
dcov$B[,1] = dcov$B[,1]*-1
dcov_irf = IRF_fast(Bcoef, dcov$B, 36) %>% trans_irf 
dcov_ip  = dcov_irf %>% plot_irf() 
dcov_ip
ggsave('../plot/ip_M1_dcov.pdf', dcov_ip, width = 12, height = 8)

library(steadyICA)
pvals = vector("numeric", length = 45)
.get_pval = function(x){
  xx = t(solve(B) %*% t(u))
  permTest(xx, group=1:3, R=199, FUN=c('gmultidcov'), alpha=1, symmetric=TRUE)
}
pvals = pblapply(pvals, .get_pval, cl = 7)
pvals = pvals %>% unlist()
pvals %>% hist(10)
