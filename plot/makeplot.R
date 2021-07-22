rm(list = ls())
library(magrittr)
library(kdensity)
library(pbapply)
setwd('../run/')
source('reduced_form.R')
Rcpp::sourceCpp('../fun/MS_ICA.cpp')

# variable(oil production, real activities, oil price)
Sig = crossprod(u)/(var3$obs - 3 * var3$p - 1)
C   = Sig %>% chol %>% t

# erg_list  = readRDS( '../out/a_grid.rds')
# erg_list  = readRDS( '../out/a_grid_fein.rds')
# erg_list  = readRDS( '../out/a_M3.rds')
erg_list  = readRDS( '../out/a_M2_ngml.rds')
# erg_vals  = erg_list %>% lapply('[[', 2) %>% unlist
erg_vals  = erg_list %>% lapply('[[', 1) %>% unlist
erg_optim = erg_list[[which.min(erg_vals)]]
# par_optim = erg_optim$par
par_optim = erg_optim$estimate

# states distribution -----------------------------------------------------
# TProb = filter_MS_ICA(par_optim, u, C, init = c(.5, .5))
# pdf("../plot/TP_filtered.pdf", width = 12, height = 6)
# TProb = filter_MS_ICA_M3(par_optim, u, C, init = c(1/3, 1/3, 1/3))
# pdf("../plot/TP_filtered_M3.pdf", width = 12, height = 6)
pdf("../plot/TP_filtered_M2_ngml.pdf", width = 12, height = 6)
TProb = filter_MS_ngml(par_optim, u, init = c(.5, .5))
par(mfrow = c(2,1))
par(mar = c(2,2,1,1))
plot(data_time[-1], TProb[-1,1], type = 'l', main = 'state1')
plot(data_time[-1], TProb[-1,2], type = 'l', main = 'state2')
#plot(data_time[-1], TProb[-1,3], type = 'l', main = 'state3')
par(mfrow = c(1,1))
dev.off()

# structural multiplier ---------------------------------------------------
theta_optim1 = par_optim[1:3]
theta_optim2 = par_optim[4:6]
(B1 = getB(theta_optim1, C))
(B2 = getB(theta_optim2, C))

# shock(oil supply, aggregate demand, oil stock demand)
B1 = B1[,c(2,3,1)]; 
B1[,1] = B1[,1]*-1
B1 


B2[,c(1,2)] = B2[,c(1,2)] * -1
B2 = B2[,c(1,3,2)]
B2


# ngml
source('../fun/ng_ica.R')
theta_optim1 = par_optim[1:9]
theta_optim2 = par_optim[10:18]
K = 3
(str_par1 = get_Bt(theta_optim1, C))
(str_par2 = get_Bt(theta_optim2, C))
B1 = str_par1$B
B2 = str_par2$B

B1[,c(1,3)] = B1[,c(1,3)]*-1
B2[,1] = B2[,1]*-1

# IRF ---------------------------------------------------------------------
Rcpp::sourceCpp('../fun/IRF.cpp')
source('../fun/trans_irf.R')
irf_1 = IRF_fast(Bcoef, B1, 36) %>% trans_irf 
irf_2 = IRF_fast(Bcoef, B2, 36) %>% trans_irf 

ip1 = irf_1 %>% plot_irf() + ggtitle('state1')
ip2 = irf_2 %>% plot_irf() + ggtitle('state2')

library(gridExtra)
margin = theme(plot.margin = unit(c(1, 1, -4, 0), "mm"))
ips    = grid.arrange(grobs = lapply(list(ip1, ip2), "+", margin))
ggsave('../plot/ips.pdf', ips, width = 12, height = 8)
