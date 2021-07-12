rm(list = ls())
#source('reduced_form.R')

library(magrittr)
library(kdensity)
Rcpp::sourceCpp('../fun/kl_fast.cpp')
# Rcpp::sourceCpp('../fun/foo.cpp')
# fa = runif(6,0,1)
# fa
# getMarkov(fa)

# fun ---------------------------------------------------------------------
# this function simulates unobserved independent component from 
# standardized t distribution and mixes them with a B matrix
# the output is the observable
make_dat = function(x, Tob, dof, B, dist = 'tdist'){
  Epsilon = matrix(NA, nrow = Tob, ncol = 2)
  for (k in 1:2) {
    Epsilon[,k] = rt(Tob, dof[k])*sqrt((dof[k] - 2)/dof[k])
  }
  u = B %*% t(Epsilon)
  t(u)
}

# this function calculates kernel log likelihood
get_kl0 = function(theta, u, C, kernel = "gaussian", bw = NULL){
  
  Q  = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, ncol = 2, byrow = T)
  xx = t(solve(C %*% Q) %*% t(u))
  out = 0
  for (i in 1:2) {
    kd  = kdensity(xx[,i])#, kernel = kernel, bw = bw)
    out = out + sum(log(kd(xx[,i])))
  }
  out
}


# DGP ---------------------------------------------------------------------
Tob   = 60
dof   = c(3,5)
rho   = 0.4 #cov(u1,u2)
theta = pi/3
C     = matrix(c(1, 0, rho, sqrt(1-rho^2)), nrow = 2, ncol = 2, byrow = T)
Q     = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, ncol = 2, byrow = T)
B     = C %*% Q

U = make_dat(0, Tob = Tob, dof = dof, B = B)


# calculate kernel log likelihood at true theta ---------------------------
get_kl0(theta, u = U, C = C)
kl_fast(theta, u = U, C = C)
kl_fast2D(theta, u = U, C = C)



# compare speed -----------------------------------------------------------
bench = microbenchmark::microbenchmark("R" = get_kl0(theta, u = U, C = C),
                            "Rcpp1" = kl_fast(theta, u = U, C = C), 
                            "Rcpp2" = kl_fast2D(theta, u = U, C = C), 
                            times = 100)
require(ggfortify)
bench %>% autoplot
bench %>% show



