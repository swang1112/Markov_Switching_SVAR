rm(list = ls())
library(magrittr)
library(kdensity)

source('reduced_form.R')
Rcpp::sourceCpp('../fun/MS_ICA.cpp')
Rcpp::sourceCpp('../fun/kl_fast.cpp')
