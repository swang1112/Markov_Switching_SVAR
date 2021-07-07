#include <RcppArmadillo.h>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// kernel density estimation
NumericVector kdensity(arma::vec r)
{
  Environment pkg = Rcpp::Environment::namespace_env("kdensity");
  Function f = pkg["kdensity"];
  Function res =f(r);
  NumericVector density_res=res(r);
  return density_res;
}

// negative Log Kernel Likelihood
// [[Rcpp::export]]
double log_eps( arma::vec& eps) 
{
  
  int K = eps.n_cols;

  double llv = 0.00;

  for (int i = 0; i < K; i++)
  {
   llv-=log(kdensity(eps.col(i)));
  }

  return llv;
}

