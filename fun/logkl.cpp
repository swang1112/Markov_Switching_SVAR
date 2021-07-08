#include <RcppArmadillo.h>
#include <omp.h>
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
double nlkl_eps( arma::mat& eps) 
{
  
  int K = eps.n_cols;
  
  double llv = 0.00;
  #pragma omp parallel for shared(llv, eps) reduction(-: llv)
  for (int k = 0; k < K; k++)
  {
    llv-=sum(log(kdensity(eps.col(k))));
  }

  return llv;
}
