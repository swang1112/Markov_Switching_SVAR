#include <RcppArmadillo.h>
#include <omp.h>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;


// Gaussian kernel density estimation at point x
// [[Rcpp::export]]
double kl_gauss(double &x, arma::vec & X, int & Tob, double &h)
{
  arma::vec x_vec(Tob);
  x_vec.fill(x);
  arma::vec t = x_vec - X;
  arma::vec out = arma::exp(arma::pow(t/h, 2) * (-0.5)) / sqrt(2*3.14159265358979323846);
  return arma::accu(out)/(Tob*h);
}

// Gaussian kernel density estimation at all points in X
// [[Rcpp::export]]
arma::vec kl_gauss_vec(arma::vec  X, int & Tob, double &h)
{
  arma::vec out(Tob);
  for (int i = 0; i < Tob; i++)
  {
    out(i) = kl_gauss(X(i), X, Tob, h);
  }
  return out;
}


// negative Log Kernel Likelihood
// [[Rcpp::export]]
double nlkl_eps( arma::mat& eps) 
{
  
  int K     = eps.n_cols;
  int Tob   = eps.u.n_rows;
  double bw = (4/(3*Tob))^(1/5);

  double llv = 0.00;
  #pragma omp parallel for shared(llv, eps) reduction(-: llv)
  for (int k = 0; k < K; k++)
  {
    llv-=sum(log(kl_gauss_vec(eps.col(k), Tob, bw)));
  }

  return llv;
}

/*
// kernel density estimation
NumericVector kdensity(arma::vec r)
{
  Environment pkg = Rcpp::Environment::namespace_env("kdensity");
  Function f = pkg["kdensity"];
  Function res =f(r);
  NumericVector density_res=res(r);
  return density_res;
}
*/