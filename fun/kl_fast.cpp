#include <RcppArmadillo.h>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// choose 2 int from N (as index)
arma::imat choose2_fast(int &N)
{
  arma::imat Out(N * (N - 1) / 2, 2);
  int c = 0;
  for (int i = 0; i < (N - 1); i++)
  {
    for (int j = i + 1; j < N; j++)
    {
      Out(c, 0) = i;
      Out(c, 1) = j;
      c++;
    }
  }
  return Out;
}

// Givens rotation
arma::mat givensQ_fast(arma::vec &thetas, int K)
{
  arma::mat Out = arma::eye(K, K);
  arma::imat Cmat = choose2_fast(K);
  for (int i = 0; i < K * (K - 1) / 2; i++)
  {
    arma::mat temp = arma::eye(K, K);
    temp(Cmat(i, 0), Cmat(i, 0)) = cos(thetas(i));
    temp(Cmat(i, 1), Cmat(i, 1)) = cos(thetas(i));
    temp(Cmat(i, 0), Cmat(i, 1)) = -sin(thetas(i));
    temp(Cmat(i, 1), Cmat(i, 0)) = sin(thetas(i));
    
    Out = Out * temp;
  }
  return Out.t();
}

NumericVector kdensity(arma::vec r){
  Environment pkg = Rcpp::Environment::namespace_env("kdensity");
  Function f = pkg["kdensity"];
  Function res =f(r);
  NumericVector density_res=res(r);
  return density_res;
}

// kernel likelihood: L(theta)
// [[Rcpp::export]]
double kl_fast(arma::vec & theta, arma::mat & u,  arma::mat & C)
{

  int K = u.n_cols;
  arma::mat Q = givensQ_fast(theta, K);
  arma::mat E = arma::inv(C * Q) * u.t();
  double out = 0;
  for (int i = 0; i < K; i++)
  {
    out = out + sum(log(kdensity(E.row(i).t())));
  }
  return out;
}

