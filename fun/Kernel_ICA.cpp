#include <RcppArmadillo.h>
#include <omp.h>

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
arma::mat givensQ_fast(arma::vec thetas, int K)
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

// get mixing matrix
// [[Rcpp::export]]
arma::mat getB(arma::vec thetas, arma::mat& C)
{
  int K = C.n_cols;
  arma::mat Q = givensQ_fast(thetas, K);
  return C * Q;
}

// Gaussian kernel density estimation at point x
double kl_gauss(double &x, arma::vec & X, int & Tob, double &h)
{
  arma::vec x_vec(Tob);
  x_vec.fill(x);
  arma::vec t = x_vec - X;
  arma::vec out = arma::exp(arma::pow(t/h, 2) * (-0.5)) / sqrt(2*3.14159265358979323846);
  return arma::accu(out)/(Tob*h);
}

// Gaussian kernel density estimation at all points in X
arma::vec kl_gauss_vec(arma::vec  X, int & Tob, double &h)
{
  arma::vec out(Tob);
  for (int i = 0; i < Tob; i++)
  {
    out(i) = kl_gauss(X(i), X, Tob, h);
  }
  return out;
}

// Negative Log-Likelihood function
// [[Rcpp::export]]
double fast_kernel_ICA(arma::vec& par,  arma::mat& u,  arma::mat& C) 
{
  
  int K     = u.n_cols;
  int Tob   = u.n_rows;
  
  arma::mat B   = getB(par, C);
  arma::mat eps = arma::inv(B) * u.t();
  eps = eps.t();
  
  double bw = pow(4.0/(3.0*Tob), 0.2);
  
  double llv = 0.00;
  #pragma omp parallel for shared(llv, eps) reduction(-: llv)
  for (int k = 0; k < K; k++)
  {
    llv-=arma::accu(arma::log(kl_gauss_vec(eps.col(k), Tob, bw)));
  }
  
  return llv;
}


