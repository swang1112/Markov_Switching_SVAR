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
  return Out;
}

NumericVector kdensity(arma::vec r){
  Environment pkg = Rcpp::Environment::namespace_env("kdensity");
  Function f = pkg["kdensity"];
  Function res =f(r);
  NumericVector density_res=res(r);
  return density_res;
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

// Multivariate Gaussian kernel density estimation at vector x (K x 1)
// variables must be orthogonal and have homogeneous 2nd moment
// data matrix X is (Tob X K)
// [[Rcpp::export]]
double Mkl_gauss(arma::rowvec x, arma::mat & X, int &K, int & Tob, double &h)
{
  arma::mat x_mat(Tob, K);
  x_mat.each_row() = x;
  arma::mat D = (x_mat - X)/h;
  
  double out = 0.0;
  for (int i = 0; i < Tob; i++)
  {
    out+= exp(arma::as_scalar(D.row(i) * D.row(i).t()) * (-0.5)) / sqrt(pow(2*3.14159265358979323846, K));
  }  
  return out/(Tob*pow(h, K));
}

// Multivariate Gaussian kernel density of observing all data in X (Tob X K)
// [[Rcpp::export]]
arma::vec Mkl_gauss_vec( arma::mat X, int &K, int & Tob, double &h)
{
  arma::vec out(Tob);
  for (int i = 0; i < Tob; i++)
  {
    out(i) = Mkl_gauss(X.row(i), X, K, Tob, h);
  }
  return out;
}


// kernel likelihood: L(theta)
// [[Rcpp::export]]
double kl_fast(arma::vec & theta, arma::mat & u,  arma::mat & C)
{

  int K     = u.n_cols;
  int Tob   = u.n_rows;
  double bw = pow(4.0/(3.0*Tob), 0.2);
  
  arma::mat Q = givensQ_fast(theta, K);
  arma::mat E = arma::inv(C * Q) * u.t();
  double out = 0;
  for (int i = 0; i < K; i++)
  {
    out = out + arma::accu(arma::log(kl_gauss_vec(E.row(i).t(), Tob, bw)));
  }
  return out;
}

// kernel likelihood: L(theta) for bivariate model
// [[Rcpp::export]]
double kl_fast2D(double & theta, arma::mat & u,  arma::mat & C)
{
  int Tob   = u.n_rows;
  double bw = pow(4.0/(3.0*Tob), 0.2);
  
  // K = 2
  arma::mat Q = arma::mat(2,2);
  Q(0, 0) = cos(theta);
  Q(1, 1) = cos(theta);
  Q(0, 1) = -sin(theta);
  Q(1, 0) = sin(theta);
  
  arma::mat E = arma::inv(C * Q) * u.t();
  double out = 0;
  for (int i = 0; i < 2; i++)
  {
    out = out + arma::accu(arma::log(kl_gauss_vec(E.row(i).t(), Tob, bw)));
  }
  return out;
}


/*
// kernel likelihood: L(theta) for trivariate model (tbd)
// [[Rcpp::export]]
double kl_fast3D(arma::vec & theta, arma::mat & u,  arma::mat & C)
{
  
  // K = 3
  arma::mat Q = arma::eye(3, 3);
  for (int i = 0; i < 3; i++)
  {
    arma::mat temp = arma::eye(3, 3);
    temp(0, 0) = cos(theta(i));
    temp(1, 1) = cos(theta(i));
    temp(0, 1) = -sin(theta(i));
    temp(1, 0) = sin(theta(i));
    Q = Q * temp;
  }
  arma::mat E = arma::inv(C * Q) * u.t();
  double out = 0;
  for (int i = 0; i < 3; i++)
  {
    out = out + sum(log(kdensity(E.row(i).t())));
  }
  return out;
}
*/
