#include <RcppArmadillo.h>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// get markov matrix
// [[Rcpp::export]]
arma::mat getMarkov(arma::vec thetas)
{
  arma::mat Out(3,3);
  Out.row(0) = {1-thetas(0)-thetas(1), thetas(0), thetas(1)};
  Out.row(1) = {thetas(2), 1-thetas(2)-thetas(3), thetas(3)};
  Out.row(2) = {thetas(4), thetas(5), 1-thetas(4)-thetas(5)};
  Rcout << arma::all(arma::all(Out > 0  && Out < 1)) << '\n';
  return Out;
}

