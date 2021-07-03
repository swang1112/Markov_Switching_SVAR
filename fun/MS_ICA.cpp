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

// Log-Likelihood function

// [[Rcpp::export]]
double loglike_MS_ICA(const arma::vec& theta, const arma::mat& r, const arma::mat& C) {
 {
  

 \\ dimensions
  int K = r.n_cols;
  
  // Length of each series
  int NoOBs = r.n_rows;
  //rotation angles

  //state 1
  double delta_1_1=theta(0);
  double delta_2_1=theta(1);
  double delta_3_1=theta(2);
  //state 2
  double delta_1_2=theta(3);
  double delta_2_2=theta(4);
  double delta_3_2=theta(5);

  //Switching Probabilities
  double p11 = theta(6);
  double p12 = 1-p11;
  double p22 = theta(7);
  double p21 = 1-p22;

  //check for identification (we might not need this)
  if(p11<=0 || p11>=1 || p22<=0 || p22>= 1){
    return -1e25;
  }

  //initial state
  double p1t = 1/2;
  double p2t = 1/2;
  
  // rotation matrices
  arma::mat Q_1 = givensQ_fast(theta(0:2),K);
  arma::mat Q_2 = givensQ_fast(theta(3:5),K);
  
  //transform series (B = CQ)
  arma::mat series_state1 = arma::inv(C * Q_1)*r.t();
  arma::mat series_state2 = arma::inv(C * Q_2)*r.t();
  
  //KDE for the first state
  arma::vec dens_1_1= kdensity(*series_state1.row(0).t());
  arma::vec dens_2_1= kdensity(*series_state1.row(1).t());
  arma::vec dens_3_1= kdensity(*series_state1.row(2).t());
  
  //KDE for the first state
  arma::vec dens_1_2= kdensity(*series_state2.row(0).t());
  arma::vec dens_2_2= kdensity(*series_state2.row(1).t());
  arma::vec dens_3_2= kdensity(*series_state2.row(2).t());

  double llv = log(arma::as_scalar(p1t*dens_1_1(0)*dens_2_1(0)*dens_3_1(0) + p2t* dens_1_2(0)*dens_2_2(0)*dens_3_2(0) ));


  for (int i = 1; i < NoOBs; i++) {
    
    
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21)*dens_1_1(i)*dens_2_1(i)*dens_3_1(i));
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22)*dens_1_2(i)*dens_2_2(i)*dens_3_2(i));

    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;

    llv+=log(llv_temp);


  }
  return llv;
}