#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// matrix exponentials
// [[Rcpp::export]]
arma::mat matexp(arma::mat X, int n)
{
  if (n == 0)
  {
    return arma::eye(X.n_cols, X.n_rows);
  } else if (n == 1)
  {
    return X;
  } else 
  {
    return X * matexp(X, n-1);
  }
}


// IRFs
// [[Rcpp::export]]
List IRF_fast(arma::mat &A_hat, arma::mat &B_hat, int &horizon)
{
 int K = A_hat.n_rows;
  int p = A_hat.n_cols / K;
  List Out(horizon);
  Out[0] = B_hat;
  if (p == 1)
  {
    for(int i = 1; i < horizon; i++)
    {
      Out[i] = matexp(A_hat, i) * B_hat;
    }
    return Out;
  } else
  {
    arma::mat Mm(K * p, K * p, arma::fill::zeros);
    Mm.submat(0, 0, (K - 1), (K * p - 1)) = A_hat;
    Mm.submat(K, 0, (K * p - 1), ((p - 1) * K - 1)) = arma::eye(K * (p - 1), K * (p - 1));
    arma::mat Mm1(K * p, K * p, arma::fill::eye);
    for (int i = 0; i < (horizon - 1); i++)
    {
      Mm1 = Mm1 * Mm;
      Out[i + 1] = Mm1.submat(0, 0, (K - 1), (K - 1)) * B_hat;
    }
    return Out;
  }
  
}


//ich passe das hier noch an
// [[Rcpp::export]]
double loglike_Normal_Gumbel_GumbelSurvival(const arma::vec& bekk, const arma::vec& theta, const arma::mat& r) {
  // Log-Likelihood function
  
  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * pow(n, 2) + n * (n + 1)/2;
  //Switching Probabilities
  double p11 = theta[0];
  double p12 = 1-p11;
  double p22 = theta[1];
  double p21 = 1-p22;
  double par_gumbel =theta[6];
  double par_gumbelSurvival =theta[7];
  //check for identification
  if(valid_bekk(C,A,G)==false || p11<=0 || p12<=0 || p11>=1 || p12>=1 || p22<=0 || p23<=0 || p22>=1 || p23>=1  || p33<=0 || p32<=0 || p33>=1 || p32>=1||par_clayton<=0 || par_clayton >=17 || par_gumbel<= 1|| par_gumbel>=17){
    return -1e25;
  }
  
  // compute inital H
  arma::mat H = (r.t() * r) / r.n_rows;
  
  arma::mat CC  = C * C.t();
  arma::mat At  = A.t();
  arma::mat Gt  = G.t();
  
  
  //compute covariance matrix
  arma::mat cor_gumbel=arma::cor(BiCopSim_cpp(10000,4,par_gumbel));
  arma::mat cor_gumbel_chol=arma::chol(cor_gumbel).t();
  
  double gumbel_det = arma::det(cor_gumbel_chol);
  
  arma::mat cor_gumbelSurvival=arma::cor(BiCopSim_cpp(10000,14,par_gumbelSurvival));
  arma::mat cor_gumbelSurvival_chol=arma::chol(cor_gumbelSurvival).t();
  
  double gumbelSurvival_det = arma::det(cor_gumbelSurvival_chol);
  
  arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
  
  
  
  double p1t = 0;
  double p2t = 1;
  double p3t = 0;
  
  double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  // double d_norm_2=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t()));
  
  arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(0).t());
  arma::vec p_norm_2=pnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(0).t());
  
  double llv = log(arma::as_scalar((p1t*p+(1-p1t)*(1-q))*arma::det(H_eigen_inv)*arma::det(cor_clayton_chol)*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),3,par_clayton)));
  
  
  for (int i = 1; i < NoOBs; i++) {
    H = CC + At * r.row(i - 1).t() * r.row(i - 1) * A + Gt * H * G;
    
    arma::mat H_eigen_inv=arma::inv( eigen_value_decomposition(H));
    
    double H_det=arma::det(H_eigen_inv);
    
    double d_norm_1=arma::prod(dnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t()));
    double d_norm_2=arma::prod(dnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t()));
    
    arma::vec p_norm_1=pnorm_cpp(cor_gumbel_chol*H_eigen_inv*r.row(i).t());
    arma::vec p_norm_2=pnorm_cpp(cor_gumbelSurvival_chol*H_eigen_inv*r.row(i).t());
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21+p3t*p31)*arma::prod(dnorm_cpp(r))*H_det );
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22+p3t*p31)*d_norm_1*BiCopPDF_cpp(p_norm_1(0),p_norm_1(1),4,par_gumbel)*H_det* gumbel_det);
    double llv_temp3 = arma::as_scalar((p1t*p13+p2t*p23+p3t*p33)*d_norm_2*BiCopPDF_cpp(p_norm_2(0),p_norm_2(1),14,par_gumbelSurvival)*H_det*gumbelSurvival_det);
      
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2+llv_temp3);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = llv_temp3/llv_temp;
    llv+=log(llv_temp);
    
    
  }
  return llv;
}


 
