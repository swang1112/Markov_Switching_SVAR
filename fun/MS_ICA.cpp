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
// [[Rcpp::export]]
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

// get markov matrix
arma::mat getMarkov(arma::vec & thetas)
{
  arma::mat Out(3,3);
  Out.row(0) = {1-thetas(0)-thetas(1), thetas(0), thetas(1)};
  Out.row(1) = {thetas(2), 1-thetas(2)-thetas(3), thetas(3)};
  Out.row(2) = {thetas(4), thetas(5), 1-thetas(4)-thetas(5)};
  return Out;
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


// Log-Likelihood function
// [[Rcpp::export]]
double loglike_MS_ICA( arma::vec& theta,  arma::mat& r,  arma::mat& C, arma::vec & init) 
{
  
  int K = r.n_cols;
  int NoOBs = r.n_rows;

  //Switching Probabilities
  double p11 = theta(6);
  double p12 = 1-p11;
  double p22 = theta(7);
  double p21 = 1-p22;

  //check for identification 
  if(p11<=0 || p11>=1 || p22<=0 || p22>= 1){
    return 1e25;
  }
  
  //initial state
  double p1t = init(0);
  double p2t = init(1);
  
  // rotation matrices
  arma::mat Q_1 = givensQ_fast(theta.subvec(0, 2), K);
  arma::mat Q_2 = givensQ_fast(theta.subvec(3, 5), K);
  
  //transform series (B = CQ)
  arma::mat series_state1 = arma::inv(C * Q_1)*r.t();
  arma::mat series_state2 = arma::inv(C * Q_2)*r.t();
  
  double bw = pow(4.0/(3.0*NoOBs), 0.2);

  //KDE for the first state
  arma::vec dens_1_1= kl_gauss_vec(series_state1.row(0).t(), NoOBs, bw);
  arma::vec dens_2_1= kl_gauss_vec(series_state1.row(1).t(), NoOBs, bw);
  arma::vec dens_3_1= kl_gauss_vec(series_state1.row(2).t(), NoOBs, bw);
  
  //KDE for the first state
  arma::vec dens_1_2= kl_gauss_vec(series_state2.row(0).t(), NoOBs, bw);
  arma::vec dens_2_2= kl_gauss_vec(series_state2.row(1).t(), NoOBs, bw);
  arma::vec dens_3_2= kl_gauss_vec(series_state2.row(2).t(), NoOBs, bw);

  double llv = -log(arma::as_scalar(p1t*dens_1_1(0)*dens_2_1(0)*dens_3_1(0) + p2t* dens_1_2(0)*dens_2_2(0)*dens_3_2(0) ));

  for (int i = 1; i < NoOBs; i++) {
    
    
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21)*dens_1_1(i)*dens_2_1(i)*dens_3_1(i));
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22)*dens_1_2(i)*dens_2_2(i)*dens_3_2(i));

    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;

    llv-=log(llv_temp);

  }
  return llv;
}

arma::vec tdens_joint(arma::mat &E, int & K, int& Tob, arma::mat & S, arma::mat &Q, arma::vec &lambda)
{
  arma::vec out(Tob);
  arma::vec tdist(K);
  for (int i = 0; i < Tob; ++i) {
    arma::vec shock = E.col(i);
    
    for (int j = 0; j < K; j++) {
      tdist(j) = R::dt(shock(j), lambda(j), 0);
    }
    out(i) = arma::sum(log(tdist)) + log(arma::det(arma::inv_sympd(S) * Q.t()));
  }
  return out;
}

// Log-Likelihood function
// [[Rcpp::export]]
double loglike_ngMS_ICA( arma::vec& theta,  arma::mat & u_st, arma::vec & init) 
{
  
  int K     = u_st.n_cols;
  int Tob   = u_st.n_rows;
  
  // Switching Probabilities
  double p11 = theta(theta.n_elem - 2);
  double p12 = 1-p11;
  double p22 = theta(theta.n_elem - 1);
  double p21 = 1-p22;
  
  // structural parameters
  arma::vec rot = theta.subvec(0, K*(K-1)-1);
  arma::vec sig = theta.subvec(K*(K-1), K*(K-1) + 2*K-1);
  arma::vec lambda = theta.subvec(K*(K-1) + 2*K, theta.n_elem - 3);
  
  //check for identification 
  if(p11<=0 || p11>=1 || p22<=0 || p22>= 1 || any(sig < 0.0)|| any(lambda < 2.0)){
    return 1e25;
  }
  
  //initial state
  double p1t = init(0);
  double p2t = init(1);
  
  // rotation matrices
  arma::mat Q_1 = givensQ_fast(rot.subvec(0, 2), K);
  arma::mat Q_2 = givensQ_fast(rot.subvec(3, 5), K);
  
  // structural covariances
  arma::mat S_1 = arma::diagmat(1/sig.subvec(0, 2));
  arma::mat S_2 = arma::diagmat(1/sig.subvec(3, 5));
  
  // dof
  arma::vec dof_1 = lambda.subvec(0,2);
  arma::vec dof_2 = lambda.subvec(3,5);
  
  // transform series (B = CQ)
  arma::mat E_1 = arma::inv_sympd(S_1) * Q_1.t() * u_st.t();
  arma::mat E_2 = arma::inv_sympd(S_2) * Q_2.t() * u_st.t();
  
  
  // Joint density
  arma::vec dens_M1 = tdens_joint(E_1, K, Tob, S_1, Q_1, dof_1);
  arma::vec dens_M2 = tdens_joint(E_2, K, Tob, S_2, Q_2, dof_2);

  double llv = -log(arma::as_scalar(p1t*exp(dens_M1(0)) + p2t*exp(dens_M2(0)) ));
  
  for (int i = 1; i < Tob; i++) {
    
    
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21)*exp(dens_M1(i)));
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22)*exp(dens_M2(i)));
    
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    
    llv-=log(llv_temp);
    
  }
  return llv;
}

// Log-Likelihood function
// [[Rcpp::export]]
double loglike_MS_ICA_M3( arma::vec& theta,  arma::mat& r,  arma::mat& C, arma::vec & init) 
{
  
  //int K = r.n_cols;
  int NoOBs = r.n_rows;
  
  //Switching Probabilities
  arma::vec lamb = theta.subvec(9, 14);
  arma::mat P    = getMarkov(lamb);
  
  //check for identification 
  if(!arma::all(arma::all(P > 0  && P < 1)) ){
    return 1e25;
  }
  
  //initial state
  double p1t = init(0);
  double p2t = init(1);
  double p3t = init(2);
  
  double bw = pow(4.0/(3.0*NoOBs), 0.2);

  //transform series 
  arma::mat series_state1 = arma::inv(getB(theta.subvec(0, 2), C))*r.t();
  arma::mat series_state2 = arma::inv(getB(theta.subvec(3, 5), C))*r.t();
  arma::mat series_state3 = arma::inv(getB(theta.subvec(6, 8), C))*r.t();
  
  //KDE for the first state
  arma::vec dens_1_1= kl_gauss_vec(series_state1.row(0).t(), NoOBs, bw);
  arma::vec dens_2_1= kl_gauss_vec(series_state1.row(1).t(), NoOBs, bw);
  arma::vec dens_3_1= kl_gauss_vec(series_state1.row(2).t(), NoOBs, bw);
  
  //KDE for the first state
  arma::vec dens_1_2= kl_gauss_vec(series_state2.row(0).t(), NoOBs, bw);
  arma::vec dens_2_2= kl_gauss_vec(series_state2.row(1).t(), NoOBs, bw);
  arma::vec dens_3_2= kl_gauss_vec(series_state2.row(2).t(), NoOBs, bw);
  
  //KDE for the third state
  arma::vec dens_1_3= kl_gauss_vec(series_state3.row(0).t(), NoOBs, bw);
  arma::vec dens_2_3= kl_gauss_vec(series_state3.row(1).t(), NoOBs, bw);
  arma::vec dens_3_3= kl_gauss_vec(series_state3.row(2).t(), NoOBs, bw);
  
  double llv = -log(arma::as_scalar(p1t*dens_1_1(0)*dens_2_1(0)*dens_3_1(0) + 
                                    p2t* dens_1_2(0)*dens_2_2(0)*dens_3_2(0) + 
                                    p3t* dens_1_3(0)*dens_2_3(0)*dens_3_3(0) ));
  
  for (int i = 1; i < NoOBs; i++) {
    
    
    double llv_temp1 = arma::as_scalar((p1t*P(0,0)+p2t*P(1,0)+p3t*P(2,0))*dens_1_1(i)*dens_2_1(i)*dens_3_1(i));
    double llv_temp2 = arma::as_scalar((p1t*P(0,1)+p2t*P(1,1)+p3t*P(2,1))*dens_1_2(i)*dens_2_2(i)*dens_3_2(i));
    double llv_temp3 = arma::as_scalar((p1t*P(0,2)+p2t*P(1,2)+p3t*P(2,2))*dens_1_3(i)*dens_2_3(i)*dens_3_3(i));
    
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2+llv_temp3);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = llv_temp3/llv_temp;
    
    llv-=log(llv_temp);
    
  }
  return llv;
}

// filter
// [[Rcpp::export]]
arma::mat filter_MS_ICA( arma::vec& theta,  arma::mat& r,  arma::mat& C, arma::vec & init) 
{
  int NoOBs = r.n_rows;
  
  //Switching Probabilities
  double p11 = theta(6);
  double p12 = 1-p11;
  double p22 = theta(7);
  double p21 = 1-p22;
  
  //initial state
  double p1t = init(0);
  double p2t = init(1);
  
  arma::mat Out(NoOBs, 2);
  Out.row(0) = {p1t, p2t};
  
  double bw = pow(4.0/(3.0*NoOBs), 0.2);
  
  
  //transform series 
  arma::mat series_state1 = arma::inv(getB(theta.subvec(0, 2), C))*r.t();
  arma::mat series_state2 = arma::inv(getB(theta.subvec(3, 5), C))*r.t();
  arma::mat series_state3 = arma::inv(getB(theta.subvec(6, 8), C))*r.t();
  
  //KDE for the first state
  arma::vec dens_1_1= kl_gauss_vec(series_state1.row(0).t(), NoOBs, bw);
  arma::vec dens_2_1= kl_gauss_vec(series_state1.row(1).t(), NoOBs, bw);
  arma::vec dens_3_1= kl_gauss_vec(series_state1.row(2).t(), NoOBs, bw);
  
  //KDE for the first state
  arma::vec dens_1_2= kl_gauss_vec(series_state2.row(0).t(), NoOBs, bw);
  arma::vec dens_2_2= kl_gauss_vec(series_state2.row(1).t(), NoOBs, bw);
  arma::vec dens_3_2= kl_gauss_vec(series_state2.row(2).t(), NoOBs, bw);
  
  //KDE for the third state
  arma::vec dens_1_3= kl_gauss_vec(series_state3.row(0).t(), NoOBs, bw);
  arma::vec dens_2_3= kl_gauss_vec(series_state3.row(1).t(), NoOBs, bw);
  arma::vec dens_3_3= kl_gauss_vec(series_state3.row(2).t(), NoOBs, bw);
  
  //double llv = -log(arma::as_scalar(p1t*dens_1_1(0)*dens_2_1(0)*dens_3_1(0) + p2t* dens_1_2(0)*dens_2_2(0)*dens_3_2(0) ));
  
  for (int i = 1; i < NoOBs; i++) {
    
    
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21)*dens_1_1(i)*dens_2_1(i)*dens_3_1(i));
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22)*dens_1_2(i)*dens_2_2(i)*dens_3_2(i));
    
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    
    Out.row(i) = {p1t, p2t};
    //llv-=log(llv_temp);
    
  }
  return Out;
}

// [[Rcpp::export]]
arma::mat filter_MS_ICA_M3(arma::vec& theta,  arma::mat& r,  arma::mat& C, arma::vec & init) 
{
  //int K = r.n_cols;
  int NoOBs = r.n_rows;
  
  //Switching Probabilities
  arma::vec lamb = theta.subvec(9, 14);
  arma::mat P    = getMarkov(lamb);
  
  //initial state
  double p1t = init(0);
  double p2t = init(1);
  double p3t = init(2);
  
  double bw = pow(4.0/(3.0*NoOBs), 0.2);
  
  arma::mat Out(NoOBs, 3);
  Out.row(0) = {p1t, p2t, p3t};
  
  //transform series 
  arma::mat series_state1 = arma::inv(getB(theta.subvec(0, 2), C))*r.t();
  arma::mat series_state2 = arma::inv(getB(theta.subvec(3, 5), C))*r.t();
  arma::mat series_state3 = arma::inv(getB(theta.subvec(6, 8), C))*r.t();
  
  //KDE for the first state
  arma::vec dens_1_1= kl_gauss_vec(series_state1.row(0).t(), NoOBs, bw);
  arma::vec dens_2_1= kl_gauss_vec(series_state1.row(1).t(), NoOBs, bw);
  arma::vec dens_3_1= kl_gauss_vec(series_state1.row(2).t(), NoOBs, bw);
  
  //KDE for the first state
  arma::vec dens_1_2= kl_gauss_vec(series_state2.row(0).t(), NoOBs, bw);
  arma::vec dens_2_2= kl_gauss_vec(series_state2.row(1).t(), NoOBs, bw);
  arma::vec dens_3_2= kl_gauss_vec(series_state2.row(2).t(), NoOBs, bw);
  
  //KDE for the third state
  arma::vec dens_1_3= kl_gauss_vec(series_state3.row(0).t(), NoOBs, bw);
  arma::vec dens_2_3= kl_gauss_vec(series_state3.row(1).t(), NoOBs, bw);
  arma::vec dens_3_3= kl_gauss_vec(series_state3.row(2).t(), NoOBs, bw);
  
  //double llv = -log(arma::as_scalar(p1t*dens_1_1(0)*dens_2_1(0)*dens_3_1(0) + p2t* dens_1_2(0)*dens_2_2(0)*dens_3_2(0) ));
  
  for (int i = 1; i < NoOBs; i++) {
    
    
    double llv_temp1 = arma::as_scalar((p1t*P(0,0)+p2t*P(1,0)+p3t*P(2,0))*dens_1_1(i)*dens_2_1(i)*dens_3_1(i));
    double llv_temp2 = arma::as_scalar((p1t*P(0,1)+p2t*P(1,1)+p3t*P(2,1))*dens_1_2(i)*dens_2_2(i)*dens_3_2(i));
    double llv_temp3 = arma::as_scalar((p1t*P(0,2)+p2t*P(1,2)+p3t*P(2,2))*dens_1_3(i)*dens_2_3(i)*dens_3_3(i));
    
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2+llv_temp3);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    p3t = llv_temp3/llv_temp;

    Out.row(i) = {p1t, p2t, p3t};
    
  }
  return Out;
}

// filter
// [[Rcpp::export]]
arma::mat filter_MS_ngml( arma::vec& theta,  arma::mat & u_st, arma::vec & init) 
{
  int K     = u_st.n_cols;
  int Tob   = u_st.n_rows;
  
  // Switching Probabilities
  double p11 = theta(theta.n_elem - 2);
  double p12 = 1-p11;
  double p22 = theta(theta.n_elem - 1);
  double p21 = 1-p22;
  
  // structural parameters
  arma::vec rot = theta.subvec(0, K*(K-1)-1);
  arma::vec sig = theta.subvec(K*(K-1), K*(K-1) + 2*K-1);
  arma::vec lambda = theta.subvec(K*(K-1) + 2*K, theta.n_elem - 3);
  

  //initial state
  double p1t = init(0);
  double p2t = init(1);
  
  // rotation matrices
  arma::mat Q_1 = givensQ_fast(rot.subvec(0, 2), K);
  arma::mat Q_2 = givensQ_fast(rot.subvec(3, 5), K);
  
  // structural covariances
  arma::mat S_1 = arma::diagmat(1/sig.subvec(0, 2));
  arma::mat S_2 = arma::diagmat(1/sig.subvec(3, 5));
  
  // dof
  arma::vec dof_1 = lambda.subvec(0,2);
  arma::vec dof_2 = lambda.subvec(3,5);
  
  // transform series (B = CQ)
  arma::mat E_1 = arma::inv_sympd(S_1) * Q_1.t() * u_st.t();
  arma::mat E_2 = arma::inv_sympd(S_2) * Q_2.t() * u_st.t();
  
  
  // Joint density
  arma::vec dens_M1 = tdens_joint(E_1, K, Tob, S_1, Q_1, dof_1);
  arma::vec dens_M2 = tdens_joint(E_2, K, Tob, S_2, Q_2, dof_2);
  
  arma::mat Out(Tob, 2);
  Out.row(0) = {p1t, p2t};
  
  for (int i = 1; i < Tob; i++) {
    
    
    double llv_temp1 = arma::as_scalar((p1t*p11+p2t*p21)*exp(dens_M1(i)));
    double llv_temp2 = arma::as_scalar((p1t*p12+p2t*p22)*exp(dens_M2(i)));
    
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;
    
    Out.row(i) = {p1t, p2t};
    
  }
  return Out;
}
