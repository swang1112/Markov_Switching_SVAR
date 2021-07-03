
//ich passe das hier noch an
// [[Rcpp::export]]
double loglike_Normal_Gumbel_GumbelSurvival(const arma::vec& theta, const arma::mat& r) {
 {
  // Log-Likelihood function

  // convert to matrices
  int n = r.n_cols;
  // Length of each series
  int NoOBs = r.n_rows;
  int numb_of_vars = 2 * pow(n, 2) + n * (n + 1)/2;
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
  double p11 = theta[6];
  double p12 = 1-p11;
  double p22 = theta[7];
  double p21 = 1-p22;

  //check for identification (we might not need this)
  if(p11<=0 || p11>=1 || p22<=0 || p22>= 1){
    return -1e25;
  }



  //initial state
  double p1t = 1/2;
  double p2t = 1/2;
  
  // rotation matrices
  arma::mat rot1_1;
  arma::mat rot2_1;
  arma::mat rot3_1;
  
  arma::mat rot1_2;
  arma::mat rot2_2;
  arma::mat rot3_2;
  
  //cholesky factors about here
  arma::mat cholesky_factor_1_1 = arma::chol(rot1_1);
  arma::mat cholesky_factor_2_1 = arma::chol(rot2_1);
  arma::mat cholesky_factor_3_1 = arma::chol(rot3_1);
  
  arma::mat cholesky_factor_1_2 = arma::chol(rot1_2);
  arma::mat cholesky_factor_2_2 = arma::chol(rot2_2);
  arma::mat cholesky_factor_3_2 = arma::chol(rot3_2); 
  
  //transform series (Hier fehlt noch die B Matrix)
  arma::mat series_state1 = cholesky_factor_1_1 * cholesky_factor_2_1 * cholesky_factor_3_1 *r;
  arma::mat series_state2 = cholesky_factor_1_2 * cholesky_factor_2_2 * cholesky_factor_3_2 *r;
  
  //KDE for the first state
  arma::vec dens_1_1= kdensity(*series_state1.col(0).t());
  arma::vec dens_2_1= kdensity(*series_state1.col(1).t());
  arma::vec dens_3_1= kdensity(*series_state1.col(2).t());
  
  //KDE for the first state
  arma::vec dens_1_2= kdensity(*series_state2.col(0).t());
  arma::vec dens_2_2= kdensity(*series_state2.col(1).t());
  arma::vec dens_3_2= kdensity(*series_state2.col(2).t());
  



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
