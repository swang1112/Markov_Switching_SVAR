
//ich passe das hier noch an
// [[Rcpp::export]]
double loglike_Normal_Gumbel_GumbelSurvival(const arma::vec& theta, const arma::mat& r) {
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
  double par_gumbel =theta[6];
  double par_gumbelSurvival =theta[7];
  //check for identification (we might not need this)
  if(p11<=0 || p11>=1 || p_22<=0 || p22>= 1){
    return -1e25;
  }
  
 
  
  //initial state
  double p1t = 1/2;
  double p2t = 1/2;
  
  //KDE for the first state
  
  
 
  
  
  
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
        
    double llv_temp = arma::as_scalar(llv_temp1+llv_temp2+llv_temp3);
    p1t = llv_temp1/llv_temp;
    p2t = llv_temp2/llv_temp;

    llv+=log(llv_temp);
    
    
  }
  return llv;
}

