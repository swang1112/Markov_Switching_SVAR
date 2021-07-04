
smoothed_probabilities<- function(Filtered_Probabilities,par){
  n=nrow(Filtered_Probabilities)
  p11=par[1]
  p22=par[2]
  smooth_prob=numeric(n)
  smooth_prob[n]=Filtered_Probabilities[n]
  for(i in (1:(n-1))){
    #updated filtered probability
    p_t_t1=p11*Filtered_Probabilities[n-i]+(1-p22)*(1-Filtered_Probabilities[n-i])
    #compute smoothed probability
    smooth_prob[n-i]=Filtered_Probabilities[n-i]*p11*smooth_prob[n-i+1]/p_t_t1+Filtered_Probabilities[n-i]*(1-p22)*(1-smooth_prob[n-i+1])/p_t_t1
    
  }
  
  return(smooth_prob)
}