get_Bt = function(theta, C){
  rot = theta[1:choose(K,2)]
  Q  = givensQ_fast(rot, K)
  sigma = theta[(choose(K,2)+1): (choose(K,2)+K)]
  lambda = theta[(choose(K,2)+1+K):length(theta)]
  sigma_sort = sort.int(sigma, decreasing = T, index.return = T)
  Q = Q[, sigma_sort$ix]
  S = diag(1/sigma_sort$x)
  D = diag(sigma_sort$x^2)
  list('B' = C %*% Q %*% S,
       'D' = D, 
       'dof' = lambda[sigma_sort$ix])
}

ng_ica = function(x){
  C = cov(x) %>% chol %>% t
  u_st = t(solve(C) %*% t(x))
  sig0 = sort(rchisq(K,3), decreasing = T)
  erg = nlm(p = c(runif(choose(K, 2), 0, pi/2), sig0, runif(K, 2, 20)), 
            f = nglike_ICA, u_st = u_st, neg = TRUE, 
            gradtol = 1e-10, iterlim = 10000)
  get_Bt(erg$estimate, C)
}
