# Elbo calculation with lamba2j and prior on lambda2j

E_log_like <- function(Y, K, p, W, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz){
  n <- length(Y)
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
  pz_long <- rep(pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)%*%(diag(1, K*p) - diag(pz_long))

  res <- -(n/2)*log(2*pi) -(n/2)*(log(delta2_q) - digamma(delta1_q)) - 0.5*E_inv_sigma2*(t(Y)%*%Y - 2*t(Y)%*%W%*%diag(pz_long)%*%mu_b_q + sum(diag(((Sigma_b_q + mu_b_q%*%t(mu_b_q))%*%((t(W)%*%W)*Omega)))))
  #cat("E_log_like_i:", res, "\n")
  
  #cat("E_log_like:", res, "\n")
  return(res)
}


diff_z <- function(p, pz, a_q, b_q){
  
  # pz_logpz <- rep(NA, p)
  # pz_c_logpz_c <- rep(NA, p)
  # for(j in 1:p){if(is.na((1-pz[j])*log(1-pz[j]))){pz_c_logpz_c[j] = 0}else{pz_c_logpz_c[j] = (1-pz[j])*log(1-pz[j])}} 
  # for(j in 1:p){if(is.na(pz[j]*log(pz[j]))){pz_logpz[j] = 0}else{pz_logpz[j] = pz[j]*log(pz[j])}}
  
  pz_logpz <- rep(NA, p)
  pz_c_logpz_c <- rep(NA, p)
  for(j in 1:p){if((1-pz[j]) == 0){pz_c_logpz_c[j] = 0}else{pz_c_logpz_c[j] = (1-pz[j])*log(1-pz[j])}}  
  for(j in 1:p){if(pz[j] == 0){pz_logpz[j] = 0}else{pz_logpz[j] = pz[j]*log(pz[j])}}
  
  
  
  res <- sum(sapply(1:p, function(j){pz[j]*(digamma(a_q[j]) - digamma(a_q[j] + b_q[j])) + (1-pz[j])*(digamma(b_q[j]) - digamma(a_q[j] + b_q[j])) - pz_logpz[j] - pz_c_logpz_c[j]}))
  
  #res <- sum(sapply(1:p, function(j){pz[j]*(digamma(a[j]) - digamma(b[j]) - log(pz[j]) + log(1-pz[j])) + digamma(b[j]) - digamma(a[j] + b[j]) - log(1-pz[j])}))
  
  #cat("E_log_z:", res, "\n")
  return(res)
}



diff_b <- function(K, p, delta1_q, delta2_q, chi_q, psi_q, mu_b_q, Sigma_b_q){
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
  E_eta <- rep(NA, K*p)
  E_log_tau2 <- rep(NA, K*p)
  for(kj in 1:(K*p)){
    E_eta[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], func = "1/x")
    E_log_tau2[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], fun = "logx")
  }
  
  res <- -(K*p/2)*(log(delta2_q) - digamma(delta1_q)) -0.5*sum(E_log_tau2) -0.5*E_inv_sigma2*sum((diag(Sigma_b_q) + mu_b_q^2)*E_eta) + 0.5*log(det(Sigma_b_q)) + 0.5*K*p
  
  #cat("E_log_b:", res, "\n")
  return(res) 
}

diff_sigma2 <- function(delta1_q, delta2_q, delta2_0, delta1_0){
  
  E_inv_sigma2 <- delta1_q/delta2_q

  res <-  delta1_0*log(delta2_0) - delta1_q*log(delta2_q) - log(gamma(delta1_0)) + log(gamma(delta1_q)) + (delta1_q - delta1_0)*(log(delta2_q) - digamma(delta1_q)) + (delta2_q - delta2_0)*E_inv_sigma2
  
  #res <-  delta2_0 - log(delta2_q) - log(gamma(delta1_0)) + (delta1_q - delta1_0)*(log(delta2_q) - digamma(delta1_q)) + (delta2_q - delta2_0)*E_inv_sigma2
  
  #cat("E_log_s2", res, "\n")
  return(res)
}

diff_tau2 <- function(K, p, shape_lambda_q, rate_q, chi_q, psi_q){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  E_eta <- rep(NA, K*p)
  E_tau2 <- rep(NA, K*p)
  E_lambda2 <- rep(NA, p)
  E_log_lambda2 <- rep(NA, p)
  E_log_tau2 <- rep(NA, K*p)
  
  for(kj in 1:(K*p)){
    E_eta[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], func = "1/x")
    E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], fun = "x")
    E_log_tau2[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], fun = "logx")
  }  
  
  for(j in 1:p){
    E_lambda2[j] <- shape_lambda_q[j]/rate_q[j]
    E_log_lambda2[j] <- digamma(shape_lambda_q[j]) - log(rate_q[j])
  }

  res <-  K*sum(E_log_lambda2) -0.5*sum(sapply(1:p, function(j){E_lambda2[j]*sum(E_tau2[ids[[j]]])})) - 
    sum(sapply(1:(K*p), function(kj){0.25*(log(psi_q[kj]) - log(chi_q[kj]))-log(besselK(sqrt(chi_q[kj]*psi_q[kj]), 0.5))})) + 0.5*sum(E_log_tau2) + 0.5*sum(sapply(1:(K*p), function(kj){psi_q[kj]*E_tau2[kj] + chi_q[kj]*E_eta[kj]}))

  #cat("E_log_tau2:", res, "\n")
  return(res)
}

diff_lambda2 <- function(p, shape_lambda_0, shape_lambda_q, rate_0, rate_q){
  E_lambda2 <- rep(NA, p)
  E_log_lambda2 <- rep(NA, p)
  
  for(j in 1:p){
    E_lambda2[j] <- shape_lambda_q[j]/rate_q[j]
    E_log_lambda2[j] <- digamma(shape_lambda_q[j]) - log(rate_q[j])
  }
  
  res <- sum(sapply(1:p, function(j){shape_lambda_0*log(rate_0) - log(gamma(shape_lambda_0)) + 
      (shape_lambda_0 - 1)*E_log_lambda2[j] -rate_0*E_lambda2[j] - shape_lambda_q[j]*log(rate_q[j]) + 
      log(gamma(shape_lambda_q[j])) - (shape_lambda_q[j] - 1)*E_log_lambda2[j] + rate_q[j]*E_lambda2[j]}))
  
  #cat("E_log_lambda2:", res, "\n")
  return(res)
}


diff_theta <- function(a_q, b_q, a0, b0){
  
  res <- sum(sapply(1:p, function(j){(a0 - a_q[j])*(digamma(a_q[j]) - digamma(a_q[j] + b_q[j])) + (b0 - b_q[j])*(digamma(b_q[j]) - digamma(a_q[j] + b_q[j])) + log(gamma(a_q[j])) + log(gamma(b_q[j])) - log(gamma(a_q[j] + b_q[j]))}))
  
  #cat("E_log_theta:", res, "\n")
  return(res)
}


elbo <- function(Y, K, p, W, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz, a_q, b_q, chi_q, psi_q, delta2_0, delta1_0, shape_lambda_0, shape_lambda_q, rate_0, rate_q, a0, b0){
  
  res <- E_log_like(Y, K, p, W, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz) + diff_z(p, pz, a_q, b_q) + diff_theta(a_q, b_q, a0, b0) + diff_sigma2(delta1_q, delta2_q, delta2_0, delta1_0) + diff_b(K, p, delta1_q, delta2_q, chi_q, psi_q, mu_b_q, Sigma_b_q) + diff_tau2(K, p, shape_lambda_q, rate_q, chi_q, psi_q) + diff_lambda2(p, shape_lambda_0, shape_lambda_q, rate_0, rate_q)
  
  return(res)
}


