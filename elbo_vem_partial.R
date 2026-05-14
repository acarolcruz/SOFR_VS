# Elbo calculation with lamba2j and lambda2l fixed
E_log_like <- function(Y, K, p, q, W, Xs_mat, delta1_q, delta2_q, 
                       Sigma_b, mu_b, pz, 
                       pu, Sigma_alpha, mu_alpha){
  n <- length(Y)
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
  pz_long <- rep(pz, each = K)
  Omega_z <- pz_long%*%t(pz_long) + diag(pz_long)%*%(diag(1, K*p) - diag(pz_long))
  
  Omega_u <- pu%*%t(pu) + diag(pu)%*%(diag(1, q) - diag(pu))
  
  
  # approx  digamma(delta1_q)
  res <- -(n/2)*log(2*pi) -(n/2)*(log(delta2_q) -(log(delta1_q) - 1/(2*delta1_q))) - 0.5*E_inv_sigma2*(t(Y)%*%Y -2*t(Y)%*%Xs_mat%*%diag(pu)%*%mu_alpha -2*t(Y)%*%W%*%diag(pz_long)%*%mu_b + sum(diag(((Sigma_alpha + mu_alpha%*%t(mu_alpha))%*%((t(Xs_mat)%*%Xs_mat)*Omega_u)))) + sum(diag(((Sigma_b + mu_b%*%t(mu_b))%*%((t(W)%*%W)*Omega_z)))) + 2*sum(diag(((t(W)%*%Xs_mat)*(pz_long%*%t(pu)))%*%(mu_alpha%*%t(mu_b)))))
  
  #cat("E_log_like:", res, "\n")
  return(res)
}

diff_u <- function(q, pu, a_u, b_u){
  
  pu_logpu <- rep(NA, q)
  pu_c_logpu_c <- rep(NA, q)
  for(j in 1:q){if((1-pu[j]) == 0){pu_c_logpu_c[j] = 0}else{pu_c_logpu_c[j] = (1-pu[j])*log(1-pu[j])}}  
  for(j in 1:q){if(pu[j] == 0){pu_logpu[j] = 0}else{pu_logpu[j] = pu[j]*log(pu[j])}}
  
  res <- sum(sapply(1:q, function(j){pu[j]*(digamma(a_u[j]) - digamma(a_u[j] + b_u[j])) + (1-pu[j])*(digamma(b_u[j]) - digamma(a_u[j] + b_u[j])) - pu_logpu[j] - pu_c_logpu_c[j]}))
  
  #cat("E_log_z:", res, "\n")
  return(res)
}


diff_z <- function(p, pz, a_z, b_z){
  
  pz_logpz <- rep(NA, p)
  pz_c_logpz_c <- rep(NA, p)
  for(j in 1:p){if((1-pz[j]) == 0){pz_c_logpz_c[j] = 0}else{pz_c_logpz_c[j] = (1-pz[j])*log(1-pz[j])}}  
  for(j in 1:p){if(pz[j] == 0){pz_logpz[j] = 0}else{pz_logpz[j] = pz[j]*log(pz[j])}}
  
  res <- sum(sapply(1:p, function(j){pz[j]*(digamma(a_z[j]) - digamma(a_z[j] + b_z[j])) + (1-pz[j])*(digamma(b_z[j]) - digamma(a_z[j] + b_z[j])) - pz_logpz[j] - pz_c_logpz_c[j]}))
  
  #cat("E_log_z:", res, "\n")
  return(res)
}

diff_b <- function(K, p, delta1_q, delta2_q, chi_b_q, psi_b_q, mu_b, Sigma_b){
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
  E_eta <- rep(NA, K*p)
  E_log_tau2 <- rep(NA, K*p)
  for(kj in 1:(K*p)){
    E_eta[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], func = "1/x")
    E_log_tau2[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], fun = "logx")
  }
  
  # using approx for digamma(delta1_q)
  res <- -(K*p/2)*(log(delta2_q) - (log(delta1_q) - 1/(2*delta1_q))) -0.5*sum(E_log_tau2) -0.5*E_inv_sigma2*sum((diag(Sigma_b) + mu_b^2)*E_eta) + 0.5*log(det(Sigma_b)) + 0.5*K*p
  
  #cat("E_log_b:", res, "\n")
  return(res) 
}

diff_alpha <- function(q, delta1_q, delta2_q, chi_alpha_q, psi_alpha_q, mu_alpha, Sigma_alpha){
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
  E_nu_inv <- rep(NA, q)
  E_log_nu2 <- rep(NA,q)
  for(l in 1:q){
    E_nu_inv[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], func = "1/x")
    E_log_nu2[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], fun = "logx")
  }
  
  # using approx for digamma(delta1_q)
  res <- -(q/2)*(log(delta2_q) - (log(delta1_q) - 1/(2*delta1_q))) -0.5*sum(E_log_nu2) -0.5*E_inv_sigma2*sum((diag(Sigma_alpha) + mu_alpha^2)*E_nu_inv) + 0.5*log(det(Sigma_alpha)) + 0.5*q
  
  #cat("E_log_q:", res, "\n")
  return(res) 
}

diff_sigma2 <- function(delta1_q, delta2_q, delta2_0, delta1_0){
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
  # original (all terms)
  #res <-  delta1_0*log(delta2_0) - delta1_q*log(delta2_q) - log(gamma(delta1_0)) + log(gamma(delta1_q)) + (delta1_q - delta1_0)*(log(delta2_q) - digamma(delta1_q)) + (delta2_q - delta2_0)*E_inv_sigma2
  
  # using approx as in xian
  res <- delta1_0*log(delta2_0) - delta1_q*log(delta2_q) - log(gamma(delta1_0)) + delta1_q*log(delta1_q) - delta1_q -0.5*log(delta1_q) + (delta1_q - delta1_0)*(log(delta2_q) - (log(delta1_q) - 1/(2*delta1_q))) + (delta2_q - delta2_0)*E_inv_sigma2
  
  #cat("E_log_s2", res, "\n")
  return(res)
}


diff_tau2 <- function(K, p, lambda2_b, chi_b_q, psi_b_q){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  E_eta <- rep(NA, K*p)
  E_tau2 <- rep(NA, K*p)
  E_log_tau2 <- rep(NA, K*p)
  
  for(kj in 1:(K*p)){
    E_eta[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], func = "1/x")
    E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], fun = "x")
    E_log_tau2[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], fun = "logx")
  }  
  
  res <-  K*sum(log(lambda2_b)) -0.5*sum(sapply(1:p, function(j){lambda2_b[j]*sum(E_tau2[ids[[j]]])})) - 
    sum(sapply(1:(K*p), function(kj){0.25*(log(psi_b_q[kj]) - log(chi_b_q[kj]))-log(besselK(sqrt(chi_b_q[kj]*psi_b_q[kj]), 0.5))})) + 0.5*sum(E_log_tau2) + 0.5*sum(sapply(1:(K*p), function(kj){psi_b_q[kj]*E_tau2[kj] + chi_b_q[kj]*E_eta[kj]}))
  
  #cat("E_log_tau2:", res, "\n")
  return(res)
}

diff_nu2 <- function(q, lambda2_alpha, chi_alpha_q, psi_alpha_q){
  
  E_nu_inv <- rep(NA, q)
  E_nu2 <- rep(NA, q)
  E_log_nu2 <- rep(NA, q)
  
  for(l in 1:q){
    E_nu_inv[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], func = "1/x")
    E_nu2[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], fun = "x")
    E_log_nu2[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], fun = "logx")
  }
  
  res <-  sum(log(lambda2_alpha)) -0.5*sum(lambda2_alpha*E_nu2) -0.25*sum((log(psi_alpha_q) - log(chi_alpha_q))) + sum(sapply(1:q, function(l){log(besselK(sqrt(chi_alpha_q[l]*psi_alpha_q[l]), 0.5))})) +0.5*sum(E_log_nu2) + 0.5*sum(psi_alpha_q*E_nu2 + chi_alpha_q*E_nu_inv)
  
  #cat("E_log_nu2:", res, "\n")
  return(res)
}

diff_theta_z <- function(p, a_z, b_z, a_z_0, b_z_0){
  
  res <- sum(sapply(1:p, function(j){(a_z_0 - a_z[j])*(digamma(a_z[j]) - digamma(a_z[j] + b_z[j])) + (b_z_0 - b_z[j])*(digamma(b_z[j]) - digamma(a_z[j] + b_z[j])) + log(gamma(a_z[j])) + log(gamma(b_z[j])) - log(gamma(a_z[j] + b_z[j]))}))
  
  #cat("E_log_theta_z:", res, "\n")
  return(res)
}

diff_theta_u <- function(q, a_u, b_u, a_u_0, b_u_0){
  
  #with no constants
  res <- sum(sapply(1:q, function(l){(a_u_0 - 1)*(digamma(a_u[l]) - digamma(a_u[l] + b_u[l])) + (b_u_0 - 1)*(digamma(b_u[l]) - digamma(a_u[l] + b_u[l])) - (a_u[l] - 1)*(digamma(a_u[l]) - digamma(a_u[l] + b_u[l])) -
      (b_u[l] - 1)*(digamma(b_u[l]) - digamma(a_u[l] + b_u[l])) - log(gamma(a_u[l]) + b_u[l]) + log(gamma(a_u[l])) + log(gamma(b_u[l]))}))
  
  #cat("E_log_theta_u:", res, "\n")
  return(res)
}




elbo_lambda2 <- function(lambda2, Y, K, p_var, q, W, Xs_mat, delta1_q, delta2_q, Sigma_b, mu_b, pz, a_z, b_z, Sigma_alpha, mu_alpha, pu, a_u, b_u, chi_b_q, psi_b_q, chi_alpha_q, psi_alpha_q, delta2_0, delta1_0, a_z_0, b_z_0, a_u_0, b_u_0){
  
  lambda2_alpha <- lambda2[c((p_var+1):(p_var+q))]
  lambda2_b <- lambda2[1:p_var]
  
  res <- E_log_like(Y, K, p_var, q, W, Xs_mat, delta1_q, delta2_q, Sigma_b, mu_b, pz, pu, Sigma_alpha, mu_alpha) + diff_z(p, pz, a_z, b_z) + diff_theta_z(p, a_z, b_z, a_z_0, b_z_0) + diff_u(q, pu, a_u, b_u) + diff_theta_u(q, a_u, b_u, a_u_0, b_u_0) + diff_sigma2(delta1_q, delta2_q, delta2_0, delta1_0) + diff_b(K, p_var, delta1_q, delta2_q, chi_b_q, psi_b_q, mu_b, Sigma_b) + diff_alpha(q, delta1_q, delta2_q, chi_alpha_q, psi_alpha_q, mu_alpha, Sigma_alpha) + diff_tau2(K, p_var, lambda2_b, chi_b_q, psi_b_q) + diff_nu2(q, lambda2_alpha, chi_alpha_q, psi_alpha_q)
  return(res)
}

elbo_partial <- function(Y, K, p, q, W, Xs_mat, delta1_q, delta2_q, Sigma_b, mu_b, pz, a_z, b_z, Sigma_alpha, mu_alpha, pu, a_u, b_u, chi_b_q, psi_b_q, chi_alpha_q, psi_alpha_q, delta2_0, delta1_0, a_z_0, b_z_0, a_u_0, b_u_0, lambda2){
  
  lambda2_alpha <- lambda2[c((p+1):(p+q))]
  lambda2_b <- lambda2[1:p]
  
  res <- E_log_like(Y, K, p, q, W, Xs_mat, delta1_q, delta2_q, Sigma_b, mu_b, pz, pu, Sigma_alpha, mu_alpha) + diff_z(p, pz, a_z, b_z) + diff_theta_z(p, a_z, b_z, a_z_0, b_z_0) + diff_u(q, pu, a_u, b_u) + diff_theta_u(q, a_u, b_u, a_u_0, b_u_0) + diff_sigma2(delta1_q, delta2_q, delta2_0, delta1_0) + diff_b(K, p, delta1_q, delta2_q, chi_b_q, psi_b_q, mu_b, Sigma_b) + diff_alpha(q, delta1_q, delta2_q, chi_alpha_q, psi_alpha_q, mu_alpha, Sigma_alpha) + diff_tau2(K, p, lambda2_b, chi_b_q, psi_b_q) + diff_nu2(q, lambda2_alpha, chi_alpha_q, psi_alpha_q)
  
  return(res)
}


# needs to have the same arguments as fr
dev_elbo <- function(lambda2, Y, K, p_var, q, W, Xs_mat, delta1_q, delta2_q, Sigma_b, mu_b, pz, a_z, b_z, Sigma_alpha, mu_alpha, pu, a_u, b_u, chi_b_q, psi_b_q, chi_alpha_q, psi_alpha_q, delta2_0, delta1_0, a_z_0, b_z_0, a_u_0, b_u_0){
  
  lambda2_alpha <- lambda2[c((p+1):(p+q))]
  lambda2_b <- lambda2[1:p]
  
  ids <- split(1:(K*p_var), rep(1:p_var, each = K))
  
  E_tau2 <- rep(NA, K*p_var)
  E_nu2 <- rep(NA, q)
  
  for(l in 1:q){
    E_nu2[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], fun = "x")
  }
  
  for(kj in 1:(K*p_var)){
    E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], fun = "x")
  }
  
  res_b <- sapply(1:p_var, function(j){K/lambda2_b[j] -0.5*sum(E_tau2[ids[[j]]])})
  res_alpha <- sapply(1:q, function(l){1/lambda2_alpha[l] - 0.5*E_nu2[l]})
  return(c(res_b, res_alpha))
}

lambda2_hat <- function(K, p_var, q, chi_b_q, psi_b_q, chi_alpha_q, psi_alpha_q){
  
  ids <- split(1:(K*p_var), rep(1:p_var, each = K))
  
  E_tau2 <- rep(NA, K*p_var)
  E_nu2 <- rep(NA, q)
  
  for(l in 1:q){
    E_nu2[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], fun = "x")
  }
  
  for(kj in 1:(K*p_var)){
    E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], fun = "x")
  }
  
  lambda2_b_hat <- sapply(1:p_var, function(j){(2*K)/sum(E_tau2[ids[[j]]])})
  lambda2_alpha_hat <- sapply(1:q, function(l){2/E_nu2[l]})
  return(c(lambda2_b_hat, lambda2_alpha_hat))
}

check_convergence <- function(elbo_c, elbo_prev, convergence_threshold) {
  if(is.null(elbo_prev) == TRUE) {
    return(FALSE)
  }
  else{
    dif <- elbo_c - elbo_prev
    if(abs(dif)  <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
}

