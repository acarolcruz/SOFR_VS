E_quad_y_mixed <- function(p, q, K, Y, Xs_mat, W_mat, pz_q, pu_q, 
                           mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q){
  pz_long <- rep(pz_q, each = K)
  Omega_z <- pz_long%*%t(pz_long) + diag(pz_long)%*%(diag(1, K*p) - diag(pz_long))
  Omega_u <- pu_q%*%t(pu_q) + diag(pu_q)%*%(diag(1, q) - diag(pu_q))
  
  res <- (t(Y)%*%Y - 2*t(Y)%*%Xs_mat%*%diag(pu_q)%*%mu_alpha_q -2*t(Y)%*%W_mat%*%diag(pz_long)%*%mu_b_q + 2*sum(diag(((t(W_mat)%*%Xs_mat)*(pz_long%*%t(pu_q)))%*%(mu_alpha_q%*%t(mu_b_q)))) + sum(diag((Sigma_b_q + mu_b_q%*%t(mu_b_q))%*%((t(W_mat)%*%W_mat)*Omega_z))) + sum(diag((Sigma_alpha_q + mu_alpha_q%*%t(mu_alpha_q))%*%((t(Xs_mat)%*%Xs_mat)*Omega_u))))
  
  #2*sum(diag(((t(Xs_mat)%*%W_mat)*(pu_q%*%t(pz_long)))%*%(mu_b_q%*%t(mu_alpha_q))))
  return(res)
}

#E_quad_y_mixed(p, q, K, Y, Xs_mat, W_mat, pz_q, pu_q, mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q)
