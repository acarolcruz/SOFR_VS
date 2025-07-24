#' Variational Bayes for variable selection in partial scalar-on-function regression (SOFR)
#'
#' @param xi_1_0
#' @param xi_2_0 
#' @param delta1_0 
#' @param delta2_0 
#' @param a_z_0 
#' @param b_z_0 
#' @param a_u_0
#' @param b_u_0
#' @param shape_lambda_0 
#' @param rate_0 
#' @param initial_values 
#' @param data
#' @param data_std 
#' @param K 
#' @param q
#' @param p 
#' @param Niter 
#' @param convergence_threshold 
#' @param n 
#' @param std 
#'
#' @returns
#' @export
#'
#' @examples
#' 
VBSOFR_VS_partial <- function(delta1_0, delta2_0, a_z_0, b_z_0, a_u_0, b_u_0, 
                              xi_1_0, xi_2_0, shape_lambda_0, rate_0, 
                              initial_values, data, 
                              data_std, n, K, p, q, 
                              Niter, convergence_threshold, std){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  Xs_mat <- data_std$Xs_std
  W_mat <- data_std$W_mat
  Y_std <- data_std$Y_std
  Y <- data_std$Y
  
  # uncomment if testing pz_q with true values
  #W_mat <- data$W_mat
  #Y <- data$W_mat
  
  # Update VB parambers not updated within VB
  delta1_q <- n/2 + delta1_0 + (K*p)/2 + q/2
  shape_lambda_q <- rep(K + shape_lambda_0, p)
  xi_1_q <- xi_1_0 + 1
  
  #Initial values
  E_lambda2_f <- initial_values$E_lambda2_f
  E_lambda2_s <- initial_values$E_lambda2_s
  
  E_inv_sigma2 <- initial_values$E_inv_sigma2 #delta1_q/delta2_q#
  delta2_q <- 1/E_inv_sigma2 #(delta1_q - 1)*var(Y)
  E_eta <- initial_values$E_eta
  E_nu_inv <- initial_values$E_nu_inv
  
  pz_q <- initial_values$pz
  pu_q <- initial_values$pu
  
  # Structure
  # Structure
  psi_b_q <- rep(NA, K*p)
  psi_alpha_q <- rep(NA, q)
  rate_q <- rep(NA, p)
  xi_2_q <- c()
  E_tau2 <- rep(NA, K*p)
  E_nu2 <- rep(NA, q)
  
  iter = 1
  elbo_prev = 0
  converged <- FALSE
  start <- proc.time()
  #while(iter < Niter & converged == FALSE){
  while(iter < Niter){
    
    # Step 1: Update variational of b
    pz_long <- rep(pz_q, each = K)
    
    Omega <- pz_long%*%t(pz_long) + diag(pz_long)%*%(diag(1, K*p) - diag(pz_long))
    
    Q <- (diag(E_eta) + ((t(W_mat)%*%W_mat)*Omega))
    if(is.singular.matrix(Q)){warning(paste('Determinant of Sigma_b is close to zero for iteration', iter))}
    Sigma_b_q <- solve(as.numeric(E_inv_sigma2)*Q)
    if(std){
      mu_b_q <- solve(Q)%*%(diag(pz_long)%*%t(W_mat)%*%Y_std)
    } else{
      mu_b_q <- solve(Q)%*%(diag(pz_long)%*%t(W_mat)%*%Y)
    }  
    
    # Step 2: Update variational of alpha
    
    Omega_u <- pu_q%*%t(pu_q) + diag(pu_q)%*%(diag(1, q) - diag(pu_q))
    M <- (diag(E_nu_inv) + ((t(Xs_mat)%*%Xs_mat)*Omega_u))
    
    Sigma_alpha_q <- solve(as.numeric(E_inv_sigma2)*M)
    if(std){
      mu_alpha_q <- solve(M)%*%(diag(pu_q)%*%t(Xs_mat)%*%Y_std - diag(pu_q)%*%t(Xs_mat)%*%W_mat%*%diag(pz_long)%*%mu_b_q)
    } else{
      mu_alpha_q <- solve(M)%*%(diag(pu_q)%*%t(Xs_mat)%*%Y - diag(pu_q)%*%t(Xs_mat)%*%W_mat%*%diag(pz_long)%*%mu_b_q)
    }
    
    # Step 2: Update variational of sigma2
    if(std){
      A <- E_quad_y_mixed(p, q, K, Y_std, Xs_mat, W_mat, pz_q, pu_q, mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q) + sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2)) + sum(E_nu_inv*(diag(Sigma_alpha_q) + mu_alpha_q^2))
    } else{
      A <- E_quad_y_mixed(p, q, K, Y, Xs_mat, W_mat, pz_q, pu_q, mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q) + sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2)) + sum(E_nu_inv*(diag(Sigma_alpha_q) + mu_alpha_q^2))
    }
    
    delta2_q <- A/2 + delta2_0
    
    E_inv_sigma2 <- delta1_q/delta2_q
    
    # Step 3: Update variational of tau2
    chi_b_q <- c((diag(Sigma_b_q) + mu_b_q^2)*as.numeric(E_inv_sigma2))
    for(j in 1:p){
      psi_b_q[ids[[j]]] <- rep(E_lambda2_f[j], K)
    }
    for(kj in 1:(K*p)){
      E_eta[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], func = "1/x")
      E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], fun = "x")
    }
    #Step 3.1 update variational for lambda2
    for(j in 1:p){
      rate_q[j] <- rate_0 + 0.5*sum(E_tau2[ids[[j]]])
      E_lambda2_f[j] <- shape_lambda_q[j]/rate_q[j]
    }
    
    # Step 4: Update variational of nu2
    chi_alpha_q <- c((diag(Sigma_alpha_q) + mu_alpha_q^2)*as.numeric(E_inv_sigma2))
    psi_alpha_q <- E_lambda2_s
    
    for(l in 1:q){
      E_nu_inv[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], func = "1/x")
      E_nu2[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], fun = "x")
    }
    #Step 4.1 update variational for lambda_alpha2
    for(l in 1:q){
      xi_2_q[l] <- xi_2_0 + 0.5*E_nu2[l]
      E_lambda2_s[l] <- xi_1_q/xi_2_q[l]
    }
    
    # Step 5: Update variational of theta_u
    a_u_q <- pu_q + a_u_0
    b_u_q <- 1 - pu_q + b_u_0
    
    # Step 6: Update variational of u
    #for each l = 1,...,q
    if(std){
      for(l in 1:q){
        pu_r <- pu_q
        su_js <- c()
        for(r in 0:1){
          pu_r[l] <- r
          su_js[r+1] <- -(n/2)*(log(delta2_q) - digamma(delta1_q)) - 0.5*as.numeric(E_inv_sigma2)*E_quad_y_mixed(p, q, K, Y_std, Xs_mat, W_mat, pz_q, pu_r, mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q) + r*(digamma(a_u_q[l]) - digamma(a_u_q[l] + b_u_q[l])) + (1-r)*(digamma(b_u_q[l]) - digamma(a_u_q[l] + b_u_q[l]))
        }  
        if(sum(exp(su_js)) == 0){
          pu_q[l] <- c(0,1)[which.max(su_js)]
        } else if (sum(exp(su_js)) == Inf) {
          pu_q[l] <- c(0,1)[which.max(su_js)]
        } else {
          pu_q[l] <- exp(su_js[2])/sum(exp(su_js))
        }
      }
    }else{
      for(l in 1:q){
        pu_r <- pu_q
        su_js <- c()
        for(r in 0:1){
          pu_r[l] <- r
          su_js[r+1] <- -(n/2)*(log(delta2_q) - digamma(delta1_q)) - 0.5*as.numeric(E_inv_sigma2)*E_quad_y_mixed(p, q, K, Y, Xs_mat, W_mat, pz_q, pu_r, mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q) + r*(digamma(a_u_q[l]) - digamma(a_u_q[l] + b_u_q[l])) + (1-r)*(digamma(b_u_q[l]) - digamma(a_u_q[l] + b_u_q[l]))
        }  
        if(sum(exp(su_js)) == 0){
          pu_q[l] <- c(0,1)[which.max(su_js)]
        } else if (sum(exp(su_js)) == Inf) {
          pu_q[l] <- c(0,1)[which.max(su_js)]
        } else {
          pu_q[l] <- exp(su_js[2])/sum(exp(su_js))
        }
      }
    }
    
    # Step 7: Update variational of theta_b
    a_z_q <- pz_q + a_z_0
    b_z_q <- 1 - pz_q + b_z_0
    
    # Step 8: Update variational of Z
    #for each j = 1, ..., p
    if(std){
      for(j in 1:p){
        pz_r <- pz_q
        sz_js <- c()
        for(r in 0:1){
          pz_r[j] <- r
          sz_js[r+1] <- -(n/2)*(log(delta2_q) - digamma(delta1_q)) - 0.5*as.numeric(E_inv_sigma2)*E_quad_y_mixed(p, q, K, Y_std, Xs_mat, W_mat, pz_r, pu_q, mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q) + r*(digamma(a_z_q[j]) - digamma(a_z_q[j] + b_z_q[j])) + (1-r)*(digamma(b_z_q[j]) - digamma(a_z_q[j] + b_z_q[j]))
        }  
        if(sum(exp(sz_js)) == 0){
          pz_q[j] <- c(0,1)[which.max(sz_js)]
        } else if (sum(exp(sz_js)) == Inf) {
          pz_q[j] <- c(0,1)[which.max(sz_js)]
        } else {
          pz_q[j] <- exp(sz_js[2])/sum(exp(sz_js))
        }
      }
    }else{
      for(j in 1:p){
        
        pz_r <- pz_q
        sz_js <- rep(NA, p)
        for(r in 0:1){
          
          pz_r[j] = r
          
          sz_js[r+1] <- -(n/2)*(log(delta2_q) - digamma(delta1_q)) -0.5*as.numeric(E_inv_sigma2)*E_quad_y_mixed(p, q, K, Y, Xs_mat, W_mat, pz_r, pu_q, mu_b_q, mu_alpha_q, Sigma_b_q, Sigma_alpha_q) + r*(digamma(a_z_q[j]) - digamma(a_z_q[j] + b_z_q[j])) + (1-r)*(digamma(b_z_q[j]) - digamma(a_z_q[j] + b_z_q[j]))
        }   
        
        if(sum(exp(sz_js)) == 0){
          pz_q[j] <- c(0,1)[which.max(sz_js)]
        } else if (sum(exp(sz_js)) == Inf) {
          pz_q[j] <- c(0,1)[which.max(sz_js)]
        } else {
          pz_q[j] <- exp(sz_js[2])/sum(exp(sz_js))
        }
      }
    }
    
    
    iter = iter + 1
    # if(std){
    #   elbo_c <- elbo_mixed(Y_std, K, p, q, Xs_mat, W_mat, delta1_q, delta2_q, xi_1_q,
    #                        xi_2_q, xi_1_0, xi_2_0, Sigma_b_q, mu_b_q, pz_q, 
    #                        Sigma_alpha_q, mu_alpha_q, pu_q, chi_b_q, psi_b_q,
    #                        a_q, b_q, chi_alpha_q, psi_alpha_q, delta2_0, delta1_0, 
    #                        shape_lambda_0, shape_lambda_q, rate_0, rate_q, 
    #                        c0, d0, a0, b0)
    # } else{
    #   elbo_c <- elbo_mixed(Y, K, p, q, Xs_mat, W_mat, delta1_q, delta2_q, xi_1_q,
    #                        xi_2_q, xi_1_0, xi_2_0, Sigma_b_q, mu_b_q, pz_q, 
    #                        Sigma_alpha_q, mu_alpha_q, pu_q, chi_b_q, psi_b_q,
    #                        a_q, b_q, chi_alpha_q, psi_alpha_q, delta2_0, delta1_0, 
    #                        shape_lambda_0, shape_lambda_q, rate_0, rate_q, 
    #                        c0, d0, a0, b0)
    # }
    # 
    # converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)
    # 
    # elbo_prev <- elbo_c
    #print(elbo_c)
  }
  runtime_VB <- proc.time() - start
  
  res <- list(mu_b = mu_b_q, Sigma_b = Sigma_b_q, pz = pz_q, pu = pu_q, 
              delta1 = delta1_q, delta2 = delta2_q, 
              chi_b = chi_b_q, psi_b = psi_b_q,
              chi_alpha = chi_alpha_q, psi_alpha = psi_alpha_q,
              mu_alpha = mu_alpha_q, Sigma_alpha = Sigma_alpha_q, 
              a_z_q = a_z_q, b_z_q = b_z_q, a_u_q = a_u_q, b_u_q = b_u_q,
              E_lambda2f = E_lambda2_f, E_lambda2s = E_lambda2_s,
              N_iter = iter, runtime = runtime_VB[[3]], elbo = NULL)
  
  return(res)
}