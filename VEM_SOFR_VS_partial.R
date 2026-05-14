#' Variational EM for variable selection in partially scalar-on-function regression (SOFR)
#'
#' @param delta1_0 
#' @param delta2_0 
#' @param a0 
#' @param b0 
#' @param initial_values 
#' @param data 
#' @param data_std 
#' @param K 
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
VEM_SOFR_VS_partial <- function(delta1_0, delta2_0, a_z_0, b_z_0, a_u_0, b_u_0, 
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
  
  # Update VB parameters not updated within VB
  delta1_q <- n/2 + delta1_0 + (K*p)/2 + q/2
  
  #Initial values
  lambda2 <- initial_values$lambda2
  lambda2_b <- lambda2[1:p]
  lambda2_alpha <- lambda2[c((p+1):(p+q))]
  mu_alpha_q <- initial_values$E_alpha
  
  E_inv_sigma2 <- initial_values$E_inv_sigma2 #delta1_q/delta2_q#
  delta2_q <- 1/E_inv_sigma2 #(delta1_q - 1)*var(Y)
  E_eta <- initial_values$E_eta
  E_nu_inv <- initial_values$E_nu_inv
  
  pz_q <- initial_values$pz
  pu_q <- initial_values$pu
  
  # Structure
  psi_b_q <- rep(NA, K*p)
  psi_alpha_q <- rep(NA, q)
  E_tau2 <- rep(NA, K*p)
  E_nu2 <- rep(NA, q)
  
  iter = 1
  elbo_prev = 0
  converged <- FALSE
  start <- proc.time()
  
  while(iter < Niter & converged == FALSE){
    
    # Step 1: Update variational of b
    pz_long <- rep(pz_q, each = K)
    
    Omega <- pz_long%*%t(pz_long) + diag(pz_long)%*%(diag(1, K*p) - diag(pz_long))
    
    Q <- (diag(E_eta) + ((t(W_mat)%*%W_mat)*Omega))
    if(is.singular.matrix(Q)){warning(paste('Determinant of Sigma_b is close to zero for iteration', iter))}
    Sigma_b_q <- solve(as.numeric(E_inv_sigma2)*Q)
    if(std){
      mu_b_q <- solve(Q)%*%((diag(pz_long)%*%t(W_mat))%*%(Y_std - Xs_mat%*%diag(pu_q)%*%mu_alpha_q))
    } else{
      mu_b_q <- solve(Q)%*%((diag(pz_long)%*%t(W_mat))%*%(Y - Xs_mat%*%diag(pu_q)%*%mu_alpha_q))
    }  
    
    # Step 2: Update variational of alpha
    Omega_u <- pu_q%*%t(pu_q) + diag(pu_q)%*%(diag(1, q) - diag(pu_q))
    M <- (diag(E_nu_inv) + ((t(Xs_mat)%*%Xs_mat)*Omega_u))
    
    Sigma_alpha_q <- solve(as.numeric(E_inv_sigma2)*M)
    if(std){
      mu_alpha_q <- solve(M)%*%((diag(pu_q)%*%t(Xs_mat))%*%(Y_std - W_mat%*%diag(pz_long)%*%mu_b_q))
    } else{
      mu_alpha_q <- solve(M)%*%((diag(pu_q)%*%t(Xs_mat))%*%(Y - W_mat%*%diag(pz_long)%*%mu_b_q))
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
      psi_b_q[ids[[j]]] <- rep(lambda2_b[j], K)
    }
    for(kj in 1:(K*p)){
      E_eta[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], func = "1/x")
      E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_b_q[kj], psi = psi_b_q[kj], fun = "x")
    }
    
    
    # Step 4: Update variational of nu2
    chi_alpha_q <- c((diag(Sigma_alpha_q) + mu_alpha_q^2)*as.numeric(E_inv_sigma2))
    psi_alpha_q <- lambda2_alpha
    
    for(l in 1:q){
      E_nu_inv[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], func = "1/x")
      E_nu2[l] <- Egig(lambda = 0.5, chi = chi_alpha_q[l], psi = psi_alpha_q[l], fun = "x")
    }

    # Step 5: Update variational of theta_b
    a_z_q <- pz_q + a_z_0
    b_z_q <- 1 - pz_q + b_z_0
    
    # Step 6: Update variational of Z
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
    
    # Step 7: Update variational of theta_u
    a_u_q <- pu_q + a_u_0
    b_u_q <- 1 - pu_q + b_u_0
    
    # Step 8: Update variational of u
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
    
    
    # if(std){
    #   lambda2_c <- tryCatch(optim(lambda2, fn = elbo_lambda2, gr = dev_elbo, Y = Y_std, K = K, p_var = p, q = q, W = W_mat, Xs_mat = Xs_mat, delta1_q = delta1_q, delta2_q = delta2_q, Sigma_b = Sigma_b_q, mu_b =  mu_b_q, pz = pz_q, a_z = a_z_q, b_z = b_z_q, Sigma_alpha = Sigma_alpha_q, mu_alpha = mu_alpha_q, pu = pu_q, a_u = a_u_q, b_u = b_u_q, chi_b_q = chi_b_q, psi_b_q = psi_b_q, chi_alpha_q = chi_alpha_q, psi_alpha_q = psi_alpha_q, delta2_0 = delta2_0, delta1_0 = delta1_0, a_z_0 = a_z_0, b_z_0 = b_z_0, a_u_0 = a_u_0, b_u_0 = b_u_q, control = list(fnscale = -1), method = 'L-BFGS-B', lower = rep(1e-10, p+q), upper = rep(1e10, p + q))$par)
    # } else{
    #   lambda2_c <- tryCatch(optim(lambda2, fn = elbo_lambda2, gr = dev_elbo, Y = Y, K = K, p_var = p, q = q, W = W_mat, Xs_mat = Xs_mat, delta1_q = delta1_q, delta2_q = delta2_q, Sigma_b = Sigma_b_q, mu_b =  mu_b_q, pz = pz_q, a_z = a_z_q, b_z = b_z_q, Sigma_alpha = Sigma_alpha_q, mu_alpha = mu_alpha_q, pu = pu_q, a_u = a_u_q, b_u = b_u_q, chi_b_q = chi_b_q, psi_b_q = psi_b_q, chi_alpha_q = chi_alpha_q, psi_alpha_q = psi_alpha_q, delta2_0 = delta2_0, delta1_0 = delta1_0, a_z_0 = a_z_0, b_z_0 = b_z_0, a_u_0 = a_u_0, b_u_0 = b_u_q, control = list(fnscale = -1), method = 'L-BFGS-B', lower = rep(1e-10, p+q), upper = rep(1e10, p + q))$par)
    # }
    # lambda2 <- lambda2_c
    
    lambda2 <- lambda2_hat(K = K, p_var = p, q = q,
                           chi_b_q = chi_b_q, psi_b_q = psi_b_q,
                           chi_alpha_q = chi_alpha_q, psi_alpha_q = psi_alpha_q)

    lambda2_b <- lambda2[1:p]
    lambda2_alpha <- lambda2[c((p+1):(p+q))]
    
    iter = iter + 1 
    if(std){
      elbo_c <- elbo_partial(Y_std, K, p, q, W_mat, Xs_mat, delta1_q, delta2_q,
                             Sigma_b_q, mu_b_q, pz_q, a_z_q, b_z_q,
                             Sigma_alpha_q, mu_alpha_q, pu_q, a_u_q, b_u_q,
                             chi_b_q, psi_b_q, chi_alpha_q, psi_alpha_q,
                             delta2_0, delta1_0, a_z_0, b_z_0, 
                             a_u_0, b_u_0, lambda2)
    } else{
      elbo_c <- elbo_partial(Y, K, p, q, W_mat, Xs_mat, delta1_q, delta2_q,
                             Sigma_b_q, mu_b_q, pz_q, a_z_q, b_z_q,
                             Sigma_alpha_q, mu_alpha_q, pu_q, a_u_q, b_u_q,
                             chi_b_q, psi_b_q, chi_alpha_q, psi_alpha_q,
                             delta2_0, delta1_0, a_z_0, b_z_0, 
                             a_u_0, b_u_0, lambda2)
    }
    
    converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)
    
    elbo_prev <- elbo_c
  }
  runtime_VEM <- proc.time() - start
  
  res <- list(mu_b = mu_b_q, Sigma_b = Sigma_b_q, pz = pz_q, pu = pu_q, 
              delta1 = delta1_q, delta2 = delta2_q, 
              chi_b = chi_b_q, psi_b = psi_b_q,
              chi_alpha = chi_alpha_q, psi_alpha = psi_alpha_q,
              mu_alpha = mu_alpha_q, Sigma_alpha = Sigma_alpha_q,
              a_z_q = a_z_q, b_z_q = b_z_q, a_u_q = a_u_q, b_u_q = b_u_q,
              lambda2_b = lambda2_b, lambda2_alpha = lambda2_alpha,
              N_iter = iter, runtime = runtime_VEM[[3]], elbo = elbo_c)
  
  return(res)
}

