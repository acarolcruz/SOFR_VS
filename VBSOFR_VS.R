#' Variational Bayes for variable selection in scalcar-on-function regression (SOFR) with only functional predictors
#'
#' @param delta1_0 
#' @param delta2_0 
#' @param a0 
#' @param b0 
#' @param shape_lambda_0 
#' @param rate_0 
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
VBSOFR_VS <- function(delta1_0 = 0.0001, delta2_0 = 0.0001, a0 = 0.5, b0 = 0.5,
                      shape_lambda_0 = 2, rate_0 = 0.001, initial_values, data, 
                      data_std, n, K, p, 
                      Niter = 500, convergence_threshold = 0.001, std = TRUE){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  W_mat <- data_std$W_mat
  #W_mat <- data$W_mat
  Y_std <- data_std$Y_std
  Y <- data_std$Y
  # Update VB parambers not updated within VB
  delta1_q <- n/2 + delta1_0 + (K*p)/2
  shape_lambda_q <- rep(K + shape_lambda_0, p)
  
  #Initial values
  E_lambda2 <- initial_values$E_lambda2
  delta2_q <- (delta1_q - 1)*var(Y)
  E_eta <- initial_values$E_eta
  E_inv_sigma2 <- delta1_q/delta2_q
  pz_q <- rep(1, p)#data$Z
  
  # Structure
  psi_q <- rep(NA, K*p)
  rate_q <- rep(NA, p)
  E_tau2<- rep(NA, K*p)
  
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
    
    # Step 2: Update variational of sigma2
    if(std){
      A <- E_quad_y(Y_std, pz_q, mu_b_q, W_mat, K, Sigma_b_q) + sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2))
    } else{
      A <- E_quad_y(Y, pz_q, mu_b_q, W_mat, K, Sigma_b_q) + sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2))
    }
    
    delta2_q <- A/2 + delta2_0
    
    E_inv_sigma2 <- delta1_q/delta2_q
    
    # Step 3: Update variational of tau2
    chi_q <- c((diag(Sigma_b_q) + mu_b_q^2)*as.numeric(E_inv_sigma2))
    for(j in 1:p){
      psi_q[ids[[j]]] <- rep(E_lambda2[j], K)
    }
    for(kj in 1:(K*p)){
      E_eta[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], func = "1/x")
      E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], fun = "x")
    }
    #Step 3.1 update variational for lambda2
    # for(j in 1:p){
    #   rate_q[j] <- rate_0 + 0.5*sum(E_tau2[ids[[j]]])
    #   E_lambda2[j] <- shape_lambda_q[j]/rate_q[j]
    # }
    # Step 4: Update variational of theta
    a_q <- pz_q + a0
    b_q <- 2 - pz_q - b0
    
    #Step 5: Update variational of Z
    #for each j = 1, ..., p
    if(std){
      pz_r <- pz_q
      uz_js <- c()
      for(j in 1:p){
        for(r in 0:1){
          pz_r[j] <- r
          uz_js[r+1] <- -(n/2)*(log(delta2_q) - digamma(delta1_q)) - 0.5*as.numeric(E_inv_sigma2)*E_quad_y(Y_std, pz_r, mu_b_q, W_mat, K, Sigma_b_q) + r*(digamma(a_q[j]) - digamma(a_q[j] + b_q[j])) + (1-r)*(digamma(b_q[j]) - digamma(a_q[j] + b_q[j]))
        }  
        if(sum(exp(uz_js)) == 0){
          pz_q[j] <- c(0,1)[which.max(uz_js)]
        } else if (sum(exp(uz_js)) == Inf) {
          pz_q[j] <- c(0,1)[which.max(uz_js)]
        } else {
          pz_q[j] <- exp(uz_js[2])/sum(exp(uz_js))
        }
      }
    }else{
      for(j in 1:p){
        pz_r <- pz_q
        uz_js <- rep(NA, p)
        for(r in 0:1){
          
          pz_r[j] = r
          
          uz_js[r+1] <- -(n/2)*(log(delta2_q) - digamma(delta1_q)) -0.5*as.numeric(E_inv_sigma2)*E_quad_y(Y, pz_r, mu_b_q, W_mat, K, Sigma_b_q) + r*(digamma(a_q[j]) - digamma(a_q[j] + b_q[j])) + (1-r)*(digamma(b_q[j]) - digamma(a_q[j] + b_q[j]))
          
          # truth
          #uz_js[r+1] <- -0.5*as.numeric(1/sigma2)*(t(Y - W_mat%*%diag(rep(pz_r, each = K))%*%c(b1,b2))%*%(Y - W_mat%*%diag(rep(pz_r, each = K))%*%c(b1,b2))) + r*(digamma(a_q[j]) - digamma(a_q[j] + b_q[j])) + (1-r)*(digamma(b_q[j]) - digamma(a_q[j] + b_q[j]))
        }   
        
        if(sum(exp(uz_js)) == 0){
          pz_q[j] <- c(0,1)[which.max(uz_js)]
        } else if (sum(exp(uz_js)) == Inf) {
          pz_q[j] <- c(0,1)[which.max(uz_js)]
        } else {
          pz_q[j] <- exp(uz_js[2])/sum(exp(uz_js))
        }
      }
    }
    # mu_b_q_res <- array(mu_b_q, c(K, 1, p))
    # beta_hat <- array(NA, c(nt, 1, p))
    # for(j in 1:p){
    #   beta_hat[,,j] <- B[[j]]%*%mu_b_q_res[,,j]
    # }
    
    iter = iter + 1 
    # if(std){
    #   elbo_c <- elbo(Y_std, K, p, W_mat, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz_q, a_q, b_q, chi_q, psi_q, delta2_0, delta1_0, shape_lambda_0, shape_lambda_q, rate_0, rate_q, a0, b0)
    # } else{
    #   elbo_c <- elbo(Y, K, p, W_mat, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz_q, a_q, b_q, chi_q, psi_q, delta2_0, delta1_0, shape_lambda_0, shape_lambda_q, rate_0, rate_q, a0, b0)
    # }
    # 
    # converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)
    # 
    # elbo_prev <- elbo_c
  }
  runtime_VB <- proc.time() - start
  
  res <- list(mu_b = mu_b_q, Sigma_b = Sigma_b_q, delta1 = delta1_q, delta2 = delta2_q, chi_q = chi_q, psi_q = psi_q, a = a_q, b = b_q, pz = pz_q, E_lambda2, N_iter = iter, runtime = runtime_VB[[3]])#, elbo = elbo_c)
  
  return(res)
}