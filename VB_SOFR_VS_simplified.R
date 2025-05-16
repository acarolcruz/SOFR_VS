# Implementation - Simplify form (b does depend on sigma2)
library(GPBayes)
library(lqr)
library(fda)

# source R files
source('elbo_v2.R')

# expectations + aux functions

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

E_quad_b_Z <- function(Sigma, mu, pz, W, K, p){
  pz_long <- rep(pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, K*p) - diag(pz_long))
  
  sum(diag((Sigma + mu%*%t(mu))%*%((t(W)%*%W)*Omega)))
}

E_quad_y <- function(Y, pz, mu, W, K, Sigma){
  pz_long <- rep(pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, K*p) - diag(pz_long))
  
  t(Y)%*%Y - 2*t(Y)%*%W%*%diag(pz_long)%*%mu + sum(diag((Sigma + mu%*%t(mu))%*%((t(W)%*%W)*Omega)))
}

#E_quad_y(Y, Z, mu_q, W_mat, K_b, Sigma_q)

Sum_zi_notzj <- function(j, p, W, Sigma, mu, ids, pz){
  notj <- (1:p)[1:p != j]
  W_j <- W[,ids[[j]]]
  mu_qj <- mu[ids[[j]]]
  sum(sapply(notj, function(i){
    mu_qi <- mu[ids[[i]]]
    W_i <- W[,ids[[i]]]
    Sigma_ij <- Sigma[ids[[j]], ids[[i]]]
    
    pz[i]*sum(diag((t(W_j)%*%W_i)%*%Sigma_ij + as.vector(mu_qj%*%(t(W_j)%*%W_i))%*%t(mu_qi)))}))
}

# Run VB for each simulated dataset
sim <- function(seed, nsim, n, sigma2, folder, Z, p, K, nt){
  
  dir.create(folder, showWarnings = FALSE)
  
  seed <- seed + nsim
  time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  data <- gen_data_vs(seed, p, n, nt, K, Z, sigma2)
  
  # use Y
  Y_std <- data$Y_std
  Y <- data$Y
  W_mat <- data$W_mat
  B <- data$B
  beta <- data$beta
  beta_std <- data$beta_std
  Z <- data$Z
  
  save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  delta1_0 <- 0.01#0.0001 
  delta2_0 <- 0.01#0.0001 
  a0 <- 0.5 #change these values
  b0 <- 0.5
  E_lambda2 <- rep(1, p)#c(5000,100)
  shape_lambda_0 <- 2 #2 #1/3 #0.001 #stan recomend 2,0
  rate_0 <- 0.001
  
  # Are not updated within VB
  delta1_q <- n/2 + delta1_0 + (K*p)/2
  shape_lambda_q <- rep(K + shape_lambda_0, p)
  
  # initial values
  Sigma0 = diag(0.01, K*p)
  mu0 = as.vector(sapply(1:p, function(j){as.vector(lm(beta_std[,,j] ~ B[[j]] - 1)$coef)}))
  #plot(time_points[,j], beta[,,j], type = "l", col = "red")
  #lines(time_points[,j], B[[j]]%*%mu0[ids[[j]]], col = "blue", lty = 2)
  
  #Initial values
  delta2_q <- (delta1_q - 1)*var(Y_std)
  E_inv_sigma2 <- delta1_q/delta2_q
  
  Sigma_b_q <- Sigma0
  mu_b_q <- mu0
  pz_q <- Z #rep(1,p)
  
  E_eta <- c()
  E_tau2 <- c()
  psi_q <- rep(NA, K*p)
  rate_q <- rep(NA, p)
  
  Niter = 1000
  iter = 1
  elbo_prev = 0
  converged <- FALSE
  convergence_threshold = 0.01
  start <- proc.time()
  while(iter < Niter & converged == FALSE){
    #while(iter < Niter){
    
    # Step 1: Update variational of tau2
    chi_q <- (diag(Sigma_b_q) + mu_b_q^2)*as.numeric(E_inv_sigma2) 
    
    for(j in 1:p){
      psi_q[ids[[j]]] <- rep(E_lambda2[j], K)
    }
    
    for(kj in 1:(K*p)){
      E_eta[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], func = "1/x")
      E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], fun = "x")
    }  
    
    # Step 1.1 update variational for lambda2
    for(j in 1:p){
      rate_q[j] <- rate_0 + 0.5*sum(E_tau2[ids[[j]]])
      E_lambda2[j] <- shape_lambda_q[j]/rate_q[j]
    }
    
    # Step 2: Update variational of sigma2
    A <- E_quad_y(Y_std, pz_q, mu_b_q, W_mat, K, Sigma_b_q)+ sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2))
    
    delta2_q <- A/2 + delta2_0
    
    E_inv_sigma2 <- delta1_q/delta2_q
    if(E_inv_sigma2 < 0){stop('variance cannot be negative!')}
    
    # Step 3: Update variational of b
    pz_long <- rep(pz_q, each = K)
    
    Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, K*p) - diag(pz_long))
    
    Q <- as.numeric(E_inv_sigma2)*diag(E_eta) + as.numeric(E_inv_sigma2)*((t(W_mat)%*%W_mat)*Omega)
    if(is.singular.matrix(Q)){warning()}#; print(det(Q))} #if(is.singular.matrix(Q)){stop(); print(det(Q))}
    Sigma_b_q <- solve(Q)
    mu_b_q <- Sigma_b_q%*%(as.numeric(E_inv_sigma2)*(diag(pz_long)%*%t(W_mat)%*%Y_std))
    
    # Step 4: Update variational of theta
    a_q <- pz_q + a0
    b_q <- 2 - pz_q - b0
    
    # Step 5: Update variational of Z
    # for each j = 1, ..., p
    
    for(j in 1:p){
      mu_qj <- mu_b_q[ids[[j]]]
      W_j <- W_mat[,ids[[j]]]
      Sigma_qj <- Sigma_b_q[ids[[j]], ids[[j]]]
      
      uzj <- digamma(a_q[j]) - digamma(b_q[j]) +
        as.numeric(E_inv_sigma2)*(t(mu_qj)%*%t(W_j)%*%Y_std -
                                    sum(diag(t(W_j)%*%W_j%*%Sigma_qj +
                                               as.vector(mu_qj%*%(t(W_j)%*%W_j))%*%t(mu_qj)))/2) -
        Sum_zi_notzj(j, p, W_mat, Sigma_b_q, mu_b_q, ids, pz_q)
      
      pz_q[j] <- if(uzj > 709){
        1
      } else{
        exp(uzj)/(1+exp(uzj))
      }
    }
    
    mu_b_q_res <- array(mu_b_q, c(K, 1, p))
    beta_hat <- array(NA, c(nt, 1, p))
    for(j in 1:p){
      beta_hat[,,j] <- B[[j]]%*%mu_b_q_res[,,j]
    }
    
    
    
    
    iter = iter + 1 
    
    elbo_c <- elbo(Y_std, K, p, W_mat, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz_q, a_q, b_q, chi_q, psi_q, delta2_0, delta1_0, shape_lambda_0, shape_lambda_q, rate_0, rate_q, a0, b0)
    
    converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)
    
    
    elbo_prev <- elbo_c
    #print(elbo_c)
    #print(iter)
  }
  runtime_VB <- proc.time() - start
  
  for(j in 1:p){
    plot(time_points[,j], beta_std[,,j], ylab = paste0('beta',j), xlab = expression(t), type = 'l', col = 'red')
    lines(time_points[,j], beta_hat[,,j], col = "blue")
  }
  
  Z_hat <- ifelse(pz_q > 0.5, 1, 0)
  yhat_std <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q_res[,,j])}))
  y_hat <- yhat_std + mean(Y)
  
  res <- list(mu_b = mu_b_q, Sigma_b = Sigma_b_q, delta1 = delta1_q, delta2 = delta2_q, a = a_q, b = b_q, pz = pz_q, E_lambda2, N_iter = iter, runtime = runtime_VB[[3]])
  
  #source('Plot_results.R')
  
  
  return(res)
}    



