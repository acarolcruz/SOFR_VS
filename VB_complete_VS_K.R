# Implementation - Simplify form (b does depend on sigma2) with different K
library(GPBayes)
library(lqr)
library(fda)

# source R files
#source('elbo_v2.R')

# expectations + aux functions
# Testing different K
source('VB_complete_VS_K.R')

# Example 4 with variable selection
# Consistent variable selection for functional regression models
library(fda)
library(abind)
library(refund)
library(matrixcalc)

gen_data_vs <- function(seed, p, n, nt, K_vector, Z, sigma2, gamma = 0){
  time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))
  
  set.seed(seed)
  a1 <- rnorm(n, -4, 3)
  a2 <- rnorm(n, 7, 1.5)
  b1 <- runif(n, 3, 7)
  b2 <- rnorm(n, 0, 1)
  c1 <- rnorm(n, -3, 1.2)
  c2 <- rnorm(n, 2, 0.5)
  c3 <- rnorm(n, -2, 1)
  d1 <- rnorm(n, -2, 1)
  d2 <- rnorm(n, 3, 1.5)
  e1 <- runif(n, 2, 7)
  e2 <- rnorm(n, 2, 0.4)
  f1 <- rnorm(n, 4, 2)
  f2 <- rnorm(n, -3, 0.5)
  f3 <- rnorm(n, 1, 1)
  
  X1 <- t(sapply(1:n, function(i){cos(2*pi*(time_points[,1] - a1[i])) + a2[i]}))
  X2 <- t(sapply(1:n, function(i){b1[i]*sin(pi*time_points[,2]) + b2[i]}))
  X3 <- t(sapply(1:n, function(i){c1[i]*time_points[,3]^3 + c2[i]*time_points[,3]^2 + c3[i]*time_points[,3]}))
  X4 <- t(sapply(1:n, function(i){sin(2*(time_points[,4] - d1[i])) + d2[i]*time_points[,4]}))
  X5 <- t(sapply(1:n, function(i){e1[i]*cos(2*time_points[,5]) + e2[i]*time_points[,5]}))
  X6 <- t(sapply(1:n, function(i){f1[i]*exp(-time_points[,6]/3) + f2[i]*time_points[,6] + f3[i]}))
  
  beta1 <- sin(time_points[,1])
  beta2 <- sin(2*time_points[,2])
  beta3 <- -gamma*time_points[,3]^2
  beta4 <- sin(2*time_points[,4])
  beta5 <- gamma*sin(pi*time_points[,5])
  beta6 <- rep(0, nt)
  
  Xt <- array(NA, c(n, nt, p))
  Xt[,,1] <- X1
  Xt[,,2] <- X2
  Xt[,,3] <- X3
  Xt[,,4] <- X4
  Xt[,,5] <- X5
  Xt[,,6] <- X6
  
  # Expand X(t)
  X_smooth <- array(NA, dim = c(n, nt, p))
  A <- lapply(1:p, function(j){matrix(NA, nrow = n, ncol = K_vector[j])
  })
  J <- lapply(1:p, function(j){matrix(NA, nrow = K_vector[j], ncol = K_vector[j])
  })
  B <- list()
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points[,j]), nbasis = K_vector[j])
      B[[j]] <- getbasismatrix(time_points[,j], basis, nderiv = 0)
      res <- smooth.basis(time_points[,j], Xt[i,,j], basis)#Xt
      A[[j]][i,] <- res$fd$coefs
      X_smooth[i,,j] <- eval.fd(time_points[,j], res$fd)
      J[[j]] <- inprod(basis, basis)
      #plotfit.fd(Xt[i,,j], time_points[,j], res$fd)
    }
  }
  
  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[[j]][i,]%*%J[[j]])
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }
  
  
  ids <- split(1:sum(K_vector), rep(1:p, K_vector))
  
  
  beta <- array(NA, c(nt, 1, p))
  beta[,,1] <- beta1
  beta[,,2] <- beta2
  beta[,,3] <- beta3
  beta[,,4] <- beta4
  beta[,,5] <- beta5
  beta[,,6] <- beta6
  #plot(time_points, beta[,,1])
  
  
  delta_t <- lapply(1:p, function(j){c(time_points[2,j], diff(time_points[,j]))})
  #g_ui <- rowSums(sapply(1:p, function(j){Z[j]*(crossprod(t(Xt[,,j]),beta[,,j]*delta_t[[j]]))}))
  
  g_ui = sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points[,j], (X_smooth[i,,j]*(Z[j]*beta[,,j])))}))})
  set.seed(seed)
  Y <- g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  data_std <- std_pred(Xt, Y, beta, K, nt, p, n)
  
  return(list(Y = Y, Xt = Xt, B = B, beta = beta, Z = Z, sigma2 = sigma2, 
              W_mat = data_std$W_mat, Y_std = data_std$Y_std, beta_std = data_std$beta_std, 
              Xbar_t = data_std$Xbar_t, sd_t = data_std$sd_t))
}  

std_pred <- function(Xt, Y, beta, K_vector, nt, p, n){
  
  X_bar_t <- array(NA, dim = c(1, nt, p))
  for(j in 1:p){
    X_bar_t[,,j] <- apply(Xt[,,j], 2, mean)
  }
  
  Var_func_t <- array(NA, dim = c(1, nt, p))
  for(j in 1:p){
    Var_func_t[,,j] <- rowSums(sapply(1:n, function(i){(Xt[i,,j] - X_bar_t[,,j])^2}))/(n-1)
  }
  
  
  X_sd_t <- array(NA, dim = c(n, nt, p))
  for(j in 1:p){
    X_sd_t[,,j] <- t(sapply(1:n, function(i){(Xt[i,,j] - X_bar_t[1,,j])/sqrt(Var_func_t[1,,j])}))
  }
  
  y_std <- Y - mean(Y)
  
  time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))
  
  # Expand X_std(t)
  X_smooth_est <- array(NA, dim = c(n, nt, p))
  A <- lapply(1:p, function(j){matrix(NA, nrow = n, ncol = K_vector[j])
  })
  J <- lapply(1:p, function(j){matrix(NA, nrow = K_vector[j], ncol = K_vector[j])
  })
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points[,j]), nbasis = K_vector[j])
      B[[j]] <- getbasismatrix(time_points[,j], basis, nderiv = 0)
      res <- smooth.basis(time_points[,j], X_sd_t[i,,j], basis)#Xt
      A[[j]][i,] <- res$fd$coefs
      X_smooth_est[i,,j] <- eval.fd(time_points[,j], res$fd)
      J[[j]] <- inprod(basis, basis)
      #plotfit.fd(Xt[i,,j], time_points[,j], res$fd)
    }
  }
  
  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[[j]][i,]%*%J[[j]])
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }
  
  beta_std <- array(NA, c(nt, 1, p))
  for(j in 1:p){
    beta_std[,,j] <- beta[,,j]*sqrt(Var_func_t[1,,j])
  }
  
  return(list(W_mat = W_mat, Y_std = y_std, beta_std = beta_std, Xbar_t = X_bar_t, sd_t = sqrt(Var_func_t)))
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

E_quad_b_Z <- function(Sigma, mu, pz, W, K_vector, p){
  pz_long <- rep(pz, K_vector)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, sum(K_vector)) - diag(pz_long))
  
  sum(diag((Sigma + mu%*%t(mu))%*%((t(W)%*%W)*Omega)))
}

E_quad_y <- function(Y, pz, mu, W, K_vector, Sigma){
  pz_long <- rep(pz, K_vector)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, sum(K_vector)) - diag(pz_long))
  
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
sim <- function(seed, nsim, n, sigma2, folder, Z, p, K_vector, nt){
  
  dir.create(folder, showWarnings = FALSE)
  
  seed <- seed + nsim
  time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))
  
  ids <- split(1:(sum(K_vector)), rep(1:p, K_vector))
  
  data <- gen_data_vs(seed, p, n, nt, K_vector, Z, sigma2)
  
  # use Y
  Y_std <- data$Y_std
  Y <- data$Y
  W_mat <- data$W_mat
  B <- data$B
  beta <- data$beta
  beta_std <- data$beta_std
  Z <- data$Z
  
  #save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  delta1_0 <- 0.01#0.0001 
  delta2_0 <- 0.0001#0.0001 
  a0 <- 0.5 #change these values
  b0 <- 0.5
  E_lambda2 <- rep(0.001, p)#c(5000,100)
  shape_lambda_0 <- 2 #2 #1/3 #0.001 #stan recomend 2,0
  rate_0 <- 0.001
  
  # Are not updated within VB
  delta1_q <- n/2 + delta1_0 + (sum(K_vector))/2
  shape_lambda_q <- K_vector + shape_lambda_0
  
  # initial values
  mu0 <- c()
  Sigma0 = diag(1, sum(K_vector))
  for(j in 1:p){
  mu0 = c(mu0, as.vector(lm(beta_std[,,j] ~ B[[j]] - 1)$coef))
  }
  plot(time_points[,1], beta_std[,,1], type = "l", col = "red")
  lines(time_points[,1], B[[1]]%*%mu0[ids[[1]]], col = "blue", lty = 2)
  
  #Initial values
  delta2_q <- (delta1_q - 1)*var(Y_std)
  E_inv_sigma2 <- delta1_q/delta2_q
  
  Sigma_b_q <- Sigma0
  mu_b_q <- mu0
  pz_q <- rep(1,p)
  
  E_eta <- c()
  E_tau2 <- c()
  psi_q <- rep(NA, sum(K_vector))
  rate_q <- rep(NA, p)
  
  Niter = 50
  iter = 1
  elbo_prev = 0
  converged <- FALSE
  convergence_threshold = 0.0001
  start <- proc.time()
  #while(iter < Niter & converged == FALSE){
  while(iter < Niter){
    
    # Step 1: Update variational of tau2
    chi_q <- (diag(Sigma_b_q) + mu_b_q^2)*as.numeric(E_inv_sigma2) 
    
    for(j in 1:p){
      psi_q[ids[[j]]] <- rep(E_lambda2[j], K_vector[j])
    }
    
    for(kj in 1:sum(K_vector)){
      E_eta[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], func = "1/x")
      E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], fun = "x")
      
      #E_eta[kj] <- Egig(lambda = 0.5, chi = psi_q[kj], psi = chi_q[kj], func = "1/x")
      #E_tau2[kj] <- Egig(lambda = 0.5, chi = psi_q[kj], psi = chi_q[kj], fun = "x")
    }  
    
    # Step 1.1 update variational for lambda2
    for(j in 1:p){
      rate_q[j] <- rate_0 + 0.5*sum(E_tau2[ids[[j]]])
      E_lambda2[j] <- shape_lambda_q[j]/rate_q[j]
    }
    
    # Step 2: Update variational of sigma2
    A <- E_quad_y(Y_std, pz_q, mu_b_q, W_mat, K_vector, Sigma_b_q)+ sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2))
    
    delta2_q <- A/2 + delta2_0
    
    E_inv_sigma2 <- delta1_q/delta2_q
    if(E_inv_sigma2 < 0){stop('variance cannot be negative!')}
    
    # Step 3: Update variational of b
    pz_long <- rep(pz_q, K_vector)
    
    Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, sum(K_vector)) - diag(pz_long))
    
    Q <- as.numeric(E_inv_sigma2)*diag(E_eta) + as.numeric(E_inv_sigma2)*((t(W_mat)%*%W_mat)*Omega)
    if(is.singular.matrix(Q)){warning()}#; print(det(Q))} #if(is.singular.matrix(Q)){stop(); print(det(Q))}
    Sigma_b_q <- solve(Q)
    mu_b_q <- Sigma_b_q%*%(as.numeric(E_inv_sigma2)*(diag(pz_long)%*%t(W_mat)%*%Y_std))
    
    # Step 4: Update variational of theta
    #a_q <- pz_q + a0
    #b_q <- 2 - pz_q - b0
    
    # Step 5: Update variational of Z
    # for each j = 1, ..., p
    
    # for(j in 1:p){
    #   mu_qj <- mu_b_q[ids[[j]]]
    #   W_j <- W_mat[,ids[[j]]]
    #   Sigma_qj <- Sigma_b_q[ids[[j]], ids[[j]]]
    #   
    #   uzj <- digamma(a_q[j]) - digamma(b_q[j]) +
    #     as.numeric(E_inv_sigma2)*(t(mu_qj)%*%t(W_j)%*%Y_std -
    #                                 sum(diag(t(W_j)%*%W_j%*%Sigma_qj +
    #                                            as.vector(mu_qj%*%(t(W_j)%*%W_j))%*%t(mu_qj)))/2) -
    #     Sum_zi_notzj(j, p, W_mat, Sigma_b_q, mu_b_q, ids, pz_q)
    #   
    #   pz_q[j] <- if(uzj > 709){
    #     1
    #   } else{
    #     exp(uzj)/(1+exp(uzj))
    #   }
    # }
    # 
    #mu_b_q_res <- lapply(1:p, function(j){matrix(NA, c(K, 1, p))
    beta_hat <- array(NA, c(nt, 1, p))
    for(j in 1:p){
      beta_hat[,,j] <- B[[j]]%*%mu_b_q[ids[[j]]]
    }
    
    
    
    
    iter = iter + 1 
    
    #elbo_c <- elbo(Y_std, K, p, W_mat, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz_q, a_q, b_q, chi_q, psi_q, delta2_0, delta1_0, shape_lambda_0, shape_lambda_q, rate_0, rate_q, a0, b0)
    
    #converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)
    
    
    #elbo_prev <- elbo_c
    #print(elbo_c)
    #print(iter)
  }
  runtime_VB <- proc.time() - start
  
  for(j in 1:p){
    plot(time_points[,j], beta_std[,,j], ylab = paste0('beta',j), xlab = expression(t), type = 'l', col = 'red')
    lines(time_points[,j], beta_hat[,,j], col = "blue", lty = 2)
  }
  
  Z_hat <- ifelse(pz_q > 0.5, 1, 0)
  yhat_std <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q_res[,,j])}))
  y_hat <- yhat_std + mean(Y)
  
  res <- list(mu_b = mu_b_q, Sigma_b = Sigma_b_q, delta1 = delta1_q, delta2 = delta2_q, a = a_q, b = b_q, pz = pz_q, E_lambda2, N_iter = iter, runtime = runtime_VB[[3]])
  
  #source('Plot_results.R')
  
  
  return(res)
}    



