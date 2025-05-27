# Oracle simulation 

library(fda)
library(pracma)
library(matrixcalc)

p <- 2
K <- 10
sigma2 <- 0.01
nt <- 50
n <- 100
Z <- c(1,0)

# same range of points of evaluations for all predictors and functional coefficients
time_points <- seq(from = 0, to = 1, length= nt)
basis <- create.bspline.basis(rangeval = range(time_points), nbasis = K)
B <- lapply(1:p, function(j){getbasismatrix(time_points, basis, nderiv = 0)})
ids <- split(1:(K*p), rep(1:p, each = K))
# Step 1: Fix/generate lambda2j
E_lambda2 <- c(1, 150)

# Step 2: Generate taus from prior
set.seed(1234)
tau2 <- c(rexp(K, rate = E_lambda2[1]/2), rexp(K, rate = E_lambda2[2]/2))

# Step 3: Generate bs from prior
set.seed(1234)
b1 <- rnorm(K, sd = sqrt(tau2[ids[[1]]]*sigma2))
b2 <- rnorm(K, sd = sqrt(tau2[ids[[2]]]*sigma2))

# Generate curves for coefficients with K b-splines
beta1 <- B[[1]]%*%b1
beta2 <- rep(0, nt)

beta <- array(NA, dim =  c(nt, 1, p))
beta[,,1] <- beta1
beta[,,2] <- beta2

plot(time_points, beta1, type = "l", col = "red")
plot(time_points, beta2, type = "l", col = "red")


# Run Vb ----

# Implementation - Simplify form (b does depend on sigma2)
library(GPBayes)
library(lqr)
library(fda)

# source R files
source('elbo_v2.R')

# expectations + aux functions
std_pred <- function(Xt, Y, beta, K , nt, p, n){
  
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
  
  # Standardized only for estimation 
  X_smooth_est <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  J <- array(NA, dim = c(K, K, p))
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points), nbasis = K)
      res <- smooth.basis(time_points, X_sd_t[i,,j], basis)#Xt
      A[i,,j] <- res$fd$coefs
      X_smooth_est[i,,j] <- eval.fd(time_points, res$fd)
      J[,,j] <- inprod(basis, basis)
      #plotfit.fd(X_smooth[i,,j], time_points, res$fd)
    }
  }
  
  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[i,,j]%*%J[,,j])
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }
  
  
  beta_std <- array(NA, c(nt, 1, p))
  for(j in 1:p){
    beta_std[,,j] <- beta[,,j]*sqrt(Var_func_t[1,,j])
  }
  
  
  
  
  return(list(W_mat = W_mat, Y_std = y_std, beta_std = beta_std, Xbar_t = X_bar_t, sd_t = sqrt(Var_func_t), Xt_std = X_smooth_est))
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
  
  res <- sum(sapply(notj, function(i){
    mu_qi <- mu[ids[[i]]]
    W_i <- W[,ids[[i]]]
    Sigma_ij <- Sigma[ids[[i]], ids[[j]]]
    
    pz[i]*(sum(diag((t(W_j)%*%W_i)%*%(Sigma_ij + mu_qi%*%t(mu_qj)))))}))
  return(res)
}

sim_oracle <- function(seed, nsim, n, sigma2, folder, Z, p, K, nt, beta, B, time_points, std = TRUE){
  
  seed = seed + nsim
  dir.create(folder, showWarnings = FALSE)
  
  #generate data
  # Generate Xs
  X <- array(NA, dim = c(n, nt, p))
  set.seed(1234)
  mean_bs_X <- c(rnorm(K, mean = 10, sd = 2), rnorm(K, mean = 2, sd = 1))
  for(i in 1:n){
    set.seed(1234 + i + nsim)
    X[i,,1] <- B[[1]]%*%mean_bs_X[1:K] + rnorm(nt, 0, 20)
    X[i,,2] <- B[[2]]%*%(mean_bs_X[(K+1):(2*K)]) + rnorm(nt, 0, 20)
  }
  
  
  plot(time_points, X[1,,1], type = "l");for(i in 2:n){lines(time_points, X[i,,1], col = "grey")}
  plot(time_points, X[1,,2], type = "l");for(i in 2:n){lines(time_points, X[i,,2], col = "grey")}
  
  
  # Generate y and smoothed X
  # Expand X(t)
  X_smooth <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  J <- array(NA, dim = c(K, K, p))
  B <- list()
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points), nbasis = K)
      B[[j]] <- getbasismatrix(time_points, basis, nderiv = 0)
      res <- smooth.basis(time_points, X[i,,j], basis)#Xt
      A[i,,j] <- res$fd$coefs
      X_smooth[i,,j] <- eval.fd(time_points, res$fd)
      J[,,j] <- inprod(basis, basis)
    }
  }
  
  plot(time_points, X_smooth[1,,1], type = "l");for(i in 2:n){lines(time_points, X_smooth[i,,1], col = i)}
  plot(time_points, X_smooth[1,,2], type = "l");for(i in 2:n){lines(time_points, X_smooth[i,,2], col = i)}
  
  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[i,,j]%*%J[,,j])
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }
  
  g_ui = sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points, (X_smooth[i,,j]*(Z[j]*beta[,,j])))}))})
  set.seed(seed)
  Y <- g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  data <- list(Y = Y, beta = beta, Xt = X,  Xt_smooth = X_smooth, sigma2 = sigma2)
  save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  # Standardized predictors to fit model
  if(std == TRUE){
    data_std <- std_pred(data$Xt_smooth, Y, beta, K, nt, p, n)
    #save(data_std, file = paste0(folder,"/datastd_", nsim, ".RData"))
    
   
    
    # Run VB for each simulated dataset
    Xt_std <- data_std$Xt_std
    Y_std <- data_std$Y_std
    W_mat <- data_std$W_mat
    beta_std <- data_std$beta_std
    mean_t <- data_std$Xbar_t
    sd_t <- data_std$sd_t
    
    
    plot(time_points, Xt_std[1,,1], type = "l", ylim = c(-2,2));for(i in 2:n){lines(time_points, Xt_std[i,,1], col = "grey")}
    plot(time_points, Xt_std[1,,2], type = "l", ylim = c(-2,2));for(i in 2:n){lines(time_points, Xt_std[i,,2], col = "grey")}
  }
  # priors
  delta1_0 <- 0.0001 
  delta2_0 <- 0.0001 
  a0 <- 0.5 #change these values
  b0 <- 0.5
  shape_lambda_0 <- 1/3#2 #1/3 #0.001 #stan recommend 2,0
  rate_0 <- 0.001
  E_lambda2 <- c(1,150)
  
  # Are not updated within VB
  delta1_q <- n/2 + delta1_0 + (K*p)/2
  shape_lambda_q <- rep(K + shape_lambda_0, p)
  
  # initial values
  Sigma0 = diag(0.1, K*p)
  if(std == TRUE){
    mu0 = as.vector(sapply(1:p, function(j){as.vector(lm(beta_std[,,j] ~ B[[j]]-1)$coef)}))
    plot(time_points, beta[,,1], type = "l", col = "red")
    lines(time_points, B[[1]]%*%mu0[ids[[1]]]/sd_t[,,1], col = "blue", lty = 2)
  } else{
    mu0 = as.vector(sapply(1:p, function(j){as.vector(lm(beta[,,j] ~ B[[j]]-1)$coef)}))
    plot(time_points, beta[,,1], type = "l", col = "red")
    lines(time_points, B[[1]]%*%mu0[ids[[1]]], col = "blue", lty = 2)
  }
  
  #Initial values
  delta2_q <- 1#(delta1_q - 1)*var(Y_std)
  E_inv_sigma2 <- delta1_q/delta2_q
  
  Sigma_b_q <- Sigma0
  mu_b_q <- mu0
  pz_q <- Z#rep(1,p)
  

  psi_q <- rep(NA, K*p)
  rate_q <- rep(NA, p)

  E_eta <- rep(NA, K*p)#1/tau2
  E_tau2<- rep(NA, K*p)#tau2
  
  Niter = 100 #100
  iter = 1
  elbo_prev = 0
  converged <- FALSE
  convergence_threshold = 0.001
  start <- proc.time()
  mu_b_q_c <- mu0
  #while(iter < Niter & converged == FALSE){
  while(iter < Niter){
    
    # Step 1: Update variational of tau2
    chi_q <- (diag(Sigma_b_q) + mu_b_q^2)*as.numeric(E_inv_sigma2)

    for(j in 1:p){
      psi_q[ids[[j]]] <- rep(E_lambda2[j], K)
    }

    for(kj in 1:(K*p)){
      E_eta[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], func = "1/x")
      E_tau2[kj] <- Egig(lambda = 0.5, chi = chi_q[kj], psi = psi_q[kj], fun = "x")

      #incorrect:
      #E_eta[kj] <- Egig(lambda = 0.5, chi = psi_q[kj], psi = chi_q[kj], func = "1/x")
      #E_tau2[kj] <- Egig(lambda = 0.5, chi = psi_q[kj], psi = chi_q[kj], fun = "x")
    }
    
    # Step 1.1 update variational for lambda2
    for(j in 1:p){
      rate_q[j] <- rate_0 + 0.5*sum(E_tau2[ids[[j]]])
      E_lambda2[j] <- shape_lambda_q[j]/rate_q[j]
    }
    
    # Step 2: Update variational of sigma2
    if(std == TRUE){
      A <- E_quad_y(Y_std, pz_q, mu_b_q, W_mat, K, Sigma_b_q) + sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2))#E_quad_b_Z(Sigma_b_q, mu_b_q, pz_q, W_mat, K, p) + sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2))
    } else{
      A <- E_quad_y(Y, pz_q, mu_b_q, W_mat, K, Sigma_b_q) + sum(E_eta*(diag(Sigma_b_q) + mu_b_q^2))
    }
    
    delta2_q <- A/2 + delta2_0
    
    E_inv_sigma2 <- delta1_q/delta2_q
    if(E_inv_sigma2 < 0){stop('variance cannot be negative!')}
    
    # Step 3: Update variational of b
    pz_long <- rep(pz_q, each = K)
    
    Omega <- pz_long%*%t(pz_long) + diag(pz_long)%*%(diag(1, K*p) - diag(pz_long))
    
    Q <- (diag(E_eta) + ((t(W_mat)%*%W_mat)*Omega))
    if(is.singular.matrix(Q)){warning()}#; print(det(Q))} #if(is.singular.matrix(Q)){stop(); print(det(Q))}
    Sigma_b_q <- solve(as.numeric(E_inv_sigma2)*Q)
    if(std == TRUE){
      mu_b_q <- solve(Q)%*%(diag(pz_long)%*%t(W_mat)%*%Y_std)
    } else{
      mu_b_q <- solve(Q)%*%(diag(pz_long)%*%t(W_mat)%*%Y)
    }
   
    plot(time_points, B[[1]]%*%mu_b_q[ids[[1]]], type = "l")
    #print(abs(sum(mu_b_q[ids[[1]]] - mu_b_q_c[ids[[1]]])))
    #mu_b_q_c <- mu_b_q
    
    
    # Step 4: Update variational of theta
    a_q <- pz_q + a0
    b_q <- 2 - pz_q - b0
    
    #Step 5: Update variational of Z
    #for each j = 1, ..., p
    if(std == TRUE){
      for(j in 1:p){
        mu_qj <- mu_b_q[ids[[j]]]
        W_j <- W_mat[,ids[[j]]]
        Sigma_qj <- Sigma_b_q[ids[[j]], ids[[j]]]

        uzj <- digamma(a_q[j]) - digamma(b_q[j]) + as.numeric(E_inv_sigma2)*(t(mu_qj)%*%t(W_j)%*%Y_std - (sum(diag((t(W_j)%*%W_j)%*%(Sigma_qj + mu_qj%*%t(mu_qj)))))/2 - Sum_zi_notzj(j, p, W_mat, Sigma_b_q, mu_b_q, ids, pz_q))

        pz_q[j] <- if(uzj > 709){
          1
        } else{
          exp(uzj)/(1+exp(uzj))
        }
      }
    } else{
      for(j in 1:p){
        mu_qj <- mu_b_q[ids[[j]]]
        W_j <- W_mat[,ids[[j]]]
        Sigma_qj <- Sigma_b_q[ids[[j]], ids[[j]]]

        uzj <- digamma(a_q[j]) - digamma(b_q[j]) + as.numeric(E_inv_sigma2)*(t(mu_qj)%*%t(W_j)%*%Y - sum(diag((t(W_j)%*%W_j)%*%(Sigma_qj + mu_qj%*%t(mu_qj))))/2 - Sum_zi_notzj(j, p, W_mat, Sigma_b_q, mu_b_q, ids, pz_q))
        
        sum(diag(t(W_j)%*%W_j%*%Sigma_qj)) + t(mu_qj)%*%(t(W_j)%*%W_j)%*%mu_qj

        pz_q[j] <- if(uzj > 709){
          1
        } else{
          exp(uzj)/(1+exp(uzj))
        }
      }
    }

    print(pz_q)
    
    mu_b_q_res <- array(mu_b_q, c(K, 1, p))
    beta_hat <- array(NA, c(nt, 1, p))
    for(j in 1:p){
      beta_hat[,,j] <- B[[j]]%*%mu_b_q_res[,,j]
    }
    
    
    
    
    iter = iter + 1 
    # if(std == TRUE){
    #   elbo_c <- elbo(Y_std, K, p, W_mat, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz_q, a_q, b_q, chi_q, psi_q, delta2_0, delta1_0, shape_lambda_0, shape_lambda_q, rate_0, rate_q, a0, b0)
    # } else{
    #   elbo_c <- elbo(Y, K, p, W_mat, delta1_q, delta2_q, Sigma_b_q, mu_b_q, pz_q, a_q, b_q, chi_q, psi_q, delta2_0, delta1_0, shape_lambda_0, shape_lambda_q, rate_0, rate_q, a0, b0)
    # }
    # 
    # 
    # converged <- check_convergence(elbo_c, elbo_prev, convergence_threshold)
    # 
    # 
    # elbo_prev <- elbo_c
    # print(elbo_c)
    print(iter)
    
  }  
  
  runtime_VB <- proc.time() - start
  
  if(std == TRUE){
    for(j in 1:p){
      plot(time_points, beta[,,j], ylab = paste0('beta',j), xlab = expression(t), type = 'l', col = 'red')
      lines(time_points, beta_hat[,,j]/sd_t[,,j], col = "blue")
    }
    
    Z_hat <- ifelse(pz_q > 0.5, 1, 0)
    yhat_std <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q_res[,,j])}))
    y_hat <- yhat_std + mean(Y)
    
    plot(Y, y_hat)
    abline(0,1)
  } else{
    for(j in 1:p){
      plot(time_points, beta[,,j], ylab = paste0('beta',j), xlab = expression(t), type = 'l', col = 'red')
      lines(time_points, beta_hat[,,j], col = "blue")
    }
    
    Z_hat <- ifelse(pz_q > 0.5, 1, 0)
    y_hat <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q_res[,,j])}))
    plot(Y, y_hat)
    abline(0,1)
    
  }
  
  res <- list(mu_b = mu_b_q, Sigma_b = Sigma_b_q, delta1 = delta1_q, delta2 = delta2_q, a = a_q, b = b_q, pz = pz_q, E_lambda2, N_iter = iter, runtime = runtime_VB[[3]], elbo = elbo_c)
  
  #source('Plot_results.R')
  
  
  return(res)
}    

res <- sim_oracle(1234,1,n,0.1, 'TESTE3/tau2j', Z, p, K, nt, beta, B, time_points, std = TRUE)


results <- lapply(1:100, function(i){sim_oracle(1234,i,n,0.01, 'TESTE3/tau2j', Z, p, K, nt, beta, B, time_points, std = FALSE)})
save(results, file = 'TESTE3/tau2j/results.RData')


ids <- split(1:(K*p), rep(1:p, each = K))

res <- c()
res2 <- c()
res3 <- c()
load(paste0('TESTE3/tau2j',"/results.RData"))




# results for selection
nsim = 100
res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
res_pz <- do.call(rbind, lapply(results, `[[`, 7))
res_lambda2 <- do.call(rbind,lapply(results, `[[`, 8))

Zhat <- ifelse(res_pz > 0.5, 1,0)
n_sel_j <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == Z)})))))
n_sel_j2 <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == 1)})))))
n_sel <- sum(sapply(1:nsim, function(sim){sum(Zhat[sim,] == Z) == p}))



# AMSE
amse <- mean(sapply(1:nsim, function(i){rss_std(i, data, results)}))

library(xtable)
options(xtable.include.rownames = FALSE, xtable.booktabs = TRUE, xtable.caption.placement = "top", xtable.sanitize.text.function = function(x){x})

# if std == TRUE
beta_hat <- array(NA, c(nt, 1, p))
beta_hat_all <- array(NA, c(nsim, nt, p))
for(j in 1:p){
  for(s in 1:nsim){
    load(paste0('TESTE3/std',"/datastd_", j, ".RData"))
    sd_t <- data_std$sd_t
    beta_hat_all[s,,j] <- (B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j]
  }
  beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
}

# if std == FALSE
beta_hat <- array(NA, c(nt, 1, p))
beta_hat_all <- array(NA, c(nsim, nt, p))
for(j in 1:p){
  for(s in 1:nsim){
    beta_hat_all[,,j] <- t(sapply(1:nsim, function(s){(B[[j]]%*%(res_b[s,ids[[j]]]))}))
  }
  beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
}


pdf(paste0(case,'/results_VB.pdf'), width = 6, height = 6)
for(j in 1:p){
  plot(time_points, beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", lty = 2)
  for(s in 1:nsim){
    lines(time_points, beta_hat_all[s,,j], col = "gray")
  }
  lines(time_points, beta_hat[,,j], col = "blue")
  lines(time_points, beta[,,j], col = 'red', lty = 2)
  #legend("topleft", c('True', 'Estimated'), col = c("red", "blue"), lty = c(2,1))
}
dev.off()

