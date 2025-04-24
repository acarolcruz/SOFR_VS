# VB for SOFR - Estimation only
library(ghyp)

# Expectations: 

E_quad_y_wb <- function(Y, mu_b_q, Sigma_b_q, W){
  res <- t(Y)%*%Y - 2*t(mu_b_q)%*%t(W)%*%Y + sum(diag(Sigma_b_q)) + t(mu_b_q)%*%(t(W)%*%W)%*%mu_b_q
 
  return(res) 
}

source('data_example3.R')

sim <- function(seed, sim){
  seed <- seed + sim
  p <- 2
  K <- 10
  n <- 100
  nt <- 100
  sigma2 <- 0.1
  time_points <- seq(0, 1, length.out = nt)
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  data <- gen_data(seed, p, n, nt, K, sigma2)
  Y <- data$Y
  W_mat <- data$W_mat
  B <- data$B
  beta <- data$beta
  b <- data$b
  
  save(data, file = paste0("Simulation Oracle/data", "_", sim, ".RData"))
  
  delta1_0 <- 0.0001
  delta2_0 <- 0.0001
  E_lambda2 <- c(1,1)#c(5000,100)
  shape_lambda_0 <- 0.001
  rate_0 <- 0.001
  Niter = 10
  delta1_q <- n/2 + delta1_0 + (K*p)/2
  
  #Initial values
  delta2_q <- (delta1_q - 1)*0.1
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
  # initial values
  Sigma0 = diag(0.1, K*p)
  mu0 = as.vector(sapply(1:p, function(j){as.vector(lm(beta[,,j] ~ B[[j]] - 1)$coef)}))
  #plot(time_points, B[[j]]%*%mu0[ids[[1]]])
  Sigma_b_q <- Sigma0
  mu_b_q <- mu0
  iter = 0
  E_eta <- c()
  E_tau2 <- c()
  a_q <- rep(NA, K*p)
  while(iter < Niter){
    
    #Step 1: tau2
    # a_q <- rep(lambda2, K*p) #only one lambda for all p
    for(j in 1:p){
      a_q[ids[[j]]] <- rep(E_lambda2[j], K)
    }
    b_q <- (mu_b_q^2 + diag(Sigma_b_q))*as.numeric(E_inv_sigma2) 
    for(kj in 1:(K*p)){
      E_eta[kj] <- Egig(lambda = 0.5, chi = a_q[kj], psi = b_q[kj], func = "1/x")
      E_tau2[kj] <- Egig(lambda = 0.5, chi = a_q[kj], psi = b_q[kj], fun = "x")
    }  
    #print(E_eta)
    if(sum(E_eta < 0) >= 1){stop('variance cannot be negative!')}
    
    for(j in 1:p){
       shape_lambda_q <- K + shape_lambda_0
       rate_q <- 0.5*rate_0*sum(E_tau2[ids[[j]]])
       E_lambda2[j] <- shape_lambda_q/rate_q
    }
    
    
    #Step 2: sigma2
    delta2_q <- E_quad_y_wb(Y, mu_b_q, Sigma_b_q, W_mat)/2 + delta2_0 + sum(E_eta*(mu_b_q^2 + diag(Sigma_b_q)))/2
    E_inv_sigma2 <- delta1_q/delta2_q
    if(E_inv_sigma2 < 0){stop('variance cannot be negative!')}
    
    #print(E_inv_sigma2)
    #Step 3: bs
    Q <- as.numeric(E_inv_sigma2)*(t(W_mat)%*%W_mat + diag(E_eta))
    Sigma_b_q <- solve(Q)
    mu_b_q <- Sigma_b_q%*%(as.numeric(E_inv_sigma2)*(t(W_mat)%*%Y))
    
    mu_b_q_res <- array(mu_b_q, c(K, 1, p))
    beta_hat <- array(NA, c(nt, 1, p))
    for(j in 1:p){
      beta_hat[,,j] <- B[[j]]%*%mu_b_q_res[,,j]
    }
    
    for(j in 1:p){
      plot(time_points, beta[,,j], type = 'l', col = 'red')
      lines(time_points, beta_hat[,,j], col = "blue")
    }
    
    iter = iter + 1
  }
  
  res <- list(mu_b_q, Sigma_b_q, delta1_q, delta2_q, a_q, b_q)
  
  #source('Plot_results.R')
  
  
  return(res)
}  

nsim <- 1
sim(seed = 1234, nsim)

results <- lapply(1:100, function(i) sim(seed = 1234, i))

load("Simulation Oracle/data_1.RData")
K <- 10
p <- 2
nt <- 100
B <- data$B
ids <- split(1:(K*p), rep(1:p, each = K))
time_points <- seq(0, 1, length.out = nt)
beta <- data$beta

res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
colMeans(res_b)

beta_hat <- array(NA, c(nt, 1, p))
for(j in 1:p){
  beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])
}

for(j in 1:p){
  plot(time_points, beta[,,j], type = 'l', col = 'red')
  lines(time_points, beta_hat[,,j], col = "blue")
}





