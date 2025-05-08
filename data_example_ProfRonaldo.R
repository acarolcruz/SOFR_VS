# Example 4
# Consistent variable selection for functional regression models
library(fda)
library(abind)
library(refund)

gen_data <- function(seed, p, n, nt, K, sigma2){
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
  X2 <- t(sapply(1:n, function(i){b1[i]*cos(pi*time_points[,2]) + b2[i]}))
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
  
  # Expand Xs(t)
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
  
  
  #Xt <- X_sd_t 
  X_smooth <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  J <- array(NA, dim = c(K, K, p))
  B <- list()
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points[,j]), nbasis = K)
      B[[j]] <- getbasismatrix(time_points[,j], basis, nderiv = 0)
      res <- smooth.basis(time_points[,j], X_sd_t[i,,j], basis)#Xt
      A[i,,j] <- res$fd$coefs
      X_smooth[i,,j] <- eval.fd(time_points[,j], res$fd)
      J[,,j] <- inprod(basis, basis)
      #plotfit.fd(Xt[i,,j], time_points, res$fd)
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
  
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  
  beta <- array(NA, c(nt, 1, p))
  beta[,,1] <- beta1
  beta[,,2] <- beta2
  beta[,,3] <- beta3
  beta[,,4] <- beta4
  beta[,,5] <- beta5
  beta[,,6] <- beta6
  #plot(time_points, beta[,,1])
  
  
  delta_t <- lapply(1:p, function(j){c(time_points[2,j], diff(time_points[,j]))})
  set.seed(seed)
  Y <- rowSums(sapply(1:p, function(j){(crossprod(t(X_sd_t[,,j]),beta[,,j]*delta_t[[j]]))})) + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  return(list(Y = Y, Xt = Xt, W_mat = W_mat, B = B, beta = beta, sigma2 = sigma2))
  
}

# Run VB
sim <- function(seed, sim){
  seed <- seed + sim
  p <- 6
  K <- 6
  n <- 100
  nt <- 50
  sigma2 <- 0.05
  time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  data <- gen_data(seed, p, n, nt, K, sigma2)
  Y <- data$Y
  W_mat <- data$W_mat
  B <- data$B
  beta <- data$beta
  
  save(data, file = paste0("Simulation Ronaldo/data", "_", sim, ".RData"))
  
  delta1_0 <- 0.0001
  delta2_0 <- 0.0001
  E_lambda2 <- rep(1, p)#c(5000,100)
  shape_lambda_0 <-3 #1/3 #0.001
  rate_0 <- 0.01
  Niter = 50
  delta1_q <- n/2 + delta1_0 + (K*p)/2
  
  # initial values
  Sigma0 = diag(0.01, K*p)
  mu0 = as.vector(sapply(1:p, function(j){as.vector(lm(beta[,,j] ~ B[[j]] - 1)$coef)}))
  plot(time_points[,j], beta[,,j], type = "l", col = "red")
  lines(time_points[,j], B[[j]]%*%mu0[ids[[j]]], col = "blue", lty = 2)
  
  #Initial values
  delta2_q <- (delta1_q - 1)*var(Y)
  
  E_inv_sigma2 <- delta1_q/delta2_q
  
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
      plot(time_points[,j], beta[,,j], ylab = paste0('beta',j), xlab = expression(t), type = 'l', col = 'red')
      lines(time_points[,j], beta_hat[,,j], col = "blue")
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

K <- 6
nt <- 50
n <- 100
gamma <- 0
p <- 6
sigma2 <- 0.05

load("Simulation Ronaldo/data_1.RData")
B <- data$B
beta <- data$beta

time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))

ids <- split(1:(K*p), rep(1:p, each = K))

res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
colMeans(res_b)

beta_hat <- array(NA, c(nt, 1, p))
for(j in 1:p){
  beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])
}
pdf('results_VB.pdf', width = 6, height = 6)
for(j in 1:p){
  plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", ylim = c(-1.5,1.5))
  lines(time_points[,j], beta_hat[,,j], col = "blue")
  legend("topleft", c('True', 'Estimated'), col = c("red", "blue"), lty = 1)
}
dev.off()
save(results, file = 'Simulation Ronaldo/results.RData')
  






