# model only estimation with tau2kj with IG as prior

#source('data_example_ProfRonaldo.R')

p = 6
K = 6
nt = 50
n = 100
sigma2 = 0.01

E_quad_y_wb <- function(Y, mu_b_q, Sigma_b_q, W){
  res <- t(Y)%*%Y - 2*t(mu_b_q)%*%t(W)%*%Y + sum(diag(Sigma_b_q)) + t(mu_b_q)%*%(t(W)%*%W)%*%mu_b_q
  
  return(res) 
}

data <- gen_data(1234, p, n, nt, K, sigma2)

time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))

ids <- split(1:(K*p), rep(1:p, each = K))

Y <- data$Y
W_mat <- data$W_mat
B <- data$B
beta <- data$beta

a0 = 2
b0 = 0.001
delta2_0 = 0.001
delta1_0 = 0.001
E_lambda2 <- rep(0.0001, K)

a_tau2_q <- a0 + 0.5
delta1_q <- n/2 + delta1_0 + K*p/2  
delta2_q <- (delta1_q - 1)*var(Y)

E_inv_sigma2 <- delta1_q/delta2_q

Sigma0 = diag(0.01, K*p)
mu0 = as.vector(sapply(1:p, function(j){as.vector(lm(beta[,,j] ~ B[[j]] - 1)$coef)}))
plot(time_points[,2], beta[,,2], type = "l", col = "red")
lines(time_points[,2], B[[2]]%*%mu0[ids[[2]]], col = "blue", lty = 2)

Sigma_b_q <- Sigma0
mu_b_q <- mu0

iter = 1
while(iter < 100){
  # Update tau2
  b_tau2_q <- (diag(Sigma_b_q) + mu_b_q^2)*as.numeric(E_inv_sigma2)/2 + b0
  E_eta <- rep(0.001, K*p)
  #a_tau2_q/b_tau2_q
  
  # Update sigma2
  delta2_q <- E_quad_y_wb(Y, mu_b_q, Sigma_b_q, W_mat)/2 + delta2_0 + sum(E_eta*(mu_b_q^2 + diag(Sigma_b_q)))/2
  E_inv_sigma2 <- delta1_q/delta2_q
  if(E_inv_sigma2 < 0){stop('variance cannot be negative!')}
  
  Q <- (t(W_mat)%*%W_mat + diag(E_eta))
  if(is.singular.matrix(Q)){warning()}
  Sigma_b_q <- solve(as.numeric(E_inv_sigma2)*Q)
  mu_b_q <- solve(Q)%*%(t(W_mat)%*%Y)
  
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

# only tau2kj no sigma2 in b


# model only estimation with tau2kj with IG as prior

#source('data_example_ProfRonaldo.R')

p = 6
K = 6
nt = 50
n = 100
sigma2 = 0.01

E_quad_y_wb <- function(Y, mu_b_q, Sigma_b_q, W){
  res <- t(Y)%*%Y - 2*t(mu_b_q)%*%t(W)%*%Y + sum(diag(Sigma_b_q)) + t(mu_b_q)%*%(t(W)%*%W)%*%mu_b_q
  
  return(res) 
}

data <- gen_data(1234, p, n, nt, K, sigma2)

time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))

ids <- split(1:(K*p), rep(1:p, each = K))

Y <- data$Y
W_mat <- data$W_mat
B <- data$B
beta <- data$beta

a0 = 2
b0 = 0.001
delta2_0 = 0.001
delta1_0 = 0.001
E_lambda2 <- rep(0.0001, K)

a_tau2_q <- a0 + 0.5
delta1_q <- n/2 + delta1_0 + K*p/2  
delta2_q <- (delta1_q - 1)*var(Y)

E_inv_sigma2 <- delta1_q/delta2_q

Sigma0 = diag(1, K*p)
mu0 = as.vector(sapply(1:p, function(j){as.vector(lm(beta[,,j] ~ B[[j]] - 1)$coef)}))
plot(time_points[,2], beta[,,2], type = "l", col = "red")
lines(time_points[,2], B[[2]]%*%mu0[ids[[2]]], col = "blue", lty = 2)

Sigma_b_q <- Sigma0
mu_b_q <- mu0

iter = 1
while(iter < 100){
  # Update tau2
  b_tau2_q <- (diag(Sigma_b_q) + mu_b_q^2)/2 + b0
  E_eta <- rep(0.0001, K*p)
  #<- a_tau2_q/b_tau2_q
  
  # Update sigma2
  delta2_q <- E_quad_y_wb(Y, mu_b_q, Sigma_b_q, W_mat)/2 + delta2_0 + sum((mu_b_q^2 + diag(Sigma_b_q)))/2
  E_inv_sigma2 <- delta1_q/delta2_q
  if(E_inv_sigma2 < 0){stop('variance cannot be negative!')}
  
  Q <- (t(W_mat)%*%W_mat) + diag(E_eta))
  if(is.singular.matrix(Q)){warning()}
  Sigma_b_q <- solve(as.numeric(E_inv_sigma2)*Q)
  mu_b_q <- solve(Q)%*%(t(W_mat)%*%Y)
  
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

