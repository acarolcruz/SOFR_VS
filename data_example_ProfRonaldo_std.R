source('VB_SOFR_VS_complete.R')

# Example 4 with variable selection
# Consistent variable selection for functional regression models
library(fda)
library(abind)
library(refund)
library(matrixcalc)
#seq(0, 1, length.out = nt)

gen_data_vs <- function(seed, p, n, nt, K, Z, sigma2, gamma = 0){
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
  
  plot(time_points[,1], Xt[1,,1], type = "l");for(i in 2:n){lines(time_points[,1], Xt[i,,1], col = "grey")}
  plot(time_points[,2], Xt[1,,2], type = "l");for(i in 2:n){lines(time_points[,2], Xt[i,,2], col = "grey")}
  plot(time_points[,3], Xt[1,,3], type = "l");for(i in 2:n){lines(time_points[,3], Xt[i,,3], col = "grey")} 
  plot(time_points[,4], Xt[1,,4], type = "l");for(i in 2:n){lines(time_points[,4], Xt[i,,4], col = "grey")}
  plot(time_points[,5], Xt[1,,5], type = "l");for(i in 2:n){lines(time_points[,5], Xt[i,,5], col = "grey")}
  plot(time_points[,6], Xt[1,,6], type = "l");for(i in 2:n){lines(time_points[,6], Xt[i,,6], col = "grey")}
  
  # Expand X(t)
  X_smooth <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  J <- array(NA, dim = c(K, K, p))
  B <- list()
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points[,j]), nbasis = K)
      B[[j]] <- getbasismatrix(time_points[,j], basis, nderiv = 0)
      res <- smooth.basis(time_points[,j], Xt[i,,j], basis)#Xt
      A[i,,j] <- res$fd$coefs
      X_smooth[i,,j] <- eval.fd(time_points[,j], res$fd)
      J[,,j] <- inprod(basis, basis)
      #plotfit.fd(Xt[i,,j], time_points[,j], res$fd)
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
  
  for(j in c(1,2,4)){
    plot(time_points[,j], beta[,,j], type = "l", ylab = paste0('beta', j))
  }
  
  delta_t <- lapply(1:p, function(j){c(time_points[2,j], diff(time_points[,j]))})
  g_ui <- rowSums(sapply(1:p, function(j){Z[j]*(crossprod(t(Xt[,,j]),beta[,,j]*delta_t[[j]]))}))
  
  #g_ui = sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points[,j], (X_smooth[i,,j]*(Z[j]*beta[,,j])))}))})
  set.seed(seed)
  Y <- g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  data_std <- std_pred(Xt, Y, beta, K, nt, p, n)
  
  return(list(Y = Y, Xt = Xt, B = B, beta = beta, Z = Z, sigma2 = sigma2, 
              W_mat = data_std$W_mat, Y_std = data_std$Y_std, beta_std = data_std$beta_std, 
              Xbar_t = data_std$Xbar_t, sd_t = data_std$sd_t))
}  

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
  
  time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))
  
  # Standardized only for estimation 
  X_smooth_est <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  J <- array(NA, dim = c(K, K, p))
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points[,j]), nbasis = K)
      res <- smooth.basis(time_points[,j], X_sd_t[i,,j], basis)#Xt
      A[i,,j] <- res$fd$coefs
      X_smooth_est[i,,j] <- eval.fd(time_points[,j], res$fd)
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
  
    
  
  
  return(list(W_mat = W_mat, Y_std = y_std, beta_std = beta_std, Xbar_t = X_bar_t, sd_t = sqrt(Var_func_t)))
}
  


# run example



Z <- c(1,1,0,1,0,0)
p <- 6
K <- 6
nt <- 50

sim(2024,80,300,0.01,'TESTE2', Z, p, K, nt)

results <- lapply(1:100, function(i){sim(1234,i,300,0.01,'TESTE2/K10', Z, p, K, nt)})
save(results, file = 'TESTE2/K10/results.RData')


nsim = 100
n = 50
K <- 10
nt <- 50
gamma <- 0
Z <- c(1,1,0,1,0,0)
p <- 6
ids <- split(1:(K*p), rep(1:p, each = K))

res <- c()
res2 <- c()
res3 <- c()
load(paste0('TESTE2/K10',"/results.RData"))
load(paste0('TESTE2/K10',"/data_1.RData"))




# results for selection
res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
res_pz <- do.call(rbind, lapply(results, `[[`, 7))
Zhat <- ifelse(res_pz > 0.5, 1,0)
n_sel_j <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == Z)})))))
n_sel_j2 <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == 1)})))))
n_sel <- sum(sapply(1:nsim, function(sim){sum(Zhat[sim,] == Z) == p}))

W_mat = data$W_mat

# AMSE
amse <- mean(sapply(1:nsim, function(i){rss_std(i, data, results)}))

library(xtable)
options(xtable.include.rownames = FALSE, xtable.booktabs = TRUE, xtable.caption.placement = "top", xtable.sanitize.text.function = function(x){x})


B <- data$B
beta <- data$beta
sd_t <- data$sd_t

time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))

ids <- split(1:(K*p), rep(1:p, each = K))

res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
res_pz <- do.call(rbind, lapply(results, `[[`, 7))



beta_hat <- array(NA, c(nt, 1, p))
beta_hat_all <- array(NA, c(nsim, nt, p))
for(j in 1:p){
  #beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])/sd_t[,,j]
  beta_hat_all[,,j] <- t(sapply(1:nsim, function(s){(B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j]}))
  beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
}
pdf(paste0(case,'/results_VB.pdf'), width = 6, height = 6)
for(j in 1:p){
  plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", ylim = c(-1.5,1.5), lty = 2)
  for(s in 1:nsim){
    lines(time_points[,j], beta_hat_all[s,,j], col = "gray")
  }
  lines(time_points[,j], beta_hat[,,j], col = "blue")
  lines(time_points[,j], beta[,,j], col = 'red', lty = 2)
  legend("topleft", c('True', 'Estimated'), col = c("red", "blue"), lty = c(2,1))
}
dev.off()

