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
  
  #delta_t <- lapply(1:p, function(j){c(time_points[2,j], diff(time_points[,j]))})
  #g_ui <- rowSums(sapply(1:p, function(j){Z[j]*(crossprod(t(Xt[,,j]),beta[,,j]*delta_t[[j]]))}))
  
  # use X_smooth to generate y
  g_ui = sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points[,j], (X_smooth[i,,j]*(Z[j]*beta[,,j])))}))})
  
  set.seed(seed)
  Y <- g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  # std the smoothed X
  data_std <- std_pred_fun(X_smooth, Y, beta, K, nt, p, n, time_points = time_points)
  
  return(list(Y = Y, Xt = Xt, B = B, beta = beta, Z = Z, sigma2 = sigma2, 
              W_mat = data_std$W_mat, Y_std = data_std$Y_std, beta_std = data_std$beta_std, 
              Xbar_t = data_std$Xbar_t, sd_t = data_std$sd_t))
}  