gen_x <- function(seed, p, nt, time_points){
  set.seed(seed + p)
  # number of basis used in paper 2025 = 50
  ui <- sapply(1:10, function(k){rnorm(1, 0, k^(-2))}) # morris also used this one in of his papers
  phi <- cbind(rep(1,nt), sapply(2:10, function(k){sqrt(2)*cos(k*pi*time_points)}))
  
  X <- 5*rowSums(sapply(1:10, function(k){ui[k]*phi[,k]}))
  return(X)
}

gen_data_ex2 <- function(seed, n, nt, p, ordem, K, Z, sigma2){
  time_points <- matrix(rep(seq(0, 1, length.out = nt), p), ncol = p, nrow = nt)
  
  Xt <- array(NA, c(n, nt, p))
  Xt[,,1] <- t(sapply(1:n, function(n){gen_x(seed + n, 1, nt, time_points[,1])}))
  Xt[,,2] <- t(sapply(1:n, function(n){gen_x(seed + n, 2, nt, time_points[,2])}))
  Xt[,,3] <- t(sapply(1:n, function(n){gen_x(seed + n, 3, nt, time_points[,3])}))
  Xt[,,4] <- t(sapply(1:n, function(n){gen_x(seed + n, 4, nt, time_points[,4])}))
  
  beta <- array(NA, c(nt, 1, p))
  beta1 <- 2*sin(pi*time_points[,1])# 2*sin(pi*2*time_points[,1])
  beta2 = beta4 = rep(0, nt)
  beta3 <- 1.25*sin(pi*3*time_points[,3]) #2*cos(pi*2*time_points[,1])
  
  # plot(time_points[,1], Xt[1,,1], type = "l");for(i in 2:n){lines(time_points[,1], Xt[i,,1], col = "grey")}
  # plot(time_points[,2], Xt[1,,2], type = "l");for(i in 2:n){lines(time_points[,2], Xt[i,,2], col = "grey")}
  # plot(time_points[,3], Xt[1,,3], type = "l");for(i in 2:n){lines(time_points[,3], Xt[i,,3], col = "grey")}
  # plot(time_points[,4], Xt[1,,4], type = "l");for(i in 2:n){lines(time_points[,4], Xt[i,,4], col = "grey")}
  
  
  # plot(time_points[,1], beta1, type = "l")
  # plot(time_points[,3], beta3, type = "l")
  
  
  # Expand X(t)
  X_smooth <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  J <- array(NA, dim = c(K, K, p))
  B <- list()
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points[,j]), norder = ordem, nbasis = K)
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
  
  # plot(time_points[,1], X_smooth[1,,1], type = "l");for(i in 2:n){lines(time_points[,1], X_smooth[i,,1], col = "grey")}
  # plot(time_points[,2], X_smooth[1,,2], type = "l");for(i in 2:n){lines(time_points[,2], X_smooth[i,,2], col = "grey")}
  # plot(time_points[,3], X_smooth[1,,3], type = "l");for(i in 2:n){lines(time_points[,3], X_smooth[i,,3], col = "grey")}
  # plot(time_points[,4], X_smooth[1,,4], type = "l");for(i in 2:n){lines(time_points[,4], X_smooth[i,,4], col = "grey")}
  
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  
  beta <- array(NA, c(nt, 1, p))
  beta[,,1] <- beta1
  beta[,,2] <- beta2
  beta[,,3] <- beta3
  beta[,,4] <- beta4
  #plot(time_points, beta[,,1])
  
  # for(j in c(1,2,4)){
  #   plot(time_points[,j], beta[,,j], type = "l", ylab = paste0('beta', j))
  # }
  
  #delta_t <- lapply(1:p, function(j){c(time_points[2,j], diff(time_points[,j]))})
  #g_ui <- rowSums(sapply(1:p, function(j){Z[j]*(crossprod(t(Xt[,,j]),beta[,,j]*delta_t[[j]]))}))
  
  # use X_smooth to generate y
  g_ui = sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points[,j], (X_smooth[i,,j]*(Z[j]*beta[,,j])))}))})
  
  set.seed(seed)
  Y <- 20 + g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  # std the smoothed X
  data_std <- std_pred_fun(X_smooth, Y, beta, K, nt, p, n, time_points = time_points, ordem = ordem)
  
  return(list(Y = Y, Xt = Xt, B = B, beta = beta, Z = Z, sigma2 = sigma2, 
              W_mat = data_std$W_mat, Y_std = data_std$Y_std, beta_std = data_std$beta_std, 
              Xbar_t = data_std$Xbar_t, sd_t = data_std$sd_t))
}  




