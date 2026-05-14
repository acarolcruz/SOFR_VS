gen_data_oracle <- function(seed = 1234, nsim, folder, n, nt, Z, K, p, sigma2){
  
  time_points <- matrix(rep(seq(from = 0, to = 1, length= nt), p), ncol = p, nrow = nt)
  basis <- create.bspline.basis(rangeval = range(time_points[,1]), nbasis = K)
  B <- lapply(1:p, function(j){getbasismatrix(time_points[,1], basis, nderiv = 0)})
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  # Step 1: Fix/generate lambda2j
  lambda2 <- c(0.001, 100)
  
  # Step 2: Generate taus from prior
  set.seed(2025)
  tau2 <- c(rexp(K, rate = lambda2[1]/2), rexp(K, rate = lambda2[2]/2))
  
  # Step 3: Generate bs from prior
  set.seed(1389)
  b1 <- rnorm(K, sd = sqrt(tau2[ids[[1]]]*sigma2))
  b2 <- rnorm(K, sd = sqrt(tau2[ids[[2]]]*sigma2))
  
  # Generate curves for coefficients with K b-splines
  beta1 <- B[[1]]%*%b1
  beta2 <- rep(0, nt)
  
  beta <- array(NA, dim =  c(nt, 1, p))
  beta[,,1] <- beta1
  beta[,,2] <- beta2
  
  #plot(time_points[,1], beta1)

  X <- array(NA, dim = c(n, nt, p))
  
  set.seed(seed)
  mean_bs_X <- c(rnorm(K, mean = 5, sd = 10), rnorm(K, mean = 2, sd = 1))
  for(i in 1:n){
    set.seed(seed + i)
    X[i,,1] <- B[[1]]%*%mean_bs_X[1:K] + rnorm(nt, 0, 10)
    #X[i,,1] <- B[[1]]%*%rnorm(K, 1, 10)
    X[i,,2] <- B[[2]]%*%(mean_bs_X[(K+1):(2*K)]) + rnorm(nt, 0, 10)
    #X[i,,1] <- B[[1]]%*%rnorm(K, 20, 10)
    #X[i,,2] <- B[[2]]%*%rnorm(K, 0, 1)# B[[2]]%*%rnorm(K, 0, 1)
    
  }
  

  #plot(time_points[,1], X[i,,1], type = "l")
  #for(i in 2:n){lines(time_points[,1], X[i,,1], type = "l")}
  
  # Generate y and smoothed X
  X_smooth <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  J <- array(NA, dim = c(K, K, p))
  B <- list()
  for(i in 1:n){
    for(j in 1:p){
      basis <- create.bspline.basis(rangeval = range(time_points[,j]), nbasis = K)
      B[[j]] <- getbasismatrix(time_points[,j], basis, nderiv = 0)
      res <- smooth.basis(time_points[,j], X[i,,j], basis)#Xt
      A[i,,j] <- res$fd$coefs
      X_smooth[i,,j] <- eval.fd(time_points[,j], res$fd)
      J[,,j] <- inprod(basis, basis)
    }
  }
  
  #plot(time_points[,1], X_smooth[i,,1], type = "l")
  #for(i in 2:n){lines(time_points[,1], X_smooth[i,,1], type = "l")}
  
  
  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[i,,j]%*%J[,,j])
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }
  
  g_ui <- sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points[,j], (X[i,,j]*(Z[j]*beta[,,j])))}))})
  
  set.seed(seed)
  Y <- 10 +  g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))
  #boxplot(Y)
  
  data <- list(Y = Y, beta = beta, Xt = X,  W_mat = W_mat, Z = Z, B = B, Xt_smooth = X_smooth,
               sigma2 = sigma2, time_points = time_points)
  
  save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  return(data)
}
