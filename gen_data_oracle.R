gen_data_oracle <- function(seed = 1234, nsim, folder, n, nt, B, time_points, Z, K, p){
  X <- array(NA, dim = c(n, nt, p))
  
  set.seed(seed)
  mean_bs_X <- c(rnorm(K, mean = 10, sd = 2), rnorm(K, mean = 2, sd = 1))
  for(i in 1:n){
    set.seed(seed + i)
    X[i,,1] <- B[[1]]%*%mean_bs_X[1:K] + rnorm(nt, 0, 20)
    X[i,,2] <- B[[2]]%*%(mean_bs_X[(K+1):(2*K)]) + rnorm(nt, 0, 20)
  }
  
  # Generate y and smoothed X
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
  
  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[i,,j]%*%J[,,j])
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }
  
  g_ui <- sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points, (X_smooth[i,,j]*(Z[j]*beta[,,j])))}))})
  
  set.seed(seed)
  Y <- g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  data <- list(Y = Y, beta = beta, Xt = X,  W_mat = W_mat, Z = Z, B = B, Xt_smooth = X_smooth,
               sigma2 = sigma2)
  save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  return(data)
}
