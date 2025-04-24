# Example 1: 
library(fda)
library(abind)
library(refund)

gen_x <- function(seed, p, nt, time_points){
  set.seed(seed + p)
  ui = sapply(1:50, function(k){rnorm(1, 0, k^(-2))})
  phi <- cbind(rep(1,nt), sapply(2:50, function(k){sqrt(2)*cos(k*pi*time_points)}))
  
  X <- 5*rowSums(sapply(1:50, function(k){ui[k]*phi[,k]}))
  return(X)
}


gen_data <- function(seed, p, n, nt, K, sigma2){

  #seed <- 1234
  #p <- 2
  #n <- 100 # number of observations
  #nt <- 50   # number of time points
  #K <- 10  # number of basis
  time_points <- seq(0, 1, length.out = nt)
  #sigma2 <- 0.1
  
  
  # generate data
  Xt <- array(NA, c(n, nt, p))
  for(j in 1:p){
    Xt[,,j] <- t(sapply(1:n, function(n){gen_x(seed + n, j, 50, time_points)}))
  }
  
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
  
  # Using FPCA
  # B <- NULL
  # X_fpca <- array(NA, dim = c(n, nt, p))
  # for(j in 1:p){
  #   X_fpca[,,j] <- fpca.sc(Y=Xt[,,j], pve=.98, var=TRUE,center = TRUE)$Yhat
  #   B[[j]] <- fpca.sc(Y=Xt[,,j], pve=.98, var=TRUE,center = TRUE)$scores
  # }
  
  # Using B-splines
  basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = K)
  
  # Expand Xs(t)
  
  #Xt <- X_sd_t 
  X_smooth <- array(NA, dim = c(n, nt, p))
  A <- array(NA, dim = c(n, K, p))
  for(i in 1:n){
    for(j in 1:p){
      res <- smooth.basis(time_points, Xt[i,,j], basis)
      A[i,,j] <- res$fd$coefs
      X_smooth[i,,j] <- eval.fd(time_points, res$fd)
      #plotfit.fd(Xt[i,,j], time_points, res$fd)
    }
  }
  
  # choose best based on RMSE or gcv
  # for (nbasis in 4:12) {
  #   basisobj <- create.bspline.basis(c(0,1), nbasis)
  #   ys <- smooth.basis(time_points, Xt[1,,1], basisobj)
  #   xfd <- ys$fd
  #   gcv <- ys$gcv
  #   RMSE <- sqrt(mean((eval.fd(time_points, xfd) - Xt[1,,1])^2))
  #   cat(paste(nbasis,round(RMSE,3),round(gcv,3),"\n"))}
  
  J <- inprod(basis, basis)
  
  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[i,,j]%*%J)
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }
  
  basis_matrix <- getbasismatrix(time_points, basis, nderiv = 0)
  B <- lapply(1:p, function(j){basis_matrix})
  
  #basis_matrix <- getbasismatrix(time_points, basis, nderiv = 0)
  #B <- lapply(1:p, function(j){basis_matrix})
  
  #basis_data <- create.fourier.basis(range(time_points), nbasis = K_b)
  #B_simulated_data <- getbasismatrix(time_points, basis_data, nderiv = 0)[,-1]
  #B <- lapply(1:p, function(j){B_simulated_data})
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  set.seed(1234)
  beta <- array(NA, c(nt, 1, p))
  b <- c(rnorm(10, 2),rnorm(10, -2))
  beta[,,1] <- B[[1]]%*%b[ids[[1]]]
  beta[,,2] <- B[[2]]%*%b[ids[[2]]]
  
  #plot(time_points, beta[,,1])
  
  set.seed(seed)
  Y <- W_mat%*%b + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  return(list(Y = Y, Xt = Xt, W_mat = W_mat, B = B, b = b, beta = beta, sigma2 = sigma2))
}



