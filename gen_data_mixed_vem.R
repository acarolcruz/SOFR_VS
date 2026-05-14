gen_data_mixed <- function(seed = 1234, nsim, folder, n, nt, Z, u, K, p, q, sigma2){
  
  time_points <- matrix(rep(seq(from = 0, to = 1, length= nt), p), ncol = p, nrow = nt)
  basis <- create.bspline.basis(rangeval = range(time_points[,1]), nbasis = K)
  B <- lapply(1:p, function(j){getbasismatrix(time_points[,1], basis, nderiv = 0)})
  ids <- split(1:(K*p), rep(1:p, each = K))

  # Step 1: Fix/generate lambda2j
  E_lambda2f <- c(0.0001, 100)
  E_lambda2s <- c(100, 0.001)
  
  #E_lambda2f <- c(0.0001, 100)
  #E_lambda2s <- c(100, 0.0001)
  
  #E_lambda2f <- c(0.00001, 100)
  #E_lambda2s <- c(100, 0.00001)

  # Step 2: Generate taus from prior
  #set.seed(2025)set.seed(5547)(2574)
  set.seed(1485)
  tau2 <- c(rexp(K, rate = E_lambda2f[1]/2), rexp(K, rate = E_lambda2f[2]/2))
  nu2 <- c(rexp(1, rate = E_lambda2s[1]/2), rexp(1, rate = E_lambda2s[2]/2))

  # Step 3: Generate bs from prior
  #set.seed(1389)
  # set.seed(1234) better with k = 6
  #set.seed(1237)  set.seed(1997) better with k = 6
  #set.seed(2025) #good example to show to pedro and prof.camila
  set.seed(2019)
  b1 <- rnorm(K, sd = sqrt(tau2[ids[[1]]]*sigma2))
  b2 <- rnorm(K, sd = sqrt(tau2[ids[[2]]]*sigma2))
  a_coef <- c(rnorm(1, sd = sqrt(nu2[1]*sigma2)), rnorm(1, sd = sqrt(nu2[2]*sigma2)))
  a_coef <- c(0, a_coef[2])

  # Generate curves for coefficients with K b-splines
  beta1 <- B[[1]]%*%b1
  beta2 <- rep(0, nt)

  beta <- array(NA, dim =  c(nt, 1, p))
  beta[,,1] <- beta1
  beta[,,2] <- beta2

  # generate scalar predictors
  Xs <- matrix(NA, ncol = q, nrow = n)

  set.seed(seed)
  Xs[,1] <- rnorm(n, 10, sd = 2)
  Xs[,2] <- rnorm(n, 20, sd = 2)

  # generate functional predictors
  X <- array(NA, dim = c(n, nt, p))

  set.seed(seed)
  mean_bs_X <- c(rnorm(K, mean = 10, sd = 2), rnorm(K, mean = 2, sd = 1))
  for(i in 1:n){
    set.seed(seed + i)
    X[i,,1] <- B[[1]]%*%mean_bs_X[1:K] + rnorm(nt, 0, 10)
    X[i,,2] <- B[[2]]%*%(mean_bs_X[(K+1):(2*K)]) + rnorm(nt, 0, 10)
  }




  # X<- array(NA, c(n, nt, p))
  # X[,,1] <- t(sapply(1:n, function(n){gen_x(seed + n, 1, nt, time_points[,1])}))
  # X[,,2] <- t(sapply(1:n, function(n){gen_x(seed + n, 2, nt, time_points[,2])}))


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

  W_mat <- c()
  resp <- c()
  for(i in 1:n){
    for(j in 1:p){
      resp <- cbind(resp,A[i,,j]%*%J[,,j])
    }
    W_mat <- rbind(W_mat, resp)
    resp <- c()
  }

  #plot(time_points[,1], X_smooth[1,,1], type = "l", ylim = c(-10,20)); for(i in 2:n){lines(time_points[,1], X_smooth[i,,1], type = "l")}


  # generate y
  g_ui <- sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points[,j], (X_smooth[i,,j]*(Z[j]*beta[,,j])))}))})

  set.seed(seed)
  Y <- 30 + Xs%*%((diag(u))%*%a_coef) + g_ui + rnorm(n, mean = 0, sd = sqrt(sigma2))

  #print(range(Y))

  #plot(Xs[,1], Y)
  #plot(Xs[,2], Y)

  data <- list(Y = Y, beta = beta, alpha = a_coef, Xt = X, W_mat = W_mat,  Z = Z, B = B, Xt_smooth = X_smooth, Xs = Xs, sigma2 = sigma2, time_points = time_points)

  save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  return(data)
}

