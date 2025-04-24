# Example 1: 
library(fda)
library(abind)

gen_x <- function(seed, p, nt, time_points){
  set.seed(seed + p)
  ui = sapply(1:50, function(k){rnorm(1, 0, k^(-2))})
  phi <- cbind(rep(1,nt), sapply(2:50, function(k){sqrt(2)*cos(k*pi*time_points)}))
  
  X <- 5*rowSums(sapply(1:50, function(k){ui[k]*phi[,k]}))
  return(X)
}


p <- 2
n <- 100 # number of observations
nt <- 50   # number of time points
K <- 10  # number of basis
time_points <- seq(0, 1, length.out = nt)
sigma2 <- 0.1


# generate data
Xt <- array(NA, c(n, nt, p))
for(j in 1:p){
  Xt[,,j] <- t(sapply(1:n, function(n){gen_x(1234 + n, j, 50, time_points)}))
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


basis_data <- create.fourier.basis(range(time_points), nbasis = K)
B_simulated_data <- getbasismatrix(time_points, basis_data, nderiv = 0)[,-1]
B <- lapply(1:p, function(j){B_simulated_data})
K <- ncol(B[[1]])

ids <- split(1:(K*p), rep(1:p, each = K))

set.seed(1234)
beta <- array(NA, c(nt, 1, p))
b <- c(rnorm(10, 2),rnorm(10, -2))
beta[,,1] <- B[[1]]%*%b[ids[[1]]]
beta[,,2] <- B[[2]]%*%b[ids[[2]]]


#delta_t <- c(time_points[2], diff(time_points))
delta_t <- rep(1/nt, nt)
set.seed(1234)
Y <- rowSums(sapply(1:p, function(j){(crossprod(t(Xt[,,j]),beta[,,j]*delta_t))})) + rnorm(n, mean = 0, sd = sqrt(sigma2))


W <- array(NA, c(n, K, p))
for(j in 1:p){
  W[,,j] <- (crossprod(t(Xt[,,j]),B[[j]]*delta_t)) 
}

W_mat <- do.call(cbind, lapply(seq_len(p), function(i) W[,,i])) 


