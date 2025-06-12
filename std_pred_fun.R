#' Standardize functional predictors
#'
#' @param Xt 
#' @param Y 
#' @param beta 
#' @param K 
#' @param nt 
#' @param p 
#' @param n 
#'
#' @returns
#' @export
#'
#' @examples
std_pred_fun <- function(Xt, Y, beta, K , nt, p, n, time_points){
  
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
  
  
  
  
  return(list(Y = Y, W_mat = W_mat, Y_std = y_std, beta_std = beta_std, Xbar_t = X_bar_t, sd_t = sqrt(Var_func_t), Xt_std = X_smooth_est))
}