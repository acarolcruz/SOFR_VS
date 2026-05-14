#' Expectations
#'
#' @param Y 
#' @param pz 
#' @param mu 
#' @param W 
#' @param K 
#' @param Sigma 
#'
#' @returns
#' @export
#'
#' @examples
E_quad_y <- function(Y, pz, mu, W, K, Sigma){
  pz_long <- rep(pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, K*p) - diag(pz_long))
  
  t(Y)%*%Y - 2*t(Y)%*%W%*%diag(pz_long)%*%mu + sum(diag((Sigma + mu%*%t(mu))%*%((t(W)%*%W)*Omega)))
}