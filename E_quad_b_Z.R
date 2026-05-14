#' Expectation of bTGammaTGammab
#'
#' @param Sigma 
#' @param mu 
#' @param pz 
#' @param W 
#' @param K 
#' @param p 
#'
#' @returns
#' @export
#'
#' @examples
E_quad_b_Z <- function(Sigma, mu, pz, W, K, p){
  pz_long <- rep(pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, K*p) - diag(pz_long))
  
  sum(diag((Sigma + mu%*%t(mu))%*%((t(W)%*%W)*Omega)))
}
