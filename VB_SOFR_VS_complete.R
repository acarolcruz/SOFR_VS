# Implementation - Simplify form (b does depend on sigma2)
library(GPBayes)
library(lqr)
library(fda)

# source R files
source('elbo.R')

# expectations + aux functions

check_convergence <- function(elbo_c, elbo_prev, convergence_threshold) {
  if(is.null(elbo_prev) == TRUE) {
    return(FALSE)
  }
  else{
    dif <- elbo_c - elbo_prev
    if(abs(dif)  <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
}  

E_quad_b_Z <- function(Sigma, mu, pz, W, K, p){
  pz_long <- rep(pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, K*p) - diag(pz_long))
  
  sum(diag((Sigma + mu%*%t(mu))%*%((t(W)%*%W)*Omega)))
}

E_quad_y <- function(Y, pz, mu, W, K, Sigma){
  pz_long <- rep(pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)*(diag(1, K*p) - diag(pz_long))
  
  t(Y)%*%Y - 2*t(Y)%*%W%*%diag(pz_long)%*%mu + 
    sum(diag((Sigma + mu%*%t(mu))%*%((t(W)%*%W)*Omega)))
}

#E_quad_y(Y, Z, mu_q, W_mat, K_b, Sigma_q)

Sum_zi_notzj <- function(j, p, W, Sigma, mu, ids, pz){
  notj <- (1:p)[1:p != j]
  W_j <- W[,ids[[j]]]
  mu_qj <- mu[ids[[j]]]
  sum(sapply(notj, function(i){
    mu_qi <- mu[ids[[i]]]
    W_i <- W[,ids[[i]]]
    Sigma_ij <- Sigma[ids[[j]], ids[[i]]]
    
    pz[i]*sum(diag((t(W_j)%*%W_i)%*%Sigma_ij + as.vector(mu_qj%*%(t(W_j)%*%W_i))%*%t(mu_qi)))}))
}




