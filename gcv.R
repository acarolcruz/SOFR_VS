# gcv
gcv_sofr <- function(out, n, Y, W_mat, K, p){
  E_eta <- c()
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  Z_hat <- ifelse(out$pz > 0.5, 1, 0)
  yhat <- mean(Y) + rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%out$mu_b[ids[[j]]])}))
  
  for(kj in 1:(K*p)){
    E_eta[kj] <- Egig(lambda = 0.5, chi = out$chi_q[kj], psi = out$psi_q[kj], func = "1/x")
  }
  
  pz_long <- rep(out$pz, each = K)
  Omega <- pz_long%*%t(pz_long) + diag(pz_long)%*%(diag(1, K*p) - diag(pz_long))
  Q <- (diag(E_eta) + ((t(W_mat)%*%W_mat)*Omega))
  
  Sk <- W_mat%*%diag(pz_long)%*%solve(Q)%*%diag(pz_long)%*%t(W_mat)
  
  res <- ((1/n)*sum((Y - yhat)**2))/((1 - (1/n)*sum(diag(Sk)))**2)
  
  return(res)
}  

