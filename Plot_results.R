# plotting results

# Check results
mu_b_q <- array(res[[1]], c(K, 1, p))

beta_hat <- array(NA, c(nt, 1, p))
for(j in 1:p){
  beta_hat[,,j] <- B[[j]]%*%mu_b_q[,,j]
}

for(j in 1:p){
  plot(time_points, beta[,,j], type = 'l', col = 'red')
  lines(time_points, beta_hat[,,j], col = "blue")
}


yhat <- rowSums(sapply(1:p, function(j){(W_mat[,ids[[j]]]%*%mu_b_q[,,j])}))

plot(Y, yhat)
abline(0,1, col = "red")

plot(Y-yhat)
abline(0, 0, col = "red")

sum((Y-yhat)^2)/n               
