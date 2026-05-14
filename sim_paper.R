# sim_paper not to be in the R Package
sim_paper <- function(seed, nsim, initial_values, folder, 
                      n, sigma2, Z, p, K, nt, 
                      delta1_0, delta2_0, a0, b0,
                      shape_lambda_0, rate_0,
                      Niter, convergence_threshold, ordem, std = TRUE){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  seed = seed + nsim
  dir.create(folder, showWarnings = FALSE)
  
  #generate data
  # Generate Xs
  #data <- gen_data_vs(seed = seed, p = p, n = n, nt = nt, K = K, Z = Z, sigma2 = sigma2, ordem = ordem)
  data <- gen_data_vs_me(seed = seed, p = p, n = n, nt = nt, K = K, Z = Z, sigma2 = sigma2, ordem = ordem)
  save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  
  vb_res <- VBSOFR_VS(initial_values = initial_values, data = data, 
                      data_std = data, n = n, K = K, p = p,
                      delta1_0 = delta1_0, delta2_0 = delta2_0, a0 = a0, b0 = b0,
                      shape_lambda_0 = shape_lambda_0, rate_0 = rate_0, 
                      Niter = Niter, convergence_threshold = convergence_threshold, std = std)
  
  
  # mu_b_q_res <- array(vb_res$mu_b, c(K, 1, p))
  # beta_hat <- array(NA, c(nt, 1, p))
  # for(j in 1:p){
  #   beta_hat[,,j] <- B[[j]]%*%mu_b_q_res[,,j]
  # }
  # 
  # 
  # if(std){
  #   for(j in 1:p){
  #     plot(time_points, beta[,,j], ylab = paste0('beta',j), xlab = expression(t), type = 'l', col = 'red')
  #     lines(time_points, beta_hat[,,j]/sd_t[,,j], col = "blue")
  #   }
  #   
  #   Z_hat <- ifelse(vb_res$pz > 0.5, 1, 0)
  #   yhat_std <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q_res[,,j])}))
  #   y_hat <- yhat_std + mean(data$Y)
  #   
  #   plot(data$Y, y_hat)
  #   abline(0,1)
  # } else{
  #   for(j in 1:p){
  #     plot(time_points, beta[,,j], ylab = paste0('beta',j), xlab = expression(t), type = 'l', col = 'red')
  #     lines(time_points, beta_hat[,,j], col = "blue")
  #   }
  #   
  #   Z_hat <- ifelse(pz_q > 0.5, 1, 0)
  #   y_hat <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q_res[,,j])}))
  #   plot(data$Y, y_hat)
  #   abline(0,1)
  #   
  # }
  
  return(vb_res)
  
}

