# sim_oracle not to be in the R Package
sim_oracle <- function(seed, nsim, initial_values, n, sigma2, folder, Z, p, K, 
                       nt, beta, B, time_points, std = TRUE, 
                       delta1_0 = 0.01, delta2_0 = 0.01, a0 = 0.5, b0 = 0.5,
                       shape_lambda_0 = 2, rate_0 = 0.001,
                       Niter = 500, convergence_threshold = 0.001){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  seed = seed + nsim
  dir.create(folder, showWarnings = FALSE)
  
  #generate data
  # Generate Xs
  data <- gen_data_oracle(seed = seed, nsim = nsim, folder = folder, n = n, nt = nt, 
                          B = B, time_points = time_points, Z = Z, K = K, p = p)
  
  # Standardized predictors to fit model
  if(std){
    data_std <- std_pred_fun(data$Xt_smooth, data$Y, beta, K, nt, p, n)
    save(data_std, file = paste0(folder,"/datastd_", nsim, ".RData"))
    
    Xt_std <- data_std$Xt_std
    Y_std <- data_std$Y_std
    W_mat <- data_std$W_mat
    beta_std <- data_std$beta_std
    mean_t <- data_std$Xbar_t
    sd_t <- data_std$sd_t
    
    # plot(time_points, Xt_std[1,,1], type = "l", ylim = c(-2,2));for(i in 2:n){lines(time_points, Xt_std[i,,1], col = "grey")}
    # plot(time_points, Xt_std[1,,2], type = "l", ylim = c(-2,2));for(i in 2:n){lines(time_points, Xt_std[i,,2], col = "grey")}
  }

  vb_res <- VBSOFR_VS(initial_values = initial_values, data = data, 
                      data_std = data_std, n = n, K = K, p = p, delta1_0 = 0.0001, delta2_0 = 0.0001, a0 = 0.1, b0 = 0.1)
  
  
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

