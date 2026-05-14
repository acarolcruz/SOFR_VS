sim_mixed_vem <- function(seed, nsim, initial_values, n, sigma2, folder, Z, u, p, q, K, 
                      nt, std = TRUE, ordem,
                      delta1_0 = 0.01, delta2_0 = 0.01,
                      a_u_0 = 0.5, b_u_0 = 0.5, a_z_0 = 0.5, b_z_0 = 0.5,
                      Niter = 500, convergence_threshold = 0.001){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  seed = seed + nsim
  dir.create(folder, showWarnings = FALSE)
  
  #generate data
  data <- gen_data_mixed(seed = seed, nsim = nsim, folder = folder, n = n, nt = nt, 
                         Z = Z, u = u, K = K, p = p, q = q, sigma = sigma2)
  
  # Standardized predictors to fit model
  if(std){
    data_std <- std_pred_partial(data$Xt_smooth, data$Xs, data$Y, data$beta, K, nt, p, q, n, ordem = ordem, time_points = data$time_points)
    save(data_std, file = paste0(folder,"/datastd_", nsim, ".RData"))
    
    Xs_mat <- data_std$Xs_std
    Xt_std <- data_std$Xt_std
    Y_std <- data_std$Y_std
    W_mat <- data_std$W_mat
    beta_std <- data_std$beta_std
    mean_t <- data_std$Xbar_t
    mean_s <- data_std$Xbar_s
    sd_s <- data_std$sd_s
    sd_t <- data_std$sd_t
    
    # plot(time_points, Xt_std[1,,1], type = "l", ylim = c(-2,2));for(i in 2:n){lines(time_points, Xt_std[i,,1], col = "grey")}
    # plot(time_points, Xt_std[1,,2], type = "l", ylim = c(-2,2));for(i in 2:n){lines(time_points, Xt_std[i,,2], col = "grey")}
  }
  vb_res <- VEM_SOFR_VS_partial(initial_values = initial_values, data = data, 
                              data_std = data_std, n = n, K = K, p = p, q = q,
                              delta1_0 = delta1_0, delta2_0 = delta2_0, 
                              a_z_0 = a_z_0, b_z_0 = b_z_0,
                              a_u_0 = a_u_0, b_u_0 = b_u_0,
                              Niter = Niter, convergence_threshold = convergence_threshold,
                              std = std)
  
  
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
