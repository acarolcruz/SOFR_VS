# Summarize simulated results

# aux function:
rss <- function(i, data, results){
  Y <- data$Y
  W_mat <- data$W_mat
  mu_b_q <- results[[i]][[1]]
  pz_q <- results[[i]][[7]]
  Z_hat <- ifelse(pz_q > 0.5, 1, 0)
  yhat <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q[ids[[j]]])}))
  
  RSS <- sum((Y - yhat)^2/n)
  return(RSS)
}

rss_std <- function(i, data, results){
  Y <- data$Y
  W_mat <- data$W_mat
  mu_b_q <- results[[i]][[1]]
  pz_q <- results[[i]][[7]]
  Z_hat <- ifelse(pz_q > 0.5, 1, 0)
  yhat <- mean(data$Y) + rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q[ids[[j]]])}))
  
  RSS <- sum((Y - yhat)^2/n)
  return(RSS)
}


folder <- 'Simulation Ronaldo VS'
sub_folders <- list.dirs(folder)[-1]

res <- c()
res2 <- c()
res3 <- c()
for(case in sub_folders){
  
  x <- list.dirs(case)
  pattern <- "n([0-9]+)_sigma2([.0-9]+)"
  s <- grep(pattern, x)
  m <- regexec(pattern, x[s])
  m <- regmatches(x[s], m)
  
  # sapply(m, \(x) as.numeric(x[2:3]))
  n <- as.numeric(m[[1]][2])
  sigma2 <- as.numeric(m[[1]][3])
  
  model <- paste0('n = ',n, ' and sigma2 = ', sigma2)
  
  K <- 6
  nt <- 50
  gamma <- 0
  Z <- c(1,1,0,1,0,0)
  p <- 6

  load(paste0(case,"/results.RData"))
  load(paste0(case,"/data_1.RData"))
  
  # results for selection
  res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
  res_pz <- do.call(rbind, lapply(results, `[[`, 7))
  Zhat <- ifelse(res_pz > 0.5, 1,0)
  n_sel_j <- t(as.matrix(colSums(t(sapply(1:100, function(sim){as.numeric(Zhat[sim,] == Z)})))))
  n_sel_j2 <- t(as.matrix(colSums(t(sapply(1:100, function(sim){as.numeric(Zhat[sim,] == 1)})))))
  n_sel <- sum(sapply(1:100, function(sim){sum(Zhat[sim,] == Z) == p}))
  
  # AMSE
  ids <- split(1:(K*p), rep(1:p, each = K))
  amse <- mean(sapply(1:100, function(i){rss(i, data, results)}))
 
  
  res <- rbind(res, data.frame(Case = model, AMSE = round(amse, 4), 'Correct models' = n_sel))
  res2 <- rbind(res2, data.frame(Case = model, n_sel_j))
  res3 <- rbind(res3, data.frame(Case = model, n_sel_j2))
  
  B <- data$B
  beta <- data$beta
  
  time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                       seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
  res_pz <- do.call(rbind, lapply(results, `[[`, 7))
  
  
  beta_hat <- array(NA, c(nt, 1, p))
  for(j in 1:p){
    beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])
  }
  pdf(paste0(case,'/results_VB.pdf'), width = 6, height = 6)
  for(j in 1:p){
    plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", ylim = c(-1.5,1.5))
    lines(time_points[,j], beta_hat[,,j], col = "blue")
    legend("topleft", c('True', 'Estimated'), col = c("red", "blue"), lty = 1)
  }
  dev.off()
  
}  
res
res2 # correctly predicted
res3 # included in the model

# Results for estimated coefficients

B <- data$B
beta <- data$beta

time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))

ids <- split(1:(K*p), rep(1:p, each = K))

res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
res_pz <- do.call(rbind, lapply(results, `[[`, 7))


beta_hat <- array(NA, c(nt, 1, p))
for(j in 1:p){
  beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])
}
pdf(paste0(folder,'/results_VB.pdf'), width = 6, height = 6)
for(j in 1:p){
  plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", ylim = c(-1.5,1.5))
  lines(time_points[,j], beta_hat[,,j], col = "blue")
  legend("topleft", c('True', 'Estimated'), col = c("red", "blue"), lty = 1)
}
dev.off()


