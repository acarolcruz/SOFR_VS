# Summarize simulated results

# aux function:
rss <- function(i, data, results, std = TRUE){
  W_mat <- data$W_mat
  mu_b_q <- results[[i]][[1]]
  pz_q <- results[[i]][[7]]
  Z_hat <- ifelse(pz_q > 0.5, 1, 0)
  if(std == FALSE){
    yhat <- rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q[ids[[j]]])}))
  } else{
    yhat <- mean(data$Y) + rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q[ids[[j]]])}))
  }
  
  RSS <- sum((data$Y - yhat)^2/n)
  return(RSS)
}

# rss_std <- function(i, data, results){
#   Y <- data$Y
#   W_mat <- data$W_mat
#   mu_b_q <- results[[i]][[1]]
#   pz_q <- results[[i]][[7]]
#   Z_hat <- ifelse(pz_q > 0.5, 1, 0)
#   yhat <- mean(data$Y) + rowSums(sapply(1:p, function(j){Z_hat[j]*(W_mat[,ids[[j]]]%*%mu_b_q[ids[[j]]])}))
#   
#   RSS <- sum((Y - yhat)^2/n)
#   return(RSS)
# }


folder <- 'Simulation SOFR VS STD'
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
  res_time <- do.call(rbind, lapply(results, `[[`, 10))
  res_iter <- do.call(rbind, lapply(results, `[[`, 9))
  
  Zhat <- ifelse(res_pz > 0.5, 1,0)
  n_sel_j <- t(as.matrix(colSums(t(sapply(1:100, function(sim){as.numeric(Zhat[sim,] == Z)})))))
  n_sel_j2 <- t(as.matrix(colSums(t(sapply(1:100, function(sim){as.numeric(Zhat[sim,] == 1)})))))
  n_sel <- sum(sapply(1:100, function(sim){sum(Zhat[sim,] == Z) == p}))
  
  # AMSE
  ids <- split(1:(K*p), rep(1:p, each = K))
  amse <- mean(sapply(1:100, function(i){rss(i, data, results)}))
 
  
  
  res <- rbind(res, data.frame(Case = model, AMSE = round(amse, 4), Mean_time = round(mean(res_time), 4), 
                               Mean_iter = round(mean(res_iter), 4)))
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
  # for non-std case:
  # for(j in 1:p){
  #   beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])
  # }
  
  # Standardized case
  beta_hat_all <- array(NA, c(nsim, nt, p))
  for(j in 1:p){
    #beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])/sd_t[,,j]
    beta_hat_all[,,j] <- t(sapply(1:nsim, function(s){(B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j]}))
    beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
  }
  
  # EMSE
  
  EMSE <- sapply(1:p, function(j){rowSums(sapply(1:nsim, function(s){(beta_hat_all[s,,j] - beta[,,j])^2}))/nsim})
  pdf(paste0(case,'/results_EMSE.pdf'), width = 6, height = 6)
  for(j in 1:p){
    plot(time_points[,j], EMSE[,j], type = 'l', col = 'red', ylab = paste0('EMSE for beta',j), xlab = "t", ylim = c(0, 0.5))
  }
  dev.off()
  
  pdf(paste0(case,'/results_VB.pdf'), width = 6, height = 6)
  # previous version:
  # for(j in 1:p){
  #   plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", ylim = c(-1.5,1.5))
  #   lines(time_points[,j], beta_hat[,,j], col = "blue")
  #   legend("topleft", c('True', 'Estimated'), col = c("red", "blue"), lty = 1)
  # }
  
  # New version with gray curves
  for(j in 1:p){
    plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", ylim = c(-1.5,1.5), lty = 2)
    for(s in 1:nsim){
      lines(time_points[,j], beta_hat_all[s,,j], col = "gray")
    }
    lines(time_points[,j], beta_hat[,,j], col = "blue")
    lines(time_points[,j], beta[,,j], col = 'red', lty = 2)
    legend("topleft", c('True', 'Estimated'), col = c("red", "blue"), lty = c(2,1))
  }
  dev.off()  
}  
res
res2 # correctly predicted
res3 # included in the model

tab_res_proc <- data.frame(n = c("", "25", "" , "", "100", "", "", "300", ""),
                      sigma_2 = rep(c("0.01", "0.05", "0.1"), 3),
                      Mean_time = res$Mean_time[c(4,5,6,1,2,3,7,8,9)],
                      Mean_iter = res$Mean_iter[c(4,5,6,1,2,3,7,8,9)])
colnames(tab_res_proc)[c(3,4)] <- c('Average processing time', 'Average number of iterations')

# Save as LaTex
tab_res <- data.frame(n = c("", "25", "" , "", "100", "", "", "300", ""),
                      sigma_2 = rep(c("0.01", "0.05", "0.1"), 3),
                        AMSE = res$AMSE[c(4,5,6,1,2,3,7,8,9)])

library(xtable)
options(xtable.include.rownames = FALSE, xtable.booktabs = TRUE, xtable.caption.placement = "top", xtable.sanitize.text.function = function(x){x})
print(xtable(tab_res, caption = "Results from the simulated study", label = "" , digits = 4), file = "output_tables_simulation.tex")
print(xtable(tab_res_proc, caption = "Processing time and number of iterations for simulations", label = "" , digits = 4), file = "output_tables_simulation_proc.tex")

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


