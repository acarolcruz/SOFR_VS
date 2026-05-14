# Save results for gscad, gmcp, glasso, Bayesian group lasso in folder

library(grpreg)
library(MBSGS)

rm(list = ls())

folder <- 'Simulation2_VEM/n400_sigma20.05'
Z <- c(1,0,1,0)
p <- 4
K <- 7
nt <- 81
ids <- split(1:(K*p), rep(1:p, each = K))
time_points <- matrix(rep(seq(0, 1, length.out = nt), p), ncol = p, nrow = nt)

nsim = 100

b_hat_glasso_all <- matrix(NA, nrow = nsim, ncol = K*p)
b_hat_gscad_all <- matrix(NA, nrow = nsim, ncol = K*p)
b_hat_gmcp_all <- matrix(NA, nrow = nsim, ncol = K*p)
b_hat_bglasso_all <- matrix(NA, nrow = nsim, ncol = K*p)

runtime_glasso <- rep(NA, nsim)
runtime_gscad <- rep(NA, nsim)
runtime_gmcp <- rep(NA, nsim)
runtime_bglasso <- rep(NA, nsim)

for(sim in 1:nsim){
  load(paste0(folder, "/data_", sim, ".RData"))
  
  start <- proc.time()
  cvfit_glasso <- cv.grpreg(X = data$W_mat, y = data$Y_std, group = rep(1:p, each = K), penalty="grLasso")
  runtime_glasso[sim] <- proc.time() - start
  
  start <- proc.time()
  cvfit_gscad <- cv.grpreg(data$W_mat, data$Y_std, group = rep(1:p, each = K), penalty="grSCAD")
  runtime_gscad[sim] <- proc.time() - start
  
  start <- proc.time()
  cvfit_gmcp <- cv.grpreg(data$W_mat, data$Y_std, group = rep(1:p, each = K), penalty="grMCP")
  runtime_gmcp[sim] <- proc.time() - start
  
  b_hat_glasso_all[sim,] <- as.numeric(coef(cvfit_glasso)[-1])
  b_hat_gscad_all[sim,] <- as.numeric(coef(cvfit_gscad)[-1])
  b_hat_gmcp_all[sim,] <- as.numeric(coef(cvfit_gmcp)[-1])
  
  # Bayesian gLASSO
  start <- proc.time()
  bglasso <- BGLSS(data$Y_std, data$W_mat, 
                   niter=10000, burnin=5000, group_size=rep(K, p),
                   num_update = 100, niter.update = 100)
  runtime_bglasso[sim] <- proc.time() - start
  
  if(sim == 97){
    CI_bg <- apply(bglasso$coef, 1, function(x){quantile(x, c(0.025,0.975))})
  }
  
  b_hat_bglasso_all[sim,] <- bglasso$pos_median
}  

# saving results 

save(CI_bg, file = paste0(folder,'/CI_bg_sim97.RData'))
save(b_hat_glasso_all, file = paste0(folder,'/b_hat_glasso_all.RData'))
save(runtime_glasso, file = paste0(folder,'/runtime_glasso.RData'))

save(b_hat_gscad_all, file = paste0(folder,'/b_hat_gscad_all.RData'))
save(runtime_gscad, file = paste0(folder,'/runtime_gscad.RData'))

save(b_hat_gmcp_all, file = paste0(folder,'/b_hat_gmcp_all.RData'))
save(runtime_gmcp, file = paste0(folder,'/runtime_gmcp.RData'))

save(b_hat_bglasso_all, file = paste0(folder,'/b_hat_bglasso_all.RData'))
save(runtime_bglasso, file = paste0(folder,'/runtime_bglasso.RData'))

# Create plots
load(paste0(folder,'/b_hat_glasso_all.RData'))
load(paste0(folder,'/b_hat_gscad_all.RData'))
load(paste0(folder,'/b_hat_gmcp_all.RData'))
load(paste0(folder,'/b_hat_bglasso_all.RData'))

# beta_hat for all competitive models
beta_hat_glasso <- array(NA, c(nt, 1, p))
beta_hat_glasso_all <- array(NA, c(nsim, nt, p))
beta_hat_gscad <- array(NA, c(nt, 1, p))
beta_hat_gscad_all <- array(NA, c(nsim, nt, p))
beta_hat_gmcp <- array(NA, c(nt, 1, p))
beta_hat_gmcp_all <- array(NA, c(nsim, nt, p))
beta_hat_bglasso <- array(NA, c(nt, 1, p))
beta_hat_bglasso_all <- array(NA, c(nsim, nt, p))

for(j in 1:p){
  for(s in 1:nsim){
    load(paste0(folder,"/data_", s, ".RData"))
    B <- data$B
    sd_t <- data$sd_t
    beta <- data$beta
    beta_hat_glasso_all[s,,j] <- (B[[j]]%*%(b_hat_glasso_all[s,ids[[j]]]))/sd_t[,,j]
    beta_hat_gscad_all[s,,j] <- (B[[j]]%*%(b_hat_gscad_all[s,ids[[j]]]))/sd_t[,,j]
    beta_hat_gmcp_all[s,,j] <- (B[[j]]%*%(b_hat_gmcp_all[s,ids[[j]]]))/sd_t[,,j]
    beta_hat_bglasso_all[s,,j] <- (B[[j]]%*%(b_hat_bglasso_all[s,ids[[j]]]))/sd_t[,,j]
  }
  beta_hat_glasso[,,j] <- colMeans(beta_hat_glasso_all[,,j])
  beta_hat_gscad[,,j] <- colMeans(beta_hat_gscad_all[,,j])
  beta_hat_gmcp[,,j] <- colMeans(beta_hat_gmcp_all[,,j])
  beta_hat_bglasso[,,j] <- colMeans(beta_hat_bglasso_all[,,j])
}





pdf(file = paste0(folder, "/Comparative_results_average.pdf"))
for(j in 1:p){
  plot(time_points[,j], beta[,,j], col = "red", type = "l", ylim = c(-2,2),ylab = paste0('beta',j), xlab = "t", lwd = 1.5)
  #plot(time_points[,j], beta_hat[,,j], type = "l", ylab = paste0('beta',j), col = "blue", ylim = c(-5, 10))
  lines(time_points[,j], beta_hat_glasso[,,j], col = "purple", lty = 2, lwd = 1.5, type = "l", ylim = c(-5, 10))
  lines(time_points[,j], beta_hat_gscad[,,j], col = "orange", lty = 3, lwd = 1.5)
  lines(time_points[,j], beta_hat_gmcp[,,j], col = "blue", lty = 4, lwd = 1.5)
  lines(time_points[,j], beta_hat_bglasso[,,j],  col = "green", lty =  5, lwd = 1.5)
  legend("topleft", c('group Lasso', 'group SCAD', 'group MCP', 'Bayesian group Lasso'), col = c("purple","orange", "blue", 'green'), lty = c(2,2,3,4,5), lwd = 1.5, cex = 0.7)
}
dev.off()

nsim = 1
for(sim in 1:nsim){
  load(paste0(folder, "/data_", sim, ".RData"))
  
  cvfit_glasso <- cv.grpreg(data$W_mat, data$Y_std, group = 1:(K*p), penalty="grLasso", returnX = TRUE)
  
  #cvfit_gscad <- cv.grpreg(data$W_mat, data$Y_std, group = as.factor(rep(1:p, each = K)), penalty="grSCAD", returnX = TRUE)
  
  #cvfit_gmcp <- cv.grpreg(data$W_mat, data$Y_std, group = as.factor(rep(1:p, each = K)), penalty="grMCP", returnX = TRUE)
  
  
  b_hat_glasso_all[sim,] <- as.numeric(coef(cvfit_glasso)[-1]/cvfit_glasso$fit$XG$scale)
  # b_hat_gscad_all[sim,] <- as.numeric(coef(cvfit_gscad)[-1]/cvfit_gscad$fit$XG$scale)
  #b_hat_gmcp_all[sim,] <- as.numeric(coef(cvfit_gmcp)[-1]/cvfit_gmcp$fit$XG$scale)
  
  
  # Bayesian gLASSO
  
  # bglasso <- BGLSS(data$Y_std, data$W_mat, 
  #                  niter=10000, burnin=5000, group_size=rep(K, p),
  #                   num_update = 100, niter.update = 100)
  
  #b_hat_bglasso_all[sim,] <- bglasso$pos_median
}  

for(j in 1:p){
  plot(time_points[,j], beta_hat_glasso[,,j], col = "grey", lty = 2, lwd = 1.5, type = "l", ylim = c(-5, 10))
  lines(time_points[,j], beta[,,j], col = "red")
}

beta_hat_glasso <- array(NA, c(nt, 1, p))
beta_hat_glasso_all <- array(NA, c(nsim, nt, p))

for(j in 1:p){
  load(paste0(folder,"/data_", s, ".RData"))
  B <- data$B
  sd_t <- data$sd_t
  beta <- data$beta
  beta_hat_glasso_all[s,,j] <- (B[[j]]%*%(b_hat_glasso_all[s,ids[[j]]]))/sd_t[,,j]
  beta_hat_glasso[,,j] <- beta_hat_glasso_all[1,,j]
}

