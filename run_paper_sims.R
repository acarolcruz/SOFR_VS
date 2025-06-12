# Oracle simulation 
library(fda)
library(pracma)
library(matrixcalc)
library(GPBayes)
library(lqr)

# source VB and data R files
source('elbo_v2.R')
source('gen_data_vs.R')
source('VBSOFR_VS.R')

# aux functions
source('std_pred_fun.R')
source('check_convergence.R')

# expectations
source('E_quad_b_Z.R')
source('E_quad_y.R')

# oracle sim
source('sim_paper.R')


Z <- c(1,1,0,1,0,0)
p <- 6
K <- 6
nt <- 50
initial_values <- list(E_lambda2 = rep(1, p), E_eta = rep(0.01, K*p))

res <- sim_paper(seed = 2024, nsim = 2, initial_values = initial_values, 
          folder = 'TESTE5', n = 300, sigma2 = 0.05, 
          Z= Z, p = p, K = K, nt = nt)
res$pz

ids <- split(1:(K*p), rep(1:p, each = K))
load(paste0('TESTE5/', 'data_2.RData'))
B = data$B
time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))

j = 4
plot(time_points[,j], (B[[j]]%*%res$mu_b[ids[[j]]])/data$sd_t[,,j], type = "l", col = 'blue', ylim = c(-1.5, 1.5))
lines(time_points[,j], data$beta[,,j], type = "l", col = 'red')


results <- lapply(1:100, function(i){sim_paper(seed = 2024, nsim = i, initial_values = initial_values, 
                                               folder = 'TESTE5', n = 300, sigma2 = 0.05, 
                                               Z= Z, p = p, K = K, nt = nt)})
save(results, file = 'TESTE5/results.RData')


nsim = 100
n = 50
K <- 10
nt <- 50
gamma <- 0
Z <- c(1,1,0,1,0,0)
p <- 6
ids <- split(1:(K*p), rep(1:p, each = K))

res <- c()
res2 <- c()
res3 <- c()
load(paste0('TESTE5',"/results.RData"))
load(paste0('TESTE5',"/data_1.RData"))




# results for selection
res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
res_pz <- do.call(rbind, lapply(results, `[[`, 7))
Zhat <- ifelse(res_pz > 0.5, 1,0)
n_sel_j <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == Z)})))))
n_sel_j2 <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == 1)})))))
n_sel <- sum(sapply(1:nsim, function(sim){sum(Zhat[sim,] == Z) == p}))

W_mat = data$W_mat

# AMSE
amse <- mean(sapply(1:nsim, function(i){rss_std(i, data, results)}))

library(xtable)
options(xtable.include.rownames = FALSE, xtable.booktabs = TRUE, xtable.caption.placement = "top", xtable.sanitize.text.function = function(x){x})


B <- data$B
beta <- data$beta
sd_t <- data$sd_t

time_points <- cbind(seq(0, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-1, 1, length.out = nt), seq(0, pi/3, length.out = nt),
                     seq(-2, 1, length.out = nt), seq(-1, 1, length.out = nt))

ids <- split(1:(K*p), rep(1:p, each = K))

res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
res_pz <- do.call(rbind, lapply(results, `[[`, 7))



beta_hat <- array(NA, c(nt, 1, p))
beta_hat_all <- array(NA, c(nsim, nt, p))
for(j in 1:p){
  #beta_hat[,,j] <- B[[j]]%*%(colMeans(res_b)[ids[[j]]])/sd_t[,,j]
  beta_hat_all[,,j] <- t(sapply(1:nsim, function(s){(B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j]}))
  beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
}
pdf(paste0(case,'/results_VB.pdf'), width = 6, height = 6)
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

