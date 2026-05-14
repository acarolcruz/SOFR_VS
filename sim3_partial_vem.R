# Simulation 3 - SOFR partial VEM
# Implementation
library(GPBayes)
library(lqr)
library(fda)
library(pracma)
library(matrixcalc)

# source R files
source('E_quad_y_mixed.R')
source('std_pred_partial.R')
source('VEM_SOFR_VS_partial.R')
source('elbo_vem_partial.R')
source('gen_data_mixed_vem.R')
source('sim_mixed_vem.R')

p <- 2
q <- 2
K <- 6
nt <- 100
Z <- c(1,0)
u <- c(0,1)
std <- TRUE

# n = 50, sigma2 = 0.1
sigma2 <- 0.1
initial_values <- list(lambda2 = rep(1, p+q), E_inv_sigma2 = 1/sigma2, 
                       E_eta = rep(1, K*p), E_nu_inv = rep(1, q), 
                       pz = c(1,1), pu = c(1,1), E_alpha = c(1,1))


results <- lapply(1:100, function(i){sim_mixed_vem(seed = 1593, nsim = i, 
                                               initial_values = initial_values, 
                                               n = 50, sigma2 = 0.1, 
                                               folder = 'Simulation3_VEM/n50_sigma20.1', 
                                               Z = Z, u = u, p = p, q = q, K = K,
                                               nt = 100, ordem = 4,
                                               Niter = 100, convergence_threshold = 0.01)})
save(results, file = 'Simulation3_VEM/n50_sigma20.1/results.RData')

# n = 50, sigma2 = 0.5
sigma2 <- 0.5
initial_values <- list(lambda2 = rep(1, p+q), E_inv_sigma2 = 1/sigma2, 
                       E_eta = rep(1, K*p), E_nu_inv = rep(1, q), 
                       pz = c(1,1), pu = c(1,1), E_alpha = c(1,1))
results <- lapply(1:100, function(i){sim_mixed_vem(seed = 1234, nsim = i, 
                                                   initial_values = initial_values, 
                                                   n = 50, sigma2 = 0.5, 
                                                   folder = 'Simulation3_VEM/n50_sigma20.5', 
                                                   Z = Z, u = u, p = p, q = q, K = K,
                                                   nt = 100, ordem = 4,
                                                   Niter = 100, convergence_threshold = 0.01)})
save(results, file = 'Simulation3_VEM/n50_sigma20.5/results.RData')


# n = 100, sigma2 = 0.1
sigma2 <- 0.1
initial_values <- list(lambda2 = rep(1, p+q), E_inv_sigma2 = 1/sigma2, 
                       E_eta = rep(1, K*p), E_nu_inv = rep(1, q), 
                       pz = c(1,1), pu = c(1,1), E_alpha = c(1,1))
results <- lapply(1:100, function(i){sim_mixed_vem(seed = 3287, nsim = i, 
                                                   initial_values = initial_values, 
                                                   n = 100, sigma2 = 0.1, 
                                                   folder = 'Simulation3_VEM/n100_sigma20.1', 
                                                   Z = Z, u = u, p = p, q = q, K = K,
                                                   nt = 100, ordem = 4,
                                                   Niter = 100, convergence_threshold = 0.01)})
save(results, file = 'Simulation3_VEM/n100_sigma20.1/results.RData')

# n = 100, sigma2 = 0.5
sigma2 <- 0.5
initial_values <- list(lambda2 = rep(1, p+q), E_inv_sigma2 = 1/sigma2, 
                       E_eta = rep(1, K*p), E_nu_inv = rep(1, q), 
                       pz = c(1,1), pu = c(1,1), E_alpha = c(1,1))
results <- lapply(1:100, function(i){sim_mixed_vem(seed = 4872, nsim = i, 
                                                   initial_values = initial_values, 
                                                   n = 100, sigma2 = 0.5, 
                                                   folder = 'Simulation3_VEM/n100_sigma20.5', 
                                                   Z = Z, u = u, p = p, q = q, K = K,
                                                   nt = 100, ordem = 4,
                                                   Niter = 100, convergence_threshold = 0.01)})
save(results, file = 'Simulation3_VEM/n100_sigma20.5/results.RData')

folder <- 'Simulation3_VEM/n100_sigma20.05' #'TESTE3_MIXED' 'Simulation3_VEM'
load('Simulation3_VEM/n100_sigma20.01/results.RData')

# results for selection
nsim = 100
res_b <- matrix(do.call(rbind, lapply(results, `[[`, 'mu_b')), ncol = K*p, byrow = TRUE)
res_alpha <- matrix(do.call(rbind, lapply(results, `[[`, 'mu_alpha')), ncol = q, byrow = TRUE)
res_lambda2_b <- matrix(do.call(rbind, lapply(results, `[[`, 'lambda2_b')), ncol = p, byrow = FALSE)
res_lambda2_alpha <- matrix(do.call(rbind, lapply(results, `[[`, 'lambda2_alpha')), ncol = q, byrow = FALSE)
res_pz <- do.call(rbind, lapply(results, `[[`, 'pz'))
res_pu <- do.call(rbind, lapply(results, `[[`, 'pu'))


Zhat <- ifelse(res_pz > 0.5, 1, 0)
Uhat <- ifelse(res_pu > 0.5, 1,0)


# how many times the covariates were selected as important
n_sel_z <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == 1)})))))
n_sel_u <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Uhat[sim,] == 1)})))))

colMeans(res_lambda2_alpha)
colMeans(res_lambda2_b)


# true beta and alpha
load(paste0(folder,"/data_", 1, ".RData"))
load(paste0(folder,"/datastd_", 1, ".RData"))

beta <- data$beta
beta_std <- data_std$beta_std
alpha <- data$alpha

time_points <- matrix(rep(seq(from = 0, to = 1, length= nt), p), ncol = p, nrow = nt)
ids <- split(1:(K*p), rep(1:p, each = K))

# Computing beta_j_hat for each predictor and simulated dataset
beta_hat <- array(NA, c(nt, 1, p))
beta_std_hat <-array(NA, c(nt, 1, p))
beta_hat_all <- array(NA, c(nsim, nt, p))
beta_std_hat_all <- array(NA, c(nsim, nt, p))
for(j in 1:p){
  for(s in 1:nsim){
    load(paste0(folder,"/data_", s, ".RData"))
    load(paste0(folder,"/datastd_", s, ".RData"))
    B <- data$B
    sd_t <- data_std$sd_t
    beta <- data$beta
    beta_hat_all[s,,j] <- Zhat[s,j]*((B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j])
    beta_std_hat_all[s,,j] <- Zhat[s,j]*((B[[j]]%*%(res_b[s,ids[[j]]])))
  }
  beta_hat[,,j] <- apply(beta_hat_all[,,j], 2, median)#colMeans(beta_hat_all[,,j])
  #beta_std_hat[,,j] <- #colMeans(beta_std_hat_all[,,j])
}

EMSE <- sapply(1:p, function(j){rowSums(sapply(1:nsim, function(s){(Zhat[s,j]*beta_hat_all[s,,j] - beta[,,j])^2}))/nsim})

boxplot(EMSE[,2])

EMISE <- apply(EMSE, 2, median)

for(j in 1:p){
  plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", lty = 2)
  for(s in 1:nsim){
    lines(time_points[,j], beta_hat_all[s,,j], col = "gray")
  }
  lines(time_points[,j], beta_hat[,,j], col = "blue")
  lines(time_points[,j], beta[,,j], col = 'red', lty = 2)
  #legend("topleft", c('True', 'Estimated curve'), col = c("red", "blue"), lty = c(2,1))
}

for(j in 1:p){
  plot(time_points[,j], data_std$beta_std[,,j], type = 'l', col = 'red', ylab = paste0('beta_std',j), xlab = "t", lty = 2)
  for(s in 1:nsim){
    lines(time_points[,j], beta_std_hat_all[s,,j], col = "gray")
  }
  lines(time_points[,j], beta_std_hat[,,j], col = "blue")
  lines(time_points[,j], data_std$beta_std[,,j], col = 'red', lty = 2)
  #legend("topleft", c('True', 'Mean estimated'), col = c("red", "blue"), lty = c(2,1))
}

alpha_hat <- c()
alpha_hat_all <- matrix(NA, ncol = q, nrow = nsim)
for(l in 1:q){
  for(s in 1:nsim){
    load(paste0(folder,"/datastd_", s, ".RData"))
    sd_s <- data_std$sd_s
    alpha_hat_all[s,l] <- res_alpha[s,l]/sd_s[l]
  }
  alpha_hat[l] <- mean(alpha_hat_all[,l])
}





boxplot(alpha_hat_all[,1], ylim = c(-1,1))
abline(a = alpha[1], b = 0)

boxplot(alpha_hat_all[,2]) # , ylim = c(1.1, 1.2) depend on sigma2 (0.01)
abline(a = alpha[2], b = 0)

alpha_hat

# check fit


var(data$Y)

# Study of how the non-smothness of beta_std affects the estimation of beta
# Computing beta_j_hat for each predictor and simulated dataset
beta_hat <- array(NA, c(nt, 1, p))
beta_std_hat <- array(NA, c(nt, 1, p))
beta_hat_all <- array(NA, c(nsim, nt, p))
beta_std_hat_all <- array(NA, c(nsim, nt, p))
beta <- array(NA, c(nsim, nt, p))
beta_std <- array(NA, c(nsim, nt, p))
for(j in 1:p){
  for(s in 1:nsim){
    load(paste0(folder,"/datastd_", s, ".RData"))
    load(paste0(folder,"/data_", s, ".RData"))
    B <- data$B
    beta_std[s,,j] <- data_std$beta_std[,,j]
    sd_t <- data_std$sd_t
    print(sd_t)
    beta_hat_all[s,,j] <- (B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j]
    beta_std_hat_all[s,,j] <- (B[[j]]%*%(res_b[s,ids[[j]]]))
    beta[s,,j] <- beta_std[s,,j]/sd_t[,,j]
  }
  beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
  beta_std_hat[,,j] <- colMeans(beta_std_hat_all[,,j])
}


for(j in 1:p){
  
  for(s in 1:nsim){
    plot(time_points[,j], beta[s,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", lty = 2)
    lines(time_points[,j], beta_hat_all[s,,j], col = "gray")
    lines(time_points[,j], beta[s,,j], col = 'red', lty = 2)
  }
  lines(time_points[,j], beta_hat[,,j], col = "blue")

}

for(j in 1:p){
  for(s in 1:nsim){
    plot(time_points[,j], beta_std[s,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", lty = 2)
    lines(time_points[,j], beta_std_hat_all[s,,j], col = "gray")
    lines(time_points[,j], beta_std[s,,j], col = 'red')
  }

}

beta0_hat <- c()
MSE <- c()
for(j in 1:p){
  for(s in 1:nsim){
    load(paste0(folder,"/data_", s, ".RData"))
    load(paste0(folder,"/datastd_", s, ".RData"))
    B <- data$B
    sd_t <- data_std$sd_t
    sd_s <- data_std$sd_s
    beta <- data$beta
    alpha <- data$alpha
    beta0_hat[s] <- mean(data$Y) - 
      sum(sapply(1:p, function(j){trapz(time_points[,j], (data_std$Xbar_t[,,j]*(Zhat[s,j]*beta_hat_all[s,,j])))})) -sum(data_std$Xbar_s*(Uhat[s,]*alpha_hat_all[s,]))
    yhat <- beta0_hat[s] + sapply(1:n, function(i){sum(sapply(1:p, function(j){trapz(time_points[,j], (data$Xt[i,,j]*(Zhat[s,j]*beta_hat_all[s,,j])))}))}) + sapply(1:n, function(i){t(data$Xs[i,])%*%(Uhat[s,]*alpha_hat_all[s,])})
    MSE[s] <- sum((data$Y - yhat)**2)/n
  }
}
mean(MSE)
