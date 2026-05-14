# Oracle simulation 
library(fda)
library(pracma)
library(matrixcalc)
library(GPBayes)
library(lqr)

# source VB and data R files
source('elbo_vem.R')
source('gen_data_oracle.R')
source('VBSOFR_VS_VEM.R')
source('sim_oracle_VEM.R')

# aux functions
source('std_pred_fun.R')
source('check_convergence.R')

# expectations
source('E_quad_b_Z.R')
source('E_quad_y.R')


Z <- c(1,0)
p <- 2
K <- 4
nt <- 100
ordem <- 4

folder <- 'Simulation1_VEM'
std <- TRUE
nsim <- 1

initial_values <- list(pz = rep(1, p), lambda2 = c(1,1), E_eta = rep(1, K*p), E_inv_sigma2 = 1/0.01)

# Running multiple simiulations
Z <- c(1,0)
p <- 2
K <- 4

# n 50 sigma 2 0.1
initial_values <- list(pz = rep(1, p), lambda2 = c(1,1), E_eta = rep(1, K*p), E_inv_sigma2 = 1/0.1)
results <- lapply(1:100, function(i){sim_oracle_VEM(seed = 1956, nsim = i, 
                                                initial_values = initial_values, n = 50,
                                                ordem = 4,
                                                sigma2 = 0.1, 
                                                folder = 'Simulation1_VEM/n50_sigma20.1', 
                                                Z = Z, p = p, K = K, nt = 100,
                                                Niter = 500, convergence_threshold = 0.01)})
save(results, file = 'Simulation1_VEM/n50_sigma20.1/results.RData')

# n 50 sigma2 0.5
initial_values <- list(pz = rep(1, p), lambda2 = c(1,1), E_eta = rep(1, K*p), E_inv_sigma2 = 1/0.5)
results <- lapply(1:100, function(i){sim_oracle_VEM(seed = 1458, nsim = i, 
                                                    initial_values = initial_values, n = 50,
                                                    ordem = 4,
                                                    sigma2 = 0.5, 
                                                    folder = 'Simulation1_VEM/n50_sigma20.5', 
                                                    Z = Z, p = p, K = K, nt = 100,
                                                    Niter = 500, convergence_threshold = 0.01)})
save(results, file = 'Simulation1_VEM/n50_sigma20.5/results.RData')

# n 100 sigma2 0.1
initial_values <- list(pz = rep(1, p), lambda2 = c(1,1), E_eta = rep(1, K*p), E_inv_sigma2 = 1/0.1)
results <- lapply(1:100, function(i){sim_oracle_VEM(seed = 2360, nsim = i, 
                                                    initial_values = initial_values, n = 100,
                                                    ordem = 4,
                                                    sigma2 = 0.1, 
                                                    folder = 'Simulation1_VEM/n100_sigma20.1', 
                                                    Z = Z, p = p, K = K, nt = 100,
                                                    Niter = 500, convergence_threshold = 0.01)})
save(results, file = 'Simulation1_VEM/n100_sigma20.1/results.RData')


# n 100 sigma 2 0.5
initial_values <- list(pz = rep(1, p), lambda2 = c(1,1), E_eta = rep(1, K*p), E_inv_sigma2 = 1/0.5)
results <- lapply(1:100, function(i){sim_oracle_VEM(seed = 3500, nsim = i, 
                                                    initial_values = initial_values, n = 100,
                                                    ordem = 4,
                                                    sigma2 = 0.5, 
                                                    folder = 'Simulation1_VEM/n100_sigma20.5', 
                                                    Z = Z, p = p, K = K, nt = 100,
                                                    Niter = 500, convergence_threshold = 0.01)})
save(results, file = 'Simulation1_VEM/n100_sigma20.5/results.RData')


# n 200 sigma2 0.1
initial_values <- list(pz = rep(1, p), lambda2 = c(1,1), E_eta = rep(1, K*p), E_inv_sigma2 = 1/0.1)
results <- lapply(1:100, function(i){sim_oracle_VEM(seed = 1234, nsim = i, 
                                                    initial_values = initial_values, n = 200,
                                                    ordem = 4,
                                                    sigma2 = 0.1, 
                                                    folder = 'Simulation1_VEM/n200_sigma20.1', 
                                                    Z = Z, p = p, K = K, nt = 100,
                                                    Niter = 500, convergence_threshold = 0.01)})
save(results, file = 'Simulation1_VEM/n200_sigma20.1/results.RData')

# n 200 sigma2 0.5
initial_values <- list(pz = rep(1, p), lambda2 = c(1,1), E_eta = rep(1, K*p), E_inv_sigma2 = 1/0.5)
results <- lapply(1:100, function(i){sim_oracle_VEM(seed = 7890, nsim = i, 
                                                    initial_values = initial_values, n = 200,
                                                    ordem = 4,
                                                    sigma2 = 0.5, 
                                                    folder = 'Simulation1_VEM/n200_sigma20.5', 
                                                    Z = Z, p = p, K = K, nt = 100,
                                                    Niter = 500, convergence_threshold = 0.01)})
save(results, file = 'Simulation1_VEM/n200_sigma20.5/results.RData')

K = 4
p = 2
Z = c(1,0)
nt = 100
#sigma2 = 0.01

time_points <- matrix(rep(seq(from = 0, to = 1, length= nt), p), ncol = p, nrow = nt)

ids <- split(1:(K*p), rep(1:p, each = K))

folder <- 'Simulation1_VEM'

res <- c()
res2 <- c()
res3 <- c()
load(paste0(folder,"/results.RData"))


# results for selection
nsim = 100
res_b <- matrix(do.call(rbind, lapply(results, `[[`, "mu_b")), ncol = K*p, byrow = TRUE)
res_pz <- do.call(rbind, lapply(results, `[[`, "pz"))
res_pz <- do.call(rbind, lapply(results, `[[`, "pz"))
Zhat <- ifelse(res_pz > 0.5, 1,0)
n_sel_j <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == Z)})))))
n_sel_j2 <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == 1)})))))
n_sel <- sum(sapply(1:nsim, function(sim){sum(Zhat[sim,] == Z) == p}))

range(res_pz[which(Zhat[,2] == 1),2])
print(summary(res_pz[which(Zhat[,2] == 1),2]))

#mean(res_delta2)/159.1


# AMSE - update function
#amse <- mean(sapply(1:nsim, function(i){rss_std(i, data, results)}))

library(xtable)
options(xtable.include.rownames = FALSE, xtable.booktabs = TRUE, xtable.caption.placement = "top", xtable.sanitize.text.function = function(x){x})


res_elbo <- do.call(rbind, lapply(results, `[[`, 13))
mean(res_elbo)



# Computing beta_j_hat for each predictor and simulated dataset
beta_hat <- array(NA, c(nt, 1, p))
beta_std_hat <-array(NA, c(nt, 1, p))
beta_hat_all <- array(NA, c(nsim, nt, p))
beta_std_hat_all <- array(NA, c(nsim, nt, p))
mse <- c()
for(j in 1:p){
  for(s in 1:nsim){
    load(paste0(folder,"/data_", s, ".RData"))
    load(paste0(folder,"/datastd_", s, ".RData"))
    B <- data$B
    #mse <- c(mse, rss_std(s, data, results))
    sd_t <- data_std$sd_t
    beta <- data$beta
    beta_hat_all[s,,j] <- (B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j]
    beta_std_hat_all[s,,j] <- (B[[j]]%*%(res_b[s,ids[[j]]]))
  }
  beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
  beta_std_hat[,,j] <- colMeans(beta_std_hat_all[,,j])
}

#amse <- mean(mse)

for(j in 1:p){
  plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", lty = 2)
  for(s in 1:nsim){
    lines(time_points[,j], beta_hat_all[s,,j], col = "gray")
  }
  lines(time_points[,j], beta_hat[,,j], col = "blue")
  lines(time_points[,j], beta[,,j], col = 'red', lty = 2)
  legend("topleft", c('True', 'Mean estimated curve'), col = c("red", "blue"), lty = c(2,1))
}

for(j in 1:p){
  plot(time_points[,j], data$beta_std[,,j], type = 'l', col = 'red', ylab = paste0('beta_std',j), xlab = "t", ylim = c(-0.5, 3), lty = 2)
  for(s in 1:nsim){
    lines(time_points[,j], beta_std_hat_all[s,,j], col = "gray")
  }
  lines(time_points[,j], beta_std_hat[,,j], col = "blue")
  lines(time_points[,j], data$beta_std[,,j], col = 'red', lty = 2)
  legend("topleft", c('True', 'Mean estimated'), col = c("red", "blue"), lty = c(2,1))
}


# when I have more cases:






for(case in cases){
  
  x <- list.dirs(case)
  pattern <- "n([0-9]+)_sigma2([.0-9]+)"
  s <- grep(pattern, x)
  m <- regexec(pattern, x[s])
  m <- regmatches(x[s], m)
  
  # sapply(m, \(x) as.numeric(x[2:3]))
  n <- as.numeric(m[[1]][2])
  sigma2 <- as.numeric(m[[1]][3])
  
  model <- paste0('n = ',n, ' and sigma2 = ', sigma2)
  
  #res <- c()
  #res2 <- c()
  #res3 <- c()
  load(paste0(case,"/results.RData"))
  
  
  
  
  # results for selection
  nsim = 100
  res_b <- matrix(do.call(rbind, lapply(results, `[[`, 1)), ncol = K*p, byrow = TRUE)
  res_pz <- do.call(rbind, lapply(results, `[[`, 9))
  res_lambda2 <- do.call(rbind,lapply(results, `[[`, 10))
  res_elbo <- do.call(rbind,lapply(results, `[[`, 13))
  
  Zhat <- ifelse(res_pz > 0.5, 1,0)
  n_sel_j <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == Z)})))))
  n_sel_j2 <- t(as.matrix(colSums(t(sapply(1:nsim, function(sim){as.numeric(Zhat[sim,] == 1)})))))
  n_sel <- sum(sapply(1:nsim, function(sim){sum(Zhat[sim,] == Z) == p}))
  
  
  # AMSE
  #amse <- mean(sapply(1:nsim, function(i){rss_std(i, data, results)}))
  
  library(xtable)
  options(xtable.include.rownames = FALSE, xtable.booktabs = TRUE, xtable.caption.placement = "top", xtable.sanitize.text.function = function(x){x})
  
  # if std == TRUE
  beta_hat <- array(NA, c(nt, 1, p))
  beta_hat_all <- array(NA, c(nsim, nt, p))
  for(j in 1:p){
    for(s in 1:nsim){
      load(paste0(case,"/datastd_", s, ".RData"))
      sd_t <- data_std$sd_t
      beta_hat_all[s,,j] <- (B[[j]]%*%(res_b[s,ids[[j]]]))/sd_t[,,j]
    }
    beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
  }
  
  # if std == FALSE
  # beta_hat <- array(NA, c(nt, 1, p))
  # beta_hat_all <- array(NA, c(nsim, nt, p))
  # for(j in 1:p){
  #   for(s in 1:nsim){
  #     beta_hat_all[,,j] <- t(sapply(1:nsim, function(s){(B[[j]]%*%(res_b[s,ids[[j]]]))}))
  #   }
  #   beta_hat[,,j] <- colMeans(beta_hat_all[,,j])
  # }
  
  
  pdf(paste0(case,'/results_VB.pdf'), width = 6, height = 6)
  for(j in 1:p){
    plot(time_points[,j], beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", lty = 2)
    for(s in 1:nsim){
      lines(time_points[,j], beta_hat_all[s,,j], col = "gray")
    }
    lines(time_points[,j], beta_hat[,,j], col = "blue")
    #lines(time_points, B[[j]]%*%mu0[ids[[j]]], col = "purple")
    lines(time_points[,j], beta[,,j], col = 'red', lty = 2)
    legend("bottomright", c('True', 'Estimated'), col = c("red", "blue"), lty = c(2,1))
  }
  dev.off()
  
  res <- rbind(res, data.frame(model, n_sel_j2))
}  
res

