# Oracle simulation 
library(fda)
library(pracma)
library(matrixcalc)
library(GPBayes)
library(lqr)

# source VB and data R files
source('elbo_v2.R')
source('gen_data_oracle.R')
source('VBSOFR_VS.R')

# aux functions
source('std_pred_fun.R')
source('check_convergence.R')

# expectations
source('E_quad_b_Z.R')
source('E_quad_y.R')
 
# oracle sim
source('sim_oracle.R')


K = 10
p = 2
Z = c(1,0)
nt = 100
sigma2 = 0.0025

time_points <- seq(from = 0, to = 1, length= nt)
basis <- create.bspline.basis(rangeval = range(time_points), nbasis = K)
B <- lapply(1:p, function(j){getbasismatrix(time_points, basis, nderiv = 0)})
ids <- split(1:(K*p), rep(1:p, each = K))

# Step 1: Fix/generate lambda2j
lambda2 <- c(1, 150)

# Step 2: Generate taus from prior
set.seed(1234)
tau2 <- c(rexp(K, rate = lambda2[1]/2), rexp(K, rate = lambda2[2]/2))

# Step 3: Generate bs from prior
set.seed(1234)
b1 <- rnorm(K, sd = sqrt(tau2[ids[[1]]]*sigma2))
b2 <- rnorm(K, sd = sqrt(tau2[ids[[2]]]*sigma2))

# Generate curves for coefficients with K b-splines
beta1 <- B[[1]]%*%b1
beta2 <- rep(0, nt)

beta <- array(NA, dim =  c(nt, 1, p))
beta[,,1] <- beta1
beta[,,2] <- beta2

plot(time_points, beta1, type = "l", col = "red")
plot(time_points, beta2, type = "l", col = "red")


initial_values <- list(E_lambda2 = c(1,1), E_eta = rep(0.01, K*p))

res <- sim_oracle(seed = 2024, nsim = 2, initial_values = initial_values, 
                  n = 100, sigma2 = 0.0025, folder = 'TESTE4/n100', 
                  Z = Z, p = p, K = K, nt = 100, beta = beta, 
                  B = B, time_points = time_points, 
                  Niter = 500, convergence_threshold = 0.001)

results <- lapply(1:100, function(i){sim_oracle(seed = 2024, nsim = i, 
                                                shape_lambda_0 = 2, rate_0 = 0.001,
                                                initial_values = initial_values, n = 100, 
                                                sigma2 = 0.0025, folder = 'TESTE4/n100', Z = Z, p = p, K = K, 
                                                nt = 100, beta = beta, B = B, time_points = time_points, 
                                                Niter = 500, convergence_threshold = 0.001)})
save(results, file = 'TESTE4/n100/results.RData')


ids <- split(1:(K*p), rep(1:p, each = K))

#cases = c('Oracle/n25_sigma20.01','Oracle/n100_sigma20.01','Oracle/n300_sigma20.01')

case <- 'TESTE4/n100'

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
  res_pz <- do.call(rbind, lapply(results, `[[`, 7))
  res_lambda2 <- do.call(rbind,lapply(results, `[[`, 8))
  
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
    plot(time_points, beta[,,j], type = 'l', col = 'red', ylab = paste0('beta',j), xlab = "t", lty = 2)
    for(s in 1:nsim){
      lines(time_points, beta_hat_all[s,,j], col = "gray")
    }
    lines(time_points, beta_hat[,,j], col = "blue")
    #lines(time_points, B[[j]]%*%mu0[ids[[j]]], col = "purple")
    lines(time_points, beta[,,j], col = 'red', lty = 2)
    legend("bottomright", c('True', 'Estimated'), col = c("red", "blue"), lty = c(2,1))
  }
  dev.off()
  
  res <- rbind(res, data.frame(model, n_sel_j2))
}  
res

