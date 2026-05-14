# sugar dataset
library(fda)
library(pracma)
library(dplyr)
library(matrixcalc)
library(GPBayes)
library(lqr)

# source VEM and data R files
source('elbo_vem.R')
source('VBSOFR_VS_VEM.R')


# aux functions
source('std_pred_fun.R')
source('check_convergence.R')
source('gcv.R')

# expectations
source('E_quad_b_Z.R')
source('E_quad_y.R')
library(R.matlab)

data <- readMat('/Users/carol/Documents/Phd STATISTICS - Western/3rd project/sugar_Process/data.mat')
str(data)

Y <- data$y[,3]
y_new <- (Y - mean(Y))
p <- 7
nt <- 571
n <- length(Y)

time_points <- data$EmAx
time_points <- matrix(rep(time_points, p), ncol = p, nrow = nt)

# Converting X to be an third dimensional array
# Xt are already quite smooth, so just standardize them
Xt <- array(data$X, dim = c(268, 571, 7))

truecol <- "#e41a1c"
bglsscol <- "#ff7f00"#"#4daf4a"
glassocol <- "#377eb8"
gscadcol <- "#984ea3"
gmcpcol <- "#4daf4a"
vbcol <- "#377eb8"
cbcol <- "black"



plot_dim <- c(4,4)

# Include in the paper
pdf('Sugar_raw.pdf', width = 7, height = 5)
par(mar = c(4, 5, 2, .1)) 
plot(time_points[,1], Xt[100,,1], type = "l", ylab = "Intensity", xlab = "Emission spectra", lwd = 1.5)
for(i in 2:7){
  lines(time_points[,7], Xt[100,,i], col = i, lwd = 1.5)
}
legend("topright", title = "Excitation wavelength", paste0(as.character(data$ExAx), 'nm'), bty = 'n', col = (1:7), lty = 1)
dev.off()
# Test different Ks (4, 6, 10, 15)

K <- 6
basis <- create.bspline.basis(rangeval = range(time_points[,1]), nbasis = K)
B <- lapply(1:p, function(j){getbasismatrix(time_points[,1], basis, nderiv = 0)})

sugardata_std <- std_pred_fun(Xt, Y, K = K, nt = nt, beta = NULL,p = p, n = n, ordem = 4, time_points = time_points, sim_yes = FALSE)

plot(time_points[,1], sugardata_std$Xt_std[100,,1], type = "l", ylim = c(-1,1))
for(i in 2:7){
  lines(time_points[,7], sugardata_std$Xt_std[100,,i], col = i)}


# 50 random initializations: 1236
pz_ini <- t(sapply(1:74, function(ini){set.seed(1236 + ini);rbinom(7,1,0.5)}))
pz_ini <- as.matrix(distinct(as.data.frame(pz_ini)))


res <- sapply(1:50, function(ini){
  initial_values <- list(pz = pz_ini[ini,], lambda2 = rep(100,p), E_eta = rep(1, K*p), 
                         E_inv_sigma2 = 1/var(Y))
  
  res_elbo <- VBSOFR_VS_VEM(initial_values = initial_values, data = sugardata_std,
              data_std = sugardata_std, n = n, K = K, p = p, 
              delta1_0 = 0.1, delta2_0 = 0.1,
              a0 = 0.5, b0 = 0.5,
              Niter = 500,
              convergence_threshold = 0.0001,
              std = TRUE)$elbo})

res_gcv <- function(K, Xt, Y, nt, p, n, time_points){
  
  sugardata_std <- std_pred_fun(Xt, Y, K = K, nt = nt, beta = NULL, p = p, 
                                n = n, ordem = 4, time_points = time_points, sim_yes = FALSE)
  
  res_elbo <- sapply(1:50, function(ini){
    initial_values <- list(pz = pz_ini[ini,], lambda2 = rep(100,p), E_eta = rep(1, K*p), 
                           E_inv_sigma2 = 1/var(Y))
    
    res <- VBSOFR_VS_VEM(initial_values = initial_values, data = sugardata_std,
                         data_std = sugardata_std, n = n, K = K, p = p, 
                         delta1_0 = 0.1, delta2_0 = 0.1,
                         a0 = 0.5, b0 = 0.5,
                         Niter = 500,
                         convergence_threshold = 0.001,
                         std = TRUE)$elbo})
  
  optimal_pz <- pz_ini[which.max(res_elbo),]
  
  initial_values <- list(pz = optimal_pz, lambda2 = rep(100,p), E_eta = rep(1, K*p), 
                         E_inv_sigma2 = 1/var(Y))
  
  res_vb_optimal <- VBSOFR_VS_VEM(initial_values = initial_values, data = sugardata_std,
                                  data_std = sugardata_std, n = n, K = K, p = p, 
                                  delta1_0 = 0.1, delta2_0 = 0.1,
                                  a0 = 0.5, b0 = 0.5,
                                  Niter = 500,
                                  convergence_threshold = 0.001,
                                  std = TRUE)
  
  gcv_res <- gcv_sofr(res_vb_optimal, n = n, sugardata_std$Y, sugardata_std$W_mat, K = K, p = p)
  
  return(list(gcv = gcv_res, pz_op = optimal_pz))
}

res_gcv_k <- lapply(c(5,6,10,12), function(k){res_gcv(K = k, Xt = Xt, Y = Y, nt = nt, p = p, n = n, time_points = time_points)})



res_gcv_ks <- do.call(rbind, lapply(res_gcv_k, `[[`, "gcv"))
plot(c(5,6,10,12), res_gcv_ks)

df <- data.frame(k = c(5,6,10,12), elbo = res_gcv_ks) 

elbow <- function(df){
  
  df <- df[order(df$k),]
  ks <- df$k
  elbos <- df$elbo
  
  nks <- length(ks)
  
  # find line that passes through the end points
  b <- (elbos[nks] - elbos[1])/(ks[nks] - ks[1])
  a <- elbos[1] - b*ks[1]
  
  line_y <- a + b*ks
  
  # compute distance between each point to the line
  dist <- abs(elbos - line_y)
  
  case <- which.max(dist)
  best_k <- ks[which.max(dist)]
  
  return(c(case,best_k))
}

elbow(df)



res_gcv_ks 

# [1,] 2.187921
# [2,] 2.176226
# [3,] 2.024173
# [4,] 1.959832

# res_gcv_k[[2]]$pz_op:  0  0  0  1  0  1  1
z_optimal <- res_gcv_k[[2]]$pz_op

K <- 6
initial_values <- list(pz = z_optimal, lambda2 = rep(100, p), E_eta = rep(1, K*p), E_inv_sigma2 = 1/var(Y))
sugardata_std <- std_pred_fun(Xt, Y, K = K, nt = nt, beta = NULL, p = p, 
                              n = n, ordem = 4, time_points = time_points, sim_yes = FALSE)
start <- proc.time()
vb_res_K8 <- VBSOFR_VS_VEM(initial_values = initial_values, data = sugardata_std,
                            data_std = sugardata_std, n = n, K = K, p = p, 
                            delta1_0 = 0.1, delta2_0 = 0.1,
                            a0 = 0.5, b0 = 0.5,
                            Niter = 500,
                            convergence_threshold = 0.001,
                            std = TRUE)
proc.time() - start
#   user  system elapsed 
# 1.140   0.049   1.195 
# 0.590   0.028   0.629 

basis <- create.bspline.basis(rangeval = range(time_points[,1]), nbasis = K)
B <- lapply(1:p, function(j){getbasismatrix(time_points[,1], basis, nderiv = 0)})
ids <- split(1:(K*p), rep(1:p, each = K))
sd_t <- sugardata_std$sd_t
W_mat <- sugardata_std$W_mat

# credible band
LL <- NULL
UL <- NULL
estimates <- NULL
out <- vb_res_K8
b_samples <- MASS::mvrnorm(200, mu = out$mu_b,  Sigma = out$Sigma_b)
z_sample <- matrix(rbinom(p*200, 1, prob = out$pz), ncol = p, byrow = TRUE)

estimates <- b_samples*z_sample[,rep(1:p, each = K)]

beta_curves <- array(NA, dim = c(nt, 200, p))
for(s in 1:200){
  for(j in 1:p){
    beta_curves[,s,j] <- (B[[j]]%*%estimates[s,ids[[j]]])
  }
}  

LL <- list()
UL <- list()
for(j in 1:p){
  LL[[j]] <- apply(beta_curves[,,j],1,function(i)quantile(i,probs = c(0.025,0.975)))[1,]
  UL[[j]] <- apply(beta_curves[,,j],1,function(i)quantile(i,probs = c(0.025,0.975)))[2,]
}

pdf('Sugar_res_VEM.pdf', width = plot_dim[1]+0.5, height = plot_dim[2]+0.5)
par(mar = c(4, 5, 2, .1))  
for(j in 1:7){
  plot(time_points[,j], (B[[j]]%*%out$mu_b[ids[[j]]])/sd_t[,,j], 
       type = "l", col = vbcol, ylab = bquote(hat(beta[.(j)])(t)),
       xlab = "Emission (nm)", ylim = c(-0.04, 0.06), lwd = 2)
  abline(a = 0, b = 0, col = truecol, lty = 2)
  lines(time_points[,j], LL[[j]]/sd_t[,,j], col = cbcol, lty = 3, lwd = 2)
  lines(time_points[,j], UL[[j]]/sd_t[,,j], col = cbcol, lty = 3, lwd = 2)
  #polygon(c(time_points[,j], rev(time_points[,j])), c(LL[[j]]/sd_t[,,j], rev(UL[[j]]/sd_t[,,j])), col = "#00000030", border = NA)
  legend("topright", c('Estimated curve', 'Credible band'), col = c(vbcol,cbcol), lty = c(1,3), bty = 'n', lwd = 2, cex = 0.8)
}  
dev.off()

# Run competitive methods
library(grpreg)
library(MBSGS)

set.seed(1234)
cvfit_glasso <- cv.grpreg(sugardata_std$W_mat, sugardata_std$Y_std, group = rep(1:p, each = K), penalty="grLasso")

set.seed(1234)
cvfit_gscad <- cv.grpreg(sugardata_std$W_mat, sugardata_std$Y_std, group = rep(1:p, each = K), penalty="grSCAD")

set.seed(1234)
cvfit_gmcp <- cv.grpreg(sugardata_std$W_mat, sugardata_std$Y_std, group = rep(1:p, each = K), penalty="grMCP")

b_hat_glasso <- coef(cvfit_glasso)[-1]
b_hat_gscad <- coef(cvfit_gscad)[-1]
b_hat_gmcp <- coef(cvfit_gmcp)[-1]

# Bayesian gLASSO
start <- proc.time()
set.seed(1234)
bglasso <- BGLSS(sugardata_std$Y_std,sugardata_std$W_mat,niter=10000,burnin=5000,group_size=rep(K, p))
proc.time() - start

#user  system elapsed 
# 6.244   0.194   7.151 

CI_bg <- apply(bglasso$coef, 1, function(x){quantile(x, c(0.025,0.975))})

LL_bg <- NULL 
UL_bg <- NULL
for(j in 1:p){
  LL_bg[[j]] <- (B[[j]]%*%CI_bg[1,ids[[j]]])/sd_t[,,j]
  UL_bg[[j]] <- (B[[j]]%*%CI_bg[2,ids[[j]]])/sd_t[,,j]
}


bglasso$pos_median!=0


b_hat_bglasso <- bglasso$pos_median


Zhat_bglasso <- c()
Zhat_glasso <- c()
Zhat_gmcp <- c()
Zhat_gscad <- c()
for(j in 1:p){
  Zhat_bglasso[j] <- ifelse(sum(ifelse(b_hat_bglasso == 0, 0, 1)[ids[[j]]] > 0) >= 1, 1, 0)
  Zhat_glasso[j] <- ifelse(sum(ifelse(b_hat_glasso == 0, 0, 1)[ids[[j]]] > 0) >= 1, 1, 0)
  Zhat_gmcp[j]  <- ifelse(sum(ifelse(b_hat_gmcp == 0, 0, 1)[ids[[j]]] > 0) >= 1, 1, 0)
  Zhat_gscad[j] <- ifelse(sum(ifelse(b_hat_gscad == 0, 0, 1)[ids[[j]]] > 0) >= 1, 1, 0)
}

Z_hat <- ifelse(out$pz > 0.5, 1, 0)

# functional predictors back in the original scale
beta_hat_vb <- array(NA, dim = c(1,nt,p))
beta_hat_bglasso <- array(NA, dim = c(1,nt,p))
beta_hat_glasso <- array(NA, dim = c(1,nt,p))
beta_hat_gmcp <- array(NA, dim = c(1,nt,p))
beta_hat_gscad <- array(NA, dim = c(1,nt,p))
for(j in 1:p){ 
  beta_hat_vb[,,j] <- (B[[j]]%*%out$mu_b[ids[[j]]])/sd_t[,,j] 
  beta_hat_bglasso[,,j] <- (B[[j]]%*%b_hat_bglasso[ids[[j]]])/sd_t[,,j] 
  beta_hat_glasso[,,j] <- (B[[j]]%*%b_hat_glasso[ids[[j]]])/sd_t[,,j] 
  beta_hat_gmcp[,,j] <- (B[[j]]%*%b_hat_gmcp[ids[[j]]])/sd_t[,,j] 
  beta_hat_gscad[,,j] <- (B[[j]]%*%b_hat_gscad[ids[[j]]])/sd_t[,,j]
}

pdf('Sugar_res_BGLSS.pdf', width = plot_dim[1]+0.5, height = plot_dim[2]+0.5)
par(mar = c(4, 5, 2, .1))  
for(j in 1:7){
  plot(time_points[,j], (B[[j]]%*%b_hat_bglasso[ids[[j]]])/sd_t[,,j], 
       type = "l", col = bglsscol, ylab = bquote(hat(beta[.(j)])(t)),
       xlab = "Emission (nm)", ylim = c(-0.04, 0.06), lwd = 2)
  abline(a = 0, b = 0, col = truecol, lty = 2)
  lines(time_points[,j], LL_bg[[j]], col = cbcol, lty = 3, lwd = 2)
  lines(time_points[,j], UL_bg[[j]], col = cbcol, lty = 3, lwd = 2)
  #polygon(c(time_points[,j], rev(time_points[,j])), c(LL[[j]]/sd_t[,,j], rev(UL[[j]]/sd_t[,,j])), col = "#00000030", border = NA)
  legend("topright", c('Estimated curve', 'Credible band'), col = c(bglsscol,cbcol), lty = c(1,3), bty = 'n', lwd = 2, cex = 0.8)
}  
dev.off()

mse <- function(sugardata_std, time_points, Z_hat, b_hat, n, p, K){ 
  #beta0_hat <- mean(sugardata_std$Y) - 
  #sum(sapply(1:p, function(j){trapz(time_points[,j], (sugardata_std$Xbar_t[,,j]*(Z_hat[j]*beta_hat[,,j])))}))
  
  yhat <- mean(sugardata_std$Y) + rowSums(sapply(1:p, function(j){Z_hat[j]*(sugardata_std$W_mat[,ids[[j]]]%*%b_hat[ids[[j]]])}))
  
  MSE <- sum((sugardata_std$Y - yhat)**2)/n
  
  ratio <- sum((sugardata_std$Y - yhat)^2)/sum((sugardata_std$Y- mean(sugardata_std$Y))^2) # same as Gertheidss 2013
  
  R2 <- 1 - (n - 1)*sum((sugardata_std$Y - yhat)^2)/((n - sum(Z_hat)*K)*(sum((sugardata_std$Y - mean(sugardata_std$Y))^2)))
  
  plot(sugardata_std$Y, yhat)
  abline(0,1)
  
  return(c(MSE, ratio, R2))
}

mse(sugardata_std, time_points, Z_hat, out$mu_b, n, p, K)
mse(sugardata_std, time_points, Zhat_bglasso, b_hat_bglasso, n, p, K)
mse(sugardata_std, time_points, Zhat_glasso, b_hat_glasso, n, p, K)
mse(sugardata_std, time_points, Zhat_gmcp, b_hat_gmcp, n, p, K)
mse(sugardata_std, time_points, Zhat_gscad, b_hat_gscad, n, p, K)

pdf(file = "Comparative_results_sugar_vem.pdf", width = plot_dim[1]+0.5, height = plot_dim[2]+0.5)
par(mar = c(4, 5, 2, .1)) 
for(j in 1:p){
  plot(time_points[,j], B[[j]]%*%b_hat_gscad[ids[[j]]]/sd_t[,,j], type = "l", col = gscadcol, lty = 2, lwd = 1.5, ylab = bquote(hat(beta[.(j)])(t)), xlab = "Emission (nm)", ylim = c(-0.04, 0.06))
  abline(a = 0, b = 0, col = truecol, lty = 2)
  lines(time_points[,j], B[[j]]%*%b_hat_gmcp[ids[[j]]]/sd_t[,,j], col = gmcpcol, lty = 3, lwd = 2)
  lines(time_points[,j], B[[j]]%*%b_hat_glasso[ids[[j]]]/sd_t[,,j], col = glassocol, lty = 4, lwd = 2)
  lines(time_points[,j], (B[[j]]%*%b_hat_bglasso[ids[[j]]])/sd_t[,,j], col = bglsscol, lwd = 2, lty = 2)
  legend("topright", c('grLASSO', 'grMCP', 'grSCAD', 'BGLSS'), col = c(glassocol,gmcpcol,gscadcol,bglsscol), lty = c(2,3,4,5,2), lwd = 2, cex = 0.8, bty = 'n')
}
dev.off()
save(out, file = "sugar_VEM.rdata")
