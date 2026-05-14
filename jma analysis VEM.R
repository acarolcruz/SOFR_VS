# JMA dataset
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


load('/Users/carol/Documents/Phd STATISTICS - Western/Research/SOFR-VS/Data/table-1.RData')
X1 <- y
rm(y)

load('/Users/carol/Documents/Phd STATISTICS - Western/Research/SOFR-VS/Data/table-2.RData')
X2 <- y
rm(y)

load('/Users/carol/Documents/Phd STATISTICS - Western/Research/SOFR-VS/Data/table-3.RData')
X3 <- y
rm(y)

load('/Users/carol/Documents/Phd STATISTICS - Western/Research/SOFR-VS/Data/table-6.RData')
X4 <- y
rm(y)

load('/Users/carol/Documents/Phd STATISTICS - Western/Research/SOFR-VS/Data/table-7.RData')
X5 <- y
rm(y)

load('/Users/carol/Documents/Phd STATISTICS - Western/Research/SOFR-VS/Data/table-12.RData')
X6 <- y
rm(y)

load('/Users/carol/Documents/Phd STATISTICS - Western/Research/SOFR-VS/Data/table-13.RData')

X <- as.matrix(cbind(X1, X2, X3, X4, X5, X6), ncol = 12*6, nrow = 157)
which(complete.cases(X))

which(complete.cases(y$y))

complete_ids <- intersect(which(complete.cases(X)), which(complete.cases(y$y)))


Xt <- array(X[complete_ids,], dim = c(length(complete_ids), 12, 6))
y_complete <- y$y[complete_ids]



p <- 6
K <- 4
nt <- 12
n <- length(y_complete)
ordem <- 4

truecol <- "#e41a1c"
bglsscol <- "#ff7f00" #"#4daf4a"
glassocol <- "#377eb8"
gscadcol <- "#984ea3"
gmcpcol <- "#4daf4a"
vbcol <- "#377eb8"
cbcol <- "black"

plot_dim <- c(4,4)

time_points <- matrix(rep(1:12, p), ncol = p, nrow = nt)

# Smooth X(t)
X_smooth <- array(NA, dim = c(n, nt, p))
A <- array(NA, dim = c(n, K, p))
J <- array(NA, dim = c(K, K, p))
B <- list()

#finer grid
nt_finer <- length(seq(1,12,0.1))
X_smooth_finer <- array(NA, dim = c(n, nt_finer, p))
time_finer <- matrix(rep(seq(1,12,0.1), p), ncol = p, nrow = nt_finer)
for(i in 1:n){
  for(j in 1:p){
    basis <- create.bspline.basis(rangeval = range(time_points[,j]), norder = ordem, nbasis = K)
    B[[j]] <- getbasismatrix(time_points[,j], basis, nderiv = 0)
    res <- smooth.basis(time_points[,j], Xt[i,,j], basis)#Xt
    A[i,,j] <- res$fd$coefs
    X_smooth[i,,j] <- eval.fd(time_points[,j], res$fd)
    X_smooth_finer[i,,j] <- eval.fd(time_finer[,j], res$fd) # for plot
    J[,,j] <- inprod(basis, basis)
    #plotfit.fd(Xt[i,,j], time_points[,j], res$fd)
  }
}

sd_t_finer <- std_pred_fun(X_smooth_finer, y_complete, K = K, nt = length(time_finer[,1]), beta = NULL, p = p, n = n, ordem = 4, time_points = time_finer, sim_yes = FALSE)$sd_t

# 1 - Temperature ylim = c(-10, 30)
# 2 - Maximum temperature ylim = c(-10, 30)
# 3 - Minimum temperature ylim = c(-10, 30)
# 4 - Pressure  ylim = c(900, 1020)
# 5 - Humidity ylim = c(40, 100)
# 6 - Daylight ylim = c(10, 300) Average of the Total amount of hours of daylight from 1991 to 2020

# to include in the paper - finer grid
ylim_grid = list(c(-10, 30),c(-10, 30),c(-10, 30),c(900,1020), c(40, 100),c(10, 300))
var_names <- c('Average monthly temperature', 'Average monthly maximum temperature', "Average monthly minimum temperature", "Average monthly pressure", 'Average monthly humidity', 'Average monthly daylight duration in hours')
pdf('JMA_pred_curves.pdf', width = 12, height = 8)
for(j in 1:p){
  plot(time_finer[,j], X_smooth_finer[83,,j], type = "l", ylim = ylim_grid[[j]], ylab = var_names[j], xlab = "Month")
  for(i in 1:82){
    lines(time_finer[,j], X_smooth_finer[i,,j], type = "l")}
}
dev.off()

jma_std <- std_pred_fun(Xt, y_complete, K = K, nt = nt, beta = NULL, p = p, n = n, ordem = 4, time_points = time_points, sim_yes = FALSE)
#c(0,0,0,1,1,0,1)

plot(time_points[,1], jma_std$Xt_std[83,,1], type = "l")
for(j in 2:6){
  lines(time_points[,j], jma_std$Xt_std[83,,j], col = j)}


# random initializations:
pz_ini <- t(sapply(1:105, function(ini){set.seed(1234 + ini);rbinom(6,1,0.5)}))
pz_ini <- as.matrix(distinct(as.data.frame(pz_ini)))


res_elbo <- sapply(1:50, function(ini){
  initial_values <- list(pz = as.numeric(pz_ini[ini,]), lambda2 = rep(1, p), E_eta = rep(1, K*p), E_inv_sigma2 = 1/400)
  
  elbos <- VBSOFR_VS_VEM(initial_values = initial_values, data = jma_std,
                         data_std = jma_std, n = n, K = K, p = p, 
                         delta1_0 = 0.1, delta2_0 = 0.1,
                         a0 = 0.5, b0 = 0.5,
                         Niter = 500,
                         convergence_threshold = 0.001,
                         std = TRUE)$elbo}
  )

# sigma = range/4
res_gcv <- function(K, Xt, Y, nt, p, n, time_points){
  
  data_std <- std_pred_fun(Xt, Y, K = K, nt = nt, beta = NULL, p = p, 
                           n = n, ordem = 4, time_points = time_points, sim_yes = FALSE)
  
  res_elbo <- sapply(1:50, function(ini){
    initial_values <- list(pz = pz_ini[ini,], lambda2 = rep(100, p), E_eta = rep(1, K*p), 
                           E_inv_sigma2 = 1/100)
    
    res <- VBSOFR_VS_VEM(initial_values = initial_values, data = data_std,
                         data_std = data_std, n = n, K = K, p = p, 
                         delta1_0 = 0.1, delta2_0 = 0.1,
                         a0 = 0.5, b0 = 0.5,
                         Niter = 500,
                         convergence_threshold = 0.001,
                         std = TRUE)$elbo})
  
  optimal_pz <- pz_ini[which.max(res_elbo),]
  
  initial_values <- list(pz = optimal_pz, lambda2 = rep(100,p), E_eta = rep(1, K*p), 
                         E_inv_sigma2 = 1/100)
  
  res_vb_optimal <- VBSOFR_VS_VEM(initial_values = initial_values, data = data_std,
                                 data_std = data_std, n = n, K = K, p = p, 
                                 delta1_0 = 0.1, delta2_0 = 0.1,
                                 a0 = 0.5, b0 = 0.5,
                                 Niter = 500,
                                 convergence_threshold = 0.001,
                                 std = TRUE)
  
  gcv_res <- gcv_sofr(res_vb_optimal, n = n, data_std$Y, data_std$W_mat, K = K, p = p)
  
  return(list(gcv = gcv_res, pz_op = optimal_pz))
}

res_gcv_k <- lapply(c(6,8,10,12), function(k){res_gcv(K = k, Xt = Xt, Y = y_complete, nt = nt, p = p, n = n, time_points = time_points)})

#c(5,6,7,8,10,12) = selects only daylight = 6
res_gcv_ks <- do.call(rbind, lapply(res_gcv_k, `[[`, "gcv"))
plot(c(6,8,10,12), res_gcv_ks)
 
res_gcv_ks 

# [1,] 98576.43
# [2,] 83222.82
# [3,] 77202.35
# [4,] 77132.81

df <- data.frame(k = c(6,8,10,12), elbo = res_gcv_ks) 

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

res_gcv_k[[2]]$pz_op #1 0  0  0  1  1 

pz_optimal <- c(1,0,0,0,1,1)

K <- 6
initial_values <- list(pz = pz_optimal, lambda2 = rep(1, p), E_eta = rep(1, K*p), E_inv_sigma2 = 1/400)
jmadata_std <- std_pred_fun(Xt, y_complete, K = K, nt = nt, beta = NULL, p = p, 
                            n = n, ordem = 4, time_points = time_points, sim_yes = FALSE)
start <- proc.time()
jma_K8 <- VBSOFR_VS_VEM(initial_values = initial_values, data = NULL,
                        data_std = jmadata_std, n = n, K = K, p = p, 
                        delta1_0 = 0.1, delta2_0 = 0.1,
                        a0 = 0.5, b0 = 0.5,
                        Niter = 500,
                        convergence_threshold = 0.001,
                        std = TRUE)
proc.time() - start

#user  system elapsed 
#9.704   0.042  10.086 
# user  system elapsed 
# 0.625   0.010   0.635 

basis <- create.bspline.basis(rangeval = range(time_points[,1]), nbasis = K)
B <- lapply(1:p, function(j){getbasismatrix(time_points[,1], basis, nderiv = 0)})
ids <- split(1:(K*p), rep(1:p, each = K))
sd_t <- jmadata_std$sd_t

# credible band
LL <- NULL
UL <- NULL
estimates <- NULL
out <- jma_K8
b_samples <- MASS::mvrnorm(200, mu = out$mu_b,  Sigma = out$Sigma_b)
z_sample <- matrix(rbinom(p*200, 1, prob = out$pz), ncol = p, byrow = TRUE)

estimates <- b_samples*z_sample[,rep(1:p, each = K)]

# finer grid
basis_finer <- create.bspline.basis(rangeval = range(time_finer[,1]), nbasis = K)
B_finer <- lapply(1:p, function(j){getbasismatrix(time_finer[,1], basis_finer, nderiv = 0)})
sd_t_finer <- std_pred_fun(X_smooth_finer, y_complete, K = K, nt = length(time_finer[,1]), beta = NULL, p = p, n = n, ordem = 4, time_points = time_finer, sim_yes = FALSE)$sd_t

nr_finer <- length(time_finer[,1])
beta_curves <- array(NA, dim = c(nt_finer, 200, p))
for(s in 1:200){
  for(j in 1:p){
    beta_curves[,s,j] <- (B_finer[[j]]%*%estimates[s,ids[[j]]])
  }
}  

LL <- list()
UL <- list()
for(j in 1:p){
  LL[[j]] <- apply(beta_curves[,,j],1,function(i)quantile(i,probs = c(0.025,0.975)))[1,]
  UL[[j]] <- apply(beta_curves[,,j],1,function(i)quantile(i,probs = c(0.025,0.975)))[2,]
}


ylim_values <- list(c(-30,40), c(-5,5), c(-5,5), c(-5,5), c(-8,15), c(-20,10))
pdf('JMA_res_VEM.pdf', width = plot_dim[1], height = plot_dim[2])
par(mar = c(4, 5, 2, .1))    
for(j in 1:6){
  plot(time_finer[,j], (B_finer[[j]]%*%out$mu_b[ids[[j]]])/sd_t_finer[,,j], 
       type = "l", col = vbcol, ylab = bquote(hat(beta[.(j)])(t)),
       xlab = "Month", ylim = ylim_values[[j]], lwd = 2)
  abline(a = 0, b = 0, col = truecol, lty = 2)
  lines(time_finer[,j], LL[[j]]/sd_t_finer[,,j], col = cbcol, lty = 3, lwd = 2)
  lines(time_finer[,j], UL[[j]]/sd_t_finer[,,j], col = cbcol, lty = 3, lwd = 2)
  #polygon(c(time_points[,j], rev(time_points[,j])), c(LL[[j]]/sd_t[,,j], rev(UL[[j]]/sd_t[,,j])), col = "#00000030", border = NA)
  legend("topright", c('Estimated curve', 'Credible band'), col = c(vbcol,cbcol), lty = c(1,3), lwd = 2, bty = "n", cex = 0.8)
}  
dev.off()

# Run competitive methods
library(grpreg)
library(MBSGS)

set.seed(1234)
cvfit_glasso <- cv.grpreg(jmadata_std$W_mat, jmadata_std$Y_std, group = rep(1:p, each = K), penalty="grLasso")

set.seed(1234)
cvfit_gscad <- cv.grpreg(jmadata_std$W_mat, jmadata_std$Y_std, group = rep(1:p, each = K), penalty="grSCAD")

set.seed(1234)
cvfit_gmcp <- cv.grpreg(jmadata_std$W_mat, jmadata_std$Y_std, group = rep(1:p, each = K), penalty="grMCP")

b_hat_glasso <- coef(cvfit_glasso)[-1]
b_hat_gscad <- coef(cvfit_gscad)[-1]
b_hat_gmcp <- coef(cvfit_gmcp)[-1]

# Bayesian gLASSO

start <- proc.time()
set.seed(1234)
bglasso <- BGLSS(jmadata_std$Y_std, jmadata_std$W_mat, niter=10000,burnin=5000,group_size=rep(K, p))
proc.time() - start

#    user  system elapsed 
# 5.096   0.118   5.564 
bglasso$pos_median!=0

b_hat_bglasso <- bglasso$pos_median

CI_bg <- apply(bglasso$coef, 1, function(x){quantile(x, c(0.025,0.975))})

LL_bg <- NULL 
UL_bg <- NULL
for(j in 1:p){
  LL_bg[[j]] <- (B_finer[[j]]%*%CI_bg[1,ids[[j]]])/sd_t_finer[,,j]
  UL_bg[[j]] <- (B_finer[[j]]%*%CI_bg[2,ids[[j]]])/sd_t_finer[,,j]
}


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

mse <- function(data_std, time_points, Z_hat, b_hat, n, p, K){ 
  #beta0_hat <- mean(sugardata_std$Y) - 
  #sum(sapply(1:p, function(j){trapz(time_points[,j], (sugardata_std$Xbar_t[,,j]*(Z_hat[j]*beta_hat[,,j])))}))
  
  yhat <- mean(data_std$Y) + rowSums(sapply(1:p, function(j){(data_std$W_mat[,ids[[j]]]%*%b_hat[ids[[j]]])}))
  
  yhat <- mean(data_std$Y) + rowSums(sapply(1:p, function(j){Z_hat[j]*(data_std$W_mat[,ids[[j]]]%*%b_hat[ids[[j]]])}))
  
  
  MSE <- sum((data_std$Y - yhat)**2)/n
  
  ratio <- sum((data_std$Y - yhat)^2)/sum((data_std$Y- mean(data_std$Y))^2) # same as Gertheidss 2013
  
  R2 <- 1 - (n - 1)*sum((data_std$Y - yhat)^2)/((n - sum(Z_hat)*K)*(sum((data_std$Y - mean(data_std$Y))^2)))
  
  plot(data_std$Y, yhat)
  abline(0,1)
  
  return(c(MSE, ratio, R2))
}

round(mse(jmadata_std, time_points, Z_hat, out$mu_b, n, p, K),4)
round(mse(jmadata_std, time_points, Zhat_bglasso, b_hat_bglasso, n, p, K),4)
round(mse(jmadata_std, time_points, Zhat_glasso, b_hat_glasso, n, p, K),4)
round(mse(jmadata_std, time_points, Zhat_gmcp, b_hat_gmcp, n, p, K),4)
round(mse(jmadata_std, time_points, Zhat_gscad, b_hat_gscad, n, p, K),4)

extra_ylim <- list(c(-20,30), c(-30,60), c(-30,45), c(0,0), c(-30,20), c(-5,5))

pdf('JMA_res_BGLSS.pdf', width = plot_dim[1], height = plot_dim[2])
par(mar = c(4, 5, 2, .1))    
for(j in 1:6){
  plot(time_finer[,j], (B_finer[[j]]%*%b_hat_bglasso[ids[[j]]])/sd_t_finer[,,j], 
       type = "l", col = bglsscol, ylab = bquote(hat(beta[.(j)])(t)),
       xlab = "Month", ylim = ylim_values[[j]] , lwd = 2)
  abline(a = 0, b = 0, col = truecol, lty = 2)
  # lines(time_finer[,j], LL_bg[[j]], col = cbcol, lty = 3, lwd = 2)
  #lines(time_finer[,j], UL_bg[[j]], col = cbcol, lty = 3, lwd = 2)
  #polygon(c(time_points[,j], rev(time_points[,j])), c(LL[[j]]/sd_t[,,j], rev(UL[[j]]/sd_t[,,j])), col = "#00000030", border = NA)
  legend("topright", c('Estimated curve', 'Credible band'), col = c(bglsscol,cbcol), lty = c(1,3), lwd = 2, bty = "n", cex = 0.8)
}  
dev.off()

ylim_other <- list(c(-600,600), c(-100,100), c(-1,3), c(-400,400), c(-100,100), c(-40,40))
pdf(file = "Comparative_results_jma_vem.pdf", width = plot_dim[1], height = plot_dim[2])
par(mar = c(4, 5, 2, .1))
for(j in 1:p){
  plot(time_finer[,j], (B_finer[[j]]%*%b_hat_bglasso[ids[[j]]])/sd_t_finer[,,j], type = "l",ylab = bquote(hat(beta[.(j)])(t)), xlab = "Month", ylim = ylim_other[[j]] ,
       col = bglsscol, lty = 2, lwd = 2)
  lines(time_finer[,j], (B_finer[[j]]%*%b_hat_glasso[ids[[j]]])/sd_t_finer[,,j], 
        col = glassocol, lwd = 2)
  abline(a = 0, b = 0, col = truecol, lty = 2)
  lines(time_finer[,j], (B_finer[[j]]%*%b_hat_gmcp[ids[[j]]])/sd_t_finer[,,j], col = gmcpcol, lty =  3, lwd = 2)
  lines(time_finer[,j], (B_finer[[j]]%*%b_hat_gscad[ids[[j]]])/sd_t_finer[,,j], col = gscadcol, lty = 4, lwd = 2)
  legend("bottomright", c('grLASSO', 'grMCP', 'grSCAD', 'BGLSS'), col = c(glassocol, gmcpcol, gscadcol, bglsscol), lty = c(1,3,4), lwd = 2, cex = 0.8, bty = 'n')
}
dev.off()

save(out, file = "jma_VEM.rdata")
