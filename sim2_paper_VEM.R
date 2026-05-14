# sim_paper_prop not to be in the R Package
sim2_paper_VEM <- function(seed, nsim, initial_values, folder, 
                          n, sigma2, Z, p, K, nt, 
                          delta1_0, delta2_0, a0, b0,
                          Niter, convergence_threshold, ordem, std = TRUE){
  
  ids <- split(1:(K*p), rep(1:p, each = K))
  
  seed = seed + nsim
  dir.create(folder, showWarnings = FALSE)
  
  #generate data
  data <- gen_data_ex2(seed = seed, n = n, nt = nt, p = p, ordem = ordem, K = K, Z = Z, sigma2 = sigma2)
  save(data, file = paste0(folder,"/data_", nsim, ".RData"))
  
  
  vb_res <- VBSOFR_VS_VEM(initial_values = initial_values, data = data, 
                      data_std = data, n = n, K = K, p = p,
                      delta1_0 = delta1_0, delta2_0 = delta2_0, a0 = a0, b0 = b0,
                      Niter = Niter, convergence_threshold = convergence_threshold, std = std)

  
  return(vb_res)
  
}

