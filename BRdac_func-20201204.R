BRdac <- function(y, covariates_1, covariates_2, covariates_3, locs){
  output <- list()
  
  p_1 <- dim(covariates_1)[2]
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  p <- p_1+p_2+p_3
  
  S <- nrow(locs)
  s_k <- 2
  K <- S/s_k
  N <- nrow(y)
  
  z_1 <- z_constructor(covariates_1, S, N, p_1)
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  indicator <- c(1:K,K:1)
  estimates <- matrix(0, K, 2+p)
  psi_list <- list()
  sensitivity_list <- list()
  
  cov.start <- fitcovariance(y, locs, "brown", marge = "emp")$param
  init_alpha <- cov.start[2]
  init_phi <- cov.start[1]
  init_beta_1 <- sapply(1:p_1, function(x) (1+min(y))/max(z_1[,,x]) )-1
  init_values <- c(init_alpha, init_phi, init_beta_1, rep(1,p_2+p_3))
  
  for(k in 1:K){
    opt_2 <-
      nlminb(start=init_values, objective=logCL_all_wrap, gradient=Cscore_all_sum_wrap, 
             y=y[,indicator==k], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
             z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
             lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
    
    estimates[k,] <- opt_2$par
    psi_list[[k]] <- do.call(cbind, Cscore_all(opt_2$par[1], opt_2$par[2], opt_2$par[3:(p_1+2)], opt_2$par[(3+p_1):(2+p_1+p_2)], 
                                                opt_2$par[(3+p_1+p_2):(2+p)], y=y[,indicator==k], 
                                                z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                                                z_3=z_3[indicator==k,,,drop=FALSE], 
                                                locs=locs[indicator==k,])) 
    sensitivity_list[[k]] <- psi_list[[k]]%*%t(psi_list[[k]])
  }
  
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  V <- psi %*% t(psi)
  V_inv <- MASS::ginv(V)
  #sensitivity <- do.call(cbind, lapply(1:(K*(K-1)/2), function(x) V[((x-1)*(2+p)+1):(x*(2+p)),((x-1)*(2+p)+1):(x*(2+p))]))
  output$coefficients <- as.vector(solve(sensitivity %*% V_inv %*% t(sensitivity)) %*% 
                                            sensitivity %*% V_inv %*% 
                                            unlist(lapply(1:K, 
                                                          function(x) sensitivity[,((x-1)*(2+p)+1):(x*(2+p))] %*% estimates[x,]))) 
  
  
}
