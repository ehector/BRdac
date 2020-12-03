BRdac <- function(y, covariates_1, covariates_2, covariates_3, locs){
  output <- list()
  
  p_1 <- dim(covariates_1)[2]
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  p <- p_1+p_2+p_3
  
  z_1 <- z_constructor(covariates_1, S, N, p_1)
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  S <- nrow(locs)
  s_k <- 2
  K <- S/s_k
  indicator <- sort(rep(1:K, s_k))
  
  estimates <- matrix(0, K*(K-1)/2, 2+p)
  psi <- list()
  it <- 1
  
  for(k in 1:(K-1)){
    for(j in (k+1):K){
      opt_2 <-
        nlminb(start=c(opt$par, rep(1,p)), objective=logCL_all_wrap, gradient=Cscore_all_sum_wrap, 
               y=y[,c(k,j)], z_1=z_1[c(k,j),,], z_2=z_2[c(k,j),,], z_3=z_3[c(k,j),,], locs=locs[c(k,j),],
               lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
      
      estimates[it,] <- opt_2$par
      psi[[it]] <- do.call(cbind, Cscore_all(opt_2$par[1], opt_2$par[2], opt_2$par[3:(p_1+2)], opt_2$par[(3+p_1):(2+p_1+p_2)], 
                                            opt_2$par[(3+p_1+p_2):(2+p)], y=y[,c(k,j)], 
                                            z_1=z_1[c(k,j),,], z_2=z_2[c(k,j),,], z_3=z_3[c(k,j),,], 
                                            locs=locs[c(k,j),])) 
      it <- it + 1
    }
  }
  
  psi_all <- do.call(rbind,psi)
  V <- psi_all %*% t(psi_all)
  V_inv <- MASS::ginv(V)
  sensitivity <- do.call(cbind, lapply(1:(K*(K-1)/2), function(x) V[((x-1)*(2+p)+1):(x*(2+p)),((x-1)*(2+p)+1):(x*(2+p))]))
  output$coefficients <- as.vector(solve(sensitivity %*% V_inv %*% t(sensitivity)) %*% 
                                            sensitivity %*% V_inv %*% 
                                            unlist(lapply(1:(K*(K-1)/2), 
                                                          function(x) sensitivity[,((x-1)*(2+p)+1):(x*(2+p))] %*% estimates[x,]))) 
  
  output$coefficients <- rep(0, 2+p)
  
}