BRdac <- function(y, covariates_1, covariates_2, covariates_3, locs, K=NULL, s_k=NULL){
  time <- proc.time()
  output <- list()
  
  p_1 <- dim(covariates_1)[2]
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  p <- p_1+p_2+p_3
  
  S <- nrow(locs)
  N <- nrow(y)
  if(!is.null(K) & !is.null(s_k)) if(s_k!=S/K) stop("Incorrect specification of s_k and K.")
  if(is.null(K) & !is.null(s_k)) K <- S/s_k
  if(!is.null(K) & is.null(s_k)) s_k <- S/K
  if(is.null(K) & is.null(s_k)) {
    K <- 5
    s_k <- S/K
  }
  
  z_1 <- z_constructor(covariates_1, S, N, p_1)
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  set.seed(1)
  indicator <- rep(1:K, s_k)[sample(1:S, S, replace=FALSE)]
  estimates <- matrix(0, K, 2+p)
  psi_list <- list()
  sensitivity_list <- list()
  
  #cov.start <- SpatialExtremes::fitcovariance(y, locs, "brown", marge = "emp")$param
  #mu <- sigma <- xi <- rep(0, S)
  #for (i in 1:S) {
  #  marg.param <- gevmle(y[, i])
  #  mu[i] <- marg.param["loc"]
  #  sigma[i] <- marg.param["scale"]
  #  xi[i] <- marg.param["shape"]
  #}
  #init_beta_1 <- lm(rep(mu,N) ~ 0+covariates_1)$coefficients
  #init_beta_2 <- lm(rep(log(sigma),N) ~ 0+covariates_2)$coefficients
  #init_beta_3 <- lm(rep(xi,N) ~ 0+covariates_3)$coefficients
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    cov.start <- SpatialExtremes::fitcovariance(y[,indicator==k], locs[indicator==k,], "brown", marge = "emp")$param
    mu <- sigma <- xi <- rep(0, s_k)
    for (i in 1:s_k) {
      marg.param <- gevmle(y[,which(indicator==k)[i]])
      mu[i] <- marg.param["loc"]
      sigma[i] <- marg.param["scale"]
      xi[i] <- marg.param["shape"]
    }
    x_1 <- covariates_1[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_2 <- covariates_2[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_3 <- covariates_3[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    
    init_values <- c(cov.start[2], cov.start[1], solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                     solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N))
    
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
    #sensitivity_list[[k]] <- numDeriv::hessian(function(x) logCL_all(x[1], x[2], x[3:(p_1+2)], x[(3+p_1):(2+p_1+p_2)], 
    #                                                                 x[(3+p_1+p_2):(2+p)], y=y[,indicator==k], z_1=z_1[indicator==k,,,drop=FALSE], 
    #                                                                 z_2=z_2[indicator==k,,,drop=FALSE], z_3=z_3[indicator==k,,,drop=FALSE], 
    #                                                                 locs=locs[indicator==k,]) , x=opt_2$par)
    sensitivity_list[[k]] <- -Chessian_all_sum(opt_2$par[1], opt_2$par[2], opt_2$par[3:(p_1+2)], opt_2$par[(3+p_1):(2+p_1+p_2)], 
                                      opt_2$par[(3+p_1+p_2):(2+p)], y=y[,indicator==k], 
                                      z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                                      z_3=z_3[indicator==k,,,drop=FALSE], 
                                      locs=locs[indicator==k,])
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  time_after <- proc.time()
  
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  V <- psi %*% t(psi)
  V_inv <- solve(V)
  output$coefficients <- as.vector(solve(sensitivity %*% V_inv %*% t(sensitivity)) %*% 
                                            sensitivity %*% V_inv %*% 
                                            unlist(lapply(1:K, 
                                                          function(x) sensitivity[,((x-1)*(2+p)+1):(x*(2+p))] %*% estimates[x,]))) 
  output$vcov <- solve(sensitivity %*% V_inv %*% t(sensitivity))*sqrt(N)/K
  time_after <- proc.time()-time_after
  output$time <- time_before + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]]
  return(output)
}