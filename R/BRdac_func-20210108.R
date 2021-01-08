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
  indicator <- rep(1:K,s_k)
  estimates <- matrix(0, K, 2+p)
  psi_list <- list()
  sensitivity_list <- list()
  
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
    
    init_values <- as.vector(c(cov.start[2], cov.start[1], solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                               solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
    
    opt_2 <-
      nlminb(start=init_values, objective=logCL_all, 
             y=y[,indicator==k], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
             z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
             lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
    
    estimates[k,] <- opt_2$par
    #psi_list[[k]] <- Cscore_all(estimates[k,], y=y[,indicator==k], 
    #                            z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
    #                            z_3=z_3[indicator==k,,,drop=FALSE], 
    #                            locs=locs[indicator==k,])
    #sensitivity_list[[k]] <- -Chessian_all_sum(estimates[k,], y=y[,indicator==k], 
    #                                           z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
    #                                           z_3=z_3[indicator==k,,,drop=FALSE], 
    #                                           locs=locs[indicator==k,])
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  time_mean <- system.time(mean_estimates <- colMeans(estimates))
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()

    psi_list[[k]] <- Cscore_all(mean_estimates, y=y[,indicator==k], 
                                z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                                z_3=z_3[indicator==k,,,drop=FALSE], 
                                locs=locs[indicator==k,])
    sensitivity_list[[k]] <- -Chessian_all_sum(mean_estimates, y=y[,indicator==k], 
                                               z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                                               z_3=z_3[indicator==k,,,drop=FALSE], 
                                               locs=locs[indicator==k,])
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  V <- psi %*% t(psi)
  V_inv <- solve(V)
  output$coefficients <- as.vector(solve(sensitivity %*% V_inv %*% t(sensitivity)) %*% 
                                            sensitivity %*% V_inv %*% 
                                            unlist(lapply(1:K, 
                                                          function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
  output$vcov <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  return(output)
}

BRdac_quadratic <- function(y, covariates_1, covariates_2, covariates_3, locs, K=NULL, s_k=NULL){
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
  indicator <- rep(1:K,s_k)
  init_values <- matrix(0, K, 2+p)
  psi_list <- list()
  sensitivity_list <- list()
  
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
    
    init_values[k,] <- as.vector(c(cov.start[2], cov.start[1], solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                               solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
    time_k[[k]] <- proc.time() - time_k[[k]]
  }

  time_after <- proc.time()
  
  init_values <- colMeans(init_values)
  opt <-
    nlminb(start=init_values, objective=quadratic_form_par, y=y[,order(indicator)], z_1=z_1[order(indicator),,,drop=FALSE], 
           z_2=z_2[order(indicator),,,drop=FALSE], z_3=z_3[order(indicator),,,drop=FALSE], locs=locs[order(indicator),], 
           K=K, lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
  output$coefficients <- opt$par
  
  psi <- matrix(0, (p+2)*K, N)
  sensitivity <- matrix(0, p+2, (p+2)*K)
  for(k in 1:K){
    result_k <- Chessian_all_raw(output$coefficients, y=y[,indicator==k], 
                                 z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                                 z_3=z_3[indicator==k,,,drop=FALSE], 
                                 locs=locs[indicator==k,])
    psi[((k-1)*(p+2)+1):(k*(p+2)),] <- result_k$score
    sensitivity[,((k-1)*(p+2)+1):(k*(p+2))] <- result_k$sensitivity %*% t(result_k$sensitivity)
  }
  output$vcov <- solve(sensitivity %*% solve(psi %*% t(psi)) %*% t(sensitivity))
  time_after <- proc.time() - time_after
  
  output$time <- time_before + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]]
  return(output)
}

BRdac_BR <- function(x, locs, K=NULL, s_k=NULL){
  time <- proc.time()
  output <- list()
  
  S <- nrow(locs)
  N <- nrow(x)
  if(!is.null(K) & !is.null(s_k)) if(s_k!=S/K) stop("Incorrect specification of s_k and K.")
  if(is.null(K) & !is.null(s_k)) K <- S/s_k
  if(!is.null(K) & is.null(s_k)) s_k <- S/K
  if(is.null(K) & is.null(s_k)) {
    K <- 5
    s_k <- S/K
  }
  
  set.seed(1)
  indicator <- rep(1:K,s_k)
  estimates <- matrix(0, K, 2)
  psi_list <- list()
  sensitivity_list <- list()
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    cov.start <- SpatialExtremes::fitcovariance(x[,indicator==k], locs[indicator==k,], "brown", marge = "emp")$param
    
    init_values <- as.vector(cov.start[2:1])
    
    opt_2 <- nlminb(start=init_values, objective=logCL_BR, x=x[,indicator==k], locs=locs[indicator==k,], lower=c(1e-5,1e-5), upper=c(2-1e-5,Inf))
    
    estimates[k,] <- opt_2$par
    psi_list[[k]] <- Cscore_BR(estimates[k,], x=x[,indicator==k], locs=locs[indicator==k,])
    #sensitivity_list[[k]] <- numDeriv::hessian(function(x) logCL_all(x, y=y[,indicator==k], z_1=z_1[indicator==k,,,drop=FALSE], 
    #                                                                 z_2=z_2[indicator==k,,,drop=FALSE], z_3=z_3[indicator==k,,,drop=FALSE], 
    #                                                                 locs=locs[indicator==k,]) , x=estimates[k,])
    sensitivity_list[[k]] <- -psi_list[[k]]%*%t(psi_list[[k]])
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
                                                   function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
  output$vcov <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  #output <- coeff_vcov_comp(sensitivity, V, unlist(lapply(1:K, function(x) sensitivity_list[[x]] %*% estimates[x,])))
  #output$coefficients <- as.vector(output$coefficients)
  time_after <- proc.time()-time_after
  output$time <- time_before + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]]
  return(output)
}
