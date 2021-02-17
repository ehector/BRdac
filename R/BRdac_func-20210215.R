BRdac <- function(y, covariates_1, covariates_2, covariates_3, locs, indicator){
  time <- proc.time()
  output <- list()
  
  p_1 <- dim(covariates_1)[2]
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  p <- p_1+p_2+p_3
  
  S <- nrow(locs)
  N <- nrow(y)
  K <- length(unique(indicator))
  
  z_1 <- z_constructor(covariates_1, S, N, p_1)
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  estimates <- matrix(0, K, 2+p)
  psi_list <- sensitivity_list <- list()
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    n.pairs <- sum(indicator==k) * (sum(indicator==k) - 1)/2
    extcoeff <- SpatialExtremes::fitextcoeff(y[,indicator==k], locs[indicator==k,], estim = "Smith", plot = FALSE, loess = FALSE, marge = "emp", std.err = FALSE)
    weights <- rep(1, n.pairs)
    extcoeff <- extcoeff[, "ext.coeff"]
    dist <- distance(locs[indicator==k,])
    funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                        as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
    start <- list(range = max(dist)/2.73, smooth = 0.75)
    names(formals(funBR)) <- c("range", "smooth")
    obj.fun <- function(p, ...) funBR(p, ...)
    body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
    cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
    
    mu <- sigma <- xi <- rep(0, sum(indicator==k))
    for (i in 1:(sum(indicator==k))) {
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
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  time_mean <- system.time(mean_estimates <- colMeans(estimates))
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_raw(mean_estimates, y=y[,indicator==k], 
                                 z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                                 z_3=z_3[indicator==k,,,drop=FALSE], 
                                 locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- result_k$sensitivity %*% t(result_k$sensitivity)
  
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  V <- as.matrix(Matrix::bdiag(lapply(psi_list, function(x) x %*% t(x))))
  #V <- psi%*%t(psi)/kronecker(abs(outer(1:K, 1:K , "-"))+1, matrix(1,p+2,p+2))
  #V <- psi%*%t(psi)
  #for(k in 1:(K-2)){
  #  V[((k-1)*(p+2)+1):(k*(p+2)),((k-1)*(p+2)+2*(p+2)+1):nrow(psi)] <- 
  #    V[((k-1)*(p+2)+2*(p+2)+1):nrow(psi), ((k-1)*(p+2)+1):(k*(p+2))] <- 0
  #}
  V_inv <- solve(V)
  weight <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(1:K, 
                                                   function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%psi%*%t(psi)%*%V_inv%*%t(sensitivity) %*% weight
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  return(output)
}

BRdacVCMmu <- function(y, covariates_2, covariates_3, locs, indicator){
  time <- proc.time()
  output <- list()
  
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  
  S <- nrow(locs)
  N <- nrow(y)
  K <- length(unique(indicator))
  
  new_z_1 <- list()
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  estimates <- list()
  psi_list <- sensitivity_list <- list()
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    pos_k <- locs[indicator==k,]
    
    n.pairs <- sum(indicator==k) * (sum(indicator==k) - 1)/2
    extcoeff <- SpatialExtremes::fitextcoeff(y[,indicator==k], pos_k, estim = "Smith", plot = FALSE, loess = FALSE, marge = "emp", std.err = FALSE)
    weights <- rep(1, n.pairs)
    extcoeff <- extcoeff[, "ext.coeff"]
    dist <- distance(pos_k)
    funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                        as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
    start <- list(range = max(dist)/2.73, smooth = 0.75)
    names(formals(funBR)) <- c("range", "smooth")
    obj.fun <- function(p, ...) funBR(p, ...)
    body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
    cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
    
    s_k <- sum(indicator==k)
    mu <- sigma <- xi <- rep(0, s_k)
    for (i in 1:s_k) {
      marg.param <- gevmle(y[,which(indicator==k)[i]])
      mu[i] <- marg.param["loc"]
      sigma[i] <- marg.param["scale"]
      xi[i] <- marg.param["shape"]
    }

    x_2 <- covariates_2[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_3 <- covariates_3[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    
    knots_mat <- pos_k[round(seq(1,s_k,length.out=10)),]
    distances <- apply(knots_mat, 1, function(x) sqrt((pos_k[,1]-x[1])^2 + (pos_k[,2]-x[2])^2))
    basis_1 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.05, 0, 1, t))))
    #basis_2 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.1, 0, 1, t))))
    new_X <- cbind(rep(1,N*s_k), do.call(rbind, lapply(1:N, function(x) pos_k)), do.call(rbind, lapply(1:N, function(x) basis_1)))
    p_1 <- ncol(new_X)
    new_z_1[[k]] <- z_constructor(new_X, s_k, N, p_1)
    
    init_values <- as.vector(c(cov.start[2], cov.start[1], solve(t(new_X)%*%new_X)%*%t(new_X)%*%rep(mu,N), 
                               solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
    
    p <- p_1+p_2+p_3
    opt_2 <-
      nlminb(start=init_values, objective=logCL_all, 
             y=y[,indicator==k], z_1=new_z_1[[k]], z_2=z_2[indicator==k,,,drop=FALSE], 
             z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k,
             lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
    
    estimates[[k]] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_mean <- proc.time()
  estimates <- do.call(rbind, estimates)
  mean_estimates <- colMeans(estimates[,c(1:2,(2+p_1+1):(2+p_1+p_2+p_3))])
  time_mean <- proc.time() - time_mean
  
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_raw(c(mean_estimates[1:2], estimates[k,-c(1:2,(2+K*p_1+1):(2+K*p_1+p_2+p_3))], mean_estimates[-c(1:2)]), y=y[,indicator==k], 
                                 z_1=new_z_1[[k]], z_2=z_2[indicator==k,,,drop=FALSE], 
                                 z_3=z_3[indicator==k,,,drop=FALSE], 
                                 locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- result_k$sensitivity %*% t(result_k$sensitivity)
    
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  sensitivity <- matrix(0, 2+K*p_1+p_2+p_3, K*(2+p_1+p_2+p_3))
  sensitivity[c(1:2,(2+K*p_1+1):(2+K*p_1+p_2+p_3)),] <- do.call(cbind, lapply(sensitivity_list, function(x) x[c(1:2,(2+p_1+1):(2+p_1+p_2+p_3)),]))
  sensitivity[-c(1:2,(2+K*p_1+1):(2+K*p_1+p_2+p_3)),] <- 
    as.matrix(Matrix::bdiag(lapply(sensitivity_list, function(x) x[-c(1:2,(2+p_1+1):(2+p_1+p_2+p_3)),])))
  V <- as.matrix(Matrix::bdiag(lapply(psi_list, function(x) x %*% t(x))))
  V_inv <- solve(V)
  weight <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(1:K, 
                                                   function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%psi%*%t(psi)%*%V_inv%*%t(sensitivity) %*% weight
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  beta_1 <- apply(locs, 1, function(x) sqrt(x[1]^2+x[2]^2))/10
  X_tilde <- as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,1,])))
  pos_cols <- do.call(rbind, lapply(new_z_1, function(x) x[,1,2:3]))
  X_tilde <- X_tilde[order(pos_cols[,2], pos_cols[,1]),]
  output$mu.fitted.values <- as.vector(X_tilde%*%output$coefficients[-c(1:2,(2+K*p_1+1):(2+K*p_1+p_2+p_3))])
  se <- sqrt(diag(X_tilde %*% output$vcov[-c(1:2,(2+K*p_1+1):(2+K*p_1+p_2+p_3)),-c(1:2,(2+K*p_1+1):(2+K*p_1+p_2+p_3))] %*%t(X_tilde)))
  output$PCP <- mean(beta_1 <= output$mu.fitted.values + 1.96*se & beta_1 >= output$mu.fitted.values-1.96*se)
  
  return(output)
}

BRdac_two <- function(y, covariates_1, covariates_2, covariates_3, locs, indicator_dep, indicator_marg){
  time <- proc.time()
  output <- list()
  
  p_1 <- dim(covariates_1)[2]
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  p <- p_1+p_2+p_3
  
  S <- nrow(locs)
  N <- nrow(y)
  K_dep <- length(unique(indicator_dep))
  K_marg <- length(unique(indicator_marg))
  
  z_1 <- z_constructor(covariates_1, S, N, p_1)
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  estimates_dep <- init_values_dep <- matrix(0, K_dep, 2+p)
  estimates_marg <- init_values_marg <- matrix(0, K_marg, 2+p)
  psi_list_dep <- psi_list_marg <- list()
  sensitivity_list_dep <- sensitivity_list_marg <- list()
  
  time_k_dep <- time_k_marg <- list()
  time_before <- proc.time() - time
  for(k in 1:K_dep){
    time_k_dep[[k]] <- proc.time()
    
    n.pairs <- sum(indicator_dep==k) * (sum(indicator_dep==k) - 1)/2
    extcoeff <- SpatialExtremes::fitextcoeff(y[,indicator_dep==k], locs[indicator_dep==k,], estim = "Smith", plot = FALSE, loess = FALSE, marge = "emp", std.err = FALSE)
    weights <- rep(1, n.pairs)
    extcoeff <- extcoeff[, "ext.coeff"]
    dist <- distance(locs[indicator_dep==k,])
    funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                        as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
    start <- list(range = max(dist)/2.73, smooth = 0.75)
    names(formals(funBR)) <- c("range", "smooth")
    obj.fun <- function(p, ...) funBR(p, ...)
    body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
    cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
    
    mu <- sigma <- xi <- rep(0, )
    for (i in 1:(sum(indicator_dep==k))) {
      marg.param <- gevmle(y[,which(indicator_dep==k)[i]])
      mu[i] <- marg.param["loc"]
      sigma[i] <- marg.param["scale"]
      xi[i] <- marg.param["shape"]
    }
    x_1 <- covariates_1[c(t(sapply(which(indicator_dep==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_2 <- covariates_2[c(t(sapply(which(indicator_dep==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_3 <- covariates_3[c(t(sapply(which(indicator_dep==k), function(x) x+seq(0,(N-1)*S,S)))),]
    
    init_values_dep[k,] <- as.vector(c(cov.start[2], cov.start[1], solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                                       solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
    
    opt_2 <-
      nlminb(start=init_values_dep[k,], objective=logCL_all, 
             y=y[,indicator_dep==k], z_1=z_1[indicator_dep==k,,,drop=FALSE], z_2=z_2[indicator_dep==k,,,drop=FALSE], 
             z_3=z_3[indicator_dep==k,,,drop=FALSE], locs=locs[indicator_dep==k,],
             lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
    
    estimates_dep[k,] <- opt_2$par
    time_k_dep[[k]] <- proc.time() - time_k_dep[[k]]
  }
  time_mean <- system.time(mean_estimates_dep <- colMeans(estimates_dep))
  for(k in 1:K_marg){
    time_k_marg[[k]] <- proc.time()
    
    n.pairs <- sum(indicator_marg==k) * (sum(indicator_marg==k) - 1)/2
    extcoeff <- SpatialExtremes::fitextcoeff(y[,indicator_marg==k], locs[indicator_marg==k,], estim = "Smith", plot = FALSE, loess = FALSE, marge = "emp", std.err = FALSE)
    weights <- rep(1, n.pairs)
    extcoeff <- extcoeff[, "ext.coeff"]
    dist <- distance(locs[indicator_marg==k,])
    funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                        as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
    start <- list(range = max(dist)/2.73, smooth = 0.75)
    names(formals(funBR)) <- c("range", "smooth")
    obj.fun <- function(p, ...) funBR(p, ...)
    body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
    cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
    
    mu <- sigma <- xi <- rep(0, )
    for (i in 1:(sum(indicator_marg==k))) {
      marg.param <- gevmle(y[,which(indicator_marg==k)[i]])
      mu[i] <- marg.param["loc"]
      sigma[i] <- marg.param["scale"]
      xi[i] <- marg.param["shape"]
    }
    x_1 <- covariates_1[c(t(sapply(which(indicator_marg==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_2 <- covariates_2[c(t(sapply(which(indicator_marg==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_3 <- covariates_3[c(t(sapply(which(indicator_marg==k), function(x) x+seq(0,(N-1)*S,S)))),]
    
    init_values_marg[k,] <- as.vector(c(cov.start[2], cov.start[1], solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                                        solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
    
    opt_2 <-
      nlminb(start=init_values_marg[k,], objective=logCL_all, 
             y=y[,indicator_marg==k], z_1=z_1[indicator_marg==k,,,drop=FALSE], z_2=z_2[indicator_marg==k,,,drop=FALSE], 
             z_3=z_3[indicator_marg==k,,,drop=FALSE], locs=locs[indicator_marg==k,],
             lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
    
    estimates_marg[k,] <- opt_2$par
    time_k_marg[[k]] <- proc.time() - time_k_marg[[k]]
  }
  time_mean_2 <- system.time(mean_estimates <- colMeans(estimates_marg))
  time_k_2_dep <- time_k_2_marg <- list()
  for(k in 1:K_dep){
    time_k_2_dep[[k]] <- proc.time()
    
    both <- Chessian_all_raw(mean_estimates_dep, y=y[,indicator_dep==k], 
                             z_1=z_1[indicator_dep==k,,,drop=FALSE], z_2=z_2[indicator_dep==k,,,drop=FALSE], 
                             z_3=z_3[indicator_dep==k,,,drop=FALSE], 
                             locs=locs[indicator_dep==k,])
    psi_list_dep[[k]] <- both$score
    sensitivity_list_dep[[k]] <- (-both$sensitivity %*% t(both$sensitivity))
    time_k_2_dep[[k]] <- proc.time() - time_k_2_dep[[k]]
  }
  for(k in 1:K_marg){
    time_k_2_marg[[k]] <- proc.time()
    
    both <- Chessian_all_raw(mean_estimates, y=y[,indicator_marg==k], 
                             z_1=z_1[indicator_marg==k,,,drop=FALSE], z_2=z_2[indicator_marg==k,,,drop=FALSE], 
                             z_3=z_3[indicator_marg==k,,,drop=FALSE], 
                             locs=locs[indicator_marg==k,])
    psi_list_marg[[k]] <- both$score
    sensitivity_list_marg[[k]] <- (-both$sensitivity %*% t(both$sensitivity))
    time_k_2_marg[[k]] <- proc.time() - time_k_2_marg[[k]]
  }
  
  time_after <- proc.time()
  s_k_dep <- sum(indicator_dep==1)
  s_k_marg <- sum(indicator_marg==1)
  psi <- rbind(do.call(rbind,psi_list_dep), do.call(rbind, psi_list_marg))
  sensitivity <- as.matrix(cbind(do.call(cbind, sensitivity_list_dep), do.call(cbind, sensitivity_list_marg)))
  norm_mat <- diag(nrow(psi))
  diag(norm_mat)[1:K_dep] <- sqrt(s_k_dep)
  diag(norm_mat)[(K_dep+1):nrow(psi)] <- sqrt(s_k_marg)
  V <- norm_mat%*% psi %*% t(psi)%*%norm_mat
  V_inv <- solve(V)
  output$coefficients <- as.vector(solve(sensitivity %*% V_inv %*% t(sensitivity)) %*% 
                                     sensitivity %*% V_inv %*% 
                                     c(unlist(lapply(1:K_dep, 
                                                     function(k) sensitivity_list_dep[[k]] %*% estimates_dep[k,])),
                                       unlist(lapply(1:K_marg, 
                                                     function(k) sensitivity_list_marg[[k]] %*% estimates_marg[k,])))) 
  output$vcov <- solve(sensitivity %*% V_inv %*% t(sensitivity))%*%sensitivity%*%V_inv%*%psi%*%t(psi)%*%V_inv%*%t(sensitivity)%*%solve(sensitivity %*% V_inv %*% t(sensitivity))
  time_after <- proc.time()-time_after
  
  time_k <- rbind(time_k_dep[[which.max(sapply(time_k_dep, function(x) x[3]))]], 
                  time_k_marg[[which.max(sapply(time_k_marg, function(x) x[3]))]])
  time_k_2 <- rbind(time_k_2_dep[[which.max(sapply(time_k_2_dep, function(x) x[3]))]], 
                    time_k_2_marg[[which.max(sapply(time_k_2_marg, function(x) x[3]))]])
  output$time <- time_before + time_mean + time_mean_2 + time_after + time_k[which.max(time_k[,3]),] + time_k_2[which.max(time_k_2[,3]),]
  return(output)
}

CL <- function(y, covariates_1, covariates_2, covariates_3, locs){
  time <- proc.time()
  output <- list()
  
  p_1 <- dim(covariates_1)[2]
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  p <- p_1+p_2+p_3
  
  S <- nrow(locs)
  N <- nrow(y)
  
  z_1 <- z_constructor(covariates_1, S, N, p_1)
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  n.pairs <- S * (S - 1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim = "Smith", plot = FALSE, loess = FALSE, marge = "emp", std.err = FALSE)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  
  mu <- sigma <- xi <- rep(0, )
  for (i in 1:S) {
    marg.param <- gevmle(y[,i])
    mu[i] <- marg.param["loc"]
    sigma[i] <- marg.param["scale"]
    xi[i] <- marg.param["shape"]
  }
  x_1 <- covariates_1[c(t(sapply(1:S, function(x) x+seq(0,(N-1)*S,S)))),]
  x_2 <- covariates_2[c(t(sapply(1:S, function(x) x+seq(0,(N-1)*S,S)))),]
  x_3 <- covariates_3[c(t(sapply(1:S, function(x) x+seq(0,(N-1)*S,S)))),]
  
  init_values <- as.vector(c(cov.start[2], cov.start[1], solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                             solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
  
  opt_2 <-
    nlminb(start=init_values, objective=logCL_all, y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs,
           lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
  
  output$coefficients <- opt_2$par
  
  result <- Chessian_all_raw(output$coefficients, y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs)
  psi <- result$score
  sensitivity <- result$sensitivity %*% t(result$sensitivity)
  output$vcov <- solve(sensitivity)%*%psi%*%t(psi)%*%solve(sensitivity)
  
  output$time <- proc.time()-time
  
  return(output)
}

BRdac_quadratic <- function(y, covariates_1, covariates_2, covariates_3, locs, indicator){
  time <- proc.time()
  output <- list()
  
  p_1 <- dim(covariates_1)[2]
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  p <- p_1+p_2+p_3
  
  S <- nrow(locs)
  N <- nrow(y)
  K <- length(unique(indicator))
  
  z_1 <- z_constructor(covariates_1, S, N, p_1)
  z_2 <- z_constructor(covariates_2, S, N, p_2)
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  init_values <- matrix(0, K, 2+p)
  psi_list <- list()
  sensitivity_list <- list()
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    cov.start <- SpatialExtremes::fitcovariance(y[,indicator==k], locs[indicator==k,], "brown", marge = "emp")$param
    mu <- sigma <- xi <- rep(0, sum(indicator==k))
    for (i in 1:(sum(indicator==k))) {
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

BRdac_BR <- function(x, locs, indicator){
  time <- proc.time()
  output <- list()
  
  S <- nrow(locs)
  N <- nrow(x)
  K <- length(unique(indicator))
  
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
  time_after <- proc.time()-time_after
  output$time <- time_before + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]]
  return(output)
}

gaussian_basis <- function(phi2, tau2, sigma2, dist){
  if(dist>0) return(sigma2*exp(-phi2*dist^2))
  else return(tau2+sigma2)
}
