BRdac <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator){
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
  
  if(quantile>=0) {
    thresholds <- apply(Y, 2, function(x) quantile(x, prob=quantile))
  } else {
    thresholds <- rep(min(Y)-1, S)
  }
  frech <- t(apply(y, 1, function(x) SpatialExtremes::gev2frech(x, emp=TRUE)))
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    s_k <- sum(indicator==k)
    
    opt <- nlminb(start=c(0.75, max(distance(locs[indicator==k,]))/2.73), objective=logCL_BR_thresh, 
                  thresholds=thresholds[indicator==k], x=frech[,indicator==k], locs=locs[indicator==k,], 
                  lower=c(1e-5,1e-5), upper=c(2-1e-5,Inf))
    
    mu <- sigma <- xi <- rep(0, s_k)
    for (i in 1:s_k) {
      y_i <- y[,which(indicator==k)[i]]
      marg.param <- gevmle(y_i[y_i>quantile(y_i, prob=0.5)])
      mu[i] <- marg.param["loc"]
      sigma[i] <- marg.param["scale"]
      xi[i] <- marg.param["shape"]
    }
    x_1 <- covariates_1[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_2 <- covariates_2[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    x_3 <- covariates_3[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    
    init_values <- as.vector(c(opt$par, 
                               solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                               solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), 
                               solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
    
    opt_2 <-
      nlminb(start=init_values, objective=logCL_all_thresh, thresholds=thresholds[indicator==k],
             y=y[,indicator==k,drop=FALSE], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
             z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
             lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
    
    estimates[k,] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  time_mean <- system.time(mean_estimates <- colMeans(estimates))
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_thresh(mean_estimates, thresholds=thresholds[indicator==k], y=y[,indicator==k], 
                                    z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                                    z_3=z_3[indicator==k,,,drop=FALSE], 
                                    locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- result_k$sensitivity %*% t(result_k$sensitivity)
    
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  order <- 1:K
  psi <- do.call(rbind,psi_list[order])
  sensitivity <- do.call(cbind, sensitivity_list[order])
  V <- (psi-rowMeans(psi))%*%t(psi-rowMeans(psi))
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  weight <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(order, 
                                                   function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%psi%*%t(psi)%*%V_inv%*%t(sensitivity) %*% weight
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  return(output)
}

BRdacVCMmu <- function(y, covariates_2, covariates_3, locs, indicator, lambda_seq){
  time <- proc.time()
  output <- list()
  
  p_2 <- dim(covariates_2)[2]
  p_3 <- dim(covariates_3)[2]
  
  S <- nrow(locs)
  N <- nrow(y)
  K <- length(unique(indicator))
  
  new_z_1 <- list()
  new_z_2 <- z_constructor(covariates_2, S, N, p_2)
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
  
  V <- (psi-rowMeans(psi))%*%t(psi-rowMeans(psi))
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  
  quad_form <- sensitivity %*% V_inv %*% t(sensitivity)
  
  pos_cols <- do.call(rbind, lapply(new_z_1, function(x) x[,1,2:3]))
  beta_1 <- apply(locs, 1, function(x) sqrt(x[1]^2+x[2]^2))/10
  X_tilde_1 <- as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,1,])))
  X_tilde_1 <- X_tilde_1[order(pos_cols[,2], pos_cols[,1]),]
  
  output$coefficients <- output$vcov <- output$mu.fitted.values <- output$mu.se <- list()
  output$mu.PCP <- output$GCV <- vector("numeric", length(lambda_seq))
  time_k_4 <- list()
  
  for(l in 1:length(lambda_seq)){
    lambda <- lambda_seq[l]
    
    penalty <- matrix(0, 2+K*p_1+p_2+p_3, 2+K*p_1+p_2+p_3)
    diag(penalty)[3:(2+K*p_1)] <- N*lambda
    
    weight <- solve(quad_form + penalty)
    output$coefficients[[l]] <- as.vector(weight %*% 
                                       sensitivity %*% V_inv %*% 
                                       unlist(lapply(1:K, 
                                                     function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
    if(sum(is.na(output$coefficients[[l]]))!=0) break
    output$vcov[[l]] <- weight %*% sensitivity%*%V_inv%*%psi%*%t(psi)%*%V_inv%*%t(sensitivity) %*% weight
    
    output$mu.fitted.values[[l]] <- as.vector(X_tilde_1%*%output$coefficients[[l]][3:(2+K*p_1)])
    output$mu.se[[l]] <- sqrt(diag(X_tilde_1 %*% output$vcov[[l]][3:(2+K*p_1), 3:(2+K*p_1)] %*%t(X_tilde_1)))
    output$mu.PCP[l] <- mean(beta_1 <= output$mu.fitted.values[[l]] + 1.96*output$mu.se[[l]] & 
                            beta_1 >= output$mu.fitted.values[[l]]-1.96*output$mu.se[[l]])
    time_k_3 <- list()
    for(k in 1:K){
      time_k_3[[k]] <- proc.time()
      
      result_k <- Chessian_all_raw(c(output$coefficients[[l]][1:2], output$coefficients[[l]][(2+(k-1)*p_1+1):(2+k*p_1)], output$coefficients[[l]][(2+K*p_1+1):(2+K*p_1+p_2+p_3)]), 
                                   y=y[,indicator==k], z_1=new_z_1[[k]], z_2=z_2[indicator==k,,,drop=FALSE], 
                                   z_3=z_3[indicator==k,,,drop=FALSE], 
                                   locs=locs[indicator==k,])
      
      psi_list[[k]] <- result_k$score
      sensitivity_list[[k]] <- result_k$sensitivity %*% t(result_k$sensitivity)
      
      time_k_3[[k]] <- proc.time() - time_k_3[[k]]
    }
    time_k_4[[l]] <- time_k_3[[which.max(sapply(time_k_3, function(x) x[3]))]]
    psi <- do.call(rbind,psi_list)
    GCV_numerator <- N^(-1)*as.numeric(rowMeans(psi)%*%V_inv%*%rowMeans(psi))
    output$GCV[l] <- GCV_numerator / ( 1 - N^(-1)*sum(diag( weight %*% quad_form )))^2
    
  }
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  return(output)
}

BRdacVCMmusigma <- function(y, covariates_3, locs, indicator, lambda_seq){
  time <- proc.time()
  output <- list()
  
  p_3 <- dim(covariates_3)[2]
  
  S <- nrow(locs)
  N <- nrow(y)
  K <- length(unique(indicator))
  
  new_z_1 <- new_z_2 <- list()
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

    x_3 <- covariates_3[c(t(sapply(which(indicator==k), function(x) x+seq(0,(N-1)*S,S)))),]
    
    knots_mat <- pos_k[round(seq(1,s_k,length.out=10)),]
    distances <- apply(knots_mat, 1, function(x) sqrt((pos_k[,1]-x[1])^2 + (pos_k[,2]-x[2])^2))
    basis_1 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.05, 0, 1, t))))
    #basis_2 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.1, 0, 1, t))))
    new_X <- cbind(rep(1,N*s_k), do.call(rbind, lapply(1:N, function(x) pos_k)), do.call(rbind, lapply(1:N, function(x) basis_1)))
    p_1 <- p_2 <- ncol(new_X)
    new_z_1[[k]] <- new_z_2[[k]] <- z_constructor(new_X, s_k, N, p_1)
    
    init_values <- as.vector(c(cov.start[2], cov.start[1], solve(t(new_X)%*%new_X)%*%t(new_X)%*%rep(mu,N), 
                               solve(t(new_X)%*%new_X)%*%t(new_X)%*%rep(log(sigma),N), solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
    
    p <- p_1+p_2+p_3
    opt_2 <-
      nlminb(start=init_values, objective=logCL_all, 
             y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]], 
             z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k,
             lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
    
    estimates[[k]] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_mean <- proc.time()
  estimates <- do.call(rbind, estimates)
  mean_estimates <- colMeans(estimates[,c(1:2,(2+p_1+p_2+1):(2+p_1+p_2+p_3))])
  time_mean <- proc.time() - time_mean
  
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_raw(c(mean_estimates[1:2], estimates[k,c(3:(2+p_1+p_2))], mean_estimates[-c(1:2)]), 
                                 y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]], z_3=z_3[indicator==k,,,drop=FALSE], 
                                 locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- result_k$sensitivity %*% t(result_k$sensitivity)
    
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  
  sensitivity <- do.call(cbind, sensitivity_list)
  sensitivity <- matrix(0, 2+K*(p_1+p_2)+p_3, K*(2+p_1+p_2+p_3))
  sensitivity[c(1:2,(2+K*(p_1+p_2)+1):(2+K*(p_1+p_2)+p_3)),] <- do.call(cbind, lapply(sensitivity_list, function(x) x[c(1:2,(2+p_1+p_2+1):(2+p_1+p_2+p_3)),]))
  sensitivity[-c(1:2,(2+K*(p_1+p_2)+1):(2+K*(p_1+p_2)+p_3)),] <- 
    as.matrix(Matrix::bdiag(lapply(sensitivity_list, function(x) x[-c(1:2,(2+p_1+p_2+1):(2+p_1+p_2+p_3)),])))
  
  V <- (psi-rowMeans(psi))%*%t(psi-rowMeans(psi))
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  
  quad_form <- sensitivity %*% V_inv %*% t(sensitivity)
  
  pos_cols <- do.call(rbind, lapply(new_z_1, function(x) x[,1,2:3]))
  beta_1 <- apply(locs, 1, function(x) sqrt(x[1]^2+x[2]^2))/10
  X_tilde_1 <- as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,1,])))
  X_tilde_1 <- X_tilde_1[order(pos_cols[,2], pos_cols[,1]),]
  beta_2 <- apply(locs, 1, function(x) sqrt(x[1]^2+x[2]^2))/10
  X_tilde_2 <- as.matrix(Matrix::bdiag(lapply(new_z_2, function(x) x[,1,])))
  X_tilde_2 <- X_tilde_2[order(pos_cols[,2], pos_cols[,1]),]
  
  output$coefficients <- output$vcov <- output$mu.fitted.values <- output$mu.se <- output$sigma.fitted.values <- output$sigma.se <- list()
  output$mu.PCP <- output$sigma.PCP <- output$GCV <- vector("numeric", length(lambda_seq))
  time_k_4 <- list()
  
  for(l in 1:length(lambda_seq)){
    lambda <- lambda_seq[l]
    
    penalty <- matrix(0, 2+K*(p_1+p_2)+p_3, 2+K*(p_1+p_2)+p_3)
    diag(penalty)[3:(2+K*p_1+K*p_2)] <- N*lambda
    
    weight <- solve(quad_form + penalty)
    output$coefficients[[l]] <- as.vector(weight %*% 
                                            sensitivity %*% V_inv %*% 
                                            unlist(lapply(1:K, 
                                                          function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
    if(sum(is.na(output$coefficients[[l]]))!=0) break
    output$vcov[[l]] <- weight %*% sensitivity%*%V_inv%*%psi%*%t(psi)%*%V_inv%*%t(sensitivity) %*% weight
    
    output$mu.fitted.values[[l]] <- as.vector(X_tilde_1%*%output$coefficients[[l]][3:(2+K*p_1)])
    output$mu.se[[l]] <- sqrt(diag(X_tilde_1 %*% output$vcov[[l]][3:(2+K*p_1), 3:(2+K*p_1)] %*%t(X_tilde_1)))
    output$mu.PCP[l] <- mean(beta_1 <= output$mu.fitted.values[[l]] + 1.96*output$mu.se[[l]] & 
                               beta_1 >= output$mu.fitted.values[[l]]-1.96*output$mu.se[[l]])
    
    output$sigma.fitted.values[[l]] <- as.vector(X_tilde_2%*%output$coefficients[[l]][(2+K*p_1+1):(2+K*p_1+K*p_2)])
    output$sigma.se[[l]] <- sqrt(diag(X_tilde_2 %*% output$vcov[[l]][(2+K*p_1+1):(2+K*p_1+K*p_2), (2+K*p_1+1):(2+K*p_1+K*p_2)] %*%t(X_tilde_2)))
    output$sigma.PCP[l] <- mean(beta_2 <= output$sigma.fitted.values[[l]] + 1.96*output$sigma.se[[l]] & beta_2 >= output$sigma.fitted.values[[l]]-1.96*output$sigma.se[[l]])
    
    time_k_3 <- list()
    for(k in 1:K){
      time_k_3[[k]] <- proc.time()
      
      result_k <- Chessian_all_raw(c(output$coefficients[[l]][1:2], output$coefficients[[l]][(2+(k-1)*p_1+1):(2+k*(p_1+p_2))], 
                                     output$coefficients[[l]][(2+K*(p_1+p_2)+1):(2+K*p_1+p_2+p_3)]), 
                                   y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]], 
                                   z_3=z_3[indicator==k,,,drop=FALSE], 
                                   locs=locs[indicator==k,])
      
      psi_list[[k]] <- result_k$score
      sensitivity_list[[k]] <- result_k$sensitivity %*% t(result_k$sensitivity)
      
      time_k_3[[k]] <- proc.time() - time_k_3[[k]]
    }
    time_k_4[[l]] <- time_k_3[[which.max(sapply(time_k_3, function(x) x[3]))]]
    psi <- do.call(rbind,psi_list)
    GCV_numerator <- N^(-1)*as.numeric(rowMeans(psi)%*%V_inv%*%rowMeans(psi))
    output$GCV[l] <- GCV_numerator / ( 1 - N^(-1)*sum(diag( weight %*% quad_form )))^2
    
  }
  
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  return(output)
}

CL <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs){
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
  
  if(quantile>=0) {
    thresholds <- apply(Y, 2, function(x) quantile(x, prob=quantile))
  } else {
    thresholds <- rep(min(Y)-1, S)
  }
  frech <- t(apply(y, 1, function(x) SpatialExtremes::gev2frech(x, emp=TRUE)))
  
  opt <- nlminb(start=c(0.75, max(distance(locs))/2.73), objective=logCL_BR_thresh, 
                thresholds=thresholds, x=frech, locs=locs, 
                lower=c(1e-5,1e-5), upper=c(2-1e-5,Inf))
  
  mu <- sigma <- xi <- rep(0, s_k)
  for (i in 1:S) {
    y_i <- y[,i]
    marg.param <- gevmle(y_i[y_i>quantile(y_i, prob=0.5)])
    mu[i] <- marg.param["loc"]
    sigma[i] <- marg.param["scale"]
    xi[i] <- marg.param["shape"]
  }
  x_1 <- covariates_1[c(t(sapply(1:S, function(x) x+seq(0,(N-1)*S,S)))),]
  x_2 <- covariates_2[c(t(sapply(1:S, function(x) x+seq(0,(N-1)*S,S)))),]
  x_3 <- covariates_3[c(t(sapply(1:S, function(x) x+seq(0,(N-1)*S,S)))),]
  
  init_values <- as.vector(c(opt$par, 
                             solve(t(x_1)%*%x_1)%*%t(x_1)%*%rep(mu,N), 
                             solve(t(x_2)%*%x_2)%*%t(x_2)%*%rep(log(sigma),N), 
                             solve(t(x_3)%*%x_3)%*%t(x_3)%*%rep(xi,N)))
  
  opt_2 <-
    nlminb(start=init_values, objective=logCL_all_thresh, thresholds=thresholds,
           y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs,
           lower=c(1e-5,1e-5,rep(-Inf,p)), upper=c(2-1e-5,rep(Inf,1+p)))
  
  output$coefficients <- opt_2$par
  
  result <- Chessian_all_thresh(output$coefficients, thresholds=thresholds, y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs)
  psi <- result$score
  sensitivity <- result$sensitivity %*% t(result$sensitivity)
  output$vcov <- solve(sensitivity)%*%psi%*%t(psi)%*%solve(sensitivity)
  
  output$time <- proc.time()-time
  
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
  
  if(quantile>=0) {
    thresholds <- apply(x, 2, function(y) quantile(y, prob=quantile))
  } else {
    thresholds <- rep(min(x)-1, S)
  }
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    cov.start <- SpatialExtremes::fitcovariance(x[,indicator==k], locs[indicator==k,], "brown", marge = "emp")$param
    
    init_values <- as.vector(cov.start[2:1])
    
    opt_2 <- nlminb(start=init_values, objective=logCL_BR_thresh, thresholds=thresholds, x=x[,indicator==k], 
                    locs=locs[indicator==k,], lower=c(1e-5,1e-5), upper=c(2-1e-5,Inf))
    
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
