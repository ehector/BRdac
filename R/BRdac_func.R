BRdacObjective <- function(output, y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator){
  time <- proc.time()
  
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
  
  psi_list <- list()
  
  if(quantile>=0) {
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y, na.rm=T)-1, S)
  }
  
  time_k_2 <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_thresh_reparm(output$coefficients, thresholds=thresholds[indicator==k], y=y[,indicator==k],
                                    z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE],
                                    z_3=z_3[indicator==k,,,drop=FALSE],
                                    locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix/N
  score <- rowSums(psi, na.rm=T)
  output$Q <- t(score) %*% V_inv %*% score
  output$df <- qr(V_inv)$rank - length(output$coefficients)
  output$p.value <- 1 - as.numeric(pchisq(output$Q, output$df))
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_after + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  return(output)
}

invBRdacObjective <- function(output, y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator){
  time <- proc.time()
  
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
  
  psi_list <- list()
  
  if(quantile>=0) {
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y, na.rm=T)-1, S)
  }
  
  time_k_2 <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- iMSP_Chessian_all_thresh_reparm(output$coefficients, thresholds=thresholds[indicator==k], y=y[,indicator==k],
                                           z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE],
                                           z_3=z_3[indicator==k,,,drop=FALSE],
                                           locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix/N
  score <- rowSums(psi, na.rm=T)
  output$Q <- t(score) %*% V_inv %*% score
  output$df <- qr(V_inv)$rank - length(output$coefficients)
  output$p.value <- 1 - as.numeric(pchisq(output$Q, output$df))
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_after + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  return(output)
}

BRdac <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator, cluster=1){
  if(cluster==1) output <- BRdac_sequential(y, covariates_1, covariates_2, covariates_3, quantile, locs, indicator)
  if(cluster>1) output <- BRdac_parallel(y, covariates_1, covariates_2, covariates_3, quantile, locs, indicator)
  return(output)
}

BRdac_sequential <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator){
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
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y, na.rm=T)-1, S)
  }
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    mu[i] <- opt_1$par[1]
    sigma[i] <- exp(opt_1$par[2])
    xi[i] <- opt_1$par[3]
  }
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(init_cov,
                             solve(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             solve(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             solve(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  
  y_long <- c(t(y))
  thresh_long <- rep(thresholds,N)
  scaled <- unlist(sapply(1:(N*S), function(i) {
    if(!is.na(y_long[i])) {
      if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu[i])*xi[i]/sigma[i])
      else return((thresh_long[i] - mu[i])*xi[i]/sigma[i])
    }
  }))
  fudge <- which(scaled < -1)
  
  if(length(fudge)>0){
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
    
    opt_fix <- optim(par=init_values[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds, 
                     y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS", control=list(maxit=5000)) 
    
    init_values <- c(init_values[1:2], opt_fix$par)  
  }
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    opt_2 <- optim(par=init_values, fn=logCL_all_thresh_reparm, gr=score_all_thresh_reparm, thresholds=thresholds[indicator==k],
                   y=y[,indicator==k,drop=FALSE], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
                   method="BFGS")
    
    estimates[k,] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(estimates)
  time_mean <- proc.time() - time_mean
  
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_thresh_reparm(mean_estimates, thresholds=thresholds[indicator==k], y=y[,indicator==k],
                                           z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE],
                                           z_3=z_3[indicator==k,,,drop=FALSE],
                                           locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity)
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  weight <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(1:K, 
                                                   function(k) {
                                                     sensitivity_list[[k]] %*% estimates[k,]
                                                   }))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%V%*%V_inv%*%t(sensitivity) %*% weight
  
  output$coefficients_reparm <- c(2*exp(output$coefficients[1])/(1+exp(output$coefficients[1])), exp(output$coefficients[2]),
                                  output$coefficients[-c(1:2)])
  output$var_reparm <- diag(output$vcov)
  gradient <- diag(c(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])), output$coefficients_reparm[2], rep(1,p) ))
  output$var_reparm[1] <- output$var_reparm[1]*(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])))^2
  output$var_reparm[2] <- output$var_reparm[2]*(output$coefficients_reparm[2])^2
  output$vcov_reparm <- gradient %*% output$vcov %*% gradient
  
  names(output$coefficients) <- rownames(output$vcov) <- colnames(output$vcov) <- 
    c("omega","xi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  names(output$coefficients_reparm) <- names(output$var_reparm) <- rownames(output$vcov_reparm) <- colnames(output$vcov_reparm) <-
    c("alpha","phi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  output$model <- "MSP"
  return(output)
}

BRdac_parallel <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator, cluster){
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
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y, na.rm=T)-1, S)
  }
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  sock <- parallel::makeCluster(rep("localhost", cluster), type = "SOCK")
  doParallel::registerDoParallel(sock)
  
  marg_pars_list <- foreach::foreach(i=1:S) %dopar% {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    return(c(opt_1$par[1],exp(opt_1$par[2]),opt_1$par[3]))
  }
  mu <- sapply(marg_pars_list, function(x) x[1])
  sigma <- sapply(marg_pars_list, function(x) x[2])
  xi <- sapply(marg_pars_list, function(x) x[3])
  
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(init_cov,
                             solve(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             solve(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             solve(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  
  y_long <- c(t(y))
  thresh_long <- rep(thresholds,N)
  scaled <- unlist(sapply(1:(N*S), function(i) {
    if(!is.na(y_long[i])) {
      if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu[i])*xi[i]/sigma[i])
      else return((thresh_long[i] - mu[i])*xi[i]/sigma[i])
    }
  }))
  fudge <- which(scaled < -1)
  
  if(length(fudge)>0){
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
    
    opt_fix <- optim(par=init_values[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds, 
                     y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS", control=list(maxit=5000)) 
    
    init_values <- c(init_values[1:2], opt_fix$par)  
  }
  
  time_before <- proc.time() - time
  
  result_k <- foreach::foreach(k=1:K) %dopar% {
    time_k <- proc.time()
    
    opt_2 <- optim(par=init_values, fn=logCL_all_thresh_reparm, gr=score_all_thresh_reparm, thresholds=thresholds[indicator==k],
                   y=y[,indicator==k,drop=FALSE], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
                   method="BFGS")
    return(c(opt_2$par, proc.time() - time_k))
  }
  estimates <- do.call(rbind, lapply(result_k, function(x) x[1:(p+2)]))
  time_k <- lapply(result_k, function(x) x[-c(1:(p+2))])
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(estimates)
  time_mean <- proc.time() - time_mean
  
  weights_list <- foreach::foreach(k=1:K) %dopar% {
    time_k_2 <- proc.time()
    
    result_k <- Chessian_all_thresh_reparm(mean_estimates, thresholds=thresholds[indicator==k], y=y[,indicator==k],
                                           z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE],
                                           z_3=z_3[indicator==k,,,drop=FALSE],
                                           locs=locs[indicator==k,])
    return(list(result_k$score, var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity), proc.time() - time_k_2))
  }
  psi_list <- lapply(weights_list, function(x) x[[1]])
  sensitivity_list <- lapply(weights_list, function(x) x[[2]])
  time_k_2 <- lapply(weights_list, function(x) x[[3]])
  
  parallel::stopCluster(sock)
  foreach::registerDoSEQ()
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  weight <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(1:K, 
                                                   function(k) {
                                                     sensitivity_list[[k]] %*% estimates[k,]
                                                   }))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%V%*%V_inv%*%t(sensitivity) %*% weight
  
  output$coefficients_reparm <- c(2*exp(output$coefficients[1])/(1+exp(output$coefficients[1])), exp(output$coefficients[2]),
                                  output$coefficients[-c(1:2)])
  output$var_reparm <- diag(output$vcov)
  gradient <- diag(c(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])), output$coefficients_reparm[2], rep(1,p) ))
  output$var_reparm[1] <- output$var_reparm[1]*(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])))^2
  output$var_reparm[2] <- output$var_reparm[2]*(output$coefficients_reparm[2])^2
  output$vcov_reparm <- gradient %*% output$vcov %*% gradient
  
  names(output$coefficients) <- rownames(output$vcov) <- colnames(output$vcov) <- 
    c("omega","xi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  names(output$coefficients_reparm) <- names(output$var_reparm) <- rownames(output$vcov_reparm) <- colnames(output$vcov_reparm) <-
    c("alpha","phi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  output$model <- "MSP"
  return(output)
}

invBRdac <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator, cluster=1){
  if(cluster==1) output <- invBRdac_sequential(y, covariates_1, covariates_2, covariates_3, quantile, locs, indicator)
  if(cluster>1) output <- invBRdac_parallel(y, covariates_1, covariates_2, covariates_3, quantile, locs, indicator)
  return(output)
}

invBRdac_sequential <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator){
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
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y, na.rm=T)-1, S)
  }
  
  unit_y <- -1 / log(apply(y, 2, rank) / (N  +1))
  uniform_y <- exp(-1/unit_y)
  exp_y <- -log(1-uniform_y)
  MSP_y <- 1/exp_y
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(MSP_y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    mu[i] <- opt_1$par[1]
    sigma[i] <- exp(opt_1$par[2])
    xi[i] <- opt_1$par[3]
  }
  
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(init_cov,
                             solve(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             solve(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             solve(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  
  y_long <- c(t(y))
  thresh_long <- rep(thresholds,N)
  scaled <- unlist(sapply(1:(N*S), function(i) {
    if(!is.na(y_long[i])) {
      if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu[i])*xi[i]/sigma[i])
      else return((thresh_long[i] - mu[i])*xi[i]/sigma[i])
    }
  }))
  fudge <- which(scaled < -1)
  
  if(length(fudge)>0){
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
    
    opt_fix <- optim(par=init_values[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds, 
                     y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS", control=list(maxit=5000)) 
    
    init_values <- c(init_values[1:2], opt_fix$par)  
  }
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    opt_2 <- optim(par=init_values, fn=iMSP_logCL_all_thresh_reparm, gr=iMSP_score_all_thresh_reparm, thresholds=thresholds[indicator==k],
                   y=y[,indicator==k,drop=FALSE], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
                   method="BFGS")
    
    estimates[k,] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(estimates)
  time_mean <- proc.time() - time_mean
  
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- iMSP_Chessian_all_thresh_reparm(mean_estimates, thresholds=thresholds[indicator==k], y=y[,indicator==k],
                                           z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE],
                                           z_3=z_3[indicator==k,,,drop=FALSE],
                                           locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity)
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  weight <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(1:K, 
                                                   function(k) {
                                                     sensitivity_list[[k]] %*% estimates[k,]
                                                   }))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%V%*%V_inv%*%t(sensitivity) %*% weight
  
  output$coefficients_reparm <- c(2*exp(output$coefficients[1])/(1+exp(output$coefficients[1])), exp(output$coefficients[2]),
                                  output$coefficients[-c(1:2)])
  output$var_reparm <- diag(output$vcov)
  gradient <- diag(c(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])), output$coefficients_reparm[2], rep(1,p) ))
  output$var_reparm[1] <- output$var_reparm[1]*(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])))^2
  output$var_reparm[2] <- output$var_reparm[2]*(output$coefficients_reparm[2])^2
  output$vcov_reparm <- gradient %*% output$vcov %*% gradient
  
  names(output$coefficients) <- rownames(output$vcov) <- colnames(output$vcov) <- 
    c("omega","xi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  names(output$coefficients_reparm) <- names(output$var_reparm) <- rownames(output$vcov_reparm) <- colnames(output$vcov_reparm) <-
    c("alpha","phi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  output$model <- "inverted MSP"
  return(output)
}

invBRdac_parallel <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator, cluster){
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
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y, na.rm=T)-1, S)
  }
  
  unit_y <- -1 / log(apply(y, 2, rank) / (N  +1))
  uniform_y <- exp(-1/unit_y)
  exp_y <- -log(1-uniform_y)
  MSP_y <- 1/exp_y
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(MSP_y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  sock <- parallel::makeCluster(rep("localhost", cluster), type = "SOCK")
  doParallel::registerDoParallel(sock)
  
  marg_pars_list <- foreach::foreach(i=1:S) %dopar% {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    return(c(opt_1$par[1],exp(opt_1$par[2]),opt_1$par[3]))
  }
  mu <- sapply(marg_pars_list, function(x) x[1])
  sigma <- sapply(marg_pars_list, function(x) x[2])
  xi <- sapply(marg_pars_list, function(x) x[3])
  
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(init_cov,
                             solve(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             solve(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             solve(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  
  y_long <- c(t(y))
  thresh_long <- rep(thresholds,N)
  scaled <- unlist(sapply(1:(N*S), function(i) {
    if(!is.na(y_long[i])) {
      if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu[i])*xi[i]/sigma[i])
      else return((thresh_long[i] - mu[i])*xi[i]/sigma[i])
    }
  }))
  fudge <- which(scaled < -1)
  
  if(length(fudge)>0){
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
    
    opt_fix <- optim(par=init_values[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds, 
                     y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS", control=list(maxit=5000)) 
    
    init_values <- c(init_values[1:2], opt_fix$par)  
  }
  
  time_before <- proc.time() - time
  
  result_k <- foreach::foreach(k=1:K) %dopar% {
    time_k <- proc.time()
    
    opt_2 <- optim(par=init_values, fn=iMSP_logCL_all_thresh_reparm, gr=iMSP_score_all_thresh_reparm, thresholds=thresholds[indicator==k],
                   y=y[,indicator==k,drop=FALSE], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
                   method="BFGS")
    return(c(opt_2$par, proc.time() - time_k))
  }
  estimates <- do.call(rbind, lapply(result_k, function(x) x[1:(p+2)]))
  time_k <- lapply(result_k, function(x) x[-c(1:(p+2))])
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(estimates)
  time_mean <- proc.time() - time_mean
  
  weights_list <- foreach::foreach(k=1:K) %dopar% {
    time_k_2 <- proc.time()
    
    result_k <- iMSP_Chessian_all_thresh_reparm(mean_estimates, thresholds=thresholds[indicator==k], y=y[,indicator==k],
                                           z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE],
                                           z_3=z_3[indicator==k,,,drop=FALSE],
                                           locs=locs[indicator==k,])
    return(list(result_k$score, var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity), proc.time() - time_k_2))
  }
  psi_list <- lapply(weights_list, function(x) x[[1]])
  sensitivity_list <- lapply(weights_list, function(x) x[[2]])
  time_k_2 <- lapply(weights_list, function(x) x[[3]])
  
  parallel::stopCluster(sock)
  foreach::registerDoSEQ()
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  sensitivity <- do.call(cbind, sensitivity_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  weight <- solve(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(1:K, 
                                                   function(k) {
                                                     sensitivity_list[[k]] %*% estimates[k,]
                                                   }))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%V%*%V_inv%*%t(sensitivity) %*% weight
  
  output$coefficients_reparm <- c(2*exp(output$coefficients[1])/(1+exp(output$coefficients[1])), exp(output$coefficients[2]),
                                  output$coefficients[-c(1:2)])
  output$var_reparm <- diag(output$vcov)
  gradient <- diag(c(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])), output$coefficients_reparm[2], rep(1,p) ))
  output$var_reparm[1] <- output$var_reparm[1]*(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])))^2
  output$var_reparm[2] <- output$var_reparm[2]*(output$coefficients_reparm[2])^2
  output$vcov_reparm <- gradient %*% output$vcov %*% gradient
  
  names(output$coefficients) <- rownames(output$vcov) <- colnames(output$vcov) <- 
    c("omega","xi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  names(output$coefficients_reparm) <- names(output$var_reparm) <- rownames(output$vcov_reparm) <- colnames(output$vcov_reparm) <-
    c("alpha","phi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  output$model <- "inverted MSP"
  return(output)
}

BRdacVCMmusigma <- function(y, covariates_3, quantile=0.95, locs, indicator, lambda_1, lambda_2, cluster){
  if(cluster==1) output <- BRdacVCMmusigma_sequential(y, covariates_3, quantile, locs, indicator, lambda_1, lambda_2)
  if(cluster>1) output <- BRdacVCMmusigma_parallel(y, covariates_3, quantile, locs, indicator, lambda_1, lambda_2, cluster)
  return(output)
}

BRdacVCMmusigma_sequential <- function(y, covariates_3, quantile=0.95, locs, indicator, lambda_1, lambda_2){
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
  
  if(quantile>=0) {
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y)-1, S)
  }
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    mu[i] <- opt_1$par[1]
    sigma[i] <- exp(opt_1$par[2])
    xi[i] <- opt_1$par[3]
  }
  
  init_values <- as.vector(c(init_cov, solve(t(covariates_3)%*%covariates_3)%*%t(covariates_3)%*%rep(xi,N)))
  p_1 <- p_2 <- p <- rep(0,K)
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    pos_k <- locs[indicator==k,]
    s_k <- sum(indicator==k)

    knots_mat <- fields::cover.design( pos_k, 10, num.nn=s_k/2, start=floor(seq(1,s_k,length.out=10)) )$design
    attr(knots_mat,"scaled:scale") <- attr(knots_mat,"scaled:center") <- NULL
    colnames(knots_mat) <- NULL
    distances <- apply(knots_mat, 1, function(x) sqrt((pos_k[,1]-x[1])^2 + (pos_k[,2]-x[2])^2))
    basis_1 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.05, 0, 1, t))))
    new_X <- cbind(rep(1,N*s_k), do.call(rbind, lapply(1:N, function(x) pos_k)), do.call(rbind, lapply(1:N, function(x) basis_1)))
    
    p_1[k] <- p_2[k] <- ncol(new_X)
    new_z_1[[k]] <- new_z_2[[k]] <- z_constructor(new_X, s_k, N, p_1[k])
    p[k] <- p_1[k] + p_2[k] + p_3
    
    init_par <- c(init_values[1:2],
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(mu[indicator==k],N),
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(log(sigma[indicator==k]),N),
                  init_values[-c(1:2)])
    
    mu_sub <- new_X %*% init_par[3:(2+p_1[k])]
    sigma_sub <- exp(new_X %*% init_par[(3+p_1[k]):(2+p_1[k]+p_2[k])])
    xi_sub <- covariates_3[rep(indicator==k,N),,drop=FALSE] %*% init_par[(3+p_1[k]+p_2[k]):(2+p[k])]
    
    y_long <- c(t(y[,indicator==k,drop=FALSE]))
    thresh_long <- rep(thresholds[indicator==k],N)
    scaled <- unlist(sapply(1:(N*s_k), function(i) {
      if(!is.na(y_long[i])) {
        if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
        else return((thresh_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
      }
    }))
    fudge <- which(scaled < -1)
    
    if(length(fudge)>0){
      c <- max(1, -scaled[fudge]) + 1e-5
      init_par[(3+p_1[k]+p_2[k]):(2+p[k])] <- init_par[(3+p_1[k]+p_2[k]):(2+p[k])]/c
      
      opt <- optim(par=init_par[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[indicator==k], 
                   y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000)) 
      init_par <- c(init_par[1:2], opt$par)  
    }
  
    opt_2 <-
      optim(par=init_par,
            fn=logCL_all_thresh_reparm, gr=score_all_thresh_reparm,
            thresholds=thresholds[indicator==k],
            y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
            z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000))  
    estimates[[k]] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(do.call(rbind, lapply(1:K, function(k) estimates[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),drop=FALSE])))
  time_mean <- proc.time() - time_mean
  
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_thresh_reparm(c(mean_estimates[1:2], estimates[[k]][-c(1,2,length(estimates[[k]]))], mean_estimates[-c(1:2)]), 
                                    thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], 
                                    z_2=new_z_2[[k]], z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity)
    
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  
  sensitivity <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
  sensitivity[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    do.call(cbind, lapply(1:K, function(k) sensitivity_list[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
  sensitivity[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
  
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  cum_p <- cumsum(c(0,p[-K]))
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
  V_inv <- matrix
  
  quad_form <- sensitivity %*% V_inv %*% t(sensitivity)
  
  pos_cols <- do.call(rbind, lapply(1:K, function(k) locs[indicator==k,]))
  order <- apply(locs, 1, function(x) which(equiv(pos_cols[,1],x[1]) & equiv(pos_cols[,2], x[2])))
  X_tilde_1 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,n,])))[order,] ))
  X_tilde_2 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_2, function(x) x[,n,])))[order,] ))
  
  lambda_seq <- as.matrix(expand.grid(lambda_1, lambda_2))
  output$coefficients <- output$coefficients_reparm <- output$vcov <- output$vcov_reparm <- output$var_reparm <- 
    output$mu.fitted.values <- output$mu.se <- output$sigma.fitted.values <- output$sigma.se <- list()
  output$GCV <- vector("numeric", nrow(lambda_seq))
  cum_p12 <- cumsum(c(0,colSums(rbind(p_1[-K],p_2[-K]))))
  time_after <- proc.time()-time_after
  
  time_k_4 <- list()
  time_GCV <- rep(0,5)
  print("Beginning tuning.", quote=FALSE)
  for(l in 1:nrow(lambda_seq)){
    time_l1 <- proc.time()
    lam_1 <- as.numeric(lambda_seq[l,1])
    lam_2 <- as.numeric(lambda_seq[l,2])
    
    penalty <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, 2+sum(p_1)+sum(p_2)+p_3)
    diag(penalty)[3:(2+sum(p_1))] <- N*lam_1
    diag(penalty)[(2+sum(p_1)):(2+sum(p_1)+sum(p_2))] <- N*lam_2
    
    weight <- solve(quad_form + penalty)
    
    output$coefficients[[l]] <- as.vector(weight %*%
                                            sensitivity %*% V_inv %*%
                                            unlist(lapply(1:K,
                                                          function(x) sensitivity_list[[x]] %*% estimates[[x]])))
    if(sum(is.na(output$coefficients[[l]]))!=0) break
    output$vcov[[l]] <- weight %*% sensitivity%*%V_inv%*% V %*%V_inv%*%t(sensitivity) %*% weight
    
    output$mu.fitted.values[[l]] <- as.vector(X_tilde_1%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))])
    output$mu.se[[l]] <- sqrt(apply(X_tilde_1, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))] %*%x
    }))
    
    output$sigma.fitted.values[[l]] <- as.vector(X_tilde_2%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))])
    output$sigma.se[[l]] <- sqrt(apply(X_tilde_2, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))] %*%x
    }))
    
    if(sum(is.na(output$sigma.se[[l]])) + sum(is.na(output$mu.se[[l]]))!=0) next
    
    output$coefficients_reparm[[l]] <- c(2*exp(output$coefficients[[l]][1])/(1+exp(output$coefficients[[l]][1])), exp(output$coefficients[[l]][2]),
                                         output$coefficients[[l]][-c(1:2)])
    output$var_reparm[[l]] <- diag(output$vcov[[l]])
    gradient <- diag(c(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])), output$coefficients_reparm[[l]][2], rep(1,length(output$coefficients[[l]])-2) ))
    output$var_reparm[[l]][1] <- output$var_reparm[[l]][1]*(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])))^2
    output$var_reparm[[l]][2] <- output$var_reparm[[l]][2]*(output$coefficients_reparm[[l]][2])^2
    output$vcov_reparm[[l]] <- gradient %*% output$vcov[[l]] %*% gradient
    
    names(output$coefficients[[l]]) <- rownames(output$vcov[[l]]) <- colnames(output$vcov[[l]]) <- 
      c("omega","xi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    names(output$coefficients_reparm[[l]]) <- names(output$var_reparm[[l]]) <- rownames(output$vcov_reparm[[l]]) <- colnames(output$vcov_reparm[[l]]) <-
      c("alpha","phi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    
    psi_list_GCV <- list()
    sensitivity_list_GCV <- list()
    time_l1 <- proc.time() - time_l1
    
    time_k_3 <- list()
    for(k in 1:K){
      print(k)
      time_k_3[[k]] <- proc.time()
      
      result_k <- Chessian_all_thresh_reparm(c(output$coefficients[[l]][1:2], output$coefficients[[l]][2+(cum_p12[k]+1):(cum_p12[k]+p_1[k])],
                                               output$coefficients[[l]][2+(cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])], 
                                               output$coefficients[[l]][length(output$coefficients[[l]])]),
                                             thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                                             z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
      
      psi_list_GCV[[k]] <- result_k$score
      sensitivity_list_GCV[[k]] <- var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity)
      
      time_k_3[[k]] <- proc.time() - time_k_3[[k]]
    }
    time_k_4[[l]] <- time_k_3[[which.max(sapply(time_k_3, function(x) x[3]))]]
    
    time_l2 <- proc.time()
    
    sensitivity_GCV <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
    sensitivity_GCV[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      do.call(cbind, lapply(1:K, function(k) sensitivity_list_GCV[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
    sensitivity_GCV[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list_GCV[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
    
    psi_GCV <- do.call(rbind,psi_list_GCV)
    V_GCV <- var(t(psi_GCV), na.rm=T)*N
    V_inv_GCV <- solve(V_GCV)
    matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv_GCV[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
    V_inv_GCV <- matrix
    quad_form_GCV <- sensitivity_GCV %*% V_inv_GCV %*% t(sensitivity_GCV)
    
    weight_GCV <- solve(quad_form_GCV + penalty)
    
    GCV_numerator <- N^(-1)*as.numeric(rowMeans(psi_GCV)%*%V_inv_GCV%*%rowMeans(psi_GCV))
    output$GCV[l] <- GCV_numerator / ( 1 - N^(-1)*sum(diag( weight_GCV %*% quad_form_GCV )))^2
    time_l2 <- proc.time() - time_l2
    time_GCV <- time_GCV + time_l1 + time_l2
    time_GCV <- time_GCV + time_l1
  }
  
  total_time_k_4 <- do.call(rbind, time_k_4)
  if(!is.null(total_time_k_4)) {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after +
      colSums(do.call(rbind, time_k_4)) + time_GCV
  } else {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after 
  }
  output$time_noGCV <- time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after
  output$model <- "MSP"
  return(output)
}

BRdacVCMmusigma_parallel <- function(y, covariates_3, quantile=0.95, locs, indicator, lambda_1, lambda_2, cluster){
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
  
  if(quantile>=0) {
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y)-1, S)
  }
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  sock <- parallel::makeCluster(rep("localhost", cluster), type = "SOCK")
  doParallel::registerDoParallel(sock)
  
  marg_pars_list <- foreach::foreach(i=1:S) %dopar% {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    return(c(opt_1$par[1],exp(opt_1$par[2]),opt_1$par[3]))
  }
  mu <- sapply(marg_pars_list, function(x) x[1])
  sigma <- sapply(marg_pars_list, function(x) x[2])
  xi <- sapply(marg_pars_list, function(x) x[3])
  
  init_values <- as.vector(c(init_cov, solve(t(covariates_3)%*%covariates_3)%*%t(covariates_3)%*%rep(xi,N)))

  time_before <- proc.time() - time
  
  result_k <- foreach::foreach(k=1:K) %dopar% {
    time_k <- proc.time()
    
    pos_k <- locs[indicator==k,]
    s_k <- sum(indicator==k)
    
    knots_mat <- fields::cover.design( pos_k, 10, num.nn=s_k/2, start=floor(seq(1,s_k,length.out=10)) )$design
    attr(knots_mat,"scaled:scale") <- attr(knots_mat,"scaled:center") <- NULL
    colnames(knots_mat) <- NULL
    distances <- apply(knots_mat, 1, function(x) sqrt((pos_k[,1]-x[1])^2 + (pos_k[,2]-x[2])^2))
    basis_1 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.05, 0, 1, t))))
    new_X <- cbind(rep(1,N*s_k), do.call(rbind, lapply(1:N, function(x) pos_k)), do.call(rbind, lapply(1:N, function(x) basis_1)))
    
    p_1 <- p_2 <- ncol(new_X)
    new_z_1[[k]] <- new_z_2[[k]] <- z_constructor(new_X, s_k, N, p_1)
    p <- p_1 + p_2 + p_3
    
    init_par <- c(init_values[1:2],
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(mu[indicator==k],N),
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(log(sigma[indicator==k]),N),
                  init_values[-c(1:2)])
    
    mu_sub <- new_X %*% init_par[3:(2+p_1)]
    sigma_sub <- exp(new_X %*% init_par[(3+p_1):(2+p_1+p_2)])
    xi_sub <- covariates_3[rep(indicator==k,N),,drop=FALSE] %*% init_par[(3+p_1+p_2):(2+p)]
    
    y_long <- c(t(y[,indicator==k,drop=FALSE]))
    thresh_long <- rep(thresholds[indicator==k],N)
    scaled <- unlist(sapply(1:(N*s_k), function(i) {
      if(!is.na(y_long[i])) {
        if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
        else return((thresh_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
      }
    }))
    fudge <- which(scaled < -1)
    
    if(length(fudge)>0){
      c <- max(1, -scaled[fudge]) + 1e-5
      init_par[(3+p_1+p_2):(2+p)] <- init_par[(3+p_1+p_2):(2+p)]/c
      
      opt <- optim(par=init_par[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[indicator==k], 
                   y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000)) 
      init_par <- c(init_par[1:2], opt$par)  
    }
    
    opt_2 <-
      optim(par=init_par,
            fn=logCL_all_thresh_reparm, gr=score_all_thresh_reparm,
            thresholds=thresholds[indicator==k],
            y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
            z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000))  
    
    return(list(opt_2$par, proc.time() - time_k, new_z_1[[k]], new_z_2[[k]], p_1, p_2, p))
  }
  estimates <- lapply(result_k, function(x) x[[1]])
  new_z_1 <- lapply(result_k, function(x) x[[3]])
  new_z_2 <- lapply(result_k, function(x) x[[4]])
  p_1 <- sapply(result_k, function(x) x[[5]])
  p_2 <- sapply(result_k, function(x) x[[6]])
  p <- sapply(result_k, function(x) x[[7]])
  time_k <- lapply(result_k, function(x) x[[2]])
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(do.call(rbind, lapply(1:K, function(k) estimates[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),drop=FALSE])))
  time_mean <- proc.time() - time_mean
  
  weights_list <- foreach::foreach(k=1:K) %dopar% {
    time_k_2 <- proc.time()
    
    result_k <- Chessian_all_thresh_reparm(c(mean_estimates[1:2], estimates[[k]][-c(1,2,length(estimates[[k]]))], mean_estimates[-c(1:2)]), 
                                           thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], 
                                           z_2=new_z_2[[k]], z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
    return(list(result_k$score, var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity), proc.time() - time_k_2))
  }
  psi_list <- lapply(weights_list, function(x) x[[1]])
  sensitivity_list <- lapply(weights_list, function(x) x[[2]])
  time_k_2 <- lapply(weights_list, function(x) x[[3]])
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  
  sensitivity <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
  sensitivity[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    do.call(cbind, lapply(1:K, function(k) sensitivity_list[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
  sensitivity[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
  
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  cum_p <- cumsum(c(0,p[-K]))
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
  V_inv <- matrix
  
  quad_form <- sensitivity %*% V_inv %*% t(sensitivity)
  
  pos_cols <- do.call(rbind, lapply(1:K, function(k) locs[indicator==k,]))
  order <- apply(locs, 1, function(x) which(equiv(pos_cols[,1],x[1]) & equiv(pos_cols[,2], x[2])))
  X_tilde_1 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,n,])))[order,] ))
  X_tilde_2 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_2, function(x) x[,n,])))[order,] ))
  
  lambda_seq <- as.matrix(expand.grid(lambda_1, lambda_2))
  output$coefficients <- output$coefficients_reparm <- output$vcov <- output$vcov_reparm <- output$var_reparm <- 
    output$mu.fitted.values <- output$mu.se <- output$sigma.fitted.values <- output$sigma.se <- list()
  output$GCV <- vector("numeric", nrow(lambda_seq))
  cum_p12 <- cumsum(c(0,colSums(rbind(p_1[-K],p_2[-K]))))
  time_after <- proc.time()-time_after
  
  time_k_4 <- list()
  time_GCV <- rep(0,5)
  print("Beginning tuning.", quote=FALSE)
  for(l in 1:nrow(lambda_seq)){
    time_l1 <- proc.time()
    lam_1 <- as.numeric(lambda_seq[l,1])
    lam_2 <- as.numeric(lambda_seq[l,2])
    
    penalty <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, 2+sum(p_1)+sum(p_2)+p_3)
    diag(penalty)[3:(2+sum(p_1))] <- N*lam_1
    diag(penalty)[(2+sum(p_1)):(2+sum(p_1)+sum(p_2))] <- N*lam_2
    
    weight <- solve(quad_form + penalty)
    
    output$coefficients[[l]] <- as.vector(weight %*%
                                            sensitivity %*% V_inv %*%
                                            unlist(lapply(1:K,
                                                          function(x) sensitivity_list[[x]] %*% estimates[[x]])))
    if(sum(is.na(output$coefficients[[l]]))!=0) break
    output$vcov[[l]] <- weight %*% sensitivity%*%V_inv%*% V %*%V_inv%*%t(sensitivity) %*% weight
    
    output$mu.fitted.values[[l]] <- as.vector(X_tilde_1%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))])
    output$mu.se[[l]] <- sqrt(apply(X_tilde_1, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))] %*%x
    }))
    
    output$sigma.fitted.values[[l]] <- as.vector(X_tilde_2%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))])
    output$sigma.se[[l]] <- sqrt(apply(X_tilde_2, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))] %*%x
    }))
    
    if(sum(is.na(output$sigma.se[[l]])) + sum(is.na(output$mu.se[[l]]))!=0) next
    
    output$coefficients_reparm[[l]] <- c(2*exp(output$coefficients[[l]][1])/(1+exp(output$coefficients[[l]][1])), exp(output$coefficients[[l]][2]),
                                         output$coefficients[[l]][-c(1:2)])
    output$var_reparm[[l]] <- diag(output$vcov[[l]])
    gradient <- diag(c(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])), output$coefficients_reparm[[l]][2], rep(1,length(output$coefficients[[l]])-2) ))
    output$var_reparm[[l]][1] <- output$var_reparm[[l]][1]*(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])))^2
    output$var_reparm[[l]][2] <- output$var_reparm[[l]][2]*(output$coefficients_reparm[[l]][2])^2
    output$vcov_reparm[[l]] <- gradient %*% output$vcov[[l]] %*% gradient
    
    names(output$coefficients[[l]]) <- rownames(output$vcov[[l]]) <- colnames(output$vcov[[l]]) <- 
      c("omega","xi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    names(output$coefficients_reparm[[l]]) <- names(output$var_reparm[[l]]) <- rownames(output$vcov_reparm[[l]]) <- colnames(output$vcov_reparm[[l]]) <-
      c("alpha","phi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    
    psi_list_GCV <- list()
    sensitivity_list_GCV <- list()
    time_l1 <- proc.time() - time_l1
    
    weights_list <- foreach::foreach(k=1:K) %dopar% {
      time_k_3 <- proc.time()
      
      result_k <- Chessian_all_thresh_reparm(c(output$coefficients[[l]][1:2], output$coefficients[[l]][2+(cum_p12[k]+1):(cum_p12[k]+p_1[k])],
                                               output$coefficients[[l]][2+(cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])], 
                                               output$coefficients[[l]][length(output$coefficients[[l]])]),
                                             thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                                             z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
      
      return(list(result_k$score, var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity), proc.time() - time_k_3))
    }
    psi_list_GCV <- lapply(weights_list, function(x) x[[1]])
    sensitivity_list_GCV <- lapply(weights_list, function(x) x[[2]])
    time_k_3 <- lapply(weights_list, function(x) x[[3]])
    
    time_k_4[[l]] <- time_k_3[[which.max(sapply(time_k_3, function(x) x[3]))]]
    
    time_l2 <- proc.time()
    
    sensitivity_GCV <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
    sensitivity_GCV[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      do.call(cbind, lapply(1:K, function(k) sensitivity_list_GCV[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
    sensitivity_GCV[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list_GCV[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
    
    psi_GCV <- do.call(rbind,psi_list_GCV)
    V_GCV <- var(t(psi_GCV), na.rm=T)*N
    V_inv_GCV <- solve(V_GCV)
    matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv_GCV[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
    V_inv_GCV <- matrix
    quad_form_GCV <- sensitivity_GCV %*% V_inv_GCV %*% t(sensitivity_GCV)
    
    weight_GCV <- solve(quad_form_GCV + penalty)
    
    GCV_numerator <- N^(-1)*as.numeric(rowMeans(psi_GCV)%*%V_inv_GCV%*%rowMeans(psi_GCV))
    output$GCV[l] <- GCV_numerator / ( 1 - N^(-1)*sum(diag( weight_GCV %*% quad_form_GCV )))^2
    time_l2 <- proc.time() - time_l2
    time_GCV <- time_GCV + time_l1 + time_l2
    time_GCV <- time_GCV + time_l1
  }
  parallel::stopCluster(sock)
  foreach::registerDoSEQ()
  total_time_k_4 <- do.call(rbind, time_k_4)
  if(!is.null(total_time_k_4)) {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after +
      colSums(do.call(rbind, time_k_4)) + time_GCV
  } else {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after 
  }
  output$time_noGCV <- time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after
  output$model <- "MSP"
  return(output)
}

invBRdacVCMmusigma <- function(y, covariates_3, quantile=0.95, locs, indicator, lambda_1, lambda_2, cluster){
  if(cluster==1) output <- invBRdacVCMmusigma_sequential(y, covariates_3, quantile, locs, indicator, lambda_1, lambda_2)
  if(cluster>1) output <- invBRdacVCMmusigma_parallel(y, covariates_3, quantile, locs, indicator, lambda_1, lambda_2, cluster)
  return(output)
}

invBRdacVCMmusigma_sequential <- function(y, covariates_3, quantile=0.95, locs, indicator, lambda_1, lambda_2){
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
  
  if(quantile>=0) {
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y)-1, S)
  }
  
  unit_y <- -1 / log(apply(y, 2, rank) / (N  +1))
  uniform_y <- exp(-1/unit_y)
  exp_y <- -log(1-uniform_y)
  MSP_y <- 1/exp_y
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(MSP_y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    mu[i] <- opt_1$par[1]
    sigma[i] <- exp(opt_1$par[2])
    xi[i] <- opt_1$par[3]
  }
  
  init_values <- as.vector(c(init_cov, solve(t(covariates_3)%*%covariates_3)%*%t(covariates_3)%*%rep(xi,N)))
  p_1 <- p_2 <- p <- rep(0,K)
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    pos_k <- locs[indicator==k,]
    s_k <- sum(indicator==k)
    
    knots_mat <- fields::cover.design( pos_k, 10, num.nn=s_k/2, start=floor(seq(1,s_k,length.out=10)) )$design
    attr(knots_mat,"scaled:scale") <- attr(knots_mat,"scaled:center") <- NULL
    colnames(knots_mat) <- NULL
    distances <- apply(knots_mat, 1, function(x) sqrt((pos_k[,1]-x[1])^2 + (pos_k[,2]-x[2])^2))
    basis_1 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.05, 0, 1, t))))
    new_X <- cbind(rep(1,N*s_k), do.call(rbind, lapply(1:N, function(x) pos_k)), do.call(rbind, lapply(1:N, function(x) basis_1)))
    
    p_1[k] <- p_2[k] <- ncol(new_X)
    new_z_1[[k]] <- new_z_2[[k]] <- z_constructor(new_X, s_k, N, p_1[k])
    p[k] <- p_1[k] + p_2[k] + p_3
    
    init_par <- c(init_values[1:2],
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(mu[indicator==k],N),
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(log(sigma[indicator==k]),N),
                  init_values[-c(1:2)])
    
    mu_sub <- new_X %*% init_par[3:(2+p_1[k])]
    sigma_sub <- exp(new_X %*% init_par[(3+p_1[k]):(2+p_1[k]+p_2[k])])
    xi_sub <- covariates_3[rep(indicator==k,N),,drop=FALSE] %*% init_par[(3+p_1[k]+p_2[k]):(2+p[k])]
    
    y_long <- c(t(y[,indicator==k,drop=FALSE]))
    thresh_long <- rep(thresholds[indicator==k],N)
    scaled <- unlist(sapply(1:(N*s_k), function(i) {
      if(!is.na(y_long[i])) {
        if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
        else return((thresh_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
      }
    }))
    fudge <- which(scaled < -1)
    
    if(length(fudge)>0){
      c <- max(1, -scaled[fudge]) + 1e-5
      init_par[(3+p_1[k]+p_2[k]):(2+p[k])] <- init_par[(3+p_1[k]+p_2[k]):(2+p[k])]/c
      
      opt <- optim(par=init_par[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[indicator==k], 
                   y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000)) 
      init_par <- c(init_par[1:2], opt$par)  
    }
    
    opt_2 <-
      optim(par=init_par,
            fn=iMSP_logCL_all_thresh_reparm, gr=iMSP_score_all_thresh_reparm,
            thresholds=thresholds[indicator==k],
            y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
            z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000))  
    estimates[[k]] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(do.call(rbind, lapply(1:K, function(k) estimates[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),drop=FALSE])))
  time_mean <- proc.time() - time_mean
  
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- iMSP_Chessian_all_thresh_reparm(c(mean_estimates[1:2], estimates[[k]][-c(1,2,length(estimates[[k]]))], mean_estimates[-c(1:2)]), 
                                           thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], 
                                           z_2=new_z_2[[k]], z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity)
    
    time_k_2[[k]] <- proc.time() - time_k_2[[k]]
  }
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  
  sensitivity <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
  sensitivity[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    do.call(cbind, lapply(1:K, function(k) sensitivity_list[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
  sensitivity[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
  
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  cum_p <- cumsum(c(0,p[-K]))
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
  V_inv <- matrix
  
  quad_form <- sensitivity %*% V_inv %*% t(sensitivity)
  
  pos_cols <- do.call(rbind, lapply(1:K, function(k) locs[indicator==k,]))
  order <- apply(locs, 1, function(x) which(equiv(pos_cols[,1],x[1]) & equiv(pos_cols[,2], x[2])))
  X_tilde_1 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,n,])))[order,] ))
  X_tilde_2 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_2, function(x) x[,n,])))[order,] ))
  
  lambda_seq <- as.matrix(expand.grid(lambda_1, lambda_2))
  output$coefficients <- output$coefficients_reparm <- output$vcov <- output$vcov_reparm <- output$var_reparm <- 
    output$mu.fitted.values <- output$mu.se <- output$sigma.fitted.values <- output$sigma.se <- list()
  output$GCV <- vector("numeric", nrow(lambda_seq))
  cum_p12 <- cumsum(c(0,colSums(rbind(p_1[-K],p_2[-K]))))
  time_after <- proc.time()-time_after
  
  time_k_4 <- list()
  time_GCV <- rep(0,5)
  print("Beginning tuning.", quote=FALSE)
  for(l in 1:nrow(lambda_seq)){
    time_l1 <- proc.time()
    lam_1 <- as.numeric(lambda_seq[l,1])
    lam_2 <- as.numeric(lambda_seq[l,2])
    
    penalty <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, 2+sum(p_1)+sum(p_2)+p_3)
    diag(penalty)[3:(2+sum(p_1))] <- N*lam_1
    diag(penalty)[(2+sum(p_1)):(2+sum(p_1)+sum(p_2))] <- N*lam_2
    
    weight <- solve(quad_form + penalty)
    
    output$coefficients[[l]] <- as.vector(weight %*%
                                            sensitivity %*% V_inv %*%
                                            unlist(lapply(1:K,
                                                          function(x) sensitivity_list[[x]] %*% estimates[[x]])))
    if(sum(is.na(output$coefficients[[l]]))!=0) break
    output$vcov[[l]] <- weight %*% sensitivity%*%V_inv%*% V %*%V_inv%*%t(sensitivity) %*% weight
    
    output$mu.fitted.values[[l]] <- as.vector(X_tilde_1%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))])
    output$mu.se[[l]] <- sqrt(apply(X_tilde_1, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))] %*%x
    }))
    
    output$sigma.fitted.values[[l]] <- as.vector(X_tilde_2%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))])
    output$sigma.se[[l]] <- sqrt(apply(X_tilde_2, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))] %*%x
    }))
    
    if(sum(is.na(output$sigma.se[[l]])) + sum(is.na(output$mu.se[[l]]))!=0) next
    
    output$coefficients_reparm[[l]] <- c(2*exp(output$coefficients[[l]][1])/(1+exp(output$coefficients[[l]][1])), exp(output$coefficients[[l]][2]),
                                         output$coefficients[[l]][-c(1:2)])
    output$var_reparm[[l]] <- diag(output$vcov[[l]])
    gradient <- diag(c(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])), output$coefficients_reparm[[l]][2], rep(1,length(output$coefficients[[l]])-2) ))
    output$var_reparm[[l]][1] <- output$var_reparm[[l]][1]*(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])))^2
    output$var_reparm[[l]][2] <- output$var_reparm[[l]][2]*(output$coefficients_reparm[[l]][2])^2
    output$vcov_reparm[[l]] <- gradient %*% output$vcov[[l]] %*% gradient
    
    names(output$coefficients[[l]]) <- rownames(output$vcov[[l]]) <- colnames(output$vcov[[l]]) <- 
      c("omega","xi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    names(output$coefficients_reparm[[l]]) <- names(output$var_reparm[[l]]) <- rownames(output$vcov_reparm[[l]]) <- colnames(output$vcov_reparm[[l]]) <-
      c("alpha","phi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    
    psi_list_GCV <- list()
    sensitivity_list_GCV <- list()
    time_l1 <- proc.time() - time_l1
    
    time_k_3 <- list()
    for(k in 1:K){
      print(k)
      time_k_3[[k]] <- proc.time()
      
      result_k <- iMSP_Chessian_all_thresh_reparm(c(output$coefficients[[l]][1:2], output$coefficients[[l]][2+(cum_p12[k]+1):(cum_p12[k]+p_1[k])],
                                               output$coefficients[[l]][2+(cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])], 
                                               output$coefficients[[l]][length(output$coefficients[[l]])]),
                                             thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                                             z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
      
      psi_list_GCV[[k]] <- result_k$score
      sensitivity_list_GCV[[k]] <- var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity)
      
      time_k_3[[k]] <- proc.time() - time_k_3[[k]]
    }
    time_k_4[[l]] <- time_k_3[[which.max(sapply(time_k_3, function(x) x[3]))]]
    
    time_l2 <- proc.time()
    
    sensitivity_GCV <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
    sensitivity_GCV[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      do.call(cbind, lapply(1:K, function(k) sensitivity_list_GCV[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
    sensitivity_GCV[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list_GCV[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
    
    psi_GCV <- do.call(rbind,psi_list_GCV)
    V_GCV <- var(t(psi_GCV), na.rm=T)*N
    V_inv_GCV <- solve(V_GCV)
    matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv_GCV[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
    V_inv_GCV <- matrix
    quad_form_GCV <- sensitivity_GCV %*% V_inv_GCV %*% t(sensitivity_GCV)
    
    weight_GCV <- solve(quad_form_GCV + penalty)
    
    GCV_numerator <- N^(-1)*as.numeric(rowMeans(psi_GCV)%*%V_inv_GCV%*%rowMeans(psi_GCV))
    output$GCV[l] <- GCV_numerator / ( 1 - N^(-1)*sum(diag( weight_GCV %*% quad_form_GCV )))^2
    time_l2 <- proc.time() - time_l2
    time_GCV <- time_GCV + time_l1 + time_l2
    time_GCV <- time_GCV + time_l1
  }
  
  total_time_k_4 <- do.call(rbind, time_k_4)
  if(!is.null(total_time_k_4)) {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after +
      colSums(do.call(rbind, time_k_4)) + time_GCV
  } else {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after 
  }
  output$time_noGCV <- time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after
  output$model <- "inverted MSP"
  return(output)
}

invBRdacVCMmusigma_parallel <- function(y, covariates_3, quantile=0.95, locs, indicator, lambda_1, lambda_2, cluster){
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
  
  if(quantile>=0) {
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y)-1, S)
  }
  
  unit_y <- -1 / log(apply(y, 2, rank) / (N  +1))
  uniform_y <- exp(-1/unit_y)
  exp_y <- -log(1-uniform_y)
  MSP_y <- 1/exp_y
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(MSP_y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  sock <- parallel::makeCluster(rep("localhost", cluster), type = "SOCK")
  doParallel::registerDoParallel(sock)
  
  marg_pars_list <- foreach::foreach(i=1:S) %dopar% {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    return(c(opt_1$par[1],exp(opt_1$par[2]),opt_1$par[3]))
  }
  mu <- sapply(marg_pars_list, function(x) x[1])
  sigma <- sapply(marg_pars_list, function(x) x[2])
  xi <- sapply(marg_pars_list, function(x) x[3])
  
  init_values <- as.vector(c(init_cov, solve(t(covariates_3)%*%covariates_3)%*%t(covariates_3)%*%rep(xi,N)))
  
  time_before <- proc.time() - time
  
  result_k <- foreach::foreach(k=1:K) %dopar% {
    time_k <- proc.time()
    
    pos_k <- locs[indicator==k,]
    s_k <- sum(indicator==k)
    
    knots_mat <- fields::cover.design( pos_k, 10, num.nn=s_k/2, start=floor(seq(1,s_k,length.out=10)) )$design
    attr(knots_mat,"scaled:scale") <- attr(knots_mat,"scaled:center") <- NULL
    colnames(knots_mat) <- NULL
    distances <- apply(knots_mat, 1, function(x) sqrt((pos_k[,1]-x[1])^2 + (pos_k[,2]-x[2])^2))
    basis_1 <- t(apply(distances, 1, function(x) sapply(x, function(t) gaussian_basis(0.05, 0, 1, t))))
    new_X <- cbind(rep(1,N*s_k), do.call(rbind, lapply(1:N, function(x) pos_k)), do.call(rbind, lapply(1:N, function(x) basis_1)))
    
    p_1 <- p_2 <- ncol(new_X)
    new_z_1[[k]] <- new_z_2[[k]] <- z_constructor(new_X, s_k, N, p_1)
    p <- p_1 + p_2 + p_3
    
    init_par <- c(init_values[1:2],
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(mu[indicator==k],N),
                  solve(t(new_X)%*%new_X) %*% t(new_X) %*%rep(log(sigma[indicator==k]),N),
                  init_values[-c(1:2)])
    
    mu_sub <- new_X %*% init_par[3:(2+p_1)]
    sigma_sub <- exp(new_X %*% init_par[(3+p_1):(2+p_1+p_2)])
    xi_sub <- covariates_3[rep(indicator==k,N),,drop=FALSE] %*% init_par[(3+p_1+p_2):(2+p)]
    
    y_long <- c(t(y[,indicator==k,drop=FALSE]))
    thresh_long <- rep(thresholds[indicator==k],N)
    scaled <- unlist(sapply(1:(N*s_k), function(i) {
      if(!is.na(y_long[i])) {
        if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
        else return((thresh_long[i] - mu_sub[i])*xi_sub[i]/sigma_sub[i])
      }
    }))
    fudge <- which(scaled < -1)
    
    if(length(fudge)>0){
      c <- max(1, -scaled[fudge]) + 1e-5
      init_par[(3+p_1+p_2):(2+p)] <- init_par[(3+p_1+p_2):(2+p)]/c
      
      opt <- optim(par=init_par[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[indicator==k], 
                   y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000)) 
      init_par <- c(init_par[1:2], opt$par)  
    }
    
    opt_2 <-
      optim(par=init_par,
            fn=iMSP_logCL_all_thresh_reparm, gr=iMSP_score_all_thresh_reparm,
            thresholds=thresholds[indicator==k],
            y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
            z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000))  
    
    return(list(opt_2$par, proc.time() - time_k, new_z_1[[k]], new_z_2[[k]], p_1, p_2, p))
  }
  estimates <- lapply(result_k, function(x) x[[1]])
  new_z_1 <- lapply(result_k, function(x) x[[3]])
  new_z_2 <- lapply(result_k, function(x) x[[4]])
  p_1 <- sapply(result_k, function(x) x[[5]])
  p_2 <- sapply(result_k, function(x) x[[6]])
  p <- sapply(result_k, function(x) x[[7]])
  time_k <- lapply(result_k, function(x) x[[2]])
  
  time_mean <- proc.time()
  mean_estimates <- colMeans(do.call(rbind, lapply(1:K, function(k) estimates[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),drop=FALSE])))
  time_mean <- proc.time() - time_mean
  
  weights_list <- foreach::foreach(k=1:K) %dopar% {
    time_k_2 <- proc.time()
    
    result_k <- iMSP_Chessian_all_thresh_reparm(c(mean_estimates[1:2], estimates[[k]][-c(1,2,length(estimates[[k]]))], mean_estimates[-c(1:2)]), 
                                           thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], 
                                           z_2=new_z_2[[k]], z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
    return(list(result_k$score, var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity), proc.time() - time_k_2))
  }
  psi_list <- lapply(weights_list, function(x) x[[1]])
  sensitivity_list <- lapply(weights_list, function(x) x[[2]])
  time_k_2 <- lapply(weights_list, function(x) x[[3]])
  
  time_after <- proc.time()
  psi <- do.call(rbind,psi_list)
  if(sum(is.na(psi))!=0) {
    warning(paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi))))
    output$message <- paste0("Missing values in the score kernels; model may not fit well. Proportion missingness: ", mean(is.na(psi)))
  } else {
    output$message <- NULL
  }
  
  sensitivity <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
  sensitivity[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    do.call(cbind, lapply(1:K, function(k) sensitivity_list[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
  sensitivity[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
    as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
  
  V <- var(t(psi), na.rm=T)*N
  V_inv <- solve(V)
  cum_p <- cumsum(c(0,p[-K]))
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
  V_inv <- matrix
  
  quad_form <- sensitivity %*% V_inv %*% t(sensitivity)
  
  pos_cols <- do.call(rbind, lapply(1:K, function(k) locs[indicator==k,]))
  order <- apply(locs, 1, function(x) which(equiv(pos_cols[,1],x[1]) & equiv(pos_cols[,2], x[2])))
  X_tilde_1 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,n,])))[order,] ))
  X_tilde_2 <- do.call(rbind, lapply(1:N, function(n) as.matrix(Matrix::bdiag(lapply(new_z_2, function(x) x[,n,])))[order,] ))
  
  lambda_seq <- as.matrix(expand.grid(lambda_1, lambda_2))
  output$coefficients <- output$coefficients_reparm <- output$vcov <- output$vcov_reparm <- output$var_reparm <- 
    output$mu.fitted.values <- output$mu.se <- output$sigma.fitted.values <- output$sigma.se <- list()
  output$GCV <- vector("numeric", nrow(lambda_seq))
  cum_p12 <- cumsum(c(0,colSums(rbind(p_1[-K],p_2[-K]))))
  time_after <- proc.time()-time_after
  
  time_k_4 <- list()
  time_GCV <- rep(0,5)
  print("Beginning tuning.", quote=FALSE)
  for(l in 1:nrow(lambda_seq)){
    time_l1 <- proc.time()
    lam_1 <- as.numeric(lambda_seq[l,1])
    lam_2 <- as.numeric(lambda_seq[l,2])
    
    penalty <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, 2+sum(p_1)+sum(p_2)+p_3)
    diag(penalty)[3:(2+sum(p_1))] <- N*lam_1
    diag(penalty)[(2+sum(p_1)):(2+sum(p_1)+sum(p_2))] <- N*lam_2
    
    weight <- solve(quad_form + penalty)
    
    output$coefficients[[l]] <- as.vector(weight %*%
                                            sensitivity %*% V_inv %*%
                                            unlist(lapply(1:K,
                                                          function(x) sensitivity_list[[x]] %*% estimates[[x]])))
    if(sum(is.na(output$coefficients[[l]]))!=0) break
    output$vcov[[l]] <- weight %*% sensitivity%*%V_inv%*% V %*%V_inv%*%t(sensitivity) %*% weight
    
    output$mu.fitted.values[[l]] <- as.vector(X_tilde_1%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))])
    output$mu.se[[l]] <- sqrt(apply(X_tilde_1, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+1):(cum_p12[k]+p_1[k]))))] %*%x
    }))
    
    output$sigma.fitted.values[[l]] <- as.vector(X_tilde_2%*%output$coefficients[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))])
    output$sigma.se[[l]] <- sqrt(apply(X_tilde_2, 1, function(x) {
      x%*%output$vcov[[l]][c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])))),
                           c(2+unlist(sapply(1:K, function(k) (cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k]))))] %*%x
    }))
    
    if(sum(is.na(output$sigma.se[[l]])) + sum(is.na(output$mu.se[[l]]))!=0) next
    
    output$coefficients_reparm[[l]] <- c(2*exp(output$coefficients[[l]][1])/(1+exp(output$coefficients[[l]][1])), exp(output$coefficients[[l]][2]),
                                         output$coefficients[[l]][-c(1:2)])
    output$var_reparm[[l]] <- diag(output$vcov[[l]])
    gradient <- diag(c(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])), output$coefficients_reparm[[l]][2], rep(1,length(output$coefficients[[l]])-2) ))
    output$var_reparm[[l]][1] <- output$var_reparm[[l]][1]*(output$coefficients_reparm[[l]][1]/(1+exp(output$coefficients[[l]][1])))^2
    output$var_reparm[[l]][2] <- output$var_reparm[[l]][2]*(output$coefficients_reparm[[l]][2])^2
    output$vcov_reparm[[l]] <- gradient %*% output$vcov[[l]] %*% gradient
    
    names(output$coefficients[[l]]) <- rownames(output$vcov[[l]]) <- colnames(output$vcov[[l]]) <- 
      c("omega","xi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    names(output$coefficients_reparm[[l]]) <- names(output$var_reparm[[l]]) <- rownames(output$vcov_reparm[[l]]) <- colnames(output$vcov_reparm[[l]]) <-
      c("alpha","phi",paste0("beta_1",1:sum(p_1)),paste0("beta_2",1:sum(p_2)),paste0("beta_3",1:p_3))
    
    psi_list_GCV <- list()
    sensitivity_list_GCV <- list()
    time_l1 <- proc.time() - time_l1
    
    weights_list <- foreach::foreach(k=1:K) %dopar% {
      time_k_3 <- proc.time()
      
      result_k <- iMSP_Chessian_all_thresh_reparm(c(output$coefficients[[l]][1:2], output$coefficients[[l]][2+(cum_p12[k]+1):(cum_p12[k]+p_1[k])],
                                               output$coefficients[[l]][2+(cum_p12[k]+p_1[k]+1):(cum_p12[k]+p_1[k]+p_2[k])], 
                                               output$coefficients[[l]][length(output$coefficients[[l]])]),
                                             thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                                             z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
      
      return(list(result_k$score, var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity), proc.time() - time_k_3))
    }
    psi_list_GCV <- lapply(weights_list, function(x) x[[1]])
    sensitivity_list_GCV <- lapply(weights_list, function(x) x[[2]])
    time_k_3 <- lapply(weights_list, function(x) x[[3]])
    
    time_k_4[[l]] <- time_k_3[[which.max(sapply(time_k_3, function(x) x[3]))]]
    
    time_l2 <- proc.time()
    
    sensitivity_GCV <- matrix(0, 2+sum(p_1)+sum(p_2)+p_3, K*(2+p_3)+sum(p_1)+sum(p_2))
    sensitivity_GCV[c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      do.call(cbind, lapply(1:K, function(k) sensitivity_list_GCV[[k]][c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),]))
    sensitivity_GCV[-c(1:2,(3+sum(p_1)+sum(p_2)):(2+sum(p_1)+sum(p_2)+p_3)),] <-
      as.matrix(Matrix::bdiag(lapply(1:K, function(k) sensitivity_list_GCV[[k]][-c(1:2,(3+p_1[k]+p_2[k]):(2+p[k])),])))
    
    psi_GCV <- do.call(rbind,psi_list_GCV)
    V_GCV <- var(t(psi_GCV), na.rm=T)*N
    V_inv_GCV <- solve(V_GCV)
    matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv_GCV[(cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2), (cum_p[k]+(k-1)*2+1):(cum_p[k]+p[k]+k*2)])))
    V_inv_GCV <- matrix
    quad_form_GCV <- sensitivity_GCV %*% V_inv_GCV %*% t(sensitivity_GCV)
    
    weight_GCV <- solve(quad_form_GCV + penalty)
    
    GCV_numerator <- N^(-1)*as.numeric(rowMeans(psi_GCV)%*%V_inv_GCV%*%rowMeans(psi_GCV))
    output$GCV[l] <- GCV_numerator / ( 1 - N^(-1)*sum(diag( weight_GCV %*% quad_form_GCV )))^2
    time_l2 <- proc.time() - time_l2
    time_GCV <- time_GCV + time_l1 + time_l2
    time_GCV <- time_GCV + time_l1
  }
  parallel::stopCluster(sock)
  foreach::registerDoSEQ()
  total_time_k_4 <- do.call(rbind, time_k_4)
  if(!is.null(total_time_k_4)) {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after +
      colSums(do.call(rbind, time_k_4)) + time_GCV
  } else {
    output$time <- 
      time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
      time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after 
  }
  output$time_noGCV <- time_before + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + time_mean + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]] + time_after
  output$model <- "inverted MSP"
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
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile))
  } else {
    thresholds <- rep(min(y)-1, S)
  }
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    mu[i] <- opt_1$par[1]
    sigma[i] <- exp(opt_1$par[2])
    xi[i] <- opt_1$par[3]
  }
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(init_cov,
                             solve(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             solve(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             solve(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  
  y_long <- c(t(y))
  thresh_long <- rep(thresholds,N)
  scaled <- unlist(sapply(1:(N*S), function(i) {
    if(!is.na(y_long[i])) {
      if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu[i])*xi[i]/sigma[i])
      else return((thresh_long[i] - mu[i])*xi[i]/sigma[i])
    }
  }))
  fudge <- which(scaled < -1)
  
  if(length(fudge)>0){
    
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
    
    opt_fix <- optim(par=init_values[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds, 
                     y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS", control=list(maxit=5000)) 
    init_values <- c(init_values[1:2], opt_fix$par)  
  }
  
  opt_2 <-
    optim(par=init_values, fn=logCL_all_thresh_reparm, gr=score_all_thresh_reparm, thresholds=thresholds,
          y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS")
  
  output$coefficients <- opt_2$par
  
  result <- Chessian_all_thresh_reparm(output$coefficients, thresholds=thresholds, y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs)
  psi <- result$score
  sensitivity <- var(t(result$sensitivity), na.rm=T)*ncol(result$sensitivity)
  output$vcov <- solve(sensitivity)%*%var(t(psi), na.rm=T)%*%solve(sensitivity)*N
  
  output$coefficients_reparm <- c(2*exp(output$coefficients[1])/(1+exp(output$coefficients[1])), exp(output$coefficients[2]),
                                  output$coefficients[-c(1:2)])
  output$var_reparm <- diag(output$vcov)
  gradient <- diag(c(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])), output$coefficients_reparm[2], rep(1,p) ))
  output$var_reparm[1] <- output$var_reparm[1]*(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])))^2
  output$var_reparm[2] <- output$var_reparm[2]*(output$coefficients_reparm[2])^2
  output$vcov_reparm <- gradient %*% output$vcov %*% gradient
  output$time <- proc.time()-time
  
  return(output)
}

tapered_CL <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, taper){
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
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile))
  } else {
    thresholds <- rep(min(y)-1, S)
  }
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    mu[i] <- opt_1$par[1]
    sigma[i] <- exp(opt_1$par[2])
    xi[i] <- opt_1$par[3]
  }
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(init_cov,
                             solve(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             solve(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             solve(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  
  y_long <- c(t(y))
  thresh_long <- rep(thresholds,N)
  scaled <- unlist(sapply(1:(N*S), function(i) {
    if(!is.na(y_long[i])) {
      if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu[i])*xi[i]/sigma[i])
      else return((thresh_long[i] - mu[i])*xi[i]/sigma[i])
    }
  }))
  fudge <- which(scaled < -1)
  
  if(length(fudge)>0){
    
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
    
    opt_fix <- optim(par=init_values[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds, 
                     y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS", control=list(maxit=5000)) 
    init_values <- c(init_values[1:2], opt_fix$par)  
  }
  
  opt_2 <-
    optim(par=init_values, fn=tapered_logCL_all_thresh_reparm, gr=tapered_score_all_thresh_reparm, thresholds=thresholds,
          y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, taper=taper, method="BFGS")
  
  output$coefficients <- opt_2$par
  
  result <- tapered_Chessian_all_thresh_reparm(output$coefficients, thresholds=thresholds, y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, taper=taper)
  psi <- result$score
  full_sensitivity <- result$sensitivity
  tapered_locations <- which( colSums(full_sensitivity==0) != nrow(full_sensitivity) )
  sensitivity <- var(t(full_sensitivity[,tapered_locations]), na.rm=T)*ncol(full_sensitivity[,tapered_locations])
  output$vcov <- solve(sensitivity)%*%var(t(psi), na.rm=T)%*%solve(sensitivity)*N
  
  output$coefficients_reparm <- c(2*exp(output$coefficients[1])/(1+exp(output$coefficients[1])), exp(output$coefficients[2]),
                                  output$coefficients[-c(1:2)])
  output$var_reparm <- diag(output$vcov)
  gradient <- diag(c(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])), output$coefficients_reparm[2], rep(1,p) ))
  output$var_reparm[1] <- output$var_reparm[1]*(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])))^2
  output$var_reparm[2] <- output$var_reparm[2]*(output$coefficients_reparm[2])^2
  output$vcov_reparm <- gradient %*% output$vcov %*% gradient
  output$time <- proc.time()-time
  output$taper <- taper
  
  return(output)
}

block_averaged_CL <- function(y, covariates_1, covariates_2, covariates_3, quantile=0.95, locs, indicator){
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
  variances <- list()
  psi_list <- sensitivity_list <- list()
  
  if(quantile>=0) {
    thresholds <- apply(y, 2, function(x) quantile(x, prob=quantile, na.rm=T))
  } else {
    thresholds <- rep(min(y, na.rm=T)-1, S)
  }
  
  n.pairs <- S * (S-1)/2
  extcoeff <- SpatialExtremes::fitextcoeff(y, locs, estim="ST", plot=FALSE, loess=FALSE, marge="emp", std.err=FALSE, prob=quantile)
  weights <- rep(1, n.pairs)
  extcoeff <- extcoeff[, "ext.coeff"]
  dist <- SpatialExtremes::distance(locs)
  funBR <- function(range, smooth) .C(SpatialExtremes:::C_fitbrcovariance, as.double(range), as.double(smooth), as.integer(n.pairs), 
                                      as.double(dist), as.double(extcoeff), as.double(weights), ans = double(1))$ans
  start <- list(range = max(dist)/2.73, smooth = 0.75)
  names(formals(funBR)) <- c("range", "smooth")
  obj.fun <- function(p, ...) funBR(p, ...)
  body(obj.fun) <- parse(text = paste("funBR(", paste("p[", 1:2, "]", collapse = ", "), ", ...)"))
  cov.start <- optim(unlist(start), obj.fun, hessian = FALSE)$par
  init_cov <- c(0,0)
  init_cov[1] <- log(cov.start[2]/(2-cov.start[2]))
  init_cov[2] <- log(cov.start[1])
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    init_par <- c(0, sqrt(6*var(y[,i], na.rm=TRUE))/pi, 1e-5)
    init_par[1] <- mean(y[,i], na.rm = TRUE) - 0.58*init_par[2]
    init_par[2] <- log(init_par[2])
    opt_1 <- optim(par=init_par, fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds[i],
                   y=y[,i,drop=FALSE], z_1=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   z_2=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1), z_3=z_constructor(matrix(1,nrow=N,ncol=1),1,N,1),
                   locs=locs[i,,drop=FALSE], method="BFGS", control=list(maxit=5000)) 
    mu[i] <- opt_1$par[1]
    sigma[i] <- exp(opt_1$par[2])
    xi[i] <- opt_1$par[3]
  }
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(init_cov,
                             solve(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             solve(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             solve(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  
  y_long <- c(t(y))
  thresh_long <- rep(thresholds,N)
  scaled <- unlist(sapply(1:(N*S), function(i) {
    if(!is.na(y_long[i])) {
      if(y_long[i] >= thresh_long[i]) return((y_long[i] - mu[i])*xi[i]/sigma[i])
      else return((thresh_long[i] - mu[i])*xi[i]/sigma[i])
    }
  }))
  fudge <- which(scaled < -1)
  
  if(length(fudge)>0){
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
    
    opt_fix <- optim(par=init_values[-c(1:2)], fn=logL_marg_thresh, gr=score_marg_thresh, thresholds=thresholds, 
                     y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS", control=list(maxit=5000)) 
    
    init_values <- c(init_values[1:2], opt_fix$par)  
  }
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    opt_2 <- optim(par=init_values, fn=logCL_all_thresh_reparm, gr=score_all_thresh_reparm, thresholds=thresholds[indicator==k],
                   y=y[,indicator==k,drop=FALSE], z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE], 
                   z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,],
                   method="BFGS")
    
    estimates[k,] <- opt_2$par
    
    result_k <- Chessian_all_thresh_reparm(estimates[k,], thresholds=thresholds[indicator==k], y=y[,indicator==k],
                                           z_1=z_1[indicator==k,,,drop=FALSE], z_2=z_2[indicator==k,,,drop=FALSE],
                                           z_3=z_3[indicator==k,,,drop=FALSE],
                                           locs=locs[indicator==k,])
    
    psi_list[[k]] <- result_k$score
    sensitivity_list[[k]] <- var(t(result_k$sensitivity), na.rm=T)*ncol(result_k$sensitivity)
    variances[[k]] <- N*solve(sensitivity_list[[k]] %*% solve(var(t(psi_list[[k]]), na.rm=T))%*%sensitivity_list[[k]])
    
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_after <- proc.time()
  
  output$coefficients <- colMeans(estimates)
  output$vcov <- Reduce("+",variances)/K
  
  output$coefficients_reparm <- c(2*exp(output$coefficients[1])/(1+exp(output$coefficients[1])), exp(output$coefficients[2]),
                                  output$coefficients[-c(1:2)])
  output$var_reparm <- diag(output$vcov)
  gradient <- diag(c(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])), output$coefficients_reparm[2], rep(1,p) ))
  output$var_reparm[1] <- output$var_reparm[1]*(output$coefficients_reparm[1]/(1+exp(output$coefficients[1])))^2
  output$var_reparm[2] <- output$var_reparm[2]*(output$coefficients_reparm[2])^2
  output$vcov_reparm <- gradient %*% output$vcov %*% gradient
  
  names(output$coefficients) <- rownames(output$vcov) <- colnames(output$vcov) <- 
    c("omega","xi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  names(output$coefficients_reparm) <- names(output$var_reparm) <- rownames(output$vcov_reparm) <- colnames(output$vcov_reparm) <-
    c("alpha","phi",paste0("beta_1",1:p_1),paste0("beta_2",1:p_2),paste0("beta_3",1:p_3))
  
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  return(output)
}

gaussian_basis <- function(phi2, tau2, sigma2, dist){
  if(dist>0) return(sigma2*exp(-phi2*dist^2))
  else return(tau2+sigma2)
}

equiv <- function(x, y, .tol = 1e-5) abs(x - y) <= .tol
