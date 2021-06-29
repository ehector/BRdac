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
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    marg.param <- SpatialExtremes::gevmle(y[,i])
    mu[i] <- marg.param["loc"]
    sigma[i] <- marg.param["scale"]
    xi[i] <- marg.param["shape"]
  }
  covariates_1_cc <- covariates_1[complete.cases(covariates_1),]
  covariates_2_cc <- covariates_2[complete.cases(covariates_2),]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  init_values <- as.vector(c(cov.start[2], cov.start[1],
                             matrix_inv(t(covariates_1_cc)%*%covariates_1_cc)%*%t(covariates_1_cc)%*%(rep(mu,N)[complete.cases(covariates_1)]), 
                             matrix_inv(t(covariates_2_cc)%*%covariates_2_cc)%*%t(covariates_2_cc)%*%(rep(log(sigma),N)[complete.cases(covariates_2)]), 
                             matrix_inv(t(covariates_3_cc)%*%covariates_3_cc)%*%t(covariates_3_cc)%*%(rep(xi,N)[complete.cases(covariates_3)])))
  mu <- covariates_1 %*% init_values[3:(2+p_1)]
  sigma <- exp(covariates_2 %*% init_values[(3+p_1):(2+p_1+p_2)])
  xi <- covariates_3 %*% init_values[(3+p_1+p_2):(2+p)]
  scaled <- ( c(t(y)) - mu )*xi/sigma
  fudge <- which(scaled < -1 & xi < 0)
  if(length(fudge)>0) {
    c <- max(1, -scaled[fudge]) + 1e-5
    init_values[(3+p_1+p_2):(2+p)] <- init_values[(3+p_1+p_2):(2+p)]/c
  }
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    time_k[[k]] <- proc.time()
    
    opt_2 <-
      optim(par=init_values, fn=logCL_all_thresh, gr=score_all_thresh, thresholds=thresholds[indicator==k],
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
    
    result_k <- Chessian_all_thresh(mean_estimates, thresholds=thresholds[indicator==k], y=y[,indicator==k],
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
  V_inv <- matrix_inv(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  weight <- matrix_inv(sensitivity %*% V_inv %*% t(sensitivity))
  output$coefficients <- as.vector(weight %*% 
                                     sensitivity %*% V_inv %*% 
                                     unlist(lapply(1:K, 
                                                   function(k) {
                                                     sensitivity_list[[k]] %*% estimates[k,]
                                                   }))) 
  output$vcov <- weight %*% sensitivity%*%V_inv%*%V%*%V_inv%*%t(sensitivity) %*% weight
  time_after <- proc.time()-time_after
  
  output$time <- 
    time_before + time_mean + time_after + time_k[[which.max(sapply(time_k, function(x) x[3]))]] + 
    time_k_2[[which.max(sapply(time_k_2, function(x) x[3]))]]
  
  output$quantile <- quantile
  output$indicator <- indicator
  return(output)
}

BRdacVCMmusigma <- function(y, covariates_3, quantile=0.95, locs, indicator, lambda_1, lambda_2){
  time <- proc.time()
  output <- list()
  
  p_3 <- dim(covariates_3)[2]
  covariates_3_cc <- covariates_3[complete.cases(covariates_3),]
  
  S <- nrow(locs)
  N <- nrow(y)
  K <- length(unique(indicator))
  
  new_z_1 <- new_z_2 <- list()
  z_3 <- z_constructor(covariates_3, S, N, p_3)
  
  estimates <- list()
  psi_list <- sensitivity_list <- list()
  
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
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    marg.param <- SpatialExtremes::gevmle(y[,i])
    mu[i] <- marg.param["loc"]
    sigma[i] <- marg.param["scale"]
    xi[i] <- marg.param["shape"]
  }
  
  init_values <- as.vector(c(cov.start[2], cov.start[1], matrix_inv(t(covariates_3)%*%covariates_3)%*%t(covariates_3)%*%rep(xi,N)))
  
  time_k <- list()
  time_before <- proc.time() - time
  for(k in 1:K){
    print(paste0("k=",k), quote=FALSE)
    time_k[[k]] <- proc.time()
    
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
    p <- p_1+p_2+p_3
    
    init_par <- c(init_values[1:2],
                  matrix_inv(t(new_X)%*%new_X) %*% t(new_X) %*%rep(mu[indicator==k],N),
                  matrix_inv(t(new_X)%*%new_X) %*% t(new_X) %*%rep(sigma[indicator==k],N),
                  init_values[-c(1:2)])
    opt_2 <- tryCatch(
      optim(par=init_par,
            fn=logCL_all_thresh,
            thresholds=thresholds[indicator==k],
            y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
            z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000)),
      error = function(c) list(message="error")
    )
    if(!is.null(opt_2$message)){
      if(opt_2$message == "error"){
        print("Optimization failed. Trying with analytic gradient.")
        opt_2 <-
          optim(par=init_par,
                fn=logCL_all_thresh, gr=score_all_thresh,
                thresholds=thresholds[indicator==k],
                y=y[,indicator==k,drop=FALSE], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                z_3=z_3[indicator==k,,,drop=FALSE], locs=pos_k, method="BFGS", control=list(maxit=1000))  
      }
    }
    print(paste0("Convergence = ",opt_2$convergence), quote=FALSE)
    estimates[[k]] <- opt_2$par
    time_k[[k]] <- proc.time() - time_k[[k]]
  }
  
  time_mean <- proc.time()
  estimates <- do.call(rbind, estimates)
  mean_estimates <- colMeans(estimates[,c(1:2,(2+p_1+p_2+1):(2+p_1+p_2+p_3)),drop=FALSE])
  time_mean <- proc.time() - time_mean
  
  time_k_2 <- list()
  for(k in 1:K){
    time_k_2[[k]] <- proc.time()
    
    result_k <- Chessian_all_thresh(c(mean_estimates[1:2], estimates[k,c(3:(2+p_1+p_2))], mean_estimates[-c(1:2)]), 
                                    thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], 
                                    z_2=new_z_2[[k]], z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
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
  V_inv <- matrix_inv(V)
  matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
  V_inv <- matrix
  
  quad_form <- sensitivity %*% V_inv %*% t(sensitivity)
  
  pos_cols <- do.call(rbind, lapply(new_z_1, function(x) x[,1,2:3]))
  order <- apply(locs, 1, function(x) which(equiv(pos_cols[,1],x[1]) & equiv(pos_cols[,2], x[2])))
  X_tilde_1 <- as.matrix(Matrix::bdiag(lapply(new_z_1, function(x) x[,1,])))
  X_tilde_1 <- X_tilde_1[order,]
  X_tilde_2 <- as.matrix(Matrix::bdiag(lapply(new_z_2, function(x) x[,1,])))
  X_tilde_2 <- X_tilde_2[order,]
  
  lambda_seq <- as.matrix(expand.grid(lambda_1, lambda_2))
  output$coefficients <- output$vcov <- output$mu.fitted.values <- output$mu.se <- output$sigma.fitted.values <- output$sigma.se <- list()
  output$GCV <- vector("numeric", nrow(lambda_seq))
  time_after <- proc.time()-time_after
  
  time_k_4 <- list()
  time_GCV <- rep(0,5)
  print("Beginning tuning.", quote=FALSE)
  for(l in 1:nrow(lambda_seq)){
    print(paste0("l=",l), quote=FALSE)
    time_l1 <- proc.time()
    lam_1 <- as.numeric(lambda_seq[l,1])
    lam_2 <- as.numeric(lambda_seq[l,2])
    
    penalty <- matrix(0, 2+K*(p_1+p_2)+p_3, 2+K*(p_1+p_2)+p_3)
    diag(penalty)[3:(2+K*p_1)] <- N*lam_1
    diag(penalty)[(2+K*p_1+1):(2+K*p_1+K*p_2)] <- N*lam_2
    
    weight <- matrix_inv(quad_form + penalty)
    
    output$coefficients[[l]] <- as.vector(weight %*% 
                                            sensitivity %*% V_inv %*% 
                                            unlist(lapply(1:K, 
                                                          function(x) sensitivity_list[[x]] %*% estimates[x,]))) 
    if(sum(is.na(output$coefficients[[l]]))!=0) break
    output$vcov[[l]] <- weight %*% sensitivity%*%V_inv%*% V %*%V_inv%*%t(sensitivity) %*% weight
    
    output$mu.fitted.values[[l]] <- as.vector(X_tilde_1%*%output$coefficients[[l]][c(2+sapply(1:K, function(k) ((k-1)*(p_1+p_2)+1):((k-1)*(p_1+p_2)+p_1)))])
    output$mu.se[[l]] <- sqrt(diag(X_tilde_1 %*% output$vcov[[l]][c(2+sapply(1:K, function(k) ((k-1)*(p_1+p_2)+1):((k-1)*(p_1+p_2)+p_1))), 
                                                                  c(2+sapply(1:K, function(k) ((k-1)*(p_1+p_2)+1):((k-1)*(p_1+p_2)+p_1)))] %*%t(X_tilde_1)))
    
    output$sigma.fitted.values[[l]] <- as.vector(X_tilde_2%*%output$coefficients[[l]][c(2+sapply(1:K, function(k) ((k-1)*(p_1+p_2)+p_1+1):((k-1)*(p_1+p_2)+p_1+p_2)))])
    output$sigma.se[[l]] <- sqrt(diag(X_tilde_2 %*% output$vcov[[l]][c(2+sapply(1:K, function(k) ((k-1)*(p_1+p_2)+p_1+1):((k-1)*(p_1+p_2)+p_1+p_2))), 
                                                                     c(2+sapply(1:K, function(k) ((k-1)*(p_1+p_2)+p_1+1):((k-1)*(p_1+p_2)+p_1+p_2)))] %*%t(X_tilde_2)))
    
    if(sum(is.na(output$sigma.se[[1]])) + sum(is.na(output$mu.se[[1]]))!=0) next
    psi_list_GCV <- list()
    sensitivity_list_GCV <- list()
    time_l1 <- time_l1 - proc.time()
    
    time_k_3 <- list()
    for(k in 1:K){
      time_k_3[[k]] <- proc.time()
      
      result_k <- Chessian_all_thresh(c(output$coefficients[[l]][1:2], output$coefficients[[l]][(2+(k-1)*(p_1+p_2)+1):(2+k*(p_1+p_2))],
                                        output$coefficients[[l]][(2+K*(p_1+p_2)+1):(2+K*(p_1+p_2)+p_3)]),
                                      thresholds=thresholds[indicator==k], y=y[,indicator==k], z_1=new_z_1[[k]], z_2=new_z_2[[k]],
                                      z_3=z_3[indicator==k,,,drop=FALSE], locs=locs[indicator==k,])
      
      psi_list_GCV[[k]] <- result_k$score
      sensitivity_list_GCV[[k]] <- result_k$sensitivity %*% t(result_k$sensitivity)
      
      time_k_3[[k]] <- proc.time() - time_k_3[[k]]
    }
    time_k_4[[l]] <- time_k_3[[which.max(sapply(time_k_3, function(x) x[3]))]]
    
    time_l2 <- proc.time()
    sensitivity_GCV <- do.call(cbind, sensitivity_list_GCV)
    sensitivity_GCV <- matrix(0, 2+K*(p_1+p_2)+p_3, K*(2+p_1+p_2+p_3))
    sensitivity_GCV[c(1:2,(2+K*(p_1+p_2)+1):(2+K*(p_1+p_2)+p_3)),] <- do.call(cbind, lapply(sensitivity_list_GCV, function(x) x[c(1:2,(2+p_1+p_2+1):(2+p_1+p_2+p_3)),]))
    sensitivity_GCV[-c(1:2,(2+K*(p_1+p_2)+1):(2+K*(p_1+p_2)+p_3)),] <-
      as.matrix(Matrix::bdiag(lapply(sensitivity_list_GCV, function(x) x[-c(1:2,(2+p_1+p_2+1):(2+p_1+p_2+p_3)),])))
    
    psi_GCV <- do.call(rbind,psi_list_GCV)
    V_GCV <- (psi_GCV-rowMeans(psi_GCV))%*%t(psi_GCV-rowMeans(psi_GCV))
    V_inv_GCV <- matrix_inv(V_GCV)
    matrix <- as.matrix(Matrix::bdiag(lapply(1:K, function(k) V_inv_GCV[((k-1)*(p+2)+1):(k*(p+2)), ((k-1)*(p+2)+1):(k*(p+2))])))
    V_inv_GCV <- matrix
    quad_form_GCV <- sensitivity_GCV %*% V_inv_GCV %*% t(sensitivity_GCV)
    
    weight_GCV <- matrix_inv(quad_form_GCV + penalty)
    
    GCV_numerator <- N^(-1)*as.numeric(rowMeans(psi_GCV)%*%V_inv_GCV%*%rowMeans(psi_GCV))
    output$GCV[l] <- GCV_numerator / ( 1 - N^(-1)*sum(diag( weight_GCV %*% quad_form_GCV )))^2
    time_l2 <- proc.time() - time_l2
    time_GCV <- time_GCV + time_l1 + time_l2
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
  
  mu <- sigma <- xi <- rep(0, S)
  for (i in 1:S) {
    marg.param <- SpatialExtremes::gevmle(y[,i])
    mu[i] <- marg.param["loc"]
    sigma[i] <- marg.param["scale"]
    xi[i] <- marg.param["shape"]
  }
  
  init_values <- as.vector(c(cov.start[2], cov.start[1], 
                             matrix_inv(t(covariates_1)%*%covariates_1)%*%t(covariates_1)%*%rep(mu,N), 
                             matrix_inv(t(covariates_2)%*%covariates_2)%*%t(covariates_2)%*%rep(sigma,N), 
                             matrix_inv(t(covariates_3)%*%covariates_3)%*%t(covariates_3)%*%rep(xi,N)))
  
  opt_2 <-
    optim(par=init_values, fn=logCL_all_thresh, gr=score_all_thresh, thresholds=thresholds,
          y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs, method="BFGS")
  
  output$coefficients <- opt_2$par
  
  result <- Chessian_all_thresh(output$coefficients, thresholds=thresholds, y=y, z_1=z_1, z_2=z_2, z_3=z_3, locs=locs)
  psi <- result$score
  sensitivity <- result$sensitivity %*% t(result$sensitivity)
  output$vcov <- matrix_inv(sensitivity)%*%psi%*%t(psi)%*%matrix_inv(sensitivity)
  
  output$time <- proc.time()-time
  
  return(output)
}

gaussian_basis <- function(phi2, tau2, sigma2, dist){
  if(dist>0) return(sigma2*exp(-phi2*dist^2))
  else return(tau2+sigma2)
}
