logCL_BR_wrap <- function(par, x, locs){
  return(-logCL_BR(alpha=par[1], phi=par[2], as.matrix(x), locs))
}

Cscore_BR_sum_wrap <- function(par, x, locs){
  return(-Cscore_BR_sum(alpha=par[1], phi=par[2], x, locs))
}

logCL_all_wrap <- function(par, y, z_1, z_2, z_3, locs){
  alpha <- par[1]
  phi <- par[2]
  p_1 <- dim(z_1)[3]
  p_2 <- dim(z_2)[3]
  p_3 <- dim(z_3)[3]
  beta_1 <- par[3:(2+p_1)]
  beta_2 <- par[(3+p_1):(2+p_1+p_2)]
  beta_3 <- par[(3+p_1+p_2):(2+p_1+p_2+p_3)]
  return(-logCL_all(alpha, phi, beta_1, beta_2, beta_3, as.matrix(y), z_1, z_2, z_3, as.matrix(locs)))
}

Cscore_all_sum_wrap <- function(par, y, z_1, z_2, z_3, locs){
  alpha <- par[1]
  phi <- par[2]
  p_1 <- dim(z_1)[3]
  p_2 <- dim(z_2)[3]
  p_3 <- dim(z_3)[3]
  beta_1 <- par[3:(2+p_1)]
  beta_2 <- par[(3+p_1):(2+p_1+p_2)]
  beta_3 <- par[(3+p_1+p_2):(2+p_1+p_2+p_3)]
  return(-Cscore_all_sum(alpha, phi, beta_1, beta_2, beta_3, y, z_1, z_2, z_3, locs))
}

logCL_margs_wrap <- function(par, alpha, phi, y, z_1, z_2, z_3, locs){
  p_1 <- dim(z_1)[3]
  p_2 <- dim(z_2)[3]
  p_3 <- dim(z_3)[3]
  beta_1 <- par[1:p_1]
  beta_2 <- par[(p_1+1):(p_1+p_2)]
  beta_3 <- par[(p_1+p_2+1):(p_1+p_2+p_3)]
  return(-logCL_all(alpha, phi, beta_1, beta_2, beta_3, y, z_1, z_2, z_3, locs))
}

Cscore_margs_sum_wrap <- function(par, alpha, phi, y, z_1, z_2, z_3, locs){
  p_1 <- dim(z_1)[3]
  p_2 <- dim(z_2)[3]
  p_3 <- dim(z_3)[3]
  beta_1 <- par[1:p_1]
  beta_2 <- par[(p_1+1):(p_1+p_2)]
  beta_3 <- par[(p_1+p_2+1):(p_1+p_2+p_3)]
  return(-Cscore_all_sum(alpha, phi, beta_1, beta_2, beta_3, y, z_1, z_2, z_3, locs)[-c(1:2)])
}

logCL_dep_wrap <- function(par, beta_1, beta_2, beta_3, y, z_1, z_2, z_3, locs){
  alpha <- par[1]
  phi <- par[2]
  return(-logCL_all(alpha, phi, beta_1, beta_2, beta_3, y, z_1, z_2, z_3, locs))
}

Cscore_dep_sum_wrap <- function(par, beta_1, beta_2, beta_3, y, z_1, z_2, z_3, locs){
  alpha <- par[1]
  phi <- par[2]
  return(-Cscore_all_sum(alpha, phi, beta_1, beta_2, beta_3, y, z_1, z_2, z_3, locs)[1:2])
}