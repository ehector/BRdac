logCL_BR_wrap <- function(par, x, loc_1, loc_2){
  return(-logCL_BR(alpha=par[1], phi=par[2], x, loc_1, loc_2))
}

logCL_all_wrap <- function(par, y, x, z_1_1, z_2_1, z_3_1 , z_1_2, z_2_2, z_3_2, loc_1, loc_2){
  alpha <- par[1]
  phi <- par[2]
  beta_1 <- par[3:(2+ncol(z_1_1))]
  beta_2 <- par[(3+ncol(z_1_1)):(2+ncol(z_1_1)+ncol(z_2_1))]
  beta_3 <- par[(3+ncol(z_1_1)+ncol(z_2_1)):(2+ncol(z_1_1)+ncol(z_2_1)+ncol(z_3_1))]
  return(-logCL_all(alpha, phi, beta_1, beta_2, beta_3, y, x, z_1_1, z_2_1, z_3_1 , z_1_2, z_2_2, z_3_2, loc_1, loc_2))
}

Cscore_BR_sum_wrap <- function(par, x, loc_1, loc_2){
  return(-Cscore_BR_sum(alpha=par[1], phi=par[2], x, loc_1, loc_2))
}

Cscore_all_sum_wrap <- function(par, y, x, z_1_1, z_2_1, z_3_1 , z_1_2, z_2_2, z_3_2, loc_1, loc_2){
  alpha <- par[1]
  phi <- par[2]
  beta_1 <- par[3:(2+ncol(z_1_1))]
  beta_2 <- par[(3+ncol(z_1_1)):(2+ncol(z_1_1)+ncol(z_2_1))]
  beta_3 <- par[(3+ncol(z_1_1)+ncol(z_2_1)):(2+ncol(z_1_1)+ncol(z_2_1)+ncol(z_3_1))]
  return(-Cscore_all_sum(alpha, phi, beta_1, beta_2, beta_3, y, x, z_1_1, z_2_1, z_3_1 , z_1_2, z_2_2, z_3_2, loc_1, loc_2))
}