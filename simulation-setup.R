#work_dir <- "/Users/ehector/Dropbox/Projects/SpatialExtremes-Brian/simulations/"
work_dir <- "/share/research1/ehector/"
sim <- as.numeric(Sys.getenv('LSB_JOBINDEX'))

library(SpatialExtremes)
library(BRdac)
set.seed(sim)

alpha <- 0.8
phi <- 10
N <- 1000
dim <- 20

s <- as.matrix(expand.grid(1:dim,1:dim))
S <- nrow(s)

beta_1 <- c(0.5,0.5)
beta_2 <- 1.5
beta_3 <- 0.2
p_1 <- length(beta_1)
p_2 <- length(beta_2)
p_3 <- length(beta_3)
p <- p_1 + p_2 + p_3

X <- rmaxstab(n=N, coord=s, cov.mod="brown", grid=FALSE, range=phi, smooth=alpha)

covariates_1 <- as.matrix(do.call(rbind, lapply(1:N, function(x) return(s))))
covariates_2 <- covariates_3 <- as.matrix(rep(1,N*S))
z_1 <- z_constructor(covariates_1, S, N, p_1)
z_2 <- z_constructor(covariates_2, S, N, p_2)
z_3 <- z_constructor(covariates_3, S, N, p_3)
mu <- apply(z_1, 1, function(x) x%*%beta_1)
sigma <- apply(z_2, 1, function(x) exp(x%*%beta_2))
xi <- apply(z_3, 1, function(x) x%*%beta_3)
Y <- mu + sigma/xi*(X^xi - 1)

print("Computing BRdac, s_K=25.")
height_block_dim <- 5
width_block_dim <- 5
indicator <- rep(0,S)
height_num <- dim/height_block_dim
width_num <- dim/width_block_dim
for(k in 1:height_num){
  for(j in 1:width_num){
    indicator[s[,1] > (k-1)*height_block_dim & s[,1]<=k*height_block_dim & 
                s[,2] > (j-1)*width_block_dim & s[,2]<=j*width_block_dim] <- (k-1)*width_num+j
  }
}
output_25 <- list()
output_25[[sim]]  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, quantile=0.9, locs=s, indicator)
save.image(paste0(work_dir, "thresholded_alpha", alpha, "phi",phi,"N", N, "S", S, "p", p, "sim",sim, ".RData"))

print("Computing BRdac, s_K=40.")
height_block_dim <- 10
width_block_dim <- 4
indicator <- rep(0,S)
height_num <- dim/height_block_dim
width_num <- dim/width_block_dim
for(k in 1:height_num){
  for(j in 1:width_num){
    indicator[s[,1] > (k-1)*height_block_dim & s[,1]<=k*height_block_dim & 
                s[,2] > (j-1)*width_block_dim & s[,2]<=j*width_block_dim] <- (k-1)*width_num+j
  }
}
output_40 <- list()
output_40[[sim]]  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, quantile=0.9, locs=s, indicator)
save.image(paste0(work_dir, "thresholded_alpha", alpha, "phi",phi,"N", N, "S", S, "p", p, "sim",sim, ".RData"))

print("Computing BRdac, s_K=50.")
height_block_dim <- 10
width_block_dim <- 5
indicator <- rep(0,S)
height_num <- dim/height_block_dim
width_num <- dim/width_block_dim
for(k in 1:height_num){
  for(j in 1:width_num){
    indicator[s[,1] > (k-1)*height_block_dim & s[,1]<=k*height_block_dim & 
                s[,2] > (j-1)*width_block_dim & s[,2]<=j*width_block_dim] <- (k-1)*width_num+j
  }
}
output_50 <- list()
output_50[[sim]]  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, quantile=0.9, locs=s, indicator)
save.image(paste0(work_dir, "thresholded_alpha", alpha, "phi",phi,"N", N, "S", S, "p", p, "sim",sim, ".RData"))

print("Computing CL.")
CL_output <- list()
CL_time <- system.time(CL_output[[sim]] <- CL(y=Y, covariates_1, covariates_2, covariates_3, quantile=0.9, locs=s))

rm(covariates_1, covariates_2, covariates_3, mu, X, xi, Y, z_1, z_2, z_3)
save.image(paste0(work_dir, "thresholded_alpha", alpha, "phi",phi,"N", N, "S", S, "p", p, "sim",sim, ".RData"))

