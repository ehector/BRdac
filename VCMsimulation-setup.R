#work_dir <- "/Users/ehector/Dropbox/Projects/SpatialExtremes-Brian/simulations/"
work_dir <- "/share/research1/ehector/"
sim <- as.numeric(Sys.getenv('LSB_JOBINDEX'))

library(SpatialExtremes)
library(BRdac)
set.seed(sim)

alpha <- 1
phi <- 3
N <- 1000
dim <- 10

s <- as.matrix(expand.grid(1:dim,1:dim))
S <- nrow(s)
square_block_dim <- 5
indicator <- rep(0,S)
num <- sqrt(S/square_block_dim^2)
for(k in 1:num){
  for(j in 1:num){
    indicator[s[,1] > (k-1)*square_block_dim & s[,1]<=k*square_block_dim & 
                s[,2] > (j-1)*square_block_dim & s[,2]<=j*square_block_dim] <- (k-1)*num+j
  }
}

beta_1 <- apply(s, 1, function(x) sqrt(x[1]^2+x[2]^2))/10
beta_2 <- 1.5
beta_3 <- 0.2
p_1 <- length(beta_1)
p_2 <- length(beta_2)
p_3 <- length(beta_3)
p <- p_1 + p_2 + p_3

X <- rmaxstab(n=N, coord=s, cov.mod="brown", grid=FALSE, range=phi, smooth=alpha)

covariates_1 <- do.call(rbind, lapply(1:N, function(x) diag(S)))
covariates_2 <- covariates_3 <- as.matrix(rep(1,N*S))
z_1 <- z_constructor(covariates_1, S, N, p_1)
z_2 <- z_constructor(covariates_2, S, N, p_2)
z_3 <- z_constructor(covariates_3, S, N, p_3)
mu <- t(matrix(rep(beta_1,N), S, N))
sigma <- apply(z_2, 1, function(x) exp(x%*%beta_2))
xi <- apply(z_3, 1, function(x) x%*%beta_3)
Y <- mu + sigma/xi*(X^xi - 1)

output <- list()
output[[sim]]  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, locs=s, indicator)

#package_output <- list()
#package_time <- system.time(package_output[[sim]] <- fitmaxstab(Y, s, cov.mod="brown", loc ~ 0 + s, scale ~ 1, shape ~ 1)) ## doesn't allow subject-level covariates, only spatial covariates

#CL_output <- list()
#CL_time <- system.time(CL_output[[sim]] <- CL(y=Y, covariates_1, covariates_2, covariates_3, locs=s))

save.image(paste0(work_dir, "BRdacVCM_N", N, "S", S, "p", p, "sk", paste(unique(table(indicator)), collapse=""), "sim",sim, ".RData"))

