#work_dir <- "/Users/ehector/Dropbox/Projects/SpatialExtremes-Brian/simulations/"
work_dir <- "/home/ehector/share/research1/BRdac"
sim <- as.numeric(Sys.getenv('LSB_JOBINDEX'))

library(SpatialExtremes)
library(BRdac)
set.seed(sim)

alpha <- 1
phi <- 3
N <- 1000
s <- as.matrix(expand.grid(1:10,1:10))
S <- nrow(s)
beta_1 <- c(0.5,0.5)
beta_2 <- 1.5
beta_3 <- 0.5
p_1 <- length(beta_1)
p_2 <- length(beta_2)
p_3 <- length(beta_3)
p <- p_1 + p_2 + p_3

X <- rmaxstab(n=N, coord=s, cov.mod="brown", grid=FALSE, range=phi, smooth=alpha)

covariates_1 <- do.call(rbind, lapply(1:N, function(x) return(s)))
covariates_2 <- covariates_3 <- as.matrix(rep(1,N*S))
z_1 <- z_constructor(covariates_1, S, N, p_1)
z_2 <- z_constructor(covariates_2, S, N, p_2)
z_3 <- z_constructor(covariates_3, S, N, p_3)
mu <- apply(z_1, 1, function(x) x%*%beta_1)
sigma <- apply(z_2, 1, function(x) exp(x%*%beta_2))
xi <- apply(z_3, 1, function(x) x%*%beta_3)
Y <- mu + sigma/xi*(X^xi - 1)

output <- list()
output[[sim]]  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, locs=s, s_k=5)

package_output <- list()
package_output[[sim]] <- fitmaxstab(Y, s, cov.mod="brown", loc ~ 0 + s, scale ~ 1, shape ~ 1) ## doesn't allow subject-level covariates, only spatial covariates

save.image(paste0(work_dir, "BRdac_N", N, "S", S, "p", p, "_sim",sim, ".RData"))

