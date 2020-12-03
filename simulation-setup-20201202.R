library(RandomFields)
library(SpatialExtremes)
library(Rcpp)

sourceCpp("/Users/ehector/Dropbox/Projects/SpatialExtremes-Brian/simulations/CL-funcs-20201202.cpp")
source("/Users/ehector/Dropbox/Projects/SpatialExtremes-Brian/simulations/R-wrappers-20201202.R")

alpha <- 1
phi <- 3
N <- 100
s <- as.matrix(expand.grid(1:10,1:10))
S <- nrow(s)
p_1 <- p_2 <- p_3 <- 2
beta_1 <- c(0.5,0.5)
beta_2 <- c(1,1)
beta_3 <- c(2,2)

set.seed(1234)
X <- rmaxstab(n=N, coord=s, cov.mod="brown", grid=FALSE, range=phi, smooth=alpha)

#suppressWarnings( init_values <- fitmaxstab(X, s, cov.mod="brown", fit.marge=FALSE)$fitted.values )
opt <- nlminb(start=c(1,1), objective=logCL_BR_wrap, gradient=Cscore_BR_sum_wrap, x=X, locs=s, 
             lower=c(1e-5,1e-5), upper=c(2-1e-5,Inf))
psi <- do.call(cbind, Cscore_BR(opt$par[1], opt$par[2], x=X, locs=s))

#covariates <- do.call(rbind, lapply(1:N, function(x) return(s))) 
covariates <- do.call(rbind, lapply(1:N, function(x) return(s))) /10
z_1 <- z_2 <- z_3 <- z_constructor(covariates, S, N, p_1)
mu <- apply(z_1, 1, function(x) x%*%beta_1)
sigma <- apply(z_2, 1, function(x) exp(x%*%beta_2))
xi <- apply(z_3, 1, function(x) x%*%beta_3)
Y <- mu + sigma/xi*(X^xi - 1)

opt_2 <-
  nlminb(start=c(opt$par, rep(1,p_1+p_2+p_3)), objective=logCL_all_wrap, gradient=Cscore_all_sum_wrap, 
         y=Y, z_1=z_1, z_2=z_2, z_3=z_3, locs=s,
         lower=c(1e-5,1e-5,rep(-Inf,p_1+p_2+p_3)), upper=c(2-1e-5,rep(Inf,1+p_1+p_2+p_3)))

#fitspatgev(Y, s, loc ~ 1, scale ~ 1, shape ~ 1) ## doesn't allow subject-level covariates, only spatial covariates
