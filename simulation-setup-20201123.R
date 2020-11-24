library(RandomFields)
library(SpatialExtremes)
library(Rcpp)

sourceCpp("/Users/ehector/Dropbox/Projects/SpatialExtremes-Brian/simulations/CL-funcs-20201123.cpp")
source("/Users/ehector/Dropbox/Projects/SpatialExtremes-Brian/simulations/R-wrappers-20201123.R")

alpha <- 0.5
phi <- 3
m <- 500
s <- as.matrix(expand.grid(1,1:2))

z_1_1 <- z_2_1 <- z_3_1 <- as.matrix(rnorm(m,0,1))
z_1_2 <- z_2_2 <- z_3_2 <- as.matrix(rnorm(m,0,1))

beta_1 <- 0.5
beta_2 <- 1
beta_3 <- 0.5
mu_1 <- beta_1*z_1_1
sigma_1 <- exp(beta_2*z_2_1)
xi_1 <- beta_3*z_3_1
mu_2 <- beta_1*z_1_2
sigma_2 <- exp(beta_2*z_2_2)
xi_2 <- beta_3*z_3_2

mod <- RPbrownresnick(xi=1, mu=1, s=1, phi=RMfbm(alpha=alpha, scale=phi))
invisible(capture.output(R <- RFsimulate(mod, n=m, x=s)))
X <- t(as.matrix(R@data))
#fitspatgev(X, s, loc ~ 1, scale ~ 1, shape ~ 1)
suppressWarnings( init_values <- fitmaxstab(X, s, cov.mod="brown", fit.marge=FALSE)$fitted.values )
opt <- optim(par=c(1,1), fn=logCL_BR_wrap, gr=Cscore_BR_sum_wrap, x=X, loc_1=s[1,], loc_2=s[2,], 
             method="L-BFGS-B", lower=c(1e-5,1e-5), upper=c(2,Inf), control=list(maxit=1000))
psi <- do.call(cbind, Cscore_BR(opt$par[1], opt$par[2], x=X, loc_1=s[1,], loc_2=s[2,]))

Y <- matrix(0, nrow(X), ncol(X))
Y[,1] <- mu_1 + sigma_1/xi_1*(X[,1]^xi_1 - 1)
Y[,2] <- mu_2 + sigma_2/xi_2*(X[,2]^xi_2 - 1)

opt <-
  optim(par=c(alpha,phi,beta_1,beta_2,beta_3), fn=logCL_all_wrap, gr=Cscore_all_sum_wrap, y=Y, x=X, z_1_1=z_1_1, z_2_1=z_2_1, z_3_1=z_3_1, 
        z_1_2=z_1_2, z_2_2=z_2_2, z_3_2=z_3_2, loc_1=s[1,], loc_2=s[2,],
        method="L-BFGS-B", lower=c(1e-5,1e-5,rep(-Inf,3)), upper=c(2,rep(Inf,4)), control=list(maxit=1000))

#fitspatgev(Y, data.frame(z_1_1=z_1_1, z_1_2=z_1_2, z_2_1=z_2_1, z_2_2=z_2_2, z_3_1=z_3_1, z_3_2=z_3_2), 
#           loc ~ z_1_1+z_1_2, scale ~ z_2_1 + z_2_2, shape ~ z_3_1 + z_3_2) ## doesn't allow subject-level covariates, only spatial covariates
