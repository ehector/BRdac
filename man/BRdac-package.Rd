\name{BRdac-package}
\alias{BRdac-package}
\alias{BRdac}
\docType{package}
\title{\packageTitle{BRdac}}
\description{\packageDescription{BRdac}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{BRdac}
  \packageIndices{BRdac}

  This package contains the functions to perform regression analysis of extreme spatial data using a divide-and-conquer approach in the Brown-Resnick max-stable process model, as well as the inverted Brown-Resnick max-stable process model. The R package's main files are:

- src/CL-funcs.cpp: this file defines the Rcpp functions that compute the pairwise censored composite likelihood and censored composite score for the Brown-Resnick model with unit Fréchet margins, and the pairwise censored composite likelihood and censored composite score functions for the transformed data with regression models for the location, scale and shape.

- R/BRdac_func.R: this file defines the R function for the divide-and-conquer estimation of model parameters.

The main functions are:

- BRdac: estimates regression parameters of the Brown-Resnick max-stable process model with non spatially-varying coefficients (Section 3)

- invBRdac: estimates regression parameters of the inverted Brown-Resnick max-stable process with non spatially-varying coefficient model (Section 3)

- BRdacVCMmusigma: estimates regression parameters of the Brown-Resnick max-stable process model with spatially-varying coefficient model (Section 4). One Gaussian basis function is used within each block with 10 knots chosen by minimizing a geometric space-filling criterion (Royle and Nychka, 1998). For more control over the spatially varying coefficient model, users should create their own covariate expansions and use the BRdac function.

- invBRdacVCMmusigma: estimates regression parameters of the inverted Brown-Resnick max-stable process model with spatially-varying coefficient model (Section 4). One Gaussian basis function is used within each block with 10 knots chosen by minimizing a geometric space-filling criterion (Royle and Nychka, 1998). For more control over the spatially varying coefficient model, users should create their own covariate expansions and use the invBRdac function.

- BRdacObjective: performs a goodness-of-fit test for the Brown-Resnick max-stable process model by evaluating the GMM objective function at the supplied model output and comparing it to the appropriate chi squared quantile.

- invBRdacObjective: performs a goodness-of-fit test for the inverted Brown-Resnick max-stable process model by evaluating the GMM objective function at the supplied model output and comparing it to the appropriate chi squared quantile.
}
\author{
  \packageAuthor{BRdac}
  
  Maintainer: \packageMaintainer{BRdac}
}
\references{
# Installation

The BRdac R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL BRdac_1.0-1.tar.gz
- from the downloaded and renamed BRdac folder as R CMD build BRdac and R CMD INSTALL BRdac_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the BRdac R package, please consider citing the relevant manuscript.

# References

Hector, E. C. and Song, P. X.-K. (2021). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, 116(534):805–818

Lars P. Hansen. (1982). Large sample properties of generalized method of moments estimators. Econometrica, 50(4), 1029-1054.

Huser, R. and Davison, A. C. (2014). Space–time modelling of extreme events. Journal of the Royal Statistical Society, Series B, 76(2):439–461.

Padoan, S., Ribatet, M., and Sisson, S. (2010). Likelihood-based inference for max-stable processes. Journal of the American Statistical Association, 105(489):263–277.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  SpatialExtremes
}
\examples{

\dontrun{
###################################################################
## Example without varying coefficient model on location and scale
###################################################################

alpha <- 1
phi <- 10
N <- 500
dim <- 10

s <- as.matrix(expand.grid(1:dim,1:dim))
S <- nrow(s)

beta_1 <- c(0.5,0.5)
beta_2 <- 1.5
beta_3 <- 0.2
p_1 <- length(beta_1)
p_2 <- length(beta_2)
p_3 <- length(beta_3)
p <- p_1 + p_2 + p_3

X <- SpatialExtremes::rmaxstab(n=N, coord=s, cov.mod="brown", grid=FALSE, 
                               range=phi, smooth=alpha)

covariates_1 <- as.matrix(do.call(rbind, lapply(1:N, function(x) return(s))))
covariates_2 <- covariates_3 <- as.matrix(rep(1,N*S))
z_1 <- z_constructor(covariates_1, S, N, p_1)
z_2 <- z_constructor(covariates_2, S, N, p_2)
z_3 <- z_constructor(covariates_3, S, N, p_3)
mu <- apply(z_1, 1, function(x) x\%*\%beta_1)
sigma <- apply(z_2, 1, function(x) exp(x\%*\%beta_2))
xi <- apply(z_3, 1, function(x) x\%*\%beta_3)
Y <- mu + sigma/xi*(X^xi - 1)

height_block_dim <- 5
width_block_dim <- 5
indicator <- rep(0,S)
height_num <- dim/height_block_dim
width_num <- dim/width_block_dim
for(k in 1:height_num){
  for(j in 1:width_num){
    indicator[s[,1] > (k-1)*height_block_dim & s[,1]<=k*height_block_dim & 
                s[,2] > (j-1)*width_block_dim & s[,2]<=j*width_block_dim] <- 
                (k-1)*width_num+j
  }
}

## This will run in parallel on K nodes
output  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, quantile=0.9, 
                 locs=s, indicator, cluster=length(unique(indicator)))

## This will run sequentially
output  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, quantile=0.9, 
                 locs=s, indicator, cluster=1)

###################################################################
## Example with varying coefficient model on location and scale
###################################################################

alpha <- 1
phi <- 10
N <- 1000
dim <- 10

s <- as.matrix(expand.grid(1:dim,1:dim))
s[,1] <- as.numeric(s[,1])
s[,2] <- as.numeric(s[,2])
s <- unname(s)
S <- nrow(s)

beta_1 <- geoR::grf(dim^2, grid=s, nsim=1, cov.model="matern", cov.pars=c(1, 1), 
                    kappa=6, nugget=0, RF=FALSE)$data
beta_2 <- geoR::grf(dim^2, grid=s, nsim=1, cov.model="matern", cov.pars=c(1, 1), 
                    kappa=6, nugget=0, RF=FALSE)$data
beta_3 <- 0.2
p_1 <- length(beta_1)
p_2 <- length(beta_2)
p_3 <- length(beta_3)
p <- p_1 + p_2 + p_3

X <- SpatialExtremes::rmaxstab(n=N, coord=s, cov.mod="brown", grid=FALSE, 
                               range=phi, smooth=alpha)

covariates_3 <- as.matrix(rep(1,N*S))
z_3 <- z_constructor(covariates_3, S, N, p_3)
mu <- t(matrix(rep(beta_1,N), S, N))
sigma <- t(matrix(rep(exp(beta_2),N), S, N))
xi <- apply(z_3, 1, function(x) x\%*\%beta_3)
Y <- mu + sigma/xi*(X^xi - 1)

height_block_dim <- 5
width_block_dim <- 5
indicator <- rep(0,S)
height_num <- dim/height_block_dim
width_num <- dim/width_block_dim
for(k in 1:height_num){
  for(j in 1:width_num){
    indicator[s[,1] > (k-1)*height_block_dim & s[,1]<=k*height_block_dim & 
                s[,2] > (j-1)*width_block_dim & s[,2]<=j*width_block_dim] <- 
                (k-1)*width_num+j
  }
}

## This will run in parallel on K nodes
output  <- BRdacVCMmusigma(y=Y, covariates_3, quantile=0.85, locs=s, indicator, 
                           lambda_1=seq(0,0.2,0.1), lambda_2=seq(0,0.2,0.1),
                           cluster=length(unique(indicator)))
                           
## This will run sequentially
output  <- BRdacVCMmusigma(y=Y, covariates_3, quantile=0.85, locs=s, indicator, 
                           lambda_1=seq(0,0.2,0.1), lambda_2=seq(0,0.2,0.1),
                           cluster=1)
}

}
