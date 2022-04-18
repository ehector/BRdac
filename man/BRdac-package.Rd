\name{BRdac-package}
\alias{BRdac-package}
\alias{BRdac}
\docType{package}
\title{\packageTitle{BRdac}}
\description{\packageDescription{BRdac}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{BRdac}
  \packageIndices{BRdac}

  This package contains the functions to perform regression analysis of extreme spatial data using a divide-and-conquer approach in the Brown-Resnick max stable process model. The main functions are:
    - BRdac: estimates regression parameters in non spatially-varying coefficient model (Section 3)
    - BRdacVCMmusigma: estimates regression parameters in spatially-varying coefficient model (Section 4). One Gaussian basis function is used within each block with 10 knots chosen by minimizing a geometric space-filling criterion (Royle and Nychka, 1998). For more control over the spatially varying coefficient model, users should create their own covariate expansions and use the BRdac function.
    - BRdacObjective: performs a goodness-of-fit test by evaluating the GMM objective function at the supplied model output and comparing it to the appropriate chi squared quantile.
}
\author{
  \packageAuthor{BRdac}
  
  Maintainer: \packageMaintainer{BRdac}
}
\references{
  Hector, E. C. and Song, P. X.-K. (2020). Doubly distributed supervised learning and infer- ence with high-dimensional correlated outcomes. Journal of Machine Learning Research, 21:1–35.
Huser, R. and Davison, A. C. (2013). Composite likelihood estimation for the Brown- Resnick process. Biometrika, 100(2):511–518.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  SpatialExtremes
}
\examples{

\dontrun{
alpha <- 1
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

X <- SpatialExtremes::rmaxstab(n=N, coord=s, cov.mod="brown", grid=FALSE, range=phi, smooth=alpha)

covariates_1 <- as.matrix(do.call(rbind, lapply(1:N, function(x) return(s))))
covariates_2 <- covariates_3 <- as.matrix(rep(1,N*S))
z_1 <- z_constructor(covariates_1, S, N, p_1)
z_2 <- z_constructor(covariates_2, S, N, p_2)
z_3 <- z_constructor(covariates_3, S, N, p_3)
mu <- apply(z_1, 1, function(x) x%*%beta_1)
sigma <- apply(z_2, 1, function(x) exp(x%*%beta_2))
xi <- apply(z_3, 1, function(x) x%*%beta_3)
Y <- mu + sigma/xi*(X^xi - 1)

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
output  <- BRdac(y=Y, covariates_1, covariates_2, covariates_3, quantile=0.9, locs=s, indicator)
}

}
