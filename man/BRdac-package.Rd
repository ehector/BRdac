\name{BRdac-package}
\alias{BRdac-package}
\alias{BRdac}
\docType{package}
\title{\packageTitle{BRdac}}
\description{\packageDescription{BRdac}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{BRdac}
  \packageIndices{BRdac}

  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
  \packageAuthor{BRdac}
  
  Maintainer: \packageMaintainer{BRdac}
}
\references{
  This optional section can contain literature or other references for
  background information.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  Optional links to other man pages
}
\examples{

\dontrun{
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
