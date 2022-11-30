# About
Divide-and-conquer approach for Brown-Resnick models

This is a repository for the R package to estimate Brown-Resnick and inverted Brown-Resnick max-stable model parameters using a divide-and-conquer procedure. The R package's main files are:
- src/CL-funcs.cpp: this file defines the Rcpp functions that compute the pairwise censored composite likelihood and censored composite score for the Brown-Resnick and inverted Brown-Resnick models with unit Fréchet margins, and the pairwise censored composite likelihood and censored composite score functions for the transformed data with regression models for the location, scale and shape.
- R/BRdac_func.R: this file defines the R function for the divide-and-conquer estimation of model parameters.

The BRdac man file contains an example for running the regression model from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The BRdac R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL BRdac_1.0-1.tar.gz
- from the downloaded and renamed BRdac folder as R CMD build BRdac and R CMD INSTALL BRdac_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the BRdac R package, please consider citing the relevant manuscript: E.C. Hector and B.J. Reich (2022+). Distributed inference for spatial extremes modeling in high dimensions. arXiv, arXiv:2204.14165.

# References

Emily C. Hector and Peter X.-K. Song (2020). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, pages 1–14. doi: 10.1080/01621459.2020.1736082.

Lars P. Hansen. (1982). Large sample properties of generalized method of moments estimators. Econometrica, 50(4), 1029-1054.

Huser, R. and Davison, A. C. (2014). Space–time modelling of extreme events. Journal of the Royal Statistical Society, Series B, 76(2):439–461.

Padoan, S., Ribatet, M., and Sisson, S. (2010). Likelihood-based inference for max-stable processes. Journal of the American Statistical Association, 105(489):263–277.
