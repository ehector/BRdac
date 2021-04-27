# BRdac
Divide-and-conquer approach for Brown-Resnick models

This is a repository for the R package to estimate Brown-Resnick model parameters using a divide-and-conquer procedure. The R package's main files are:
- src/CL-funcs.cpp: this file defines the Rcpp functions that compute the pairwise composite likelihood and composite score for the Brown-Resnick model with unit Fréchet margins, and the pairwise composite likelihood and composite score functions for the transformed data with regression models for the location, scale and shape.
- R/BRdac_func.R: this file defines the R function for the divide-and-conquer estimation of model parameters.

The BRdac man file contains an example for running the regression model from the paper.
