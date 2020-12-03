# BRdac
Divide-and-conquer approach for Brown-Resnick models

This is a repository for the code to simulate spatial data from the Brown-Resnick model and estimate model parameters using a divide-and-conquer procedure. The repository contains three files:
- simulation-setup.R: this file defines the parameters for a simulation study. Currently, this is where the exploration is happening.
- R-wrappers.R: this file defines the R wrapper functions that call the Rcpp composite likelihood and composite score functions.
- CL-funcs.cpp: this file defines the Rcpp functions that compute the pairwise composite likelihood and composite score for the Brown-Resnick model with unit Fréchet margins, and the pairwise composite likelihood and composite score functions for the transformed data with models for the location, scale and shape.
- BRdac_func.R: this file defines the R function for the divide-and-conquer estimation of model parameters.
