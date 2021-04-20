// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <sstream>
using namespace Rcpp; 
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::export]]
arma::cube z_constructor(const arma::mat& covariates, const int& S, const int& N, const int& p){
  arma::cube z(S, N, p);
  for(int q=0; q<p; q++){
    z.slice(q) = reshape(covariates.col(q), S, N);
  }
  return z;
}

// [[Rcpp::export]]
double logPCL_BR_thresh(const double& alpha, const double& phi, const double& u, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  
  double SUM = 0;
  double f;
  int n = x.n_rows;
  double distance = 0;
  for(int d=0; d<loc_1.size(); d++){
    distance += pow(loc_1(d) - loc_2(d), 2.0);
  }
  distance = pow(distance, 0.5);
  
  double a_12 = pow(2*pow(distance/phi, alpha), 0.5);
  double x_1, x_1_sq, x_2, x_2_sq, A_k1, A_k2, V, V_x1, V_x2, V_x1_x2;
  
  for(int i=0; i<n; i++){
    x_1 = x(i,0);
    x_1_sq = pow(x_1, 2.0) ;
    x_2 = x(i,1);
    x_2_sq = pow(x_2, 2.0);
    
    if(x_1 <= u & x_2 <= u){
      V = (1/u) * R::pnorm(a_12/2, 0, 1, true, false) + (1/u) * R::pnorm(a_12/2, 0, 1, true, false) ;
      SUM -= V; 
    }
    if(x_1 >u & x_2 <= u){
      V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/u), 0, 1, true, false) + 
        (1/u) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/u), 0, 1, true, false);
      
      V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/u), 0, 1, true, false) -
        1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/u), 2.0)/2) +
        1/(a_12*x_1*u*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/u), 2.0)/2);
      
      SUM += log(-V_x1) - V; 
    }
    if(x_1 <=u & x_2 > u){
      V = (1/u) * R::pnorm(a_12/2 - (1/a_12)*log(u/x_2), 0, 1, true, false) + 
        (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(u/x_2), 0, 1, true, false);
      
      V_x2 = 1/(a_12*u*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(u/x_2), 2.0)/2) -
        1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(u/x_2), 0, 1, true, false) -
        1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(u/x_2), 2.0)/2);
      
      SUM += log(-V_x2) - V; 
    }
    if(x_1 > u & x_2 > u){
      A_k1 = a_12/2 - (1/a_12)*log(x_1/x_2);
      A_k2 = a_12/2 + (1/a_12)*log(x_1/x_2);
      
      V = (1/x_1) * R::pnorm(A_k1, 0, 1, true, false) + (1/x_2) * R::pnorm(A_k2, 0, 1, true, false) ;
      
      V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
        1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
        1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
      
      V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
        1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
        1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
      
      V_x1_x2 = 
        1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
        1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
        1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
        1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
      
      f = exp(-V) * (V_x1*V_x2 - V_x1_x2);
      
      SUM += log(f); 
    }
  }
  
  // return the SUM
  return SUM;
}

// [[Rcpp::export]]
double logCL_BR_thresh(const arma::vec& par, const double& u, const arma::mat& x, const arma::mat& locs){
  const double& alpha = par(0);
  const double& phi = par(1);
  int S = locs.n_rows;
  double SUM = 0;
  for(int s_1=0; s_1<S-1; s_1++){
    for(int s_2=s_1+1; s_2<S; s_2++){
      SUM += logPCL_BR_thresh(alpha, phi, u, join_rows(x.col(s_1), x.col(s_2)), locs.row(s_1).t(), locs.row(s_2).t()); 
    }
  }
  return -SUM;
}

// [[Rcpp::export]]
arma::mat Cscore_BR(const arma::vec& par, const arma::mat& x, const arma::mat& locs){
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  const double& alpha = par(0);
  const double& phi = par(1);
  int S = locs.n_rows;
  int n = x.n_rows;
  arma::mat output = zeros<mat>(2, n);
  double distance, a_12, a_12_alpha, a_12_phi, 
  x_1, x_2, x_1_sq, x_2_sq, A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V, 
  V_x1, V_x2, V_x1_x2, scale;
  arma::vec score_i, loc_1, loc_2;
  for(int i=0; i<n; i++){
    score_i = score_i = zeros<vec>(2);
    for(int s_1=0; s_1<S-1; s_1++){
      for(int s_2=s_1+1; s_2<S; s_2++){
        loc_1 = locs.row(s_1).t();
        loc_2 = locs.row(s_2).t();
        
        distance = 0;
        for(int d=0; d<loc_1.size(); d++){
          distance += pow(loc_1(d) - loc_2(d), 2.0);
        }
        distance = pow(distance, 0.5);
        a_12 = pow(2 * pow(distance / phi, alpha), 0.5);
        a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
        a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
        
        x_1 = x(i,s_1);
        x_2 = x(i,s_2);
        x_1_sq = pow(x_1, 2.0);
        x_2_sq = pow(x_2, 2.0);
        A_k1 = a_12/2 - (1/a_12)*log(x_1/x_2);
        A_k2 = a_12/2 + (1/a_12)*log(x_1/x_2);
        A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / x_2))/ 2; 
        A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / x_2))/ 2; 
        A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / x_2)) / 2;
        A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / x_2)) / 2;
        
        norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
        norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
        
        V = (1/x_1) * R::pnorm(A_k1, 0, 1, true, false) + (1/x_2) * R::pnorm(A_k2, 0, 1, true, false);
        
        V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
          1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
          1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2);
        
        V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
          1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
          1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2);
        
        V_x1_x2 = 
          1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
          1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
          1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
          1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
        
        scale = 1/(V_x1 * V_x2 - V_x1_x2);
        
        score_i(0) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
          scale * (
              V_x2*(norm_A_k1*(a_12_alpha/pow(a_12, 2.0) + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
                norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*x_2) ) +
                V_x1*(norm_A_k2*(a_12_alpha/pow(a_12, 2.0) + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
                norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*x_1*x_2) ) -
                
                (
                    norm_A_k1*( a_12_alpha/(pow(a_12, 2.0)*x_1_sq*x_2) + A_k1*A_k1_alpha/(a_12*x_1_sq*x_2) -
                      2*a_12_alpha*A_k1/(pow(a_12, 3.0)*x_1_sq*x_2) + A_k1_alpha/(pow(a_12, 2.0)*x_1_sq*x_2) -
                      pow(A_k1, 2.0)*A_k1_alpha/(pow(a_12, 2.0)*x_1_sq*x_2) ) +
                      norm_A_k2*( a_12_alpha/(pow(a_12, 2.0)*x_1*x_2_sq) + A_k2*A_k2_alpha/(a_12*x_1*x_2_sq) -
                      2*a_12_alpha*A_k2/(pow(a_12, 3.0)*x_1*x_2_sq) + A_k2_alpha/(pow(a_12, 2.0)*x_1*x_2_sq) -
                      pow(A_k2, 2.0)*A_k2_alpha/(pow(a_12, 2.0)*x_1*x_2_sq) )
                )
          );
        
        score_i(1) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
          scale * (
              V_x2*(norm_A_k1*(a_12_phi/pow(a_12, 2.0) + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
                norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*x_2) ) +
                V_x1*(norm_A_k2*(a_12_phi/pow(a_12, 2.0) + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
                norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*x_1*x_2) ) -
                
                (
                    norm_A_k1*( a_12_phi/(pow(a_12, 2.0)*x_1_sq*x_2) + A_k1*A_k1_phi/(a_12*x_1_sq*x_2) -
                      2*a_12_phi*A_k1/(pow(a_12, 3.0)*x_1_sq*x_2) + A_k1_phi/(pow(a_12, 2.0)*x_1_sq*x_2) -
                      pow(A_k1, 2.0)*A_k1_phi/(pow(a_12, 2.0)*x_1_sq*x_2) ) +
                      norm_A_k2*( a_12_phi/(pow(a_12, 2.0)*x_1*x_2_sq) + A_k2*A_k2_phi/(a_12*x_1*x_2_sq) -
                      2*a_12_phi*A_k2/(pow(a_12, 3.0)*x_1*x_2_sq) + A_k2_phi/(pow(a_12, 2.0)*x_1*x_2_sq) -
                      pow(A_k2, 2.0)*A_k2_phi/(pow(a_12, 2.0)*x_1*x_2_sq) )
                )
          );
      }
    }
    
    output.col(i) = score_i ;
  }
  
  return output;
}

// [[Rcpp::export]]
double logPCL_all_thresh(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, const double& u,
                  const arma::mat& y, const arma::mat& z_1_1, const arma::mat& z_2_1, const arma::mat& z_3_1 , 
                  const arma::mat& z_1_2, const arma::mat& z_2_2, const arma::mat& z_3_2, 
                  const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  
  double SUM = 0;
  double f, J_1, J_2;
  int n = y.n_rows;
  double distance = 0;
  for(int d=0; d<loc_1.size(); d++){
    distance += pow(loc_1(d) - loc_2(d), 2.0);
  }
  distance = pow(distance, 0.5);
  
  double a_12 = pow(2*pow(distance/phi, alpha), 0.5);
  double y_1, y_2, x_1, x_1_sq, x_2, x_2_sq, V, V_x1, V_x2, V_x1_x2, mu_1, sigma_1, xi_1, mu_2, sigma_2, xi_2, A_k1, A_k2;
  
  for(int i=0; i<n; i++){
    mu_1 = as_scalar(z_1_1.row(i)*beta_1);
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    mu_2 = as_scalar(z_1_2.row(i)*beta_1);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);
    
    y_1 = y(i,0);
    y_2 = y(i,1);
    x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
    x_1_sq = pow(x_1, 2.0) ;
    x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
    x_2_sq = pow(x_2, 2.0);
    
    J_1 = 1/sigma_1 * pow(1 + xi_1 * (y_1 - mu_1) / sigma_1, 1 / xi_1-1);
    J_2 = 1/sigma_2 * pow(1 + xi_2 * (y_2 - mu_2) / sigma_2, 1 / xi_2-1);
    
    SUM += log(J_1) + log(J_2);
    
    if(y_1 <= u & y_2 <= u){
      V = (1/u) * R::pnorm(a_12/2, 0, 1, true, false) + (1/u) * R::pnorm(a_12/2, 0, 1, true, false) ;
      SUM += - V; 
    }
    if(y_1 >u & y_2 <= u){
      V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/u), 0, 1, true, false) + 
        (1/u) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/u), 0, 1, true, false);
      
      V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/u), 0, 1, true, false) -
        1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/u), 2.0)/2) +
        1/(a_12*x_1*u*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/u), 2.0)/2);
      
      SUM += log(-V_x1) - V ; 
    }
    if(y_1 <=u & y_2 > u){
      V = (1/u) * R::pnorm(a_12/2 - (1/a_12)*log(u/x_2), 0, 1, true, false) + 
        (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(u/x_2), 0, 1, true, false);
      
      V_x2 = 1/(a_12*u*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(u/x_2), 2.0)/2) -
        1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(u/x_2), 0, 1, true, false) -
        1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(u/x_2), 2.0)/2);
      
      SUM += log(-V_x2) - V ; 
    }
    if(y_1 > u & y_2 > u){
      A_k1 = a_12/2 - (1/a_12)*log(x_1/x_2);
      A_k2 = a_12/2 + (1/a_12)*log(x_1/x_2);
      
      V = (1/x_1) * R::pnorm(A_k1, 0, 1, true, false) + (1/x_2) * R::pnorm(A_k2, 0, 1, true, false) ;
      
      V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
        1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
        1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
      
      V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
        1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
        1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
      
      V_x1_x2 = 
        1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
        1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
        1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
        1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
      
      f = exp(-V) * (V_x1*V_x2 - V_x1_x2);
      
      SUM += log(f) ; 
    }
  }
  // return the SUM
  return SUM;
}

// [[Rcpp::export]]
double logCL_all_thresh(const arma::vec& par, const double& u, const arma::mat& y, 
                 const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, const arma::mat& locs){
  int S = locs.n_rows;
  double SUM = 0;
  
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const double& alpha = par(0);
  const double& phi = par(1);
  const arma::vec& beta_1 = par.subvec(2, 1+p_1);
  const arma::vec& beta_2 = par.subvec(2+p_1, 1+p_1+p_2);
  const arma::vec& beta_3 = par.subvec(2+p_1+p_2, 1+p_1+p_2+p_3);
  
  //There is an RcppArmadillo bug in subsetting a cube, we need to take the transpose to preserve dimensions
  //when there is only one slice
  //check here for updates: https://github.com/RcppCore/RcppArmadillo/issues/299
  bool bug_1 = FALSE;
  bool bug_2 = FALSE;
  bool bug_3 = FALSE;
  if(z_1.n_slices == 1) bug_1 = TRUE;
  if(z_2.n_slices == 1) bug_2 = TRUE;
  if(z_3.n_slices == 1) bug_3 = TRUE;
  
  arma::mat z_1_1, z_2_1, z_3_1, z_1_2, z_2_2, z_3_2;
  for(int s_1=0; s_1<S-1; s_1++){
    for(int s_2=s_1+1; s_2<S; s_2++){
      z_1_1 = z_1.row(s_1);
      z_2_1 = z_2.row(s_1);
      z_3_1 = z_3.row(s_1);
      z_1_2 = z_1.row(s_2);
      z_2_2 = z_2.row(s_2);
      z_3_2 = z_3.row(s_2);
      
      if(bug_1) {
        z_1_1 = z_1_1.t();
        z_1_2 = z_1_2.t();
      }
      if(bug_2) {
        z_2_1 = z_2_1.t();
        z_2_2 = z_2_2.t();
      }
      if(bug_3) {
        z_3_1 = z_3_1.t();
        z_3_2 = z_3_2.t();
      }
      
      SUM += logPCL_all_thresh(alpha, phi, beta_1, beta_2, beta_3, u, join_rows(y.col(s_1), y.col(s_2)), 
                               z_1_1, z_2_1, z_3_1, z_1_2, z_2_2, z_3_2, 
                               locs.row(s_1).t(), locs.row(s_2).t());
    }
  }
  
  if(std::isnan(SUM)) SUM = R_NegInf;
  return -SUM;
}

// [[Rcpp::export]]
List Chessian_all_thresh(const arma::vec& par, const double& u, const arma::mat& y, 
                         const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, 
                         const arma::mat& locs){
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int S = locs.n_rows;
  int n = y.n_rows;
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const double& alpha = par(0);
  const double& phi = par(1);
  const arma::vec& beta_1 = par.subvec(2, 1+p_1);
  const arma::vec& beta_2 = par.subvec(2+p_1, 1+p_1+p_2);
  const arma::vec& beta_3 = par.subvec(2+p_1+p_2, 1+p_1+p_2+p_3);
  
  arma::mat score = zeros<mat>(2 + p_1 + p_2 + p_3, n);
  arma::mat sensitivity = zeros<mat>(2 + p_1 + p_2 + p_3, n*S*(S-1)/2);
  
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi,
  f, y_1, y_2, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
  A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V,
  V_x1, V_x2, V_x1_x1,V_x2_x2, V_x1_x2, V_x1_x2_x2, V_x1_x1_x2, deriv_1, deriv_2, scale,
  mu_1, sigma_1, xi_1, mu_2, sigma_2, xi_2;
  arma::vec score_i, loc_1, loc_2;
  
  //There is an RcppArmadillo bug in subsetting a cube, we need to take the transpose to preserve dimensions
  //when there is only one slice
  
  //check here for updates: https://github.com/RcppCore/RcppArmadillo/issues/299
  bool bug_1 = FALSE;
  bool bug_2 = FALSE;
  bool bug_3 = FALSE;
  if(z_1.n_slices == 1) bug_1 = TRUE;
  if(z_2.n_slices == 1) bug_2 = TRUE;
  if(z_3.n_slices == 1) bug_3 = TRUE;
  
  arma::mat z_1_1, z_2_1, z_3_1, z_1_2, z_2_2, z_3_2;
  
  int it = 0;
  for(int s_1=0; s_1<S-1; s_1++){
    for(int s_2=s_1+1; s_2<S; s_2++){
      for(int i=0; i<n; i++){
        score_i = zeros<vec>(2 + p_1 + p_2 + p_3);
        
        loc_1 = locs.row(s_1).t();
        loc_2 = locs.row(s_2).t();
        
        distance = 0;
        for(int d=0; d<loc_1.size(); d++){
          distance += pow(loc_1(d) - loc_2(d), 2.0);
        }
        distance = pow(distance, 0.5);
        
        a_12 = pow(2 * pow(distance / phi, alpha), 0.5);
        a_12_sq = pow(a_12, 2.0);
        a_12_cu = pow(a_12, 3.0);
        a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
        a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
        
        z_1_1 = z_1.row(s_1); 
        z_2_1 = z_2.row(s_1);
        z_3_1 = z_3.row(s_1);
        z_1_2 = z_1.row(s_2);
        z_2_2 = z_2.row(s_2);
        z_3_2 = z_3.row(s_2);
        
        if(bug_1) {
          z_1_1 = z_1_1.t();
          z_1_2 = z_1_2.t();
        }
        if(bug_2) {
          z_2_1 = z_2_1.t();
          z_2_2 = z_2_2.t();
        }
        if(bug_3) {
          z_3_1 = z_3_1.t();
          z_3_2 = z_3_2.t();
        }
        
        mu_1 = as_scalar(z_1_1.row(i)*beta_1);
        sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
        xi_1 = as_scalar(z_3_1.row(i)*beta_3);
        mu_2 = as_scalar(z_1_2.row(i)*beta_1);
        sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
        xi_2 = as_scalar(z_3_2.row(i)*beta_3);
        
        y_1 = y(i,s_1);
        y_2 = y(i,s_2);
        
        if(y_1 <= u & y_2 <= u){
          double cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          double cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          double A_k = a_12/2;
          double A_k_alpha = a_12_alpha/ 2; 
          double A_k_phi = a_12_phi/ 2;
          
          score_i(0) = -A_k_alpha * exp(-pow(A_k, 2.0)/2) /(u * sqrt_2_pi) -A_k_alpha * exp(-pow(A_k, 2.0)/2) /(u * sqrt_2_pi) ;
          
          score_i(1) = -A_k_phi * exp(-pow(A_k, 2.0)/2) /(u * sqrt_2_pi) -A_k_phi * exp(-pow(A_k, 2.0)/2) /(u * sqrt_2_pi);
          
          score_i.rows(2,p_1+1) = 
            -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
            -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
            ;
          
          score_i.rows(p_1+2, p_1+p_2+1) = 
            -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
            -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
            ;
          
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
            1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
            + 1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
            ;
        }
        if(y_1 >u & y_2 <= u){
          x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
          
          x_1_sq = pow(x_1, 2.0);
          x_1_cu = pow(x_1, 3.0);
          double u_sq = pow(u, 2.0);
          
          cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(x_1/u);
          A_k2 = a_12/2 + (1/a_12)*log(x_1/u);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / u))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / u))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / u)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / u)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*x_1*u*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*x_1*u*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/u_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*u_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x1_x1 = 
            2/x_1_cu * R::pnorm(A_k1, 0, 1, true, false) 
            + norm_A_k1*( 3/(a_12*x_1_cu) - A_k1/( a_12_sq*x_1_cu ) )
            - norm_A_k2*( 1/(a_12*x_1_sq*u) + A_k2/(a_12_sq*x_1_sq*u)) ;
              
          score_i(0) = (1/V_x1)*
            (norm_A_k1*(a_12_alpha/a_12_sq + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
            norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*u) ) - 
            A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(u * sqrt_2_pi);
              
          score_i(1) = (1/V_x1)*
            (norm_A_k1*(a_12_phi/a_12_sq + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
            norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*u) ) - 
            A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(u * sqrt_2_pi);
              
          score_i.rows(2,p_1+1) = 
                -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
                -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
                -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_1_1.row(i).t();
              
          score_i.rows(p_1+2, p_1+p_2+1) = 
                -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
                -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
                -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_2_1.row(i).t();
              
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
                1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                + 1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x1/V_x1 - V_x1) * z_3_1.row(i).t();
        }
        if(y_1 <= u & y_2 > u){
          x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
          
          double u_sq = pow(u, 2.0);
          x_2_sq = pow(x_2, 2.0);
          x_2_cu = pow(x_2, 3.0);
          
          cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(u/x_2);
          A_k2 = a_12/2 + (1/a_12)*log(u/x_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(u / x_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(u / x_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(u / x_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(u / x_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/u_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*u_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*u*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*u*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2_x2 =
            2/x_2_cu * R::pnorm(A_k2, 0, 1, true, false)
              - norm_A_k1*( 1/(a_12*u*x_2_sq) + A_k1/( a_12_sq*u*x_2_sq))
              + norm_A_k2*( 3/(a_12*x_2_cu) - A_k2/( a_12_sq*x_2_cu ) ) ;
              
          score_i(0) = (1/V_x2)*
          (norm_A_k2*(a_12_alpha/a_12_sq + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
          norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*u*x_2) ) - 
                A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(u * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
              
          score_i(1) = (1/V_x2)*
            (norm_A_k2*(a_12_phi/a_12_sq + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
            norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*u*x_2) ) - 
                A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(u * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
              
          score_i.rows(2,p_1+1) = 
                -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
                -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
                -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_1_2.row(i).t();
                
          score_i.rows(p_1+2, p_1+p_2+1) = 
                -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
                -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
                -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_2_2.row(i).t();
                
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
                1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                  + 1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                  +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x2_x2/V_x2 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 > u & y_2 > u){
          x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
          x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
          
          x_1_sq = pow(x_1, 2.0);
          x_1_cu = pow(x_1, 3.0);
          x_2_sq = pow(x_2, 2.0);
          x_2_cu = pow(x_2, 3.0);
          
          cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(x_1/x_2);
          A_k2 = a_12/2 + (1/a_12)*log(x_1/x_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / x_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / x_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / x_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / x_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V = (1/x_1) * R::pnorm(A_k1, 0, 1, true, false) + (1/x_2) * R::pnorm(A_k2, 0, 1, true, false) ;
          
          V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x1_x2 = 
            1/(a_12_sq*x_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*x_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
          
          V_x1_x1 = 
            2/x_1_cu * R::pnorm(A_k1, 0, 1, true, false) 
            + norm_A_k1*( 3/(a_12*x_1_cu) - A_k1/( a_12_sq*x_1_cu ) )
            - norm_A_k2*( 1/(a_12*x_1_sq*x_2) + A_k2/(a_12_sq*x_1_sq*x_2)) ;
            
          V_x2_x2 =
            2/x_2_cu * R::pnorm(A_k2, 0, 1, true, false)
              - norm_A_k1*( 1/(a_12*x_1*x_2_sq) + A_k1/( a_12_sq*x_1*x_2_sq))
              + norm_A_k2*( 3/(a_12*x_2_cu) - A_k2/( a_12_sq*x_2_cu ) ) ;
              
          V_x1_x2_x2 =
              norm_A_k1*( (1-pow(A_k1, 2.0))/(a_12_cu*x_1_sq*x_2_sq) + 1/(a_12*x_1_sq*x_2_sq) )
                + norm_A_k2* ( 2/(a_12*x_1*x_2_cu) - 3*A_k2/(a_12_sq*x_1*x_2_cu) + (pow(A_k2, 2.0)-1)/(a_12_cu*x_1*x_2_cu) ) ;
              
          V_x1_x1_x2 = 
              norm_A_k1* ( 2/(a_12*x_1_cu*x_2) - 3*A_k1/(a_12_sq*x_1_cu*x_2) + (pow(A_k1, 2.0)-1)/(a_12_cu*x_1_cu*x_2) )
                - norm_A_k2*( (pow(A_k2, 2.0)-1)/(a_12_cu*x_1_sq*x_2_sq) - 1/(a_12*x_1_sq*x_2_sq) ) ;
              
          f = exp(-V) * (V_x1*V_x2 - V_x1_x2);
              
          deriv_1 = exp(-V)*(-pow(V_x1, 2.0)*V_x2 + 2*V_x1*V_x1_x2 + V_x1_x1*V_x2 - V_x1_x1_x2)/f;
          deriv_2 = exp(-V)*(-V_x1*pow(V_x2, 2.0) + 2*V_x1_x2*V_x2 + V_x1*V_x2_x2 - V_x1_x2_x2)/f;
              
          scale = 1/(V_x1 * V_x2 - V_x1_x2);
              
          score_i(0) = -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
                scale * (
                    V_x2*(norm_A_k1*(a_12_alpha/a_12_sq + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
                      norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*x_2) ) +
                      V_x1*(norm_A_k2*(a_12_alpha/a_12_sq + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
                      norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*x_1*x_2) ) -
                      
                      (
                          norm_A_k1*( a_12_alpha/(a_12_sq*x_1_sq*x_2) + A_k1*A_k1_alpha/(a_12*x_1_sq*x_2) -
                            2*a_12_alpha*A_k1/(a_12_cu*x_1_sq*x_2) + A_k1_alpha/(a_12_sq*x_1_sq*x_2) -
                            pow(A_k1, 2.0)*A_k1_alpha/(a_12_sq*x_1_sq*x_2) ) +
                            norm_A_k2*( a_12_alpha/(a_12_sq*x_1*x_2_sq) + A_k2*A_k2_alpha/(a_12*x_1*x_2_sq) -
                            2*a_12_alpha*A_k2/(a_12_cu*x_1*x_2_sq) + A_k2_alpha/(a_12_sq*x_1*x_2_sq) -
                            pow(A_k2, 2.0)*A_k2_alpha/(a_12_sq*x_1*x_2_sq) )
                      )
                );
              
          score_i(1) = -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
                scale * (
                    V_x2*(norm_A_k1*(a_12_phi/a_12_sq + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
                      norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*x_2) ) +
                      V_x1*(norm_A_k2*(a_12_phi/a_12_sq + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
                      norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*x_1*x_2) ) -
                      
                      (
                          norm_A_k1*( a_12_phi/(a_12_sq*x_1_sq*x_2) + A_k1*A_k1_phi/(a_12*x_1_sq*x_2) -
                            2*a_12_phi*A_k1/(a_12_cu*x_1_sq*x_2) + A_k1_phi/(a_12_sq*x_1_sq*x_2) -
                            pow(A_k1, 2.0)*A_k1_phi/(a_12_sq*x_1_sq*x_2) ) +
                            norm_A_k2*( a_12_phi/(a_12_sq*x_1*x_2_sq) + A_k2*A_k2_phi/(a_12*x_1*x_2_sq) -
                            2*a_12_phi*A_k2/(a_12_cu*x_1*x_2_sq) + A_k2_phi/(a_12_sq*x_1*x_2_sq) -
                            pow(A_k2, 2.0)*A_k2_phi/(a_12_sq*x_1*x_2_sq) )
                      )
                );
              
          score_i.rows(2,p_1+1) = 
                -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
                -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
                -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_1_1.row(i).t()
                -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_1_2.row(i).t()
                ;
              
          score_i.rows(p_1+2, p_1+p_2+1) = 
                -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
                -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
                -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
                -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
                ;
              
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
                1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                + 1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
                ;
        }
         
        sensitivity.col(it) = score_i;
        score.col(i) += score_i;
        it += 1;
      }
    }
  }
  
  return List::create(Named("score") = score,
                      Named("sensitivity") = sensitivity);
}