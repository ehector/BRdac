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
arma::mat matrix_inv(const arma::mat& X){
  arma::mat X_inv = inv(X);
  return(X_inv);
}

// [[Rcpp::export]]
arma::cube z_constructor(const arma::mat& covariates, const int& S, const int& N, const int& p){
  arma::cube z(S, N, p);
  for(int q=0; q<p; q++){
    z.slice(q) = reshape(covariates.col(q), S, N);
  }
  return z;
}

// [[Rcpp::export]]
double logPCL_BR_thresh(const double& alpha, const double& phi, const double& u_1, const double u_2, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
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
    x_2 = x(i,1);
    if(std::isnan(x_1) | std::isnan(x_2)) continue;
    x_1_sq = pow(x_1, 2.0) ;
    x_2_sq = pow(x_2, 2.0);
    
    if(x_1 <= u_1 & x_2 <= u_2){
      V = (1/u_1) * R::pnorm(a_12/2 - (1/a_12)*log(u_1/u_2), 0, 1, true, false) + 
        (1/u_2) * R::pnorm(a_12/2 + (1/a_12)*log(u_1/u_2), 0, 1, true, false);
      SUM -= V; 
    }
    if(x_1 > u_1 & x_2 <= u_2){
      V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/u_2), 0, 1, true, false) + 
        (1/u_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/u_2), 0, 1, true, false);
      
      V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/u_2), 0, 1, true, false) -
        1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/u_2), 2.0)/2) +
        1/(a_12*x_1*u_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/u_2), 2.0)/2);
      
      SUM += log(-V_x1) - V; 
    }
    if(x_1 <= u_1 & x_2 > u_2){
      V = (1/u_1) * R::pnorm(a_12/2 - (1/a_12)*log(u_1/x_2), 0, 1, true, false) + 
        (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(u_1/x_2), 0, 1, true, false);
      
      V_x2 = 1/(a_12*u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(u_1/x_2), 2.0)/2) -
        1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(u_1/x_2), 0, 1, true, false) -
        1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(u_1/x_2), 2.0)/2);
      
      SUM += log(-V_x2) - V; 
    }
    if(x_1 > u_1 & x_2 > u_2){
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
double logCL_BR_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& x, const arma::mat& locs){
  const double& alpha = par(0);
  const double& phi = par(1);
  int S = locs.n_rows;
  double SUM = 0;
  for(int s_1=0; s_1<S-1; s_1++){
    for(int s_2=s_1+1; s_2<S; s_2++){
      SUM += logPCL_BR_thresh(alpha, phi, thresholds(s_1), thresholds(s_2), join_rows(x.col(s_1), x.col(s_2)), locs.row(s_1).t(), locs.row(s_2).t()); 
    }
  }
  if(alpha > 2 | alpha < 0 ) SUM = R_NegInf;
  if(phi < 0) SUM = R_NegInf;
  return -SUM;
}

// [[Rcpp::export]]
arma::mat Cscore_BR_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& x, const arma::mat& locs){
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  const double& alpha = par(0);
  const double& phi = par(1);
  int S = locs.n_rows;
  int n = x.n_rows;
  
  arma::mat score = zeros<mat>(2, n);
  
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi, 
  x_1, x_2, u_1, u_2, x_1_sq, x_2_sq, A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2,
  V_x1, V_x2, V_x1_x2, scale;
  arma::vec score_i, loc_1, loc_2;
  for(int i=0; i<n; i++){
    for(int s_1=0; s_1<S-1; s_1++){
      for(int s_2=s_1+1; s_2<S; s_2++){
        x_1 = x(i,s_1);
        x_2 = x(i,s_2);
        if(std::isnan(x_1) | std::isnan(x_2)) continue;
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
        
        u_1 = thresholds(s_1);
        u_2 = thresholds(s_2);
        if(x_1 <= u_1 & x_2 <= u_2){
          A_k1 = a_12/2 - (1/a_12)*log(u_1/u_2);
          A_k2 = a_12/2 + (1/a_12)*log(u_1/u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(u_1 / u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(u_1 / u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(u_1 / u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(u_1 / u_2)) / 2;
          
          score(0,i) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(u_2 * sqrt_2_pi);
          
          score(1,i) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(u_2 * sqrt_2_pi);
        }
        if(x_1 > u_1 & x_2 <= u_2){
          x_1_sq = pow(x_1, 2.0);
          
          A_k1 = a_12/2 - (1/a_12)*log(x_1/u_2);
          A_k2 = a_12/2 + (1/a_12)*log(x_1/u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / u_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*x_1*u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;

          score(0,i) += (1/V_x1)*
              (norm_A_k1*(a_12_alpha/a_12_sq + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
              norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*u_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(u_2 * sqrt_2_pi);
            
          score(1,i) += (1/V_x1)*
              (norm_A_k1*(a_12_phi/a_12_sq + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
              norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*u_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(u_2 * sqrt_2_pi);
        }
        if(x_1 <= u_1 & x_2 > u_2){
          x_2_sq = pow(x_2, 2.0);
          
          A_k1 = a_12/2 - (1/a_12)*log(u_1/x_2);
          A_k2 = a_12/2 + (1/a_12)*log(u_1/x_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(u_1 / x_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(u_1 / x_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(u_1 / x_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(u_1 / x_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x2 = 1/(a_12*u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
            
          score(0,i) += (1/V_x2)*
              (norm_A_k2*(a_12_alpha/a_12_sq + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
              norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*u_1*x_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
            
          score(1,i) += (1/V_x2)*
              (norm_A_k2*(a_12_phi/a_12_sq + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
              norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*u_1*x_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
        }
        if(x_1 > u_1 & x_2 > u_2){
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
              
          scale = 1/(V_x1 * V_x2 - V_x1_x2);
              
          score(0,i) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
              
          score(1,i) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
        }
      }
    }
  }
  
  return -score;
}

// [[Rcpp::export]]
double logPCL_all_thresh(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, 
                         const double& u_1, const double& u_2, const arma::mat& y, 
                         const arma::mat& z_1_1, const arma::mat& z_2_1, const arma::mat& z_3_1, 
                         const arma::mat& z_1_2, const arma::mat& z_2_2, const arma::mat& z_3_2, 
                         const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  
  double SUM = 0;
  double f=0;
  double J_1=0;
  double J_2=0;
  int n = y.n_rows;
  double distance = 0;
  for(int d=0; d<loc_1.size(); d++){
    distance += pow(loc_1(d) - loc_2(d), 2.0);
  }
  distance = pow(distance, 0.5);
  
  double a_12 = pow(2*pow(distance/phi, alpha), 0.5);
  double y_1=0;
  double y_2=0;
  double x_1=0;
  double x_1_sq=0;
  double x_2=0;
  double x_2_sq=0;
  double V=0;
  double V_x1=0;
  double V_x2=0;
  double V_x1_x2=0;
  double mu_1=0;
  double sigma_1=0;
  double xi_1=0;
  double mu_2=0;
  double sigma_2=0;
  double xi_2=0;
  double A_k1=0;
  double A_k2=0;
  double new_u_1=0;
  double new_u_2=0;
  
  for(int i=0; i<n; i++){
    mu_1 = as_scalar(z_1_1.row(i)*beta_1);
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    mu_2 = as_scalar(z_1_2.row(i)*beta_1);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);
    
    y_1 = y(i,0);
    y_2 = y(i,1);
    if(std::isnan(y_1) | std::isnan(y_2) | std::isnan(mu_1) | std::isnan(mu_2) | 
       std::isnan(sigma_1) | std::isnan(sigma_2) | std::isnan(xi_1) | std::isnan(xi_2)) continue;
    x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
    x_1_sq = pow(x_1, 2.0) ;
    x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
    x_2_sq = pow(x_2, 2.0);
    if(std::isnan(x_1) | std::isnan(x_2)) {
      SUM = R_NaN;
      break;
    }
    
    new_u_1 = pow(1+xi_1/sigma_1*(u_1-mu_1), 1/xi_1);
    new_u_2 = pow(1+xi_2/sigma_2*(u_2-mu_2), 1/xi_2);
    
    if(y_1 <= u_1 & y_2 <= u_2){
      A_k1 = a_12/2 - (1/a_12)*log(new_u_1/new_u_2);
      A_k2 = a_12/2 + (1/a_12)*log(new_u_1/new_u_2);
      
      V = (1/new_u_1) * R::pnorm(A_k1, 0, 1, true, false) + (1/new_u_2) * R::pnorm(A_k2, 0, 1, true, false) ;
      SUM += - V ; 
    }
    if(y_1 > u_1 & y_2 <= u_2){
      J_1 = 1/sigma_1 * pow(1 + xi_1 * (y_1 - mu_1) / sigma_1, 1 / xi_1-1);
      
      V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/new_u_2), 0, 1, true, false) + 
        (1/new_u_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/new_u_2), 0, 1, true, false);
      
      V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/new_u_2), 0, 1, true, false) -
        1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/new_u_2), 2.0)/2) +
        1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/new_u_2), 2.0)/2);
      
      SUM += log(-V_x1) - V + log(J_1); 
    }
    if(y_1 <= u_1 & y_2 > u_2){
      J_2 = 1/sigma_2 * pow(1 + xi_2 * (y_2 - mu_2) / sigma_2, 1 / xi_2-1);
      
      V = (1/new_u_1) * R::pnorm(a_12/2 - (1/a_12)*log(new_u_1/x_2), 0, 1, true, false) + 
        (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(new_u_1/x_2), 0, 1, true, false);
      
      V_x2 = 1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(new_u_1/x_2), 2.0)/2) -
        1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(new_u_1/x_2), 0, 1, true, false) -
        1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(new_u_1/x_2), 2.0)/2);
      
      SUM += log(-V_x2) - V + log(J_2); 
    }
    if(y_1 > u_1 & y_2 > u_2){
      J_1 = 1/sigma_1 * pow(1 + xi_1 * (y_1 - mu_1) / sigma_1, 1 / xi_1-1);
      J_2 = 1/sigma_2 * pow(1 + xi_2 * (y_2 - mu_2) / sigma_2, 1 / xi_2-1);
      
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
      
      SUM += log(f) + log(J_1) + log(J_2); 
    }
    if(std::isnan(SUM)) break;
  }
  // return the SUM
  return SUM;
}

// [[Rcpp::export]]
double logCL_all_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                 const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, const arma::mat& locs){
  int S = locs.n_rows;
  double SUM = 0;
  
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const double& alpha = 2*exp(par(0))/(1-exp(par(0)));
  const double& phi = exp(par(1));
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
      double result = logPCL_all_thresh(alpha, phi, beta_1, beta_2, beta_3, thresholds(s_1), thresholds(s_2), join_rows(y.col(s_1), y.col(s_2)), 
                                        z_1_1, z_2_1, z_3_1, z_1_2, z_2_2, z_3_2, 
                                        locs.row(s_1).t(), locs.row(s_2).t());
      if(std::isnan(SUM)) break;
      SUM += result;
    }
  }
  
  if(std::isnan(SUM)) SUM = R_NegInf;
  if(alpha > 2 | alpha < 0 ) SUM = R_NegInf;
  if(phi < 0) SUM = R_NegInf;
  return -SUM;
}

// [[Rcpp::export]]
double logCL_all_thresh_reparm(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                               const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, const arma::mat& locs){
  int S = locs.n_rows;
  double SUM = 0;
  
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const double& alpha = 2*exp(par(0))/(1+exp(par(0)));
  const double& phi = exp(par(1));
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
      double result = logPCL_all_thresh(alpha, phi, beta_1, beta_2, beta_3, thresholds(s_1), thresholds(s_2), join_rows(y.col(s_1), y.col(s_2)), 
                                        z_1_1, z_2_1, z_3_1, z_1_2, z_2_2, z_3_2, 
                                        locs.row(s_1).t(), locs.row(s_2).t());
      if(std::isnan(SUM)) break;
      SUM += result;
    }
  }
  
  if(std::isnan(SUM)) SUM = R_NegInf;
  if(alpha > 2 | alpha < 0 ) SUM = R_NegInf;
  if(phi < 0) SUM = R_NegInf;
  return -SUM;
}

// [[Rcpp::export]]
double logL_marg_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                        const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, const arma::mat& locs){
  int S = locs.n_rows;
  double SUM = 0;
  int n = y.n_rows;
  
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const arma::vec& beta_1 = par.subvec(0, p_1-1);
  const arma::vec& beta_2 = par.subvec(p_1, p_1+p_2-1);
  const arma::vec& beta_3 = par.subvec(p_1+p_2, p_1+p_2+p_3-1);
  
  //There is an RcppArmadillo bug in subsetting a cube, we need to take the transpose to preserve dimensions
  //when there is only one slice
  //check here for updates: https://github.com/RcppCore/RcppArmadillo/issues/299
  bool bug_1 = FALSE;
  bool bug_2 = FALSE;
  bool bug_3 = FALSE;
  if(z_1.n_slices == 1) bug_1 = TRUE;
  if(z_2.n_slices == 1) bug_2 = TRUE;
  if(z_3.n_slices == 1) bug_3 = TRUE;
  
  arma::mat z_1_s, z_2_s, z_3_s;
  double mu, sigma, xi, y_s, u_s, x_s;
  for(int i=0; i<n; i++){
    for(int s=0; s<S; s++){
      z_1_s = z_1.row(s);
      z_2_s = z_2.row(s);
      z_3_s = z_3.row(s);
      
      if(bug_1) {
        z_1_s = z_1_s.t();
      }
      if(bug_2) {
        z_2_s = z_2_s.t();
      }
      if(bug_3) {
        z_3_s = z_3_s.t();
      }
      
      mu = as_scalar(z_1_s.row(i)*beta_1);
      sigma = exp(as_scalar(z_2_s.row(i)*beta_2));
      xi = as_scalar(z_3_s.row(i)*beta_3);
      
      y_s = y(i,s);
      u_s = thresholds(s);
      if(std::isnan(y_s) | std::isnan(mu) | std::isnan(sigma) | std::isnan(xi)) continue;
      
      if(y_s <= u_s){
        if(abs(xi) <= 1e-6) x_s = exp(-(u_s-mu)/sigma);
        else x_s = pow(1+xi*(u_s - mu)/sigma, -1/xi);
        SUM += -x_s;
      }
      if(y_s > u_s){
        if(abs(xi) <= 1e-6) SUM += -log(sigma) - (y_s-mu)/sigma - exp(-(y_s-mu)/sigma); //x_s = exp(-(y_s-mu)/sigma);
        else SUM += -log(sigma) - (1+1/xi)*log(1+xi*(y_s - mu)/sigma) - pow(1+xi*(y_s - mu)/sigma, -1/xi); //x_s = pow(1+xi*(y_s - mu)/sigma, -1/xi);
      }
    }
  }
  
  if(std::isnan(SUM)) SUM = R_NegInf;
  return -SUM;
}

// [[Rcpp::export]]
List Chessian_all_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
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
  f, y_1, y_2, u_1, u_2, new_u_1, new_u_2, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
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
        if(std::isnan(y_1) | std::isnan(y_2) | std::isnan(mu_1) | std::isnan(mu_2) | 
           std::isnan(sigma_1) | std::isnan(sigma_2) | std::isnan(xi_1) | std::isnan(xi_2)) continue;
        u_1 = thresholds(s_1);
        u_2 = thresholds(s_2);
        
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
        
        new_u_1 = pow(1+xi_1/sigma_1*(u_1-mu_1), 1/xi_1);
        new_u_2 = pow(1+xi_2/sigma_2*(u_2-mu_2), 1/xi_2);
        
        double new_u_1_sq = pow(new_u_1, 2.0);
        double new_u_2_sq = pow(new_u_2, 2.0);
        
        if(y_1 <= u_1 & y_2 <= u_2){
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          score_i(0) = -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score_i(1) = -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score_i.rows(2,p_1+1) = V_x1* 1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_1_1.row(i).t() + 
            V_x2*1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_1_2.row(i).t();
          
          score_i.rows(p_1+2, p_1+p_2+1) = V_x1*(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_2_1.row(i).t() +
            V_x2*(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_2_2.row(i).t();
          
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
            - V_x1 * (-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * z_3_1.row(i).t()
            - V_x2 * (-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 <= u_2){
          x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
          
          x_1_sq = pow(x_1, 2.0);
          x_1_cu = pow(x_1, 3.0);
          double new_u_2_sq = pow(new_u_2, 2.0);
          
          cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(x_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(x_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x1_x1 = 
            2/x_1_cu * R::pnorm(A_k1, 0, 1, true, false) 
            + norm_A_k1*( 3/(a_12*x_1_cu) - A_k1/( a_12_sq*x_1_cu ) )
            - norm_A_k2*( 1/(a_12*x_1_sq*new_u_2) + A_k2/(a_12_sq*x_1_sq*new_u_2)) ;
            
          V_x1_x2 = 
            1/(a_12_sq*x_1_sq*new_u_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*x_1*new_u_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1_sq*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
              
          score_i(0) = (1/V_x1)*
            (norm_A_k1*(a_12_alpha/a_12_sq + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
            norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*new_u_2) ) - 
            A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
              
          score_i(1) = (1/V_x1)*
            (norm_A_k1*(a_12_phi/a_12_sq + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
            norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*new_u_2) ) - 
            A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
              
          score_i.rows(2,p_1+1) = 
                -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
                -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_1_1.row(i).t()
                -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_1_2.row(i).t();
              
          score_i.rows(p_1+2, p_1+p_2+1) = 
                -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
                -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_2_1.row(i).t()
                -(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_2_2.row(i).t();
              
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
                1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x1/V_x1 - V_x1) * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x1_x2/V_x1 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 <= u_1 & y_2 > u_2){
          x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
          
          double new_u_1_sq = pow(new_u_1, 2.0);
          x_2_sq = pow(x_2, 2.0);
          x_2_cu = pow(x_2, 3.0);
          
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/x_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/x_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2_x2 =
            2/x_2_cu * R::pnorm(A_k2, 0, 1, true, false)
              - norm_A_k1*( 1/(a_12*new_u_1*x_2_sq) + A_k1/( a_12_sq*new_u_1*x_2_sq))
              + norm_A_k2*( 3/(a_12*x_2_cu) - A_k2/( a_12_sq*x_2_cu ) ) ;
              
          V_x1_x2 = 
              1/(a_12_sq*new_u_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
              1/(a_12_sq*new_u_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
              1/(a_12*new_u_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
              1/(a_12*new_u_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;    
              
          score_i(0) = (1/V_x2)*
          (norm_A_k2*(a_12_alpha/a_12_sq + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
          norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*new_u_1*x_2) ) - 
                A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
              
          score_i(1) = (1/V_x2)*
            (norm_A_k2*(a_12_phi/a_12_sq + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
            norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*new_u_1*x_2) ) - 
                A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
              
          score_i.rows(2,p_1+1) = 
                -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
                -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_1_1.row(i).t()
                -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_1_2.row(i).t();
                
          score_i.rows(p_1+2, p_1+p_2+1) = 
                -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
                -(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_2_1.row(i).t()
                -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_2_2.row(i).t();
                
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
                  1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                  +(-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x2/V_x2 - V_x1) * z_3_1.row(i).t()
                  +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x2_x2/V_x2 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 > u_2){
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

// [[Rcpp::export]]
List Chessian_all_thresh_reparm(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                         const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, 
                         const arma::mat& locs){
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int S = locs.n_rows;
  int n = y.n_rows;
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const double& alpha = 2*exp(par(0))/(1+exp(par(0)));
  const double& phi = exp(par(1));
  const arma::vec& beta_1 = par.subvec(2, 1+p_1);
  const arma::vec& beta_2 = par.subvec(2+p_1, 1+p_1+p_2);
  const arma::vec& beta_3 = par.subvec(2+p_1+p_2, 1+p_1+p_2+p_3);
  
  arma::mat score = zeros<mat>(2 + p_1 + p_2 + p_3, n);
  arma::mat sensitivity = zeros<mat>(2 + p_1 + p_2 + p_3, n*S*(S-1)/2);
  
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi,
  f, y_1, y_2, u_1, u_2, new_u_1, new_u_2, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
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
        if(std::isnan(y_1) | std::isnan(y_2) | std::isnan(mu_1) | std::isnan(mu_2) | 
           std::isnan(sigma_1) | std::isnan(sigma_2) | std::isnan(xi_1) | std::isnan(xi_2)) continue;
        u_1 = thresholds(s_1);
        u_2 = thresholds(s_2);
        
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
        
        new_u_1 = pow(1+xi_1/sigma_1*(u_1-mu_1), 1/xi_1);
        new_u_2 = pow(1+xi_2/sigma_2*(u_2-mu_2), 1/xi_2);
        
        double new_u_1_sq = pow(new_u_1, 2.0);
        double new_u_2_sq = pow(new_u_2, 2.0);
        
        if(y_1 <= u_1 & y_2 <= u_2){
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          score_i(0) = -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score_i(1) = -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score_i.rows(2,p_1+1) = V_x1* 1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_1_1.row(i).t() + 
            V_x2*1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_1_2.row(i).t();
          
          score_i.rows(p_1+2, p_1+p_2+1) = V_x1*(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_2_1.row(i).t() +
            V_x2*(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_2_2.row(i).t();
          
          score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
            - V_x1 * (-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * z_3_1.row(i).t()
            - V_x2 * (-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 <= u_2){
          x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
          
          x_1_sq = pow(x_1, 2.0);
          x_1_cu = pow(x_1, 3.0);
          double new_u_2_sq = pow(new_u_2, 2.0);
          
          cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(x_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(x_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x1_x1 = 
            2/x_1_cu * R::pnorm(A_k1, 0, 1, true, false) 
            + norm_A_k1*( 3/(a_12*x_1_cu) - A_k1/( a_12_sq*x_1_cu ) )
            - norm_A_k2*( 1/(a_12*x_1_sq*new_u_2) + A_k2/(a_12_sq*x_1_sq*new_u_2)) ;
            
            V_x1_x2 = 
            1/(a_12_sq*x_1_sq*new_u_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*x_1*new_u_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1_sq*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
            
            score_i(0) = (1/V_x1)*
              (norm_A_k1*(a_12_alpha/a_12_sq + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
              norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*new_u_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
            
            score_i(1) = (1/V_x1)*
              (norm_A_k1*(a_12_phi/a_12_sq + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
              norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*new_u_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
            
            score_i.rows(2,p_1+1) = 
              -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
              -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_1_1.row(i).t()
              -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_1_2.row(i).t();
              
              score_i.rows(p_1+2, p_1+p_2+1) = 
              -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
              -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_2_1.row(i).t()
              -(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_2_2.row(i).t();
              
              score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
              1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x1/V_x1 - V_x1) * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x1_x2/V_x1 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 <= u_1 & y_2 > u_2){
          x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
          
          double new_u_1_sq = pow(new_u_1, 2.0);
          x_2_sq = pow(x_2, 2.0);
          x_2_cu = pow(x_2, 3.0);
          
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/x_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/x_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2_x2 =
            2/x_2_cu * R::pnorm(A_k2, 0, 1, true, false)
            - norm_A_k1*( 1/(a_12*new_u_1*x_2_sq) + A_k1/( a_12_sq*new_u_1*x_2_sq))
            + norm_A_k2*( 3/(a_12*x_2_cu) - A_k2/( a_12_sq*x_2_cu ) ) ;
            
            V_x1_x2 = 
            1/(a_12_sq*new_u_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*new_u_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*new_u_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*new_u_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;    
            
            score_i(0) = (1/V_x2)*
              (norm_A_k2*(a_12_alpha/a_12_sq + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
              norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*new_u_1*x_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
            
            score_i(1) = (1/V_x2)*
              (norm_A_k2*(a_12_phi/a_12_sq + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
              norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*new_u_1*x_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
            
            score_i.rows(2,p_1+1) = 
              -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
              -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_1_1.row(i).t()
              -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_1_2.row(i).t();
              
              score_i.rows(p_1+2, p_1+p_2+1) = 
              -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
              -(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_2_1.row(i).t()
              -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_2_2.row(i).t();
              
              score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
              1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                +(-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x2/V_x2 - V_x1) * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x2_x2/V_x2 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 > u_2){
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
  
  score.row(0) = score.row(0) * alpha/(1+exp(par(0)));
  score.row(1) = score.row(1) * phi;
  
  sensitivity.row(0) = sensitivity.row(0) * alpha/(1+exp(par(0)));
  sensitivity.row(1) = sensitivity.row(1) * phi;
  
  return List::create(Named("score") = score,
                      Named("sensitivity") = sensitivity);
}

// [[Rcpp::export]]
arma::vec score_all_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
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
  
  arma::mat score = zeros<vec>(2 + p_1 + p_2 + p_3);
  
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi,
  f, y_1, y_2, u_1, u_2, new_u_1, new_u_2, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
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
  
  for(int s_1=0; s_1<S-1; s_1++){
    for(int s_2=s_1+1; s_2<S; s_2++){
      for(int i=0; i<n; i++){
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
        if(std::isnan(y_1) | std::isnan(y_2) | std::isnan(mu_1) | std::isnan(mu_2) | 
           std::isnan(sigma_1) | std::isnan(sigma_2) | std::isnan(xi_1) | std::isnan(xi_2)) continue;
        u_1 = thresholds(s_1);
        u_2 = thresholds(s_2);
        
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
        
        new_u_1 = pow(1+xi_1/sigma_1*(u_1-mu_1), 1/xi_1);
        new_u_2 = pow(1+xi_2/sigma_2*(u_2-mu_2), 1/xi_2);
        
        double new_u_1_sq = pow(new_u_1, 2.0);
        double new_u_2_sq = pow(new_u_2, 2.0);
        
        if(y_1 <= u_1 & y_2 <= u_2){
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          score(0) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score(1) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score.rows(2,p_1+1) += V_x1* 1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_1_1.row(i).t() + 
            V_x2*1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_1_2.row(i).t();
          
          score.rows(p_1+2, p_1+p_2+1) += V_x1*(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_2_1.row(i).t() +
            V_x2*(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_2_2.row(i).t();
          
          score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
            - V_x1 * (-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * z_3_1.row(i).t()
            - V_x2 * (-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 <= u_2){
          x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
          
          x_1_sq = pow(x_1, 2.0);
          x_1_cu = pow(x_1, 3.0);
          double new_u_2_sq = pow(new_u_2, 2.0);
          
          cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(x_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(x_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x1_x1 = 
            2/x_1_cu * R::pnorm(A_k1, 0, 1, true, false) 
            + norm_A_k1*( 3/(a_12*x_1_cu) - A_k1/( a_12_sq*x_1_cu ) )
            - norm_A_k2*( 1/(a_12*x_1_sq*new_u_2) + A_k2/(a_12_sq*x_1_sq*new_u_2)) ;
            
          V_x1_x2 = 
            1/(a_12_sq*x_1_sq*new_u_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*x_1*new_u_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1_sq*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
            
          score(0) += (1/V_x1)*
              (norm_A_k1*(a_12_alpha/a_12_sq + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
              norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*new_u_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
            
          score(1) += (1/V_x1)*
              (norm_A_k1*(a_12_phi/a_12_sq + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
              norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*new_u_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
            
          score.rows(2,p_1+1) += 
              -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
              -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_1_1.row(i).t()
              -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_1_2.row(i).t();
              
          score.rows(p_1+2, p_1+p_2+1) += 
              -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
              -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_2_1.row(i).t()
              -(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_2_2.row(i).t();
              
          score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
              1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x1/V_x1 - V_x1) * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x1_x2/V_x1 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 <= u_1 & y_2 > u_2){
          x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
          
          double new_u_1_sq = pow(new_u_1, 2.0);
          x_2_sq = pow(x_2, 2.0);
          x_2_cu = pow(x_2, 3.0);
          
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/x_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/x_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2_x2 =
            2/x_2_cu * R::pnorm(A_k2, 0, 1, true, false)
            - norm_A_k1*( 1/(a_12*new_u_1*x_2_sq) + A_k1/( a_12_sq*new_u_1*x_2_sq))
            + norm_A_k2*( 3/(a_12*x_2_cu) - A_k2/( a_12_sq*x_2_cu ) ) ;
            
          V_x1_x2 = 
            1/(a_12_sq*new_u_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*new_u_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*new_u_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*new_u_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;    
            
          score(0) += (1/V_x2)*
              (norm_A_k2*(a_12_alpha/a_12_sq + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
              norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*new_u_1*x_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
            
          score(1) += (1/V_x2)*
              (norm_A_k2*(a_12_phi/a_12_sq + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
              norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*new_u_1*x_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
            
          score.rows(2,p_1+1) += 
              -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
              -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_1_1.row(i).t()
              -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_1_2.row(i).t();
              
          score.rows(p_1+2, p_1+p_2+1) += 
              -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
              -(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_2_1.row(i).t()
              -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_2_2.row(i).t();
              
          score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
              1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                +(-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x2/V_x2 - V_x1) * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x2_x2/V_x2 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 > u_2){
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
              
          score(0) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
              
          score(1) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
              
            score.rows(2,p_1+1) += 
                -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
                -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
                -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_1_1.row(i).t()
                -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_1_2.row(i).t()
                ;
              
            score.rows(p_1+2, p_1+p_2+1) += 
                -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
                -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
                -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
                -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
                ;
              
            score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
                1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                + 1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
                ;
        }
      }
    }
  }
  
  return -score;
}

// [[Rcpp::export]]
arma::vec score_all_thresh_reparm(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                           const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, 
                           const arma::mat& locs){
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int S = locs.n_rows;
  int n = y.n_rows;
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const double& alpha = 2*exp(par(0))/(1+exp(par(0)));
  const double& phi = exp(par(1));
  const arma::vec& beta_1 = par.subvec(2, 1+p_1);
  const arma::vec& beta_2 = par.subvec(2+p_1, 1+p_1+p_2);
  const arma::vec& beta_3 = par.subvec(2+p_1+p_2, 1+p_1+p_2+p_3);
  
  arma::mat score = zeros<vec>(2 + p_1 + p_2 + p_3);
  
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi,
  f, y_1, y_2, u_1, u_2, new_u_1, new_u_2, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
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
  
  for(int s_1=0; s_1<S-1; s_1++){
    for(int s_2=s_1+1; s_2<S; s_2++){
      for(int i=0; i<n; i++){
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
        if(std::isnan(y_1) | std::isnan(y_2) | std::isnan(mu_1) | std::isnan(mu_2) | 
           std::isnan(sigma_1) | std::isnan(sigma_2) | std::isnan(xi_1) | std::isnan(xi_2)) continue;
        u_1 = thresholds(s_1);
        u_2 = thresholds(s_2);
        
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
        
        new_u_1 = pow(1+xi_1/sigma_1*(u_1-mu_1), 1/xi_1);
        new_u_2 = pow(1+xi_2/sigma_2*(u_2-mu_2), 1/xi_2);
        
        double new_u_1_sq = pow(new_u_1, 2.0);
        double new_u_2_sq = pow(new_u_2, 2.0);
        
        if(y_1 <= u_1 & y_2 <= u_2){
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / new_u_2)) / 2;
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          score(0) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score(1) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
          
          score.rows(2,p_1+1) += V_x1* 1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_1_1.row(i).t() + 
            V_x2*1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_1_2.row(i).t();
          
          score.rows(p_1+2, p_1+p_2+1) += V_x1*(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * z_2_1.row(i).t() +
            V_x2*(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * z_2_2.row(i).t();
          
          score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
            - V_x1 * (-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * z_3_1.row(i).t()
            - V_x2 * (-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 <= u_2){
          x_1 = pow( 1 + xi_1*(y_1 - mu_1)/sigma_1, 1/xi_1);
          
          x_1_sq = pow(x_1, 2.0);
          x_1_cu = pow(x_1, 3.0);
          double new_u_2_sq = pow(new_u_2, 2.0);
          
          cent_1 = 1 + xi_1 * (y_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (u_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(x_1/new_u_2);
          A_k2 = a_12/2 + (1/a_12)*log(x_1/new_u_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(x_1 / new_u_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/x_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*x_1*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/new_u_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x1_x1 = 
            2/x_1_cu * R::pnorm(A_k1, 0, 1, true, false) 
            + norm_A_k1*( 3/(a_12*x_1_cu) - A_k1/( a_12_sq*x_1_cu ) )
            - norm_A_k2*( 1/(a_12*x_1_sq*new_u_2) + A_k2/(a_12_sq*x_1_sq*new_u_2)) ;
            
            V_x1_x2 = 
            1/(a_12_sq*x_1_sq*new_u_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*x_1*new_u_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1*new_u_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*x_1_sq*new_u_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;
            
            score(0) += (1/V_x1)*
              (norm_A_k1*(a_12_alpha/a_12_sq + A_k1*A_k1_alpha/a_12-A_k1_alpha)/x_1_sq -
              norm_A_k2*(a_12_alpha/a_12 + A_k2*A_k2_alpha)/(a_12*x_1*new_u_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
            
            score(1) += (1/V_x1)*
              (norm_A_k1*(a_12_phi/a_12_sq + A_k1*A_k1_phi/a_12-A_k1_phi)/x_1_sq -
              norm_A_k2*(a_12_phi/a_12 + A_k2*A_k2_phi)/(a_12*x_1*new_u_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(new_u_2 * sqrt_2_pi);
            
            score.rows(2,p_1+1) += 
              -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
              -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_1_1.row(i).t()
              -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_1_2.row(i).t();
              
              score.rows(p_1+2, p_1+p_2+1) += 
              -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
              -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x1/V_x1 - V_x1) * z_2_1.row(i).t()
              -(u_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x1_x2/V_x1 - V_x2) * z_2_2.row(i).t();
              
              score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
              1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x1/V_x1 - V_x1) * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (u_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x1_x2/V_x1 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 <= u_1 & y_2 > u_2){
          x_2 = pow( 1 + xi_2*(y_2 - mu_2)/sigma_2, 1/xi_2);
          
          double new_u_1_sq = pow(new_u_1, 2.0);
          x_2_sq = pow(x_2, 2.0);
          x_2_cu = pow(x_2, 3.0);
          
          cent_1 = 1 + xi_1 * (u_1 - mu_1) / sigma_1;
          cent_2 = 1 + xi_2 * (y_2 - mu_2) / sigma_2;
          
          A_k1 = a_12/2 - (1/a_12)*log(new_u_1/x_2);
          A_k2 = a_12/2 + (1/a_12)*log(new_u_1/x_2);
          A_k1_alpha = a_12_alpha * (1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k2_alpha = a_12_alpha * (1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2))/ 2; 
          A_k1_phi = a_12_phi * ( 1 + pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          A_k2_phi = a_12_phi * ( 1 - pow(distance / phi, -alpha) * log(new_u_1 / x_2)) / 2;
          
          norm_A_k1 = exp(-pow(A_k1, 2.0)/2)/pow(2*pi, 0.5);
          norm_A_k2 = exp(-pow(A_k2, 2.0)/2)/pow(2*pi, 0.5);
          
          V_x1 = -1/new_u_1_sq * R::pnorm(A_k1, 0, 1, true, false) -
            1/(a_12*new_u_1_sq*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2 = 1/(a_12*new_u_1*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) -
            1/x_2_sq * R::pnorm(A_k2, 0, 1, true, false) -
            1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) ;
          
          V_x2_x2 =
            2/x_2_cu * R::pnorm(A_k2, 0, 1, true, false)
            - norm_A_k1*( 1/(a_12*new_u_1*x_2_sq) + A_k1/( a_12_sq*new_u_1*x_2_sq))
            + norm_A_k2*( 3/(a_12*x_2_cu) - A_k2/( a_12_sq*x_2_cu ) ) ;
            
            V_x1_x2 = 
            1/(a_12_sq*new_u_1_sq*x_2*pow(2*pi, 0.5)) * (A_k1) * exp(-pow(A_k1, 2.0)/2) +
            1/(a_12_sq*new_u_1*x_2_sq*pow(2*pi, 0.5)) * (A_k2) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*new_u_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(A_k2, 2.0)/2) -
            1/(a_12*new_u_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(A_k1, 2.0)/2) ;    
            
            score(0) += (1/V_x2)*
              (norm_A_k2*(a_12_alpha/a_12_sq + A_k2*A_k2_alpha/a_12- A_k2_alpha)/x_2_sq -
              norm_A_k1*(a_12_alpha/a_12 + A_k1*A_k1_alpha)/(a_12*new_u_1*x_2) ) - 
              A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
            
            score(1) += (1/V_x2)*
              (norm_A_k2*(a_12_phi/a_12_sq + A_k2*A_k2_phi/a_12- A_k2_phi)/x_2_sq -
              norm_A_k1*(a_12_phi/a_12 + A_k1*A_k1_phi)/(a_12*new_u_1*x_2) ) - 
              A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(new_u_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi);
            
            score.rows(2,p_1+1) += 
              -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
              -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_1_1.row(i).t()
              -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_1_2.row(i).t();
              
              score.rows(p_1+2, p_1+p_2+1) += 
              -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
              -(u_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * (V_x1_x2/V_x2 - V_x1) * z_2_1.row(i).t()
              -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * (V_x2_x2/V_x2 - V_x2) * z_2_2.row(i).t();
              
              score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
              1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                +(-log(cent_1) / pow(xi_1, 2.0) + (u_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * (V_x1_x2/V_x2 - V_x1) * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * (V_x2_x2/V_x2 - V_x2) * z_3_2.row(i).t();
        }
        if(y_1 > u_1 & y_2 > u_2){
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
              
              score(0) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
              
              score(1) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
              
              score.rows(2,p_1+1) += 
                -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
                -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
                -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_1_1.row(i).t()
                -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_1_2.row(i).t()
                ;
              
              score.rows(p_1+2, p_1+p_2+1) += 
                -( 1+ (y_1 - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
                -( 1+ (y_2 - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
                -(y_1-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
                -(y_2-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
                ;
              
              score.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
                1/xi_1*( (y_1 - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
                + 1/xi_2*( (y_2 - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
                +(-log(cent_1) / pow(xi_1, 2.0) + (y_1-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
                +(-log(cent_2) / pow(xi_2, 2.0) + (y_2-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
                ;
        }
      }
    }
  }
  
  score(0) = score(0) * alpha/(1+exp(par(0)));
  score(1) = score(1) * phi;
  return -score;
}

// [[Rcpp::export]]
arma::vec score_marg_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                           const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, 
                           const arma::mat& locs){
  int S = locs.n_rows;
  int n = y.n_rows;
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const arma::vec& beta_1 = par.subvec(0, p_1-1);
  const arma::vec& beta_2 = par.subvec(p_1, p_1+p_2-1);
  const arma::vec& beta_3 = par.subvec(p_1+p_2, p_1+p_2+p_3-1);
  
  arma::mat score = zeros<vec>(p_1 + p_2 + p_3);
  
  //There is an RcppArmadillo bug in subsetting a cube, we need to take the transpose to preserve dimensions
  //when there is only one slice
  
  //check here for updates: https://github.com/RcppCore/RcppArmadillo/issues/299
  bool bug_1 = FALSE;
  bool bug_2 = FALSE;
  bool bug_3 = FALSE;
  if(z_1.n_slices == 1) bug_1 = TRUE;
  if(z_2.n_slices == 1) bug_2 = TRUE;
  if(z_3.n_slices == 1) bug_3 = TRUE;
  
  arma::mat z_1_s, z_2_s, z_3_s;
  double mu, sigma, xi, y_s, u_s;
  for(int i=0; i<n; i++){
    for(int s=0; s<S; s++){
      z_1_s = z_1.row(s);
      z_2_s = z_2.row(s);
      z_3_s = z_3.row(s);
      
      if(bug_1) {
        z_1_s = z_1_s.t();
      }
      if(bug_2) {
        z_2_s = z_2_s.t();
      }
      if(bug_3) {
        z_3_s = z_3_s.t();
      }
      
      mu = as_scalar(z_1_s.row(i)*beta_1);
      sigma = exp(as_scalar(z_2_s.row(i)*beta_2));
      xi = as_scalar(z_3_s.row(i)*beta_3);
      
      y_s = y(i,s);
      u_s = thresholds(s);
      if(std::isnan(y_s) | std::isnan(mu) | std::isnan(sigma) | std::isnan(xi)) continue;
      
      if(y_s <= u_s){
        if(abs(xi) <= 1e-6){
          score.rows(0,p_1-1) += 1/sigma * exp(-(u_s-mu)/sigma) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1) += (u_s-mu)/pow(sigma, 2.0)*exp(-(u_s-mu)/sigma)* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1) += 0;
        } else {
          score.rows(0,p_1-1) += 1/sigma * pow(1+xi*(u_s-mu)/sigma, -1-1/xi) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1) += (u_s-mu)/pow(sigma, 2.0)*pow(1+xi*(u_s-mu)/sigma, -1-1/xi)* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1) += 
            pow(1+xi*(u_s-mu)/sigma, -1/xi) * ( log(xi*(u_s-mu)/sigma+1)/pow(xi, 2.0) - (u_s-mu)/(sigma*xi*(xi*(u_s-mu)/sigma+1)) ) * z_3_s.row(i).t();
        }
        
      }
      if(y_s > u_s){
        
        if(abs(xi) <= 1e-6){
          score.rows(0,p_1-1) += ((1+xi)/exp(-(y_s-mu)/sigma)-1)*(1/sigma * exp(-(y_s-mu)/sigma) ) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1) += ((1+xi)/exp(-(y_s-mu)/sigma)-1)*((y_s-mu)/pow(sigma, 2.0)*exp(-(y_s-mu)/sigma))* sigma*z_2_s.row(i).t() ;
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1) += 0;
        } else{
          score.rows(0,p_1-1) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*(1/sigma * pow(1+xi*(y_s-mu)/sigma, -1-1/xi)) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*((y_s-mu)/pow(sigma, 2.0)*pow(1+xi*(y_s-mu)/sigma, -1-1/xi))* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*
            (pow(1+xi*(y_s-mu)/sigma, -1/xi) * ( log(xi*(y_s-mu)/sigma+1)/pow(xi, 2.0) - (y_s-mu)/(sigma*xi*(xi*(y_s-mu)/sigma+1)) )) * z_3_s.row(i).t();
        }
      }
    }
  }
  
  return -score;
}

// [[Rcpp::export]]
arma::mat scoremat_marg_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                            const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, 
                            const arma::mat& locs){
  int S = locs.n_rows;
  int n = y.n_rows;
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const arma::vec& beta_1 = par.subvec(0, p_1-1);
  const arma::vec& beta_2 = par.subvec(p_1, p_1+p_2-1);
  const arma::vec& beta_3 = par.subvec(p_1+p_2, p_1+p_2+p_3-1);
  
  arma::mat score = zeros<mat>(p_1 + p_2 + p_3, n);
  
  //There is an RcppArmadillo bug in subsetting a cube, we need to take the transpose to preserve dimensions
  //when there is only one slice
  
  //check here for updates: https://github.com/RcppCore/RcppArmadillo/issues/299
  bool bug_1 = FALSE;
  bool bug_2 = FALSE;
  bool bug_3 = FALSE;
  if(z_1.n_slices == 1) bug_1 = TRUE;
  if(z_2.n_slices == 1) bug_2 = TRUE;
  if(z_3.n_slices == 1) bug_3 = TRUE;
  
  arma::mat z_1_s, z_2_s, z_3_s;
  double mu, sigma, xi, y_s, u_s;
  for(int i=0; i<n; i++){
    for(int s=0; s<S; s++){
      z_1_s = z_1.row(s);
      z_2_s = z_2.row(s);
      z_3_s = z_3.row(s);
      
      if(bug_1) {
        z_1_s = z_1_s.t();
      }
      if(bug_2) {
        z_2_s = z_2_s.t();
      }
      if(bug_3) {
        z_3_s = z_3_s.t();
      }
      
      mu = as_scalar(z_1_s.row(i)*beta_1);
      sigma = exp(as_scalar(z_2_s.row(i)*beta_2));
      xi = as_scalar(z_3_s.row(i)*beta_3);
      
      y_s = y(i,s);
      u_s = thresholds(s);
      if(std::isnan(y_s) | std::isnan(mu) | std::isnan(sigma) | std::isnan(xi)) continue;
      
      if(y_s <= u_s){
        if(abs(xi) <= 1e-6){
          score.rows(0,p_1-1).col(i) += 1/sigma * exp(-(u_s-mu)/sigma) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(i) += (u_s-mu)/pow(sigma, 2.0)*exp(-(u_s-mu)/sigma)* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(i) += 0;
        } else {
          score.rows(0,p_1-1).col(i) += 1/sigma * pow(1+xi*(u_s-mu)/sigma, -1-1/xi) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(i) += (u_s-mu)/pow(sigma, 2.0)*pow(1+xi*(u_s-mu)/sigma, -1-1/xi)* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(i) += 
            pow(1+xi*(u_s-mu)/sigma, -1/xi) * ( log(xi*(u_s-mu)/sigma+1)/pow(xi, 2.0) - (u_s-mu)/(sigma*xi*(xi*(u_s-mu)/sigma+1)) ) * z_3_s.row(i).t();
        }
        
      }
      if(y_s > u_s){
        
        if(abs(xi) <= 1e-6){
          score.rows(0,p_1-1).col(i) += ((1+xi)/exp(-(y_s-mu)/sigma)-1)*(1/sigma * exp(-(y_s-mu)/sigma) ) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(i) += ((1+xi)/exp(-(y_s-mu)/sigma)-1)*((y_s-mu)/pow(sigma, 2.0)*exp(-(y_s-mu)/sigma))* sigma*z_2_s.row(i).t() ;
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(i) += 0;
        } else{
          score.rows(0,p_1-1).col(i) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*(1/sigma * pow(1+xi*(y_s-mu)/sigma, -1-1/xi)) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(i) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*((y_s-mu)/pow(sigma, 2.0)*pow(1+xi*(y_s-mu)/sigma, -1-1/xi))* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(i) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*
            (pow(1+xi*(y_s-mu)/sigma, -1/xi) * ( log(xi*(y_s-mu)/sigma+1)/pow(xi, 2.0) - (y_s-mu)/(sigma*xi*(xi*(y_s-mu)/sigma+1)) )) * z_3_s.row(i).t();
        }
      }
    }
  }
  
  return -score;
}

// [[Rcpp::export]]
arma::mat sensitivity_marg_thresh(const arma::vec& par, const arma::vec& thresholds, const arma::mat& y, 
                               const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, 
                               const arma::mat& locs){
  int S = locs.n_rows;
  int n = y.n_rows;
  int p_1 = z_1.n_slices;
  int p_2 = z_2.n_slices;
  int p_3 = z_3.n_slices;
  
  const arma::vec& beta_1 = par.subvec(0, p_1-1);
  const arma::vec& beta_2 = par.subvec(p_1, p_1+p_2-1);
  const arma::vec& beta_3 = par.subvec(p_1+p_2, p_1+p_2+p_3-1);
  
  arma::mat score = zeros<mat>(p_1 + p_2 + p_3, n*S);
  
  //There is an RcppArmadillo bug in subsetting a cube, we need to take the transpose to preserve dimensions
  //when there is only one slice
  
  //check here for updates: https://github.com/RcppCore/RcppArmadillo/issues/299
  bool bug_1 = FALSE;
  bool bug_2 = FALSE;
  bool bug_3 = FALSE;
  if(z_1.n_slices == 1) bug_1 = TRUE;
  if(z_2.n_slices == 1) bug_2 = TRUE;
  if(z_3.n_slices == 1) bug_3 = TRUE;
  
  arma::mat z_1_s, z_2_s, z_3_s;
  double mu, sigma, xi, y_s, u_s;
  int it = 0;
  for(int i=0; i<n; i++){
    for(int s=0; s<S; s++){
      z_1_s = z_1.row(s);
      z_2_s = z_2.row(s);
      z_3_s = z_3.row(s);
      
      if(bug_1) {
        z_1_s = z_1_s.t();
      }
      if(bug_2) {
        z_2_s = z_2_s.t();
      }
      if(bug_3) {
        z_3_s = z_3_s.t();
      }
      
      mu = as_scalar(z_1_s.row(i)*beta_1);
      sigma = exp(as_scalar(z_2_s.row(i)*beta_2));
      xi = as_scalar(z_3_s.row(i)*beta_3);
      
      y_s = y(i,s);
      u_s = thresholds(s);
      if(std::isnan(y_s) | std::isnan(mu) | std::isnan(sigma) | std::isnan(xi)) continue;
      
      if(y_s <= u_s){
        if(abs(xi) <= 1e-6){
          score.rows(0,p_1-1).col(it) += 1/sigma * exp(-(u_s-mu)/sigma) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(it) += (u_s-mu)/pow(sigma, 2.0)*exp(-(u_s-mu)/sigma)* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(it) += 0;
        } else {
          score.rows(0,p_1-1).col(it) += 1/sigma * pow(1+xi*(u_s-mu)/sigma, -1-1/xi) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(it) += (u_s-mu)/pow(sigma, 2.0)*pow(1+xi*(u_s-mu)/sigma, -1-1/xi)* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(it) += 
            pow(1+xi*(u_s-mu)/sigma, -1/xi) * ( log(xi*(u_s-mu)/sigma+1)/pow(xi, 2.0) - (u_s-mu)/(sigma*xi*(xi*(u_s-mu)/sigma+1)) ) * z_3_s.row(i).t();
        }
        
      }
      if(y_s > u_s){
        
        if(abs(xi) <= 1e-6){
          score.rows(0,p_1-1).col(it) += ((1+xi)/exp(-(y_s-mu)/sigma)-1)*(1/sigma * exp(-(y_s-mu)/sigma) ) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(it) += ((1+xi)/exp(-(y_s-mu)/sigma)-1)*((y_s-mu)/pow(sigma, 2.0)*exp(-(y_s-mu)/sigma))* sigma*z_2_s.row(i).t() ;
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(it) += 0;
        } else{
          score.rows(0,p_1-1).col(it) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*(1/sigma * pow(1+xi*(y_s-mu)/sigma, -1-1/xi)) * z_1_s.row(i).t();
          
          score.rows(p_1, p_1+p_2-1).col(it) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*((y_s-mu)/pow(sigma, 2.0)*pow(1+xi*(y_s-mu)/sigma, -1-1/xi))* sigma*z_2_s.row(i).t();
          
          score.rows(p_1+p_2, p_1+p_2+p_3-1).col(it) += ((1+xi)/pow(1+xi*(y_s - mu)/sigma, -1/xi)-1)*
            (pow(1+xi*(y_s-mu)/sigma, -1/xi) * ( log(xi*(y_s-mu)/sigma+1)/pow(xi, 2.0) - (y_s-mu)/(sigma*xi*(xi*(y_s-mu)/sigma+1)) )) * z_3_s.row(i).t();
        }
      }
    }
  }
  
  return score;
}

