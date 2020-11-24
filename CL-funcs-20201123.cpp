// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <sstream>
using namespace Rcpp; 
using namespace arma;

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
double logCL_BR(const double& alpha, const double& phi, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
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
  double x_1, x_1_sq, x_2, x_2_sq, V, V_x1, V_x2, V_x1_x2;
  
  for(int i=0; i<n; i++){
    x_1 = x(i,0);
    x_1_sq = pow(x_1, 2.0) ;
    x_2 = x(i,1);
    x_2_sq = pow(x_2, 2.0);

    V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) + 
      (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false);
    
    V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x1_x2 = 
      1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (a_12/2 - (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (a_12/2 + (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) ;
    
    f = exp(-V) * (V_x1*V_x2 - V_x1_x2);
    
    SUM += log(f);
  }
  
  // return the SUM
  return SUM;
}

// [[Rcpp::export]]
double logCL_all(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, 
                 arma::mat& y, const arma::mat& x, const arma::mat& z_1_1, const arma::mat& z_2_1, const arma::mat& z_3_1 , 
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
  double x_1, x_1_sq, x_2, x_2_sq, V, V_x1, V_x2, V_x1_x2, sigma_1, xi_1, sigma_2, xi_2;
  
  for(int i=0; i<n; i++){
    x_1 = x(i,0);
    x_1_sq = pow(x_1, 2.0) ;
    x_2 = x(i,1);
    x_2_sq = pow(x_2, 2.0);
    
    V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) + 
      (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false);
    
    V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x1_x2 = 
      1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (a_12/2 - (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (a_12/2 + (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) ;
    
    f = exp(-V) * (V_x1*V_x2 - V_x1_x2);
    
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);
    J_1 =1/sigma_1 * pow(1 + xi_1 * (y(i,0) - as_scalar(z_1_1.row(i)*beta_1)) / sigma_1, 1 / xi_1);
    J_2 =1/sigma_2 * pow(1 + xi_2 * (y(i,1) - as_scalar(z_1_2.row(i)*beta_1)) / sigma_2, 1 / xi_2);
    
    SUM += log(f) + log(J_1) + log(J_2);
  }
  
  // return the SUM
  return SUM;;
}

// [[Rcpp::export]]
List Cscore_BR(const double& alpha, const double& phi, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int n = x.n_rows;
  List output(n);
  double distance = 0;
  for(int d=0; d<loc_1.size(); d++){
    distance += pow(loc_1(d) - loc_2(d), 2.0);
  }
  distance = pow(distance, 0.5);
  
  double a_12 = pow(2 * pow(distance / phi, alpha), 0.5);
  double a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
  double a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
  double x_1, x_2, x_1_sq, x_2_sq, A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V, V_x1, V_x2, V_x1_x2, scale;
  arma::vec score_i = zeros<vec>(2);
  
  for(int i=0; i<n; i++){
    x_1 = x(i,0);
    x_2 = x(i,1);
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
    
    V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) + 
      (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false);
    
    V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x1_x2 = 
      1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (a_12/2 - (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (a_12/2 + (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) ;
    
    scale = 1/(V_x1 * V_x2 - V_x1_x2);
    
    score_i(0) = -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
    
    score_i(1) = -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
  
    output(i) = score_i ;
  }
  
  return output;
}

// [[Rcpp::export]]
arma::vec Cscore_BR_sum(const double& alpha, const double& phi, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int n = x.n_rows;
  arma::vec score = zeros<vec>(2);
  List output(n);
  double distance = 0;
  for(int d=0; d<loc_1.size(); d++){
    distance += pow(loc_1(d) - loc_2(d), 2.0);
  }
  distance = pow(distance, 0.5);
  
  double a_12 = pow(2 * pow(distance / phi, alpha), 0.5);
  double a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
  double a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
  double x_1, x_2, x_1_sq, x_2_sq, A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V, V_x1, V_x2, V_x1_x2, scale;
  arma::vec score_i = zeros<vec>(2);
  
  for(int i=0; i<n; i++){
    x_1 = x(i,0);
    x_2 = x(i,1);
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
    
    V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) + 
      (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false);
    
    V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x1_x2 = 
      1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (a_12/2 - (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (a_12/2 + (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) ;
    
    scale = 1/(V_x1 * V_x2 - V_x1_x2);
    
    score_i(0) = -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
    
    score_i(1) = -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
    
    score += score_i ;
  }
  
  return score;
}

// [[Rcpp::export]]
List Cscore_all(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, 
             arma::mat& y, const arma::mat& x, const arma::mat& z_1_1, const arma::mat& z_2_1, const arma::mat& z_3_1 , 
             const arma::mat& z_1_2, const arma::mat& z_2_2, const arma::mat& z_3_2, 
             const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int n = y.n_rows;
  int p_1 = beta_1.size();
  int p_2 = beta_2.size();
  int p_3 = beta_3.size();
  List output(n);
  double distance = 0;
  for(int d=0; d<loc_1.size(); d++){
    distance += pow(loc_1(d) - loc_2(d), 2.0);
  }
  distance = pow(distance, 0.5);
  
  double a_12 = pow(2 * pow(distance / phi, alpha), 0.5);
  double a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
  double a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
  double x_1, x_2, x_1_sq, x_2_sq, A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V, V_x1, V_x2, V_x1_x2, scale, 
  mu_1, sigma_1, xi_1, mu_2, sigma_2, xi_2;
  arma::vec score_i = zeros<vec>(2 + p_1 + p_2 + p_3);
  
  for(int i=0; i<n; i++){
    x_1 = x(i,0);
    x_2 = x(i,1);
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
    
    V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) + 
      (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false);
    
    V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x1_x2 = 
      1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (a_12/2 - (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (a_12/2 + (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) ;
    
    scale = 1/(V_x1 * V_x2 - V_x1_x2);
    
    score_i(0) = -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
    
    score_i(1) = -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
    
    mu_1 = as_scalar(z_1_1.row(i)*beta_1);
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    mu_2 = as_scalar(z_1_2.row(i)*beta_1);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);

    score_i.rows(2,p_1+1) = 
      -1/sigma_1 * pow(1 + xi_1 * (y(i,0) - mu_1) / sigma_1, -1) * z_1_1.row(i) 
      -1/sigma_2 * pow(1 + xi_2 * (y(i,1) - mu_2) / sigma_2, -1) * z_1_2.row(i)
      ;
    
    score_i.rows(p_1+2, p_1+p_2+1) = 
      -( 1+ (y(i,0) - mu_1)*pow(1 + xi_1 * (y(i,0) - mu_1) / sigma_1, -1)/sigma_1 ) * z_2_1.row(i)
      -( 1+ (y(i,1) - mu_2)*pow(1 + xi_2 * (y(i,1) - mu_2) / sigma_2, -1)/sigma_2 ) * z_2_2.row(i)
      ;
    
    score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
      1/xi_1*( (y(i,0) - mu_1)/sigma_1 /( 1 + xi_1 * (y(i,0) - mu_1) / sigma_1 ) -
              1/xi_1 * log(1 + xi_1 * (y(i,0) - mu_1) / sigma_1 )
             ) * z_3_1.row(i) +
      1/xi_2*( (y(i,1) - mu_2)/sigma_2 /( 1 + xi_2 * (y(i,1) - mu_2) / sigma_2 ) -
              1/xi_2 * log(1 + xi_2 * (y(i,1) - mu_2) / sigma_2 )
             ) * z_3_2.row(i)
      ;
    
    output(i) = score_i ;
  }

  return output;
}

// [[Rcpp::export]]
arma::vec Cscore_all_sum(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, 
                arma::mat& y, const arma::mat& x, const arma::mat& z_1_1, const arma::mat& z_2_1, const arma::mat& z_3_1 , 
                const arma::mat& z_1_2, const arma::mat& z_2_2, const arma::mat& z_3_2, 
                const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int n = y.n_rows;
  int p_1 = beta_1.size();
  int p_2 = beta_2.size();
  int p_3 = beta_3.size();
  arma::vec score = zeros<vec>(2 + p_1 + p_2 + p_3);
  double distance = 0;
  for(int d=0; d<loc_1.size(); d++){
    distance += pow(loc_1(d) - loc_2(d), 2.0);
  }
  distance = pow(distance, 0.5);
  
  double a_12 = pow(2 * pow(distance / phi, alpha), 0.5);
  double a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
  double a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
  double x_1, x_2, x_1_sq, x_2_sq, A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V, V_x1, V_x2, V_x1_x2, scale, 
  mu_1, sigma_1, xi_1, mu_2, sigma_2, xi_2;
  arma::vec score_i = zeros<vec>(2 + p_1 + p_2 + p_3);
  
  for(int i=0; i<n; i++){
    x_1 = x(i,0);
    x_2 = x(i,1);
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
    
    V = (1/x_1) * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) + 
      (1/x_2) * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false);
    
    V_x1 = -1/x_1_sq * R::pnorm(a_12/2 - (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_1_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x2 = 1/(a_12*x_1*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/x_2_sq * R::pnorm(a_12/2 + (1/a_12)*log(x_1/x_2), 0, 1, true, false) -
      1/(a_12*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2);
    
    V_x1_x2 = 
      1/(pow(a_12, 2.0)*x_1_sq*x_2*pow(2*pi, 0.5)) * (a_12/2 - (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) +
      1/(pow(a_12, 2.0)*x_1*x_2_sq*pow(2*pi, 0.5)) * (a_12/2 + (1/a_12)*log(x_1/x_2)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1*x_2_sq*pow(2*pi, 0.5)) * exp(-pow(a_12/2 + (1/a_12)*log(x_1/x_2), 2.0)/2) -
      1/(a_12*x_1_sq*x_2*pow(2*pi, 0.5)) * exp(-pow(a_12/2 - (1/a_12)*log(x_1/x_2), 2.0)/2) ;
    
    scale = 1/(V_x1 * V_x2 - V_x1_x2);
    
    score_i(0) = -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
    
    score_i(1) = -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
    
    mu_1 = as_scalar(z_1_1.row(i)*beta_1);
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    mu_2 = as_scalar(z_1_2.row(i)*beta_1);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);
    
    score_i.rows(2,p_1+1) = 
      -1/sigma_1 * pow(1 + xi_1 * (y(i,0) - mu_1) / sigma_1, -1) * z_1_1.row(i) 
      -1/sigma_2 * pow(1 + xi_2 * (y(i,1) - mu_2) / sigma_2, -1) * z_1_2.row(i)
      ;
    
    score_i.rows(p_1+2, p_1+p_2+1) = 
      -( 1+ (y(i,0) - mu_1)*pow(1 + xi_1 * (y(i,0) - mu_1) / sigma_1, -1)/sigma_1 ) * z_2_1.row(i)
      -( 1+ (y(i,1) - mu_2)*pow(1 + xi_2 * (y(i,1) - mu_2) / sigma_2, -1)/sigma_2 ) * z_2_2.row(i)
      ;
    
    score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
      1/xi_1*( (y(i,0) - mu_1)/sigma_1 /( 1 + xi_1 * (y(i,0) - mu_1) / sigma_1 ) -
      1/xi_1 * log(1 + xi_1 * (y(i,0) - mu_1) / sigma_1 )
      ) * z_3_1.row(i) +
        1/xi_2*( (y(i,1) - mu_2)/sigma_2 /( 1 + xi_2 * (y(i,1) - mu_2) / sigma_2 ) -
        1/xi_2 * log(1 + xi_2 * (y(i,1) - mu_2) / sigma_2 )
        ) * z_3_2.row(i)
      ;
    
    score += score_i ;
  }
  
  return score;
}
