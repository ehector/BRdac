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
arma::cube z_constructor(const arma::mat& covariates, const int& S, const int& N, const int& p){
  arma::cube z(S, N, p);
  for(int q=0; q<p; q++){
    z.slice(q) = reshape(covariates.col(q), S, N);
  }
  return z;
}

// [[Rcpp::export]]
double logPCL_BR(const double& alpha, const double& phi, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
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
  
  // return the SUM
  return SUM;
}

// [[Rcpp::export]]
List PCscore_BR(const double& alpha, const double& phi, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
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
  double x_1, x_2, x_1_sq, x_2_sq, A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V, 
  V_x1, V_x2, V_x1_x2, scale;
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
arma::vec PCscore_BR_sum(const double& alpha, const double& phi, const arma::mat& x, const arma::vec& loc_1, const arma::vec& loc_2)
{
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int n = x.n_rows;
  arma::vec score = zeros<vec>(2);
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
double logCL_BR(const double& alpha, const double& phi, const arma::mat& x, const arma::mat& locs){
  int S = locs.n_rows;
  double SUM = 0;
  for(int s_1=0; s_1<S-1; s_1++){
    for(int s_2=s_1+1; s_2<S; s_2++){
      SUM += logPCL_BR(alpha, phi, join_rows(x.col(s_1), x.col(s_2)), locs.row(s_1).t(), locs.row(s_2).t()); 
    }
  }
  return SUM;
}

// [[Rcpp::export]]
List Cscore_BR(const double& alpha, const double& phi, const arma::mat& x, const arma::mat& locs){
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int S = locs.n_rows;
  int n = x.n_rows;
  List output(n);
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
    
    output(i) = score_i ;
  }
  
  return output;
}

// [[Rcpp::export]]
arma::vec Cscore_BR_sum(const double& alpha, const double& phi, const arma::mat& x, const arma::mat& locs){
  static const double pi = 3.14159265359; 
  static const double sqrt_2_pi = pow(2*pi, 0.5);
  
  int S = locs.n_rows;
  int n = x.n_rows;
  arma::vec score = zeros<vec>(2);
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
    
    score += score_i ;
  }
  
  return score;
}

// [[Rcpp::export]]
double logPCL_all(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, 
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
  double x_1, x_1_sq, x_2, x_2_sq, V, V_x1, V_x2, V_x1_x2, mu_1, sigma_1, xi_1, mu_2, sigma_2, xi_2;
  
  for(int i=0; i<n; i++){
    mu_1 = as_scalar(z_1_1.row(i)*beta_1);
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    mu_2 = as_scalar(z_1_2.row(i)*beta_1);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);
    
    x_1 = pow( 1 + xi_1*(y(i,0) - mu_1)/sigma_1, 1/xi_1);
    x_1_sq = pow(x_1, 2.0) ;
    x_2 = pow( 1 + xi_2*(y(i,1) - mu_2)/sigma_2, 1/xi_2);
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
    
    J_1 = 1/sigma_1 * pow(1 + xi_1 * (y(i,0) - mu_1) / sigma_1, 1 / xi_1-1);
    J_2 = 1/sigma_2 * pow(1 + xi_2 * (y(i,1) - mu_2) / sigma_2, 1 / xi_2-1);
    
    SUM += log(f) + log(J_1) + log(J_2);
  }
  
  // return the SUM
  return SUM;
}

// [[Rcpp::export]]
List PCscore_all(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, 
             const arma::mat& y, const arma::mat& z_1_1, const arma::mat& z_2_1, const arma::mat& z_3_1 , 
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
  double a_12_sq = pow(a_12, 2.0);
  double a_12_cu = pow(a_12, 3.0);
  double a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
  double a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
  double f, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
         A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V,
         V_x1, V_x2, V_x1_x1,V_x2_x2, V_x1_x2, V_x1_x2_x2, V_x1_x1_x2, deriv_1, deriv_2, scale,
         mu_1, sigma_1, xi_1, mu_2, sigma_2, xi_2;
  arma::vec score_i = zeros<vec>(2 + p_1 + p_2 + p_3);
  
  for(int i=0; i<n; i++){
    mu_1 = as_scalar(z_1_1.row(i)*beta_1);
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    mu_2 = as_scalar(z_1_2.row(i)*beta_1);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);
    
    x_1 = pow( 1 + xi_1*(y(i,0) - mu_1)/sigma_1, 1/xi_1);
    x_1_sq = pow(x_1, 2.0);
    x_1_cu = pow(x_1, 3.0);
    x_2 = pow( 1 + xi_2*(y(i,1) - mu_2)/sigma_2, 1/xi_2);
    x_2_sq = pow(x_2, 2.0);
    x_2_cu = pow(x_2, 3.0);
    
    cent_1 = 1 + xi_1 * (y(i,0) - mu_1) / sigma_1;
    cent_2 = 1 + xi_2 * (y(i,1) - mu_2) / sigma_2;
    
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
        
    //V_x1_x1_x2_x2 = 
    //  norm_A_k1*( (3*A_k1- pow(A_k1, 3.0))/(pow(a_12, 4.0)*x_1_cu*x_2_sq) + 2*(pow(A_k1, 2.0)-1)/(a_12_cu*x_1_cu*x_2_sq) +
    //  A_k1/(a_12_sq*x_1_cu*x_2_sq) - 2/(a_12*x_1_cu*x_2_sq) ) +
    //  norm_A_k2*( (3*A_k2- pow(A_k2, 3.0))/(pow(a_12, 4.0)*x_1_sq*x_2_cu) + 2*(pow(A_k2, 2.0)-1)/(a_12_cu*x_1_sq*x_2_cu) +
    //  A_k2/(a_12_sq*x_1_sq*x_2_cu) - 2/(a_12*x_1_sq*x_2_cu)) ;
    
    f = exp(-V) * (V_x1*V_x2 - V_x1_x2);
    
    //deriv = exp (-V)*( pow(V_x1*V_x2, 2.0) + 2*pow(V_x1_x2, 2.0) - pow(V_x1, 2.0)*V_x2_x2 - pow(V_x2, 2.0)*V_x1_x1 
    //                     + 2*V_x1*V_x1_x2_x2 + 2*V_x2*V_x1_x1_x2 + V_x1_x1*V_x2_x2 - V_x1_x1_x2_x2 - 4*V_x1*V_x2*V_x1_x2);
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
      -( 1+ (y(i,0) - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
      -( 1+ (y(i,1) - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
      -(y(i,0)-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
      -(y(i,1)-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
      ;
    
    score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
      1/xi_1*( (y(i,0) - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() +
      1/xi_2*( (y(i,1) - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
      +(-log(cent_1) / pow(xi_1, 2.0) + (y(i,0)-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
      +(-log(cent_2) / pow(xi_2, 2.0) + (y(i,1)-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
      ;
    
    output(i) = score_i ;
  }

  return output;
}

// [[Rcpp::export]]
arma::vec PCscore_all_sum(const double& alpha, const double& phi, const arma::vec& beta_1, const arma::vec& beta_2, const arma::vec& beta_3, 
                const arma::mat& y, const arma::mat& z_1_1, const arma::mat& z_2_1, const arma::mat& z_3_1 , 
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
  double a_12_sq = pow(a_12, 2.0);
  double a_12_cu = pow(a_12, 3.0);
  double a_12_alpha = pow(distance / phi, alpha - alpha/2) * log(distance / phi) / pow(2, 0.5);
  double a_12_phi = -alpha * pow(distance / phi, alpha - alpha / 2) / (phi * pow(2, 0.5));
  double f, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
         A_k1, A_k2, A_k1_alpha, A_k2_alpha, A_k1_phi, A_k2_phi, norm_A_k1, norm_A_k2, V,
         V_x1, V_x2, V_x1_x1,V_x2_x2, V_x1_x2, V_x1_x2_x2, V_x1_x1_x2, deriv_1, deriv_2, scale,
         mu_1, sigma_1, xi_1, mu_2, sigma_2, xi_2;
  arma::vec score_i = zeros<vec>(2 + p_1 + p_2 + p_3);
  
  for(int i=0; i<n; i++){
    mu_1 = as_scalar(z_1_1.row(i)*beta_1);
    sigma_1 = exp(as_scalar(z_2_1.row(i)*beta_2));
    xi_1 = as_scalar(z_3_1.row(i)*beta_3);
    mu_2 = as_scalar(z_1_2.row(i)*beta_1);
    sigma_2 = exp(as_scalar(z_2_2.row(i)*beta_2));
    xi_2 = as_scalar(z_3_2.row(i)*beta_3);
    
    x_1 = pow( 1 + xi_1*(y(i,0) - mu_1)/sigma_1, 1/xi_1);
    x_1_sq = pow(x_1, 2.0);
    x_1_cu = pow(x_1, 3.0);
    x_2 = pow( 1 + xi_2*(y(i,1) - mu_2)/sigma_2, 1/xi_2);
    x_2_sq = pow(x_2, 2.0);
    x_2_cu = pow(x_2, 3.0);
    
    cent_1 = 1 + xi_1 * (y(i,0) - mu_1) / sigma_1;
    cent_2 = 1 + xi_2 * (y(i,1) - mu_2) / sigma_2;
    
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
        
    //V_x1_x1_x2_x2 = 
    //  norm_A_k1*( (3*A_k1- pow(A_k1, 3.0))/(pow(a_12, 4.0)*x_1_cu*x_2_sq) + 2*(pow(A_k1, 2.0)-1)/(a_12_cu*x_1_cu*x_2_sq) +
    //  A_k1/(a_12_sq*x_1_cu*x_2_sq) - 2/(a_12*x_1_cu*x_2_sq) ) +
    //  norm_A_k2*( (3*A_k2- pow(A_k2, 3.0))/(pow(a_12, 4.0)*x_1_sq*x_2_cu) + 2*(pow(A_k2, 2.0)-1)/(a_12_cu*x_1_sq*x_2_cu) +
    //  A_k2/(a_12_sq*x_1_sq*x_2_cu) - 2/(a_12*x_1_sq*x_2_cu)) ;
        
    f = exp(-V) * (V_x1*V_x2 - V_x1_x2);
        
    //deriv = exp (-V)*( pow(V_x1*V_x2, 2.0) + 2*pow(V_x1_x2, 2.0) - pow(V_x1, 2.0)*V_x2_x2 - pow(V_x2, 2.0)*V_x1_x1 
    //                     + 2*V_x1*V_x1_x2_x2 + 2*V_x2*V_x1_x1_x2 + V_x1_x1*V_x2_x2 - V_x1_x1_x2_x2 - 4*V_x1*V_x2*V_x1_x2);
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
           )
    ;
        
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
           )
      ;
        
    score_i.rows(2,p_1+1) = 
      -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
      -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
      -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_1_1.row(i).t()
      -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_1_2.row(i).t()
      ;
        
    score_i.rows(p_1+2, p_1+p_2+1) = 
      -( 1+ (y(i,0) - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
      -( 1+ (y(i,1) - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
      -(y(i,0)-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
      -(y(i,1)-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
      ;
        
    score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
      1/xi_1*( (y(i,0) - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
      +1/xi_2*( (y(i,1) - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
      +(-log(cent_1) / pow(xi_1, 2.0) + (y(i,0)-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
      +(-log(cent_2) / pow(xi_2, 2.0) + (y(i,1)-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
      ;
    
    score += score_i ;
  }
  
  return score;
}

// [[Rcpp::export]]
double logCL_all(const arma::vec& par, const arma::mat& y, const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, const arma::mat& locs){
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
      
      SUM += logPCL_all(alpha, phi, beta_1, beta_2, beta_3, join_rows(y.col(s_1), y.col(s_2)), 
                        z_1_1, z_2_1, z_3_1, z_1_2, z_2_2, z_3_2, 
                        locs.row(s_1).t(), locs.row(s_2).t());
    }
  }
  
  if(std::isnan(SUM)) SUM = R_NegInf;
  return -SUM;
}

// [[Rcpp::export]]
arma::mat Cscore_all(const arma::vec& par, const arma::mat& y, const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, const arma::mat& locs){
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
  arma::mat output = zeros<mat>(2+p_1+p_2+p_3, n*S*(S-1)/2);
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi,
         f, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
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
  
  int it=0;
  for(int i=0; i<n; i++){
    for(int s_1=0; s_1<S-1; s_1++){
      for(int s_2=s_1+1; s_2<S; s_2++){
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
        
        x_1 = pow( 1 + xi_1*(y(i,s_1) - mu_1)/sigma_1, 1/xi_1);
        x_1_sq = pow(x_1, 2.0);
        x_1_cu = pow(x_1, 3.0);
        x_2 = pow( 1 + xi_2*(y(i,s_2) - mu_2)/sigma_2, 1/xi_2);
        x_2_sq = pow(x_2, 2.0);
        x_2_cu = pow(x_2, 3.0);
        
        cent_1 = 1 + xi_1 * (y(i,s_1) - mu_1) / sigma_1;
        cent_2 = 1 + xi_2 * (y(i,s_2) - mu_2) / sigma_2;
        
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
              -( 1+ (y(i,s_1) - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
              -( 1+ (y(i,s_2) - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
              -(y(i,s_1)-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
              -(y(i,s_2)-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
              ;
            
        score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
              1/xi_1*( (y(i,s_1) - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
              + 1/xi_2*( (y(i,s_2) - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
              +(-log(cent_1) / pow(xi_1, 2.0) + (y(i,s_1)-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
              +(-log(cent_2) / pow(xi_2, 2.0) + (y(i,s_2)-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
              ;
        output.col(it) = score_i;
        it += 1;
      }
    }
  }
  
  return output;
}

// [[Rcpp::export]]
arma::vec Cscore_all_sum(const arma::vec& par, const arma::mat& y, const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, const arma::mat& locs){
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
  arma::vec score = zeros<vec>(2 + p_1 + p_2 + p_3);
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi,
  f, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
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
  
  for(int i=0; i<n; i++){
    score_i = zeros<vec>(2 + p_1 + p_2 + p_3);
    
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
        
        x_1 = pow( 1 + xi_1*(y(i,s_1) - mu_1)/sigma_1, 1/xi_1);
        x_1_sq = pow(x_1, 2.0);
        x_1_cu = pow(x_1, 3.0);
        x_2 = pow( 1 + xi_2*(y(i,s_2) - mu_2)/sigma_2, 1/xi_2);
        x_2_sq = pow(x_2, 2.0);
        x_2_cu = pow(x_2, 3.0);
        
        cent_1 = 1 + xi_1 * (y(i,s_1) - mu_1) / sigma_1;
        cent_2 = 1 + xi_2 * (y(i,s_2) - mu_2) / sigma_2;
        
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
            
     score_i(0) += -A_k1_alpha * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_alpha * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
            
     score_i(1) += -A_k1_phi * exp(-pow(A_k1, 2.0)/2) /(x_1 * sqrt_2_pi) -A_k2_phi * exp(-pow(A_k2, 2.0)/2) /(x_2 * sqrt_2_pi) +
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
          
     score_i.rows(2,p_1+1) += 
        -(1-xi_1)/sigma_1 * pow(cent_1, -1) * z_1_1.row(i).t() 
        -(1-xi_2)/sigma_2 * pow(cent_2, -1) * z_1_2.row(i).t()
        -1/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_1_1.row(i).t()
        -1/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_1_2.row(i).t()
        ;
            
     score_i.rows(p_1+2, p_1+p_2+1) += 
        -( 1+ (y(i,s_1) - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
        -( 1+ (y(i,s_2) - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
        -(y(i,s_1)-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
        -(y(i,s_2)-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
        ;
            
      score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) += 
        1/xi_1*( (y(i,s_1) - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
        + 1/xi_2*( (y(i,s_2) - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
        +(-log(cent_1) / pow(xi_1, 2.0) + (y(i,s_1)-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
        +(-log(cent_2) / pow(xi_2, 2.0) + (y(i,s_2)-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
        ;
      }
    }
    
    score += score_i ;
  }
  
  for(int l=0; l < 2+p_1+p_2+p_3; l++){
    if(std::isnan(score(l))) score(l) = R_NegInf; 
  }
  return -score;
}

// [[Rcpp::export]]
arma::mat Chessian_all_sum(const arma::vec& par, const arma::mat& y, const arma::cube& z_1, const arma::cube& z_2, const arma::cube& z_3, 
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
  arma::mat sensitivity = zeros<mat>(2 + p_1 + p_2 + p_3, n*S*(S-1)/2);
  
  double distance, a_12, a_12_sq, a_12_cu, a_12_alpha, a_12_phi,
  f, x_1, x_2, x_1_sq, x_2_sq, x_1_cu, x_2_cu, cent_1, cent_2,
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
        
        x_1 = pow( 1 + xi_1*(y(i,s_1) - mu_1)/sigma_1, 1/xi_1);
        x_1_sq = pow(x_1, 2.0);
        x_1_cu = pow(x_1, 3.0);
        x_2 = pow( 1 + xi_2*(y(i,s_2) - mu_2)/sigma_2, 1/xi_2);
        x_2_sq = pow(x_2, 2.0);
        x_2_cu = pow(x_2, 3.0);
        
        cent_1 = 1 + xi_1 * (y(i,s_1) - mu_1) / sigma_1;
        cent_2 = 1 + xi_2 * (y(i,s_2) - mu_2) / sigma_2;
        
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
              -( 1+ (y(i,s_1) - mu_1)*(1-xi_1)*pow(cent_1, -1)/sigma_1 ) * z_2_1.row(i).t()
              -( 1+ (y(i,s_2) - mu_2)*(1-xi_2)*pow(cent_2, -1)/sigma_2 ) * z_2_2.row(i).t()
              -(y(i,s_1)-mu_1)/sigma_1 * pow(cent_1, 1/xi_1 - 1) * deriv_1 * z_2_1.row(i).t()
              -(y(i,s_2)-mu_2)/sigma_2 * pow(cent_2, 1/xi_2 - 1) * deriv_2 * z_2_2.row(i).t()
              ;
            
        score_i.rows(p_1+p_2+2, p_1+p_2+p_3+1) = 
              1/xi_1*( (y(i,s_1) - mu_1)*(1-xi_1)/sigma_1 /( cent_1 ) - 1/xi_1 * log(cent_1 ) ) * z_3_1.row(i).t() 
              + 1/xi_2*( (y(i,s_2) - mu_2)*(1-xi_2)/sigma_2 /( cent_2 ) - 1/xi_2 * log(cent_2 ) ) * z_3_2.row(i).t()
              +(-log(cent_1) / pow(xi_1, 2.0) + (y(i,s_1)-mu_1) / ( sigma_1*xi_1*cent_1 )) * pow(cent_1, 1/xi_1) * deriv_1 * z_3_1.row(i).t()
              +(-log(cent_2) / pow(xi_2, 2.0) + (y(i,s_2)-mu_2) / ( sigma_2*xi_2*cent_2 )) * pow(cent_2, 1/xi_2) * deriv_2 * z_3_2.row(i).t()
              ;
        
        sensitivity.col(it) = score_i;
        it += 1;
      }
    }
  }
  
  return sensitivity;
}

// [[Rcpp::export]]
List coeff_vcov_comp(const arma::mat& sensitivity, const arma::mat& V, const arma::vec& weighted_estimates){
  arma::mat V_inv = inv(V);
  arma::mat vcov = inv(sensitivity * V_inv * sensitivity.t());
  arma::vec coefficients = vcov * sensitivity * V_inv *weighted_estimates;
  
  return List::create(Named("coefficients") = coefficients,
                      Named("vcov") = vcov);
}

// [[Rcpp::export]]
arma::mat V_inv_eval(const arma::mat& psi){
  const int& m = psi.n_rows;
  const int& n = psi.n_cols;
  arma::mat V = zeros<mat>(m, m);;
  for(int i=0; i<n; i++){
    V += psi.col(i)*psi.col(i).t();
  }
  return inv(V);
}
