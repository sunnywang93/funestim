// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// estimate_curve.cpp
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec epaKernelSmoothingCurve(
    const arma::vec & U, // Estimation points in U 
    const arma::vec & T, // Sampling points
    const arma::vec & Y, // Curves points
    const arma::vec & b, // Smoothing bandwidths
    const double & n_obs_min // Minimal number of obs for smoothing
){
  // Get parameters
  arma::uword M = T.n_elem; // Number of sampling points
  arma::uword L = U.n_elem; // Number of estimation points
  
  // Define output
  arma::vec Y_hat(L);
  Y_hat.fill(0);
  
  // Loop over the points to estimate
  for(arma::uword i=0; i<L; i++){
    double cpt = 0;
    double weight = 0;
    // Loop over the known points
    for (arma::uword k=0; k<M; k++){
      if (std::abs(T(k) - U(i)) <= b(i)){
        Y_hat(i) += (Y(k) * (1 - std::pow((T(k) - U(i))/b(i), 2)));
        weight += (1 - std::pow((T(k) - U(i))/b(i), 2));
        cpt += 1;
      }
    }
    if(cpt >= n_obs_min){
      Y_hat(i) /= weight;
    } else{
      Y_hat(i) = R_NaN;
    }
  }
  return Y_hat;
}

// [[Rcpp::export]]
arma::vec uniKernelSmoothingCurve(
    const arma::vec & U, // Estimation points in U 
    const arma::vec & T, // Sampling points
    const arma::vec & Y, // Curves points
    const arma::vec & b, // Smoothing bandwidths
    const double & n_obs_min // Minimal number of obs for smoothing
){
  // Get parameters
  arma::uword M = T.n_elem; // Number of sampling points
  arma::uword L = U.n_elem; // Number of estimation points
  
  // Define output
  arma::vec Y_hat(L);
  Y_hat.fill(0);
  
  // Loop over the points to estimate
  for(arma::uword i=0; i<L; i++){
    double cpt = 0;
    // Loop over the known points
    for (arma::uword k=0; k<M; k++){
      if (std::abs(T(k) - U(i)) <= b(i)){
        Y_hat(i) += Y(k);
        cpt += 1;
      }
    }
    if(cpt >= n_obs_min){
      Y_hat(i) /= cpt;
    } else{
      Y_hat(i) = R_NaN;
    }
  }
  return Y_hat;
}

// [[Rcpp::export]]
arma::vec biweightKernelSmoothingCurve(
    const arma::vec & U, // Estimation points in U 
    const arma::vec & T, // Sampling points
    const arma::vec & Y, // Curves points
    const arma::vec & b, // Smoothing bandwidths
    const double & n_obs_min // Minimal number of obs for smoothing
){
  // Get parameters
  arma::uword M = T.n_elem; // Number of sampling points
  arma::uword L = U.n_elem; // Number of estimation points
  
  // Define output
  arma::vec Y_hat(L);
  Y_hat.fill(0);
  
  // Loop over the points to estimate
  for(arma::uword i=0; i<L; i++){
    double cpt = 0;
    double weight = 0;
    // Loop over the known points
    for (arma::uword k=0; k<M; k++){
      if (std::abs(T(k) - U(i)) <= b(i)){
        Y_hat(i) += (Y(k) * std::pow(1 - std::pow((T(k) - U(i))/b(i), 2), 2));
        weight += std::pow(1 - std::pow((T(k) - U(i))/b(i), 2), 2);
        cpt += 1;
      }
    }
    if(cpt >= n_obs_min){
      Y_hat(i) /= weight;
    } else{
      Y_hat(i) = R_NaN;
    }
  }
  return Y_hat;
}
