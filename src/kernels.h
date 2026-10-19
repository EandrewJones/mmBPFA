// Declarations of the sampler kernels, so the chain driver can call them directly.
#ifndef MMBPFA_KERNELS_H
#define MMBPFA_KERNELS_H

#include <RcppArmadillo.h>

Rcpp::NumericVector sample_x_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                                    const arma::mat& lambda, const arma::mat& omega, double seed, double iter);
arma::mat sample_omega_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                              const arma::mat& lambda, double seed, double iter);
Rcpp::List sample_loadings(const Rcpp::List& data, const Rcpp::NumericVector& x, const arma::mat& omega,
                           const Rcpp::List& state, const Rcpp::List& prior, double seed, double iter);
Rcpp::List item_scale_move(const Rcpp::List& data, const Rcpp::NumericVector& x, const arma::mat& lambda,
                           const arma::imat& z, const arma::mat& omega,
                           const arma::vec& theta, double seed, double iter);
Rcpp::List column_scale_move(const arma::mat& lambda, const arma::imat& z, const arma::mat& omega,
                             const arma::vec& theta, double seed, double iter);

#endif
