// Declarations of the sampler kernels, so the chain driver can call them directly.
#ifndef MMBPFA_KERNELS_H
#define MMBPFA_KERNELS_H

#include <RcppArmadillo.h>

Rcpp::NumericVector sample_x_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                                    const arma::mat& lambda, const arma::mat& omega,
                                    const arma::vec& alpha, double seed, double iter);
arma::vec sample_alpha_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                              const arma::mat& lambda, const arma::mat& omega,
                              double seed, double iter);
arma::mat sample_omega_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                              const arma::mat& lambda, const arma::vec& alpha,
                              double seed, double iter);
Rcpp::List sample_loadings(const Rcpp::List& data, const Rcpp::NumericVector& x, const arma::mat& omega,
                           const arma::vec& alpha, const Rcpp::List& state, const Rcpp::List& prior,
                           double seed, double iter);

#endif
