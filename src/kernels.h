// Declarations of the sampler kernels, so the chain driver can call them directly.
#ifndef MOSAICFA_KERNELS_H
#define MOSAICFA_KERNELS_H

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

#include <vector>
#include "data_structure.h"
namespace mosaicfa {
struct Discrete {
    std::vector<int> is_disc, g_off;
    int n_g = 0;
    Discrete(const Data& d, int max_levels);
};
void set_cut_sd(double sd);
std::vector<double> cell_means(const Data& d, const arma::mat& lambda, const arma::mat& omega);
std::vector<double> initial_cutpoints(const Data& d, const Discrete& ds, const Rcpp::NumericVector& x);
void draw_cutpoints(const Data& d, const Discrete& ds, const std::vector<double>& mu, std::vector<double>& g,
                    const arma::mat& lambda, double seed, double iter);
void draw_x_discrete(const Data& d, const Discrete& ds, const std::vector<double>& mu, const std::vector<double>& g,
                     Rcpp::NumericVector& x, double seed, double iter);
void collapsed_loadings(const Data& d, const Discrete& ds, const std::vector<double>& g, arma::mat& lambda,
                        arma::imat& z, const arma::mat& omega, const arma::vec& rho, const arma::vec& theta,
                        const std::vector<int>& cols, std::vector<double>& mu, double seed, double iter,
                        double& accepted, double& tried);
} // namespace mosaicfa

#endif
