// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

namespace alpha {
    arma::vec colSums(const arma::mat& x) {
        // get dimensions
        int p = x.n_cols;

        // sum across rows
        arma::vec sum_out(p);
        for (int j = 0; j < p; j++) {
            sum_out(j) = sum(x.col(j));
        }
        return sum_out;
    }
}

//' Sampling step for alpha (item-intercepts)
//'
//' @keywords internal
//'
//' @param x A multivariate Gaussian copula matrix of dimenion N x P.
//' @param lambda A standard multivariate normal matrix of factor loadings of dimension P x K.
//' @param omega A standard multivariate normal matrix of factor scores of dimenion N x K.
//'
//' @return alpha Numeric vector of item-intercepts of length P.
//'
// [[Rcpp::export(name = "sample_alpha")]]
arma::vec sample_alpha(
    const arma::mat& x,
    const arma::mat& lambda,
    const arma::mat& omega
    ) {
        // get dimensions
        int n = x.n_rows;
        int p = x.n_cols;

        // sample
        arma::mat pred_vals = trans(lambda * omega.t());
        arma::mat resids = pred_vals - x;
        // flat prior and unit noise variance: alpha_j | -- ~ N(mean_resids_j, 1 / n)
        arma::vec mean_resids = alpha::colSums(resids) / n;
        const double sd_alpha = sqrt(1.0 / n);
        arma::vec alpha_out(p);
        for (int j = 0; j < p; j++) {
            alpha_out(j) = R::rnorm(mean_resids(j), sd_alpha);
        }
        return alpha_out;
}