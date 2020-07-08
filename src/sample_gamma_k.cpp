// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

namespace local {
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


//' Sampling step for gamma_k (factor precisions)
//'
//' @keywords internal
//'
//' @param zeros A binary, sparsity inducing matrix of dimenion P x K.
//' @param lambda A standard multivariate normal matrix of factor loadings of dimension P x K.
//' @param c A tunable hyperparameter for the shape parameter of the gamma draw. Default is 0.
//' @param d A tunable hyperparameter for the inverse-scale (rate) parameter of the gamma draw. Default is 100
//'
//' @return gamma_k Numeric vector of factor precisions of length K.
//'
// [[Rcpp::export(name = "sample_gamma_k")]]
arma::vec sample_gamma_k(
    const arma::mat& zeros,
    const arma::mat& lambda,
    const double c = 0,
    const double d = 100
    ) {
        // get dimensions
        const int K = zeros.n_cols;

        // sum over active K
        const arma::vec m_k = local::colSums(zeros);

        // sample
        const arma::vec sum_l2 = local::colSums(pow(lambda, 2.0));
        arma::vec gamma_k_out(K);
        for (int k = 0; k < K; k++) {
            gamma_k_out(k) = R::rgamma(
                c + (m_k(k) / 2),
                1 / (sum_l2(k) + d));
        }

        return gamma_k_out;
}

//' Sampling step for d (hyperparameter for sample_gamma_k)
//'
//' @keywords internal
//'
//' @param zeros A binary, sparsity inducing matrix of dimenion P x K.
//' @param gamma_k A gamma-distributed vector of factor precisions of length K.
//' @param c A tunable hyperparameter for the shape parameter of the gamma draw. Default is 1.
//' @param c0 A tunable hyperparameter for the shape parameter of the gamma draw.
//' @param d0 A tunable hyperparameter for the inverse-scale (rate) parameter of the gamma draw.
//'
//' @return scalar double
//'
// [[Rcpp::export(name = "sample_d")]]
double sample_d(
    const arma::mat& zeros,
    const arma::vec& gamma_k,
    const double c,
    const double c0,
    const double d0
    ) {
        // get dimensions
        const int K = zeros.n_cols;

        // sample
        double d = R::rgamma(
            c0 + c * K,
            1 / (d0 + sum(gamma_k)));

        return d;
}