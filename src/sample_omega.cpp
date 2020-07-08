// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Sample N x P observations from a Standard
// Multivariate Normal given N observations, a
// vector  of P means, and a P x P cov matrix
namespace omega {
    arma::mat rmvnorm(
    int n,
    const arma::vec& mu,
    const arma::mat& Sigma
    ) {
        unsigned int p = Sigma.n_cols;

        // First draw N x P values from a N(0,1)
        Rcpp::NumericVector draw = Rcpp::rnorm(n*p);

        // Instantiate an Armadillo matrix with the
        // drawn values using advanced constructor
        // to reuse allocated memory
        arma::mat Z = arma::mat(draw.begin(), n, p, false, true);

        // Generate a sample from the Transformed
        // Multivariate Normal
        arma::mat Y = arma::repmat(mu, 1, n).t() + Z * arma::chol(Sigma);

        return Y;
    }
}

//' Sampling step for omega matrix (Factor Scores) of dimension N x K.
//'
//' @keywords internal
//'
//' @param lambda A standard multivariate normal matrix of factor loadings of dimension P x K.
//' @param x A multivariate Gaussian copula matrix of dimenion N x P.
//' @param alpha A numeric vector of item-intercepts of length P.
//'
//' @return omega A multivariate Gaussian matrix of factor scores of dimension N x K.
//'
// [[Rcpp::export(name = "sample_omega")]]
arma::mat sample_omega(
    arma::mat lambda,
    arma::mat x,
    arma::vec alpha
    ) {
        // get dimensions
        int n = x.n_rows, p = x.n_cols, K = lambda.n_cols;

        // initialize alpha matrix
        arma::mat alpha_mat(n, p);
        alpha_mat.each_row() = alpha.t();

        // sample omega
        arma::mat I(K, K, arma::fill::eye);
        arma::mat omega_var = (lambda.t() * lambda + I).i();
        arma::mat omega_mu = omega_var * lambda.t() * trans(x + alpha_mat);
        arma::mat omega_out(n, K);
        for (int i = 0; i < n; i++) {
            omega_out.row(i) = omega::rmvnorm(1, omega_mu.col(i), omega_var);
        }

        // normalize omega
        for (int k = 0; k < K; k++) {
            omega_out.col(k) = (omega_out.col(k) - arma::mean(omega_out.col(k))) / arma::stddev(omega_out.col(k));
        }

        return omega_out;
}