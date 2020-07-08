// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

namespace ibp {
    NumericVector colSums(const arma::mat& x) {
        // get dimensions
        int p = x.n_cols;

        // sum across rows
        NumericVector sum_out(p);
        for (int j = 0; j < p; j++) {
            sum_out(j) = sum(x.col(j));
        }
        return sum_out;
    }

    arma::vec ibp_lbeta(
        NumericVector m_k,
        const int p,
        const double b
        ) {
        return lbeta(m_k, p - m_k + b);
    }

    double Hpb(
        const int p,
        const double b
        ) {
        // generate sequence from 1 to P
        arma::vec p_k = arma::linspace<arma::vec>(1, p, p);
        double Hpb = sum(b / (b + p_k - 1));
        return Hpb;
    }
}

// Sampling step for Indian Buffet Process Alpha Parameter
//
// @param ibp_b A scalar double, Indian Buffet Process beta parameter.
// @param zeros A binary, sparsity inducing matrix of dimenion P x K.
// @param e A scalar double, tunable hyperparameter for the shape parameter of the gamma draw. Default is 1.
// @param f A scalar double, tunable hyperparameter for the rate parameter of the gamma draw. Default is 1.
//
// @return ibp_a Scalar double, Indian Buffet Process alpha parameter.
//
// [[Rcpp::export(name = "sample_IBP_a")]]
double sample_IBP_a(
    const double ibp_b,
    const arma::mat& zeros,
    const double e = 1,
    const double f = 1
    ) {
        // get dimensions
        const int p = zeros.n_rows;
        const int K = zeros.n_cols;

        // sample a
        double ibp_a = R::rgamma(e + K, 1 / (f + ibp::Hpb(p, ibp_b)));
        return ibp_a;
}

//' Metropolis-Hastings sampling step for Indian Buffet Process Beta Parameter
//'
//' @param ibp_a A scalar double, Indian Buffet Process alpha parameter.
//' @param ibp_b A scalar double, prior Indian Buffet Process Beta parameter .
//' @param zeros A binary, sparsity inducing matrix of dimenion P x K.
//'
//' @return ibp_b Scalar double, Indian Buffet Process beta parameter.
//'
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(name = "sample_IBP_b")]]
double sample_IBP_b(
    const double ibp_a,
    const double ibp_b,
    const arma::mat& zeros
    ) {
        // get dimensions
        const int p = zeros.n_rows;
        const int K = zeros.n_cols;

        // number of features loaded on Ks
        NumericVector m_k = ibp::colSums(zeros);

        // sample new IBP propasal
        const double new_ibp_b = R::rgamma(2, 1);

        // Metropolis-Hastings update
        const arma::vec l_beta = ibp::ibp_lbeta(m_k, p, new_ibp_b);
        const arma::vec r_beta = ibp::ibp_lbeta(m_k, p, ibp_b);
        const double prop_a = sum(l_beta - r_beta);
        const double prop_b = K * (log(new_ibp_b) - log(ibp_b)) -
                              ibp_a * (ibp::Hpb(p, new_ibp_b) - ibp::Hpb(p, ibp_b));
        const double prop = exp(prop_a + prop_b);
        if (R::runif(0, 1) < prop) {
            return new_ibp_b;
        } else {
            return ibp_b;
        }
}