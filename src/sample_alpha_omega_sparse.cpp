// Item intercepts and factor scores, conditioning on observed cells only.
#include <RcppArmadillo.h>
#include "data_structure.h"
#include "rng.h"

//' Sample item intercepts from observed cells
//'
//' alpha_j ~ N(mean over observed i of (lambda_j' omega_i - x_ij), 1 / n_j) under
//' a flat prior. Items with no observations get alpha_j = 0.
//'
//' @keywords internal
//'
//' @param data output of \code{prepare_data}.
//' @param x latent values, one per observed cell.
//' @param lambda P x K loadings.
//' @param omega N x K factor scores.
//' @param seed,iter seed and iteration; item j uses its own alpha-kernel stream.
//'
//' @return Numeric vector of length P.
// [[Rcpp::export(name = ".sample_alpha_sparse")]]
arma::vec sample_alpha_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                              const arma::mat& lambda, const arma::mat& omega,
                              double seed, double iter) {
    const mmbpfa::Data d(data);
    const arma::mat omega_t = omega.t();
    arma::vec alpha(d.p, arma::fill::zeros);
    for (int j = 0; j < d.p; j++) {
        const int begin = d.col_ptr[j], end = d.col_ptr[j + 1], n_j = end - begin;
        if (n_j == 0) continue;
        const arma::rowvec lam_j = lambda.row(j);
        double total = 0;
        for (int c = begin; c < end; c++) total += arma::dot(lam_j, omega_t.col(d.row_idx[c])) - x[c];
        mmbpfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), mmbpfa::stream_id(mmbpfa::KERNEL_ALPHA, j));
        alpha[j] = total / n_j + rng.norm() / std::sqrt(static_cast<double>(n_j));
    }
    return alpha;
}

//' Sample factor scores from observed cells
//'
//' omega_i ~ N(S_i sum_{j in obs(i)} lambda_j (x_ij + alpha_j), S_i) with
//' S_i = (sum_{j in obs(i)} lambda_j lambda_j' + I)^-1, drawn through the
//' Cholesky factor of S_i^-1. Rows with no observations draw from N(0, I).
//'
//' @keywords internal
//'
//' @param data output of \code{prepare_data}.
//' @param x latent values, one per observed cell.
//' @param lambda P x K loadings.
//' @param alpha item intercepts, length P.
//' @param seed,iter seed and iteration; row i uses its own omega-kernel stream.
//'
//' @return N x K matrix.
// [[Rcpp::export(name = ".sample_omega_sparse")]]
arma::mat sample_omega_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                              const arma::mat& lambda, const arma::vec& alpha,
                              double seed, double iter) {
    const mmbpfa::Data d(data);
    const int K = lambda.n_cols;
    const arma::mat lambda_t = lambda.t();       // K x P: one contiguous column per item
    arma::mat omega_t(K, d.n);
    arma::mat prec(K, K);
    arma::vec rhs(K), z(K);
    for (int i = 0; i < d.n; i++) {
        prec.eye();
        rhs.zeros();
        for (int r = d.row_ptr[i]; r < d.row_ptr[i + 1]; r++) {
            const int c = d.cell_of_row[r], j = d.col_idx[c];
            const auto lam_j = lambda_t.col(j);
            prec += lam_j * lam_j.t();
            rhs += lam_j * (x[c] + alpha[j]);
        }
        mmbpfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), mmbpfa::stream_id(mmbpfa::KERNEL_OMEGA, i));
        for (int k = 0; k < K; k++) z[k] = rng.norm();
        // prec = R'R; mean = prec^-1 rhs; draw = mean + R^-1 z
        const arma::mat R = arma::chol(prec);
        const arma::vec mean = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), rhs));
        omega_t.col(i) = mean + arma::solve(arma::trimatu(R), z);
    }
    return omega_t.t();
}
