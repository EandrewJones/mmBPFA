// Extended rank likelihood step on observed cells only (Hoff 2007).
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include "data_structure.h"
#include "rng.h"
#include "threads.h"

//' Sample latent x for every observed cell
//'
//' Within each column, levels are drawn in increasing order: each draw is
//' bounded below by the maximum of the level just drawn and above by the
//' minimum of the next level's current values. Each column costs O(n_j).
//'
//' @keywords internal
//'
//' @param data output of \code{prepare_data}.
//' @param x current latent values, one per observed cell.
//' @param lambda P x K loadings (zero where inactive).
//' @param omega N x K factor scores.
//' @param seed,iter seed and iteration; column j uses its own x-kernel stream.
//'
//' @return New latent values, one per observed cell.
// [[Rcpp::export(name = ".sample_x_sparse")]]
Rcpp::NumericVector sample_x_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                                    const arma::mat& lambda, const arma::mat& omega, double seed, double iter) {
    const mosaicfa::Data d(data);
    const arma::mat omega_t = omega.t();     // K x N: one contiguous column per row
    Rcpp::NumericVector out = Rcpp::clone(x);

    #pragma omp parallel for schedule(dynamic, 8) num_threads(mosaicfa::num_threads())
    for (int j = 0; j < d.p; j++) {
        const int L = d.n_levels[j];
        if (L == 0) continue;
        mosaicfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), mosaicfa::stream_id(mosaicfa::KERNEL_X, j));
        const arma::rowvec lam_j = lambda.row(j);

        std::vector<double> level_min(L, R_PosInf);
        for (int v = 0; v < L; v++)
            for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++)
                level_min[v] = std::min(level_min[v], out[c]);

        double prev_max = R_NegInf;
        for (int v = 0; v < L; v++) {
            const double hi = (v + 1 < L) ? level_min[v + 1] : R_PosInf;
            double cur_max = R_NegInf;
            for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) {
                const double mu = arma::dot(lam_j, omega_t.col(d.row_idx[c]));
                out[c] = rng.tnorm(mu, prev_max, hi);
                cur_max = std::max(cur_max, out[c]);
            }
            prev_max = cur_max;
        }
    }
    return out;
}
