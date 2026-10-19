// Factor scores, conditioning on observed cells only.
#include <RcppArmadillo.h>
#include "data_structure.h"
#include "rng.h"
#include "threads.h"

namespace {

// In-place Cholesky of the upper triangle of a K x K column-major matrix: A = R'R.
// Returns false if A is not positive definite.
bool chol_upper(double* A, int K) {
    for (int j = 0; j < K; j++) {
        double* cj = A + j * K;
        for (int i = 0; i < j; i++) {
            const double* ci = A + i * K;
            double s = cj[i];
            for (int k = 0; k < i; k++) s -= ci[k] * cj[k];
            cj[i] = s / ci[i];
        }
        double s = cj[j];
        for (int k = 0; k < j; k++) s -= cj[k] * cj[k];
        if (!(s > 0)) return false;
        cj[j] = std::sqrt(s);
    }
    return true;
}

// Solve R'y = b in place (R upper triangular, column-major).
void solve_rt(const double* R, double* b, int K) {
    for (int i = 0; i < K; i++) {
        const double* ci = R + i * K;
        double s = b[i];
        for (int k = 0; k < i; k++) s -= ci[k] * b[k];
        b[i] = s / ci[i];
    }
}

// Solve R y = b in place.
void solve_r(const double* R, double* b, int K) {
    for (int i = K - 1; i >= 0; i--) {
        b[i] /= R[i * K + i];
        const double* ci = R + i * K;
        for (int k = 0; k < i; k++) b[k] -= ci[k] * b[i];
    }
}

// add (sign = +1) or remove (sign = -1) lambda_j lambda_j' in the upper triangle
void accumulate_outer(arma::mat& prec, const arma::mat& lambda_t, int j, double sign) {
    const int K = prec.n_rows;
    const double* l = lambda_t.colptr(j);
    for (int b = 0; b < K; b++) {
        const double lb = sign * l[b];
        if (lb == 0) continue;
        double* col = prec.colptr(b);
        for (int a = 0; a <= b; a++) col[a] += l[a] * lb;
    }
}

} // namespace


//' Sample factor scores from observed cells
//'
//' omega_i ~ N(S_i sum_{j in obs(i)} lambda_j x_ij, S_i) with
//' S_i = (sum_{j in obs(i)} lambda_j lambda_j' + I)^-1, drawn through the
//' Cholesky factor of S_i^-1. Fully observed rows share one precision matrix
//' and are drawn together; rows missing fewer than half their items start from
//' the shared matrix and subtract the missing items; rows observed on fewer than
//' K/2 items use the n_i x n_i sampler of Bhattacharya, Chakraborty & Mallick (2016);
//' other rows accumulate their observed items. Rows with no observations draw from N(0, I).
//'
//' @keywords internal
//'
//' @param data output of \code{prepare_data}.
//' @param x latent values, one per observed cell.
//' @param lambda P x K loadings.
//' @param seed,iter seed and iteration; row i uses its own omega-kernel stream.
//'
//' @return N x K matrix.
// [[Rcpp::export(name = ".sample_omega_sparse")]]
arma::mat sample_omega_sparse(const Rcpp::List& data, const Rcpp::NumericVector& x,
                              const arma::mat& lambda,
                              double seed, double iter) {
    const mmbpfa::Data d(data);
    const int K = lambda.n_cols;
    const arma::mat lambda_t = lambda.t();       // K x P: one contiguous column per item

    arma::mat G_full = lambda.t() * lambda;      // shared precision of a fully observed row
    G_full.diag() += 1.0;

    // fully observed rows are drawn together; give each a slot first
    std::vector<int> full_slot(d.n, -1), full_rows;
    for (int i = 0; i < d.n; i++)
        if (d.row_ptr[i + 1] - d.row_ptr[i] == d.p) { full_slot[i] = full_rows.size(); full_rows.push_back(i); }
    const int nf = full_rows.size();

    arma::mat omega_t(K, d.n);
    arma::mat rhs_full(K, nf), z_full(K, nf);
    bool failed = false;

    #pragma omp parallel for schedule(dynamic, 16) num_threads(mmbpfa::num_threads())
    for (int i = 0; i < d.n; i++) {
        const int begin = d.row_ptr[i], end = d.row_ptr[i + 1], n_i = end - begin;
        arma::vec rhs(K, arma::fill::zeros), z(K);
        for (int r = begin; r < end; r++) {
            const int c = d.cell_of_row[r], j = d.col_idx[c];
            const double w = x[c];
            const double* l = lambda_t.colptr(j);
            for (int k = 0; k < K; k++) rhs[k] += l[k] * w;
        }
        mmbpfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), mmbpfa::stream_id(mmbpfa::KERNEL_OMEGA, i));
        for (int k = 0; k < K; k++) z[k] = rng.norm();

        if (n_i > 0 && 2 * n_i < K) {
            // Bhattacharya, Chakraborty & Mallick (2016): with Phi = lambda_obs (n_i x K) and
            // a = x_obs, draw u ~ N(0, I_K), delta ~ N(0, I_n), solve
            // (Phi Phi' + I_n) w = a - (Phi u + delta), return u + Phi' w. Exact, O(n_i^2 K).
            std::vector<int> cols(n_i);
            std::vector<double> a(n_i), w(n_i), M(n_i * n_i);
            for (int r = 0; r < n_i; r++) {
                const int c = d.cell_of_row[begin + r];
                cols[r] = d.col_idx[c];
                a[r] = x[c];
            }
            for (int r = 0; r < n_i; r++) {
                const double* lr = lambda_t.colptr(cols[r]);
                double phi_u = 0;
                for (int k = 0; k < K; k++) phi_u += lr[k] * z[k];
                w[r] = a[r] - phi_u - rng.norm();
                for (int q = 0; q <= r; q++) {
                    const double* lq = lambda_t.colptr(cols[q]);
                    double dot = 0;
                    for (int k = 0; k < K; k++) dot += lr[k] * lq[k];
                    M[r * n_i + q] = dot + (q == r ? 1.0 : 0.0);   // upper triangle, column r
                }
            }
            if (!chol_upper(M.data(), n_i)) { failed = true; continue; }
            solve_rt(M.data(), w.data(), n_i);
            solve_r(M.data(), w.data(), n_i);
            for (int r = 0; r < n_i; r++) {
                const double* lr = lambda_t.colptr(cols[r]);
                for (int k = 0; k < K; k++) z[k] += lr[k] * w[r];
            }
            omega_t.col(i) = z;
            continue;
        }

        if (full_slot[i] >= 0) {
            rhs_full.col(full_slot[i]) = rhs;
            z_full.col(full_slot[i]) = z;
            continue;
        }
        arma::mat prec(K, K);
        if (2 * n_i > d.p) {
            prec = G_full;                      // subtract the missing items
            int r = begin;
            for (int j = 0; j < d.p; j++) {
                if (r < end && d.col_idx[d.cell_of_row[r]] == j) { r++; continue; }
                accumulate_outer(prec, lambda_t, j, -1.0);
            }
        } else {
            prec.eye();                         // add the observed items
            for (int r = begin; r < end; r++) accumulate_outer(prec, lambda_t, d.col_idx[d.cell_of_row[r]], 1.0);
        }
        // prec = R'R; draw = R^-1 (R'^-1 rhs + z), i.e. mean prec^-1 rhs plus N(0, prec^-1)
        if (!chol_upper(prec.memptr(), K)) { failed = true; continue; }
        solve_rt(prec.memptr(), rhs.memptr(), K);
        rhs += z;
        solve_r(prec.memptr(), rhs.memptr(), K);
        omega_t.col(i) = rhs;
    }
    if (failed) Rcpp::stop("factor-score precision is not positive definite");

    if (nf > 0) {
        const arma::mat R = arma::chol(arma::symmatu(G_full));
        const arma::mat draws = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), rhs_full) + z_full);
        for (int f = 0; f < nf; f++) omega_t.col(full_rows[f]) = draws.col(f);
    }
    return omega_t.t();
}
