// Scale moves that leave the posterior invariant (generalized Gibbs: Liu & Wu
// 1999; Liu & Sabatti 2000). A move picks a group action g_s, s > 0, and draws s
// from p(s) ∝ pi(g_s(state)) |J_s| / s (1/s is the Haar measure), so it moves the
// chain along directions the one-block-at-a-time Gibbs updates explore slowly.
#include <RcppArmadillo.h>
#include <cmath>
#include "data_structure.h"
#include "rng.h"
#include "threads.h"

namespace {

// Exact draw from GIG(p, chi, psi), density ∝ u^(p - 1) exp(-(psi u + chi / u) / 2).
// t = log u has the log-concave density g(t) = p t - (psi e^t + chi e^-t) / 2, sampled
// by ratio-of-uniforms around its mode; the acceptance rate of this bounding box is
// at least e / 4 for any log-concave density.
double draw_gig(double p, double chi, double psi, mosaicfa::Rng& rng) {
    if (chi <= 0 || (psi <= 0 && p >= 0)) return 1.0;     // improper: leave the state where it is
    const double y = psi > 0 ? (p + std::sqrt(p * p + psi * chi)) / psi : -chi / (2 * p);
    const double m = std::log(y);
    auto g = [&](double t) { return p * t - 0.5 * (psi * std::exp(t) + chi * std::exp(-t)); };
    auto dg = [&](double t) { return p - 0.5 * (psi * std::exp(t) - chi * std::exp(-t)); };
    const double gm = g(m);
    // v bounds: extrema of (t - m) sqrt(h(t)), where d/dt log|t - m| + g'(t) / 2 = 0
    auto root = [&](double sign) {
        double near = 0, far = 1;
        while (1.0 / far + 0.5 * sign * dg(m + sign * far) > 0 && far < 1e6) { near = far; far *= 2; }
        for (int it = 0; it < 100; it++) {
            const double mid = 0.5 * (near + far);
            if (1.0 / mid + 0.5 * sign * dg(m + sign * mid) > 0) near = mid; else far = mid;
        }
        const double d = 0.5 * (near + far);
        return d * std::exp(0.5 * (g(m + sign * d) - gm));
    };
    const double v_hi = root(1.0), v_lo = -root(-1.0);
    while (true) {
        const double u = rng.unif(), v = v_lo + (v_hi - v_lo) * rng.unif();
        const double t = m + v / u;
        if (2 * std::log(u) <= g(t) - gm) return std::exp(t);
    }
}

} // namespace

//' Item scale move
//'
//' For each item j, (x_j, lambda_j) -> s_j (x_j, lambda_j). The ranks of x_j are unchanged,
//' and s_j^2 ~ Gamma((n_j + m_j) / 2, A_j / 2) with A_j = sum_obs (x_ij - lambda_j' omega_i)^2
//' + sum_k z_jk lambda_jk^2 / theta_k: the parameter-expanded data augmentation of
//' Liu & Wu (1999) for probit models, per item.
//'
//' @keywords internal
//'
//' @return list with x, lambda and the scales s.
// [[Rcpp::export(name = ".item_scale_move")]]
Rcpp::List item_scale_move(const Rcpp::List& data, const Rcpp::NumericVector& x, const arma::mat& lambda,
                           const arma::imat& z, const arma::mat& omega,
                           const arma::vec& theta, double seed, double iter) {
    const mosaicfa::Data d(data);
    const arma::mat omega_t = omega.t();
    Rcpp::NumericVector x_out = Rcpp::clone(x);
    arma::mat lambda_out = lambda;
    arma::vec s(d.p, arma::fill::ones);
    #pragma omp parallel for schedule(dynamic, 8) num_threads(mosaicfa::num_threads())
    for (int j = 0; j < d.p; j++) {
        const int begin = d.col_ptr[j], end = d.col_ptr[j + 1];
        const arma::rowvec lam_j = lambda.row(j);
        double A = 0, m = 0;
        for (int c = begin; c < end; c++) {
            const double r = x[c] - arma::dot(lam_j, omega_t.col(d.row_idx[c]));
            A += r * r;
        }
        for (arma::uword k = 0; k < lambda.n_cols; k++)
            if (z(j, k)) { A += lam_j[k] * lam_j[k] / theta[k]; m++; }
        if (A <= 0 || end - begin + m == 0) continue;
        mosaicfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), mosaicfa::stream_id(mosaicfa::KERNEL_MOVES, j));
        const double sj = std::sqrt(rng.gamma(0.5 * (end - begin + m)) / (0.5 * A));
        s[j] = sj;
        for (int c = begin; c < end; c++) x_out[c] *= sj;
        lambda_out.row(j) *= sj;
    }
    return Rcpp::List::create(Rcpp::Named("x") = x_out, Rcpp::Named("lambda") = lambda_out,
                              Rcpp::Named("s") = Rcpp::NumericVector(s.begin(), s.end()));
}

//' Column scale move
//'
//' For each column k, (lambda_k, omega_k) -> (s_k lambda_k, omega_k / s_k), which leaves
//' Lambda Omega' unchanged; s_k^2 ~ GIG((m_k - N) / 2, sum_i omega_ik^2, sum_j lambda_jk^2 / theta_k).
//'
//' @keywords internal
//'
//' @return list with lambda, omega and the scales s.
// [[Rcpp::export(name = ".column_scale_move")]]
Rcpp::List column_scale_move(const arma::mat& lambda, const arma::imat& z, const arma::mat& omega,
                             const arma::vec& theta, double seed, double iter) {
    const int H = lambda.n_cols, N = omega.n_rows;
    arma::mat lambda_out = lambda, omega_out = omega;
    arma::vec s(H, arma::fill::ones);
    for (int k = 0; k < H; k++) {
        const double m = arma::accu(z.col(k));
        const double chi = arma::dot(omega.col(k), omega.col(k));
        const double psi = arma::dot(lambda.col(k), lambda.col(k)) / theta[k];
        mosaicfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter),
                        mosaicfa::stream_id(mosaicfa::KERNEL_MOVES, (uint64_t(1) << 40) + k));
        s[k] = std::sqrt(draw_gig(0.5 * (m - N), chi, psi, rng));
        lambda_out.col(k) *= s[k];
        omega_out.col(k) /= s[k];
    }
    return Rcpp::List::create(Rcpp::Named("lambda") = lambda_out, Rcpp::Named("omega") = omega_out,
                              Rcpp::Named("s") = Rcpp::NumericVector(s.begin(), s.end()));
}
