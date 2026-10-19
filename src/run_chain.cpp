// One MCMC chain, run entirely in C++.
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include "data_structure.h"
#include "kernels.h"
#include "threads.h"

namespace {

// log P(a < Z < b) for standard normal Z, accurate in both tails
double log_interval_prob(double a, double b) {
    double prob;
    if (a > 0) prob = R::pnorm(-a, 0, 1, 1, 0) - R::pnorm(-b, 0, 1, 1, 0);
    else prob = R::pnorm(b, 0, 1, 1, 0) - R::pnorm(a, 0, 1, 1, 0);
    return std::log(std::max(prob, 1e-300));
}

// Thresholds between adjacent levels: midpoint of the gap between the largest
// x of level v and the smallest x of level v + 1. Stored flat, column by
// column, with n_levels[j] - 1 entries per column.
std::vector<double> level_thresholds(const mmbpfa::Data& d, const Rcpp::NumericVector& x) {
    std::vector<double> th;
    for (int j = 0; j < d.p; j++) {
        for (int v = 0; v + 1 < d.n_levels[j]; v++) {
            double hi_v = R_NegInf, lo_next = R_PosInf;
            for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) hi_v = std::max(hi_v, x[c]);
            for (int c = d.level_begin(j, v + 1); c < d.level_begin(j, v + 2); c++) lo_next = std::min(lo_next, x[c]);
            th.push_back(0.5 * (hi_v + lo_next));
        }
    }
    return th;
}

double log_likelihood(const mmbpfa::Data& d, const std::vector<double>& th, const arma::mat& lambda,
                      const arma::mat& omega, const arma::vec& alpha) {
    const arma::mat omega_t = omega.t();
    std::vector<int> offset(d.p + 1, 0);     // start of column j's thresholds in th
    for (int j = 0; j < d.p; j++) offset[j + 1] = offset[j] + std::max(d.n_levels[j] - 1, 0);
    std::vector<double> ll_col(d.p, 0.0);
    #pragma omp parallel for schedule(dynamic, 8) num_threads(mmbpfa::num_threads())
    for (int j = 0; j < d.p; j++) {
        const int L = d.n_levels[j];
        double ll = 0;
        if (L >= 2) {
            const arma::rowvec lam_j = lambda.row(j);
            for (int v = 0; v < L; v++) {
                const double lo = v > 0 ? th[offset[j] + v - 1] : R_NegInf;
                const double hi = v + 1 < L ? th[offset[j] + v] : R_PosInf;
                for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) {
                    const double mu = arma::dot(lam_j, omega_t.col(d.row_idx[c])) - alpha[j];
                    ll += log_interval_prob(lo - mu, hi - mu);
                }
            }
        }
        ll_col[j] = ll;
    }
    double ll = 0;
    for (int j = 0; j < d.p; j++) ll += ll_col[j];   // fixed order: independent of threads
    return ll;
}

} // namespace

//' Run one chain
//'
//' @keywords internal
//'
//' @param data output of \code{prepare_data}.
//' @param x0 starting latent values (rank consistent, see \code{initial_x}).
//' @param state0 starting loadings state (see \code{.sample_loadings}).
//' @param omega0,alpha0 starting factor scores and intercepts.
//' @param prior resolved priors.
//' @param warmup,iter,thin numbers of warmup and sampling iterations, and thinning.
//' @param seed chain seed.
//' @param threads OpenMP threads used inside the chain.
//'
//' @return List of kept draws and summaries.
// [[Rcpp::export(name = ".run_chain")]]
Rcpp::List run_chain(const Rcpp::List& data, Rcpp::NumericVector x0, Rcpp::List state0,
                     arma::mat omega0, arma::vec alpha0, const Rcpp::List& prior,
                     int warmup, int iter, int thin, double seed, int threads) {
    mmbpfa::num_threads() = threads < 1 ? 1 : threads;
    const mmbpfa::Data d(data);
    Rcpp::NumericVector x = Rcpp::clone(x0);
    Rcpp::List state = state0;
    arma::mat omega = omega0;
    arma::vec alpha = alpha0;
    const int n_keep = iter / thin;

    Rcpp::List keep_lambda(n_keep), keep_z(n_keep), keep_c(n_keep), keep_omega(n_keep),
        keep_alpha(n_keep), keep_theta(n_keep);
    Rcpp::NumericVector keep_a(n_keep), keep_ll(n_keep);
    Rcpp::IntegerVector keep_K(n_keep);
    std::vector<double> th_sum;

    int kept = 0;
    for (int t = 1; t <= warmup + iter; t++) {
        arma::mat lambda = Rcpp::as<arma::mat>(state["lambda"]);
        x = sample_x_sparse(data, x, lambda, omega, alpha, seed, t);
        state = sample_loadings(data, x, omega, alpha, state, prior, seed, t);
        lambda = Rcpp::as<arma::mat>(state["lambda"]);
        omega = sample_omega_sparse(data, x, lambda, alpha, seed, t);
        alpha = sample_alpha_sparse(data, x, lambda, omega, seed, t);

        if (t > warmup && (t - warmup) % thin == 0 && kept < n_keep) {
            const arma::imat z = Rcpp::as<arma::imat>(state["z"]);
            const arma::ivec c = Rcpp::as<arma::ivec>(state["c"]);
            int K = 0;
            for (arma::uword h = 0; h < c.n_elem; h++)
                if (c[h] > static_cast<int>(h) + 1 && arma::accu(z.col(h)) >= 3) K++;
            const std::vector<double> th = level_thresholds(d, x);
            if (th_sum.empty()) th_sum.assign(th.size(), 0.0);
            for (size_t i = 0; i < th.size(); i++) th_sum[i] += th[i];

            keep_lambda[kept] = lambda;
            keep_z[kept] = z;
            keep_c[kept] = c;
            keep_omega[kept] = omega;
            keep_alpha[kept] = alpha;
            keep_theta[kept] = Rcpp::as<arma::vec>(state["theta"]);
            keep_a[kept] = Rcpp::as<double>(state["a_cusp"]);
            keep_K[kept] = K;
            keep_ll[kept] = log_likelihood(d, th, lambda, omega, alpha);
            kept++;
        }
        if (t % 50 == 0) Rcpp::checkUserInterrupt();
    }
    for (double& v : th_sum) v /= std::max(kept, 1);

    return Rcpp::List::create(
        Rcpp::Named("lambda") = keep_lambda, Rcpp::Named("z") = keep_z, Rcpp::Named("c") = keep_c,
        Rcpp::Named("omega") = keep_omega, Rcpp::Named("alpha") = keep_alpha,
        Rcpp::Named("theta") = keep_theta, Rcpp::Named("a_cusp") = keep_a,
        Rcpp::Named("K") = keep_K, Rcpp::Named("log_lik") = keep_ll,
        Rcpp::Named("thresholds") = th_sum);
}
