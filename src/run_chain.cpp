// One MCMC chain, run entirely in C++.
#include <RcppArmadillo.h>
#include <chrono>
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
std::vector<double> level_thresholds(const mosaicfa::Data& d, const Rcpp::NumericVector& x) {
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

double log_likelihood(const mosaicfa::Data& d, const std::vector<double>& th, const arma::mat& lambda,
                      const arma::mat& omega) {
    const arma::mat omega_t = omega.t();
    std::vector<int> offset(d.p + 1, 0);     // start of column j's thresholds in th
    for (int j = 0; j < d.p; j++) offset[j + 1] = offset[j] + std::max(d.n_levels[j] - 1, 0);
    std::vector<double> ll_col(d.p, 0.0);
    #pragma omp parallel for schedule(dynamic, 8) num_threads(mosaicfa::num_threads())
    for (int j = 0; j < d.p; j++) {
        const int L = d.n_levels[j];
        double ll = 0;
        if (L >= 2) {
            const arma::rowvec lam_j = lambda.row(j);
            for (int v = 0; v < L; v++) {
                const double lo = v > 0 ? th[offset[j] + v - 1] : R_NegInf;
                const double hi = v + 1 < L ? th[offset[j] + v] : R_PosInf;
                for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) {
                    const double mu = arma::dot(lam_j, omega_t.col(d.row_idx[c]));
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

//' Log-likelihood of the observed data given latent x (test entry point)
//'
//' Thresholds are the midpoints between adjacent levels of x, as in the sampler.
//'
//' @keywords internal
// [[Rcpp::export(name = ".log_likelihood_from_x")]]
double log_likelihood_from_x(const Rcpp::List& data, const Rcpp::NumericVector& x, const arma::mat& lambda,
                             const arma::mat& omega) {
    const mosaicfa::Data d(data);
    return log_likelihood(d, level_thresholds(d, x), lambda, omega);
}

//' Run one chain
//'
//' @keywords internal
//'
//' @param data output of \code{prepare_data}.
//' @param x0 starting latent values (rank consistent, see \code{initial_x}).
//' @param state0 starting loadings state (see \code{.sample_loadings}).
//' @param omega0 starting factor scores.
//' @param prior resolved priors.
//' @param warmup,iter,thin numbers of warmup and sampling iterations, and thinning.
//' @param seed chain seed.
//' @param threads OpenMP threads used inside the chain.
//' @param moves bit mask of extra scale moves: 1 item scale, 2 column scale.
//' @param max_levels SPIKE: items with 2..max_levels levels get explicit cutpoints (0: off).
//'
//' @return List of kept draws and summaries. \code{threshold_draws} holds, for every kept draw,
//'   the thresholds of items with at most 10 levels (rows in item order).
// [[Rcpp::export(name = ".run_chain")]]
Rcpp::List run_chain(const Rcpp::List& data, Rcpp::NumericVector x0, Rcpp::List state0,
                     arma::mat omega0, const Rcpp::List& prior,
                     int warmup, int iter, int thin, double seed, int threads, int moves, int max_levels) {
    mosaicfa::num_threads() = threads < 1 ? 1 : threads;
    const mosaicfa::Data d(data);
    Rcpp::NumericVector x = Rcpp::clone(x0);
    Rcpp::List state = state0;
    arma::mat omega = omega0;
    const int n_keep = iter / thin;

    Rcpp::List keep_lambda(n_keep), keep_z(n_keep), keep_c(n_keep), keep_omega(n_keep),
        keep_theta(n_keep);
    Rcpp::NumericVector keep_a(n_keep), keep_ll(n_keep);
    Rcpp::IntegerVector keep_K(n_keep);
    std::vector<double> th_sum;
    // thresholds of items with at most max_stored_levels levels are kept for every draw
    const int max_stored_levels = 10;
    std::vector<int> stored;                        // positions in the flat threshold vector
    for (int j = 0, pos = 0; j < d.p; j++) {
        const int n_th = std::max(d.n_levels[j] - 1, 0);
        if (d.n_levels[j] <= max_stored_levels)
            for (int v = 0; v < n_th; v++) stored.push_back(pos + v);
        pos += n_th;
    }
    Rcpp::NumericMatrix keep_th(stored.size(), n_keep);

    const mosaicfa::Discrete ds(d, max_levels);
    mosaicfa::set_cut_sd(prior.containsElementNamed("cut_sd") ? Rcpp::as<double>(prior["cut_sd"]) : 1e300);
    std::vector<double> g = mosaicfa::initial_cutpoints(d, ds, x);
    std::vector<double> g_sum(ds.n_g, 0.0);
    double coll_acc = 0, coll_try = 0;

    int kept = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int t = 1; t <= warmup + iter; t++) {
        arma::mat lambda = Rcpp::as<arma::mat>(state["lambda"]);
        if (ds.n_g > 0) {
            // x integrated out: cutpoints, then loadings of discrete items on active columns
            std::vector<double> mu = mosaicfa::cell_means(d, lambda, omega);
            mosaicfa::draw_cutpoints(d, ds, mu, g, lambda, seed, t);
            arma::imat z = Rcpp::as<arma::imat>(state["z"]);
            const arma::ivec c = Rcpp::as<arma::ivec>(state["c"]);
            std::vector<int> cols;
            bool first_off = true;
            for (int k = 0; k < (int) c.n_elem; k++) {
                if (c[k] > k + 1) cols.push_back(k);
                else if (first_off) { cols.push_back(k); first_off = false; }
            }
            mosaicfa::collapsed_loadings(d, ds, g, lambda, z, omega, Rcpp::as<arma::vec>(state["rho"]),
                                         Rcpp::as<arma::vec>(state["theta"]), cols, mu, seed, t, coll_acc, coll_try);
            state["lambda"] = lambda;
            state["z"] = z;
            if (Rcpp::as<double>(prior["cut_sd"]) <= 0) {      // Dirichlet cutpoint prior: tell the DA step
                Rcpp::NumericVector G(d.p, -1.0);
                Rcpp::IntegerVector m(d.p, 0);
                for (int j = 0; j < d.p; j++) {
                    if (!ds.is_disc[j]) continue;
                    m[j] = d.n_levels[j] - 1; G[j] = 0;
                    for (int v = 0; v < m[j]; v++) G[j] += g[ds.g_off[j] + v] * g[ds.g_off[j] + v];
                }
                state["cut_G"] = G; state["cut_m"] = m;
            }
            x = sample_x_sparse(data, x, lambda, omega, seed, t);
            mosaicfa::draw_x_discrete(d, ds, mu, g, x, seed, t);
        } else {
            x = sample_x_sparse(data, x, lambda, omega, seed, t);
        }
        state = sample_loadings(data, x, omega, state, prior, seed, t);
        lambda = Rcpp::as<arma::mat>(state["lambda"]);
        omega = sample_omega_sparse(data, x, lambda, seed, t);
        if (moves) {
            const arma::imat z = Rcpp::as<arma::imat>(state["z"]);
            const arma::vec theta = Rcpp::as<arma::vec>(state["theta"]);
            if (moves & 1) {
                Rcpp::List mv = item_scale_move(data, x, lambda, z, omega, theta, seed, t);
                x = mv["x"]; lambda = Rcpp::as<arma::mat>(mv["lambda"]);
            }
            if (moves & 2) {
                Rcpp::List mv = column_scale_move(lambda, z, omega, theta, seed, t);
                lambda = Rcpp::as<arma::mat>(mv["lambda"]); omega = Rcpp::as<arma::mat>(mv["omega"]);
            }
            state["lambda"] = lambda;
        }

        if (t > warmup && (t - warmup) % thin == 0 && kept < n_keep) {
            const arma::imat z = Rcpp::as<arma::imat>(state["z"]);
            const arma::ivec c = Rcpp::as<arma::ivec>(state["c"]);
            int K = 0;
            for (arma::uword h = 0; h < c.n_elem; h++)
                if (c[h] > static_cast<int>(h) + 1 && arma::accu(z.col(h)) >= 3) K++;
            const std::vector<double> th = level_thresholds(d, x);
            if (th_sum.empty()) th_sum.assign(th.size(), 0.0);
            for (size_t i = 0; i < th.size(); i++) th_sum[i] += th[i];
            for (int i = 0; i < ds.n_g; i++) g_sum[i] += g[i];
            for (size_t r = 0; r < stored.size(); r++) keep_th(r, kept) = th[stored[r]];

            keep_lambda[kept] = lambda;
            keep_z[kept] = z;
            keep_c[kept] = c;
            keep_omega[kept] = omega;
            keep_theta[kept] = Rcpp::as<arma::vec>(state["theta"]);
            keep_a[kept] = Rcpp::as<double>(state["a_cusp"]);
            keep_K[kept] = K;
            keep_ll[kept] = log_likelihood(d, th, lambda, omega);
            kept++;
        }
        if (t % 50 == 0) Rcpp::checkUserInterrupt();
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    for (double& v : th_sum) v /= std::max(kept, 1);
    for (double& v : g_sum) v /= std::max(kept, 1);

    return Rcpp::List::create(
        Rcpp::Named("lambda") = keep_lambda, Rcpp::Named("z") = keep_z, Rcpp::Named("c") = keep_c,
        Rcpp::Named("omega") = keep_omega,
        Rcpp::Named("theta") = keep_theta, Rcpp::Named("a_cusp") = keep_a,
        Rcpp::Named("K") = keep_K, Rcpp::Named("log_lik") = keep_ll,
        Rcpp::Named("thresholds") = th_sum, Rcpp::Named("threshold_draws") = keep_th,
        Rcpp::Named("cutpoints") = g_sum, Rcpp::Named("collapsed_accept") = coll_try > 0 ? coll_acc / coll_try : NA_REAL,
        Rcpp::Named("seconds_per_sweep") = seconds / (warmup + iter));
}
