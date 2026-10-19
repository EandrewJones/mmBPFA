// Loadings prior: cumulative shrinkage over columns (Legramanti, Durante &
// Dunson 2020) with item-level inclusion within columns (Schiavon, Canale &
// Dunson 2022). Column k is "on" (slab) when c_k > k, using 1-based k.
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include "data_structure.h"
#include "rng.h"

namespace {

struct Prior {
    double slab_shape, slab_rate, spike_var, a_cusp_shape, a_cusp_rate, rho_a, rho_b;
    explicit Prior(const Rcpp::List& p)
        : slab_shape(p["slab_shape"]), slab_rate(p["slab_rate"]), spike_var(p["spike_var"]),
          a_cusp_shape(p["a_cusp_shape"]), a_cusp_rate(p["a_cusp_rate"]),
          rho_a(p["rho_a"]), rho_b(p["rho_b"]) {}
};

double inv_logit(double lo) {
    if (lo >= 0) return 1.0 / (1.0 + std::exp(-lo));
    const double e = std::exp(lo);
    return e / (1.0 + e);
}

// log density of the active loadings of one column under the spike and the
// slab (theta integrated out of N(0, theta), theta ~ InvGamma(a, b): multivariate t)
void column_loglik(const arma::mat& lambda, const arma::imat& z, int h, const Prior& pr,
                   double& ll_spike, double& ll_slab) {
    double ss = 0;
    int m = 0;
    for (arma::uword j = 0; j < lambda.n_rows; j++)
        if (z(j, h)) { ss += lambda(j, h) * lambda(j, h); m++; }
    ll_spike = -0.5 * m * std::log(2 * M_PI * pr.spike_var) - 0.5 * ss / pr.spike_var;
    const double a = pr.slab_shape, b = pr.slab_rate;
    ll_slab = std::lgamma(a + 0.5 * m) - std::lgamma(a) - 0.5 * m * std::log(2 * M_PI * b) -
              (a + 0.5 * m) * std::log1p(ss / (2 * b));
}

arma::ivec draw_column_indicators(const arma::mat& lambda, const arma::imat& z,
                                  const arma::vec& v, const Prior& pr, mmbpfa::Rng& rng) {
    const int H = lambda.n_cols;
    arma::vec logw(H);
    double cum = 0;
    for (int l = 0; l < H; l++) {
        logw[l] = std::log(v[l]) + cum;
        cum += std::log1p(-v[l]);
    }
    arma::ivec c(H);
    arma::vec lp(H);
    for (int h = 0; h < H; h++) {
        double ll_spike, ll_slab;
        column_loglik(lambda, z, h, pr, ll_spike, ll_slab);
        for (int l = 0; l < H; l++) lp[l] = logw[l] + (l <= h ? ll_spike : ll_slab);  // l <= h: c = l + 1 <= h + 1
        const double mx = lp.max();
        double total = 0;
        for (int l = 0; l < H; l++) { lp[l] = std::exp(lp[l] - mx); total += lp[l]; }
        double u = rng.unif() * total;
        int pick = H - 1;
        for (int l = 0; l < H; l++) { u -= lp[l]; if (u <= 0) { pick = l; break; } }
        c[h] = pick + 1;
    }
    return c;
}

arma::vec draw_sticks(const arma::ivec& c, double a_cusp, mmbpfa::Rng& rng) {
    const int H = c.n_elem;
    arma::vec v(H);
    for (int l = 1; l < H; l++) {                      // 1-based l = 1..H-1
        double eq = 0, gt = 0;
        for (int h = 0; h < H; h++) { eq += (c[h] == l); gt += (c[h] > l); }
        v[l - 1] = rng.beta(1 + eq, a_cusp + gt);
    }
    v[H - 1] = 1;
    return v;
}

// log p(log a | c) with the sticks integrated out:
// Gamma(a; shape, rate) * prod_{l < H} B(1 + n_l, a + g_l) / B(1, a), times the Jacobian a
double log_a_posterior(double u, const std::vector<double>& n, const std::vector<double>& g, const Prior& pr) {
    const double a = std::exp(u);
    double lp = pr.a_cusp_shape * u - pr.a_cusp_rate * a;          // (shape - 1) log a + log a (Jacobian)
    for (size_t l = 0; l < n.size(); l++)
        lp += std::lgamma(1 + n[l]) + std::lgamma(a + g[l]) - std::lgamma(1 + n[l] + a + g[l]) + std::log(a);
    return lp;
}

// Slice sampling (Neal 2003) on log a, given the column indicators only. Drawing a
// given the sticks instead gets trapped at a = 0 once every column is off.
double draw_a_cusp(const arma::ivec& c, double a_current, const Prior& pr, mmbpfa::Rng& rng) {
    const int H = c.n_elem;
    std::vector<double> n(H - 1, 0), g(H - 1, 0);
    for (int l = 1; l < H; l++)
        for (int h = 0; h < H; h++) { n[l - 1] += (c[h] == l); g[l - 1] += (c[h] > l); }
    const double u0 = std::log(a_current);
    const double level = log_a_posterior(u0, n, g, pr) + std::log(rng.unif());
    const double width = 1.0;
    double lo = u0 - width * rng.unif(), hi = lo + width;
    for (int step = 0; step < 50 && log_a_posterior(lo, n, g, pr) > level; step++) lo -= width;
    for (int step = 0; step < 50 && log_a_posterior(hi, n, g, pr) > level; step++) hi += width;
    while (true) {
        const double u = lo + (hi - lo) * rng.unif();
        if (log_a_posterior(u, n, g, pr) > level) return std::exp(u);
        if (u < u0) lo = u; else hi = u;
    }
}

arma::vec draw_theta(const arma::mat& lambda, const arma::imat& z, const arma::ivec& c,
                     const Prior& pr, mmbpfa::Rng& rng) {
    const int H = lambda.n_cols;
    arma::vec theta(H);
    for (int h = 0; h < H; h++) {
        if (c[h] > h + 1) {
            double ss = 0; int m = 0;
            for (arma::uword j = 0; j < lambda.n_rows; j++)
                if (z(j, h)) { ss += lambda(j, h) * lambda(j, h); m++; }
            theta[h] = (pr.slab_rate + 0.5 * ss) / rng.gamma(pr.slab_shape + 0.5 * m);
        } else {
            theta[h] = pr.spike_var;
        }
    }
    return theta;
}

mmbpfa::Rng column_rng(double seed, double iter, int which) {
    return mmbpfa::Rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter),
                       mmbpfa::stream_id(mmbpfa::KERNEL_COLUMNS, which));
}

} // namespace

//' Sample column indicators of the cumulative shrinkage prior (test entry point)
//' @keywords internal
// [[Rcpp::export(name = ".sample_column_indicators")]]
arma::ivec sample_column_indicators(const arma::mat& lambda, const arma::imat& z, const arma::vec& v,
                                    const Rcpp::List& prior, double seed, double iter) {
    mmbpfa::Rng rng = column_rng(seed, iter, 0);
    return draw_column_indicators(lambda, z, v, Prior(prior), rng);
}

//' Sample stick-breaking weights (test entry point)
//' @keywords internal
// [[Rcpp::export(name = ".sample_sticks")]]
Rcpp::List sample_sticks(const arma::ivec& c, double a_cusp, const Rcpp::List& prior, double seed, double iter) {
    mmbpfa::Rng rng = column_rng(seed, iter, 1);
    return Rcpp::List::create(Rcpp::Named("v") = draw_sticks(c, a_cusp, rng));
}

//' Sample the cumulative shrinkage concentration given the indicators (test entry point)
//' @keywords internal
// [[Rcpp::export(name = ".sample_a_cusp")]]
double sample_a_cusp(const arma::ivec& c, double a_current, const Rcpp::List& prior, double seed, double iter) {
    mmbpfa::Rng rng = column_rng(seed, iter, 2);
    return draw_a_cusp(c, a_current, Prior(prior), rng);
}

//' Sample column slab variances (test entry point)
//' @keywords internal
// [[Rcpp::export(name = ".sample_theta")]]
arma::vec sample_theta(const arma::mat& lambda, const arma::imat& z, const arma::ivec& c,
                       const Rcpp::List& prior, double seed, double iter) {
    mmbpfa::Rng rng = column_rng(seed, iter, 3);
    return draw_theta(lambda, z, c, Prior(prior), rng);
}

//' One Gibbs update of the loadings block
//'
//' Sweeps the columns, updating each item's inclusion and loading with a
//' running residual over observed cells only, then the column inclusion
//' probabilities, column indicators, the concentration (sticks integrated out),
//' stick weights and slab variances.
//'
//' @keywords internal
//'
//' @param data output of \code{prepare_data}.
//' @param x latent values, one per observed cell.
//' @param omega N x H factor scores.
//' @param alpha item intercepts.
//' @param state list with lambda (P x H), z (P x H integer), rho, c (1-based), v, a_cusp, theta.
//' @param prior list with slab_shape, slab_rate, spike_var, a_cusp_shape, a_cusp_rate, rho_a, rho_b.
//' @param seed,iter seed and iteration.
//'
//' @return Updated state.
// [[Rcpp::export(name = ".sample_loadings")]]
Rcpp::List sample_loadings(const Rcpp::List& data, const Rcpp::NumericVector& x, const arma::mat& omega,
                           const arma::vec& alpha, const Rcpp::List& state, const Rcpp::List& prior,
                           double seed, double iter) {
    const mmbpfa::Data d(data);
    const Prior pr(prior);
    arma::mat lambda = Rcpp::as<arma::mat>(state["lambda"]);
    arma::imat z = Rcpp::as<arma::imat>(state["z"]);
    arma::vec rho = Rcpp::as<arma::vec>(state["rho"]);
    arma::vec v = Rcpp::as<arma::vec>(state["v"]);
    const arma::vec theta_in = Rcpp::as<arma::vec>(state["theta"]);
    const int H = lambda.n_cols;
    const arma::mat omega_t = omega.t();

    // residual e_c = x_c + alpha_j - lambda_j' omega_i for each observed cell
    const int nnz = x.size();
    std::vector<double> e(nnz);
    for (int c = 0; c < nnz; c++) {
        const int j = d.col_idx[c];
        e[c] = x[c] + alpha[j] - arma::dot(lambda.row(j), omega_t.col(d.row_idx[c]));
    }

    for (int k = 0; k < H; k++) {
        const double prior_prec = 1.0 / theta_in[k];
        const double logit_rho = std::log(rho[k]) - std::log1p(-rho[k]);
        for (int j = 0; j < d.p; j++) {
            const int begin = d.col_ptr[j], end = d.col_ptr[j + 1];
            const double old = lambda(j, k);
            double ww = 0, we = 0;
            for (int c = begin; c < end; c++) {
                const double w = omega(d.row_idx[c], k);
                e[c] += old * w;
                ww += w * w;
                we += w * e[c];
            }
            const double prec = ww + prior_prec, mu = we / prec;
            const double lo = logit_rho + 0.5 * (std::log(prior_prec) - std::log(prec)) + 0.5 * prec * mu * mu;
            mmbpfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter),
                            mmbpfa::stream_id(mmbpfa::KERNEL_LOADINGS, static_cast<uint64_t>(k) * d.p + j));
            const int on = rng.unif() < inv_logit(lo);
            const double lam = on ? mu + rng.norm() / std::sqrt(prec) : 0.0;
            z(j, k) = on;
            lambda(j, k) = lam;
            if (lam != 0)
                for (int c = begin; c < end; c++) e[c] -= lam * omega(d.row_idx[c], k);
        }
    }

    mmbpfa::Rng crng = column_rng(seed, iter, 4);
    for (int k = 0; k < H; k++) {
        const double m = arma::accu(z.col(k));
        rho[k] = crng.beta(pr.rho_a + m, pr.rho_b + d.p - m);
    }
    const arma::ivec c = draw_column_indicators(lambda, z, v, pr, crng);
    const double a_cusp = draw_a_cusp(c, Rcpp::as<double>(state["a_cusp"]), pr, crng);
    v = draw_sticks(c, a_cusp, crng);
    const arma::vec theta = draw_theta(lambda, z, c, pr, crng);

    return Rcpp::List::create(Rcpp::Named("lambda") = lambda, Rcpp::Named("z") = z,
                              Rcpp::Named("rho") = rho, Rcpp::Named("c") = c, Rcpp::Named("v") = v,
                              Rcpp::Named("a_cusp") = a_cusp, Rcpp::Named("theta") = theta);
}
