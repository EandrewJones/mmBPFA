// SPIKE: explicit cutpoints for discrete items, with updates that integrate the
// latent data out. Item j with 2..max_levels levels has ordered cutpoints g_j, and a
// cell at level v has probability Phi(g_v - mu) - Phi(g_{v-1} - mu), mu = lambda_j' omega_i
// (g_0 = -Inf, g_V = Inf); cutpoints have a flat prior on the ordered set.
//
// Partially collapsed Gibbs (van Dyk & Park 2008): the cutpoints and the loadings of
// discrete items on active columns are updated with x integrated out, and x is redrawn
// from its conditional before any step that conditions on it.
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "data_structure.h"
#include "kernels.h"
#include "rng.h"
#include "threads.h"

namespace {

const double LOG_SQRT_2PI = 0.918938533204672742;
double cut_sd = 1e300;   // SPIKE: prior sd of cutpoints (set from run_chain)

double log_phi(double u) { return std::isfinite(u) ? -0.5 * u * u - LOG_SQRT_2PI : R_NegInf; }

// log(Phi(b) - Phi(a)) for a < b, accurate in both tails
double log_interval(double a, double b) {
    if (a > 0) {           // both in the upper tail: use Phi(-a) - Phi(-b)
        const double la = R::pnorm(-a, 0, 1, 1, 1), lb = R::pnorm(-b, 0, 1, 1, 1);
        return la + std::log1p(-std::exp(lb - la));
    }
    const double lb = R::pnorm(b, 0, 1, 1, 1), la = R::pnorm(a, 0, 1, 1, 1);
    return lb + std::log1p(-std::exp(la - lb));
}

// log P, d log P / d mu and d^2 log P / d mu^2 for the interval [lo, hi] and mean mu
struct CellTerms { double lp, d1, d2; };
CellTerms cell_terms(double lo, double hi, double mu) {
    const double a = lo - mu, b = hi - mu;
    CellTerms t;
    t.lp = log_interval(a, b);
    const double pa = std::exp(log_phi(a) - t.lp), pb = std::exp(log_phi(b) - t.lp);
    t.d1 = pa - pb;
    const double apa = std::isfinite(a) ? a * pa : 0.0, bpb = std::isfinite(b) ? b * pb : 0.0;
    t.d2 = apa - bpb - t.d1 * t.d1;
    return t;
}

// Cutpoint prior when cut_sd <= 0: Dirichlet(1, ..., 1) on the item's marginal category
// probabilities under x ~ N(0, s^2), s^2 = 1 + |lambda_j|^2, i.e. density prod phi(g_v / s) / s on
// ordered g. Scale-coupled, so scaling lambda and g together is not rewarded.
double dir_logprior(double G, int m, double s2) { return -0.5 * G / s2 - 0.5 * m * std::log(s2); }

double cut_lo(const std::vector<double>& g, int off, int v) { return v == 0 ? R_NegInf : g[off + v - 1]; }
double cut_hi(const std::vector<double>& g, int off, int v, int L) { return v == L - 1 ? R_PosInf : g[off + v]; }

} // namespace

namespace mosaicfa {

// per item: whether it is discrete, and the offset of its cutpoints in g
Discrete::Discrete(const Data& d, int max_levels) : is_disc(d.p, 0), g_off(d.p, 0) {
    for (int j = 0; j < d.p; j++) {
        g_off[j] = n_g;
        if (d.n_levels[j] >= 2 && d.n_levels[j] <= max_levels) { is_disc[j] = 1; n_g += d.n_levels[j] - 1; }
    }
}

// mu_c = lambda_j' omega_i for every observed cell
std::vector<double> cell_means(const Data& d, const arma::mat& lambda, const arma::mat& omega) {
    const arma::mat omega_t = omega.t();
    std::vector<double> mu(d.row_idx.size());
    #pragma omp parallel for schedule(static) num_threads(num_threads())
    for (int c = 0; c < (int) mu.size(); c++) mu[c] = arma::dot(lambda.row(d.col_idx[c]), omega_t.col(d.row_idx[c]));
    return mu;
}

// cutpoints at the midpoints between adjacent levels of x (a rank-consistent start)
std::vector<double> initial_cutpoints(const Data& d, const Discrete& ds, const Rcpp::NumericVector& x) {
    std::vector<double> g(ds.n_g);
    for (int j = 0; j < d.p; j++) {
        if (!ds.is_disc[j]) continue;
        for (int v = 0; v + 1 < d.n_levels[j]; v++) {
            double hi_v = R_NegInf, lo_next = R_PosInf;
            for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) hi_v = std::max(hi_v, x[c]);
            for (int c = d.level_begin(j, v + 1); c < d.level_begin(j, v + 2); c++) lo_next = std::min(lo_next, x[c]);
            g[ds.g_off[j] + v] = 0.5 * (hi_v + lo_next);
        }
    }
    return g;
}

double item_loglik(const Data& d, int j, const std::vector<double>& g, int off, const std::vector<double>& mu) {
    const int L = d.n_levels[j];
    double ll = 0;
    for (int v = 0; v < L; v++) {
        const double lo = cut_lo(g, off, v), hi = cut_hi(g, off, v, L);
        for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) ll += log_interval(lo - mu[c], hi - mu[c]);
    }
    return ll;
}

// Random-walk Metropolis on (g_1, log gaps), flat prior on ordered g (Jacobian = product of gaps).
void set_cut_sd(double sd) { cut_sd = sd; }

void draw_cutpoints(const Data& d, const Discrete& ds, const std::vector<double>& mu, std::vector<double>& g,
                    const arma::mat& lambda, double seed, double iter) {
    #pragma omp parallel for schedule(dynamic, 8) num_threads(num_threads())
    for (int j = 0; j < d.p; j++) {
        if (!ds.is_disc[j]) continue;
        const int L = d.n_levels[j], m = L - 1, off = ds.g_off[j];
        int min_count = d.col_ptr[j + 1] - d.col_ptr[j];
        for (int v = 0; v < L; v++) min_count = std::min(min_count, d.level_begin(j, v + 1) - d.level_begin(j, v));
        const double step = 1.5 / std::sqrt(static_cast<double>(std::max(min_count, 1)) * m);
        Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), stream_id(KERNEL_CUTPOINTS, j));
        std::vector<double> cur(g.begin() + off, g.begin() + off + m), prop(m);
        // u = (g_1, log(g_2 - g_1), ...); proposal in u
        prop[0] = cur[0] + step * rng.norm();
        double log_jac_cur = 0, log_jac_prop = 0;
        for (int v = 1; v < m; v++) {
            const double lgap = std::log(cur[v] - cur[v - 1]);
            const double lgap_new = lgap + step * rng.norm();
            prop[v] = prop[v - 1] + std::exp(lgap_new);
            log_jac_cur += lgap; log_jac_prop += lgap_new;
        }
        const double ll_cur = item_loglik(d, j, cur, 0, mu);
        const double ll_prop = item_loglik(d, j, prop, 0, mu);
        double lp_cur = 0, lp_prop = 0;
        if (cut_sd > 0) {                        // SPIKE: N(0, cut_sd^2) prior on each cutpoint
            for (int v = 0; v < m; v++) { lp_cur -= 0.5 * cur[v] * cur[v] / (cut_sd * cut_sd); lp_prop -= 0.5 * prop[v] * prop[v] / (cut_sd * cut_sd); }
        } else {
            const double s2 = 1 + arma::dot(lambda.row(j), lambda.row(j));
            double Gc = 0, Gp = 0;
            for (int v = 0; v < m; v++) { Gc += cur[v] * cur[v]; Gp += prop[v] * prop[v]; }
            lp_cur = dir_logprior(Gc, m, s2); lp_prop = dir_logprior(Gp, m, s2);
        }
        if (std::log(rng.unif()) < ll_prop + log_jac_prop + lp_prop - ll_cur - log_jac_cur - lp_cur)
            for (int v = 0; v < m; v++) g[off + v] = prop[v];
    }
}

// x_c ~ N(mu_c, 1) truncated to the cell's interval, independently
void draw_x_discrete(const Data& d, const Discrete& ds, const std::vector<double>& mu, const std::vector<double>& g,
                     Rcpp::NumericVector& x, double seed, double iter) {
    #pragma omp parallel for schedule(dynamic, 8) num_threads(num_threads())
    for (int j = 0; j < d.p; j++) {
        if (!ds.is_disc[j]) continue;
        const int L = d.n_levels[j], off = ds.g_off[j];
        Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), stream_id(KERNEL_X, (uint64_t(1) << 40) + j));
        for (int v = 0; v < L; v++) {
            const double lo = cut_lo(g, off, v), hi = cut_hi(g, off, v, L);
            for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) x[c] = rng.tnorm(mu[c], lo, hi);
        }
    }
}

// For each column k in cols and each discrete item j: independence Metropolis on (z_jk, lambda_jk)
// with x integrated out. The proposal is the Laplace approximation of the exact conditional:
// q(z = 1) from Laplace evidences, lambda | z = 1 ~ N(mode, 1 / curvature).
// mu is kept current (mu_c includes column k's contribution).
void collapsed_loadings(const Data& d, const Discrete& ds, const std::vector<double>& g, arma::mat& lambda,
                        arma::imat& z, const arma::mat& omega, const arma::vec& rho, const arma::vec& theta,
                        const std::vector<int>& cols, std::vector<double>& mu, double seed, double iter,
                        double& accepted, double& tried) {
    double acc_tot = 0, try_tot = 0;
    for (int k : cols) {
        const double th = theta[k], log_rho = std::log(rho[k]), log_1m_rho = std::log1p(-rho[k]);
        #pragma omp parallel for schedule(dynamic, 4) num_threads(num_threads()) reduction(+:acc_tot, try_tot)
        for (int j = 0; j < d.p; j++) {
            if (!ds.is_disc[j]) continue;
            const int L = d.n_levels[j], off = ds.g_off[j], m = L - 1;
            const double lam_cur = lambda(j, k);
            const bool dir = cut_sd <= 0;
            double G = 0;
            for (int v = 0; v < m; v++) G += g[off + v] * g[off + v];
            const double c_rest = 1 + arma::dot(lambda.row(j), lambda.row(j)) - lam_cur * lam_cur;
            // f(l) = sum log P(mu^{-k} + l w) (+ the cutpoint prior's dependence on l) and its derivatives
            auto eval = [&](double l, double& f, double& f1, double& f2) {
                f = 0; f1 = 0; f2 = 0;
                for (int v = 0; v < L; v++) {
                    const double lo = cut_lo(g, off, v), hi = cut_hi(g, off, v, L);
                    for (int c = d.level_begin(j, v); c < d.level_begin(j, v + 1); c++) {
                        const double w = omega(d.row_idx[c], k);
                        const CellTerms t = cell_terms(lo, hi, mu[c] + (l - lam_cur) * w);
                        f += t.lp; f1 += t.d1 * w; f2 += t.d2 * w * w;
                    }
                }
                if (dir) {
                    const double s2 = c_rest + l * l, a1 = G / (s2 * s2) - m / s2;
                    f += dir_logprior(G, m, s2);
                    f1 += l * a1;
                    f2 += a1 + l * (-4 * l * G / (s2 * s2 * s2) + 2 * l * m / (s2 * s2));
                }
            };
            double f, f1, f2;
            // Newton for the mode of f(l) - l^2 / (2 theta), from the current value
            double l = lam_cur;
            for (int it = 0; it < 4; it++) {
                eval(l, f, f1, f2);
                const double g1 = f1 - l / th, g2 = f2 - 1.0 / th;
                if (!(g2 < 0)) break;
                l -= g1 / g2;
            }
            eval(l, f, f1, f2);
            const double mode = l, prec = 1.0 / th - f2;
            const double log_m1 = f - 0.5 * std::log(2 * M_PI * th) - mode * mode / (2 * th) +
                                  0.5 * std::log(2 * M_PI / prec);
            double f0, d0, e0;
            eval(0.0, f0, d0, e0);
            const double lq1 = log_rho + log_m1, lq0 = log_1m_rho + f0;
            const double mx = std::max(lq1, lq0);
            const double log_q1 = lq1 - (mx + std::log(std::exp(lq1 - mx) + std::exp(lq0 - mx)));
            const double log_q0 = lq0 - (mx + std::log(std::exp(lq1 - mx) + std::exp(lq0 - mx)));
            Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter),
                    stream_id(KERNEL_COLLAPSED, static_cast<uint64_t>(k) * d.p + j));
            const int z_new = rng.unif() < std::exp(log_q1);
            const double lam_new = z_new ? mode + rng.norm() / std::sqrt(prec) : 0.0;
            // log target and log proposal of a state (zz, ll)
            auto log_target = [&](int zz, double ll) {
                double fv, a, b;
                eval(zz ? ll : 0.0, fv, a, b);
                return fv + (zz ? log_rho - 0.5 * std::log(2 * M_PI * th) - ll * ll / (2 * th) : log_1m_rho);
            };
            auto log_prop = [&](int zz, double ll) {
                return zz ? log_q1 - 0.5 * std::log(2 * M_PI / prec) - 0.5 * prec * (ll - mode) * (ll - mode) : log_q0;
            };
            const int z_cur = z(j, k);
            const double log_a = log_target(z_new, lam_new) - log_target(z_cur, lam_cur) +
                                 log_prop(z_cur, lam_cur) - log_prop(z_new, lam_new);
            try_tot += 1;
            if (std::log(rng.unif()) < log_a) {
                acc_tot += 1;
                for (int c = d.col_ptr[j]; c < d.col_ptr[j + 1]; c++) mu[c] += (lam_new - lam_cur) * omega(d.row_idx[c], k);
                lambda(j, k) = lam_new;
                z(j, k) = z_new;
            }
        }
    }
    accepted += acc_tot; tried += try_tot;
}

} // namespace mosaicfa

//' SPIKE test entry: apply the collapsed (z, lambda) update for column k (0-based) n_iter times
//' with everything else fixed; returns the lambda draws of every item (P x n_iter).
//' @keywords internal
// [[Rcpp::export(name = ".collapsed_chain")]]
arma::mat collapsed_chain(const Rcpp::List& data, arma::mat lambda, arma::imat z, const arma::mat& omega,
                          const arma::vec& rho, const arma::vec& theta, Rcpp::NumericVector g_in, int k,
                          int max_levels, int n_iter, double seed, double cut_sd_in) {
    mosaicfa::set_cut_sd(cut_sd_in);
    const mosaicfa::Data d(data);
    const mosaicfa::Discrete ds(d, max_levels);
    std::vector<double> g(g_in.begin(), g_in.end());
    std::vector<double> mu = mosaicfa::cell_means(d, lambda, omega);
    arma::mat out(d.p, n_iter);
    double a = 0, t = 0;
    for (int it = 0; it < n_iter; it++) {
        mosaicfa::collapsed_loadings(d, ds, g, lambda, z, omega, rho, theta, std::vector<int>{k}, mu, seed, it + 1, a, t);
        out.col(it) = lambda.col(k);
    }
    return out;
}
