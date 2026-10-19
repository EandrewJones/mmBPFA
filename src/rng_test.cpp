// Test entry points for the RNG in rng.h (not part of the user API).
#include <Rcpp.h>
#include <string>
#include "rng.h"

// [[Rcpp::export(name = ".rng_draws")]]
Rcpp::NumericVector rng_draws(std::string kind, int n, double seed, Rcpp::NumericVector params) {
    mosaicfa::Rng rng(static_cast<uint64_t>(seed), 0, 0);
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; i++) {
        if (kind == "norm") out[i] = rng.norm();
        else if (kind == "unif") out[i] = rng.unif();
        else if (kind == "gamma") out[i] = rng.gamma(params[0]);
        else if (kind == "beta") out[i] = rng.beta(params[0], params[1]);
        else if (kind == "tnorm") out[i] = rng.tnorm(params[0], params[1], params[2]);
        else Rcpp::stop("unknown kind: " + kind);
    }
    return out;
}

// [[Rcpp::export(name = ".rng_stream")]]
Rcpp::NumericVector rng_stream(double seed, double iter, double stream, int n) {
    mosaicfa::Rng rng(static_cast<uint64_t>(seed), static_cast<uint64_t>(iter), static_cast<uint64_t>(stream));
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; i++) out[i] = rng.norm();
    return out;
}
