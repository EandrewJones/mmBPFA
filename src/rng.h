// Random number generation for the sampler kernels.
//
// Each parallel unit of work (a column, a row) gets its own stream, seeded
// from (seed, iteration, stream id), so draws do not depend on the number of
// threads. The generator is xoshiro256++ seeded through splitmix64.
#ifndef MMBPFA_RNG_H
#define MMBPFA_RNG_H

#include <cmath>
#include <cstdint>
#include <limits>

namespace mmbpfa {

// Stream ids carry the kernel in their high bits, so two kernels drawing for the
// same column or row in the same iteration never share random numbers.
enum Kernel : uint64_t { KERNEL_X = 0, KERNEL_OMEGA = 2, KERNEL_LOADINGS = 3, KERNEL_COLUMNS = 4, KERNEL_MOVES = 5 };
inline uint64_t stream_id(Kernel kernel, uint64_t unit) { return (static_cast<uint64_t>(kernel) << 48) ^ unit; }

inline uint64_t splitmix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

class Rng {
public:
    Rng(uint64_t seed, uint64_t iter, uint64_t stream) : has_spare_(false), spare_(0) {
        uint64_t k = splitmix64(splitmix64(splitmix64(seed) ^ iter) ^ (stream * 0xD1B54A32D192ED03ULL));
        for (int i = 0; i < 4; i++) { k = splitmix64(k); s_[i] = k; }
    }

    uint64_t next() {
        const uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
        const uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0]; s_[3] ^= s_[1]; s_[1] ^= s_[2]; s_[0] ^= s_[3];
        s_[2] ^= t; s_[3] = rotl(s_[3], 45);
        return result;
    }

    // uniform on the open interval (0, 1)
    double unif() { return ((next() >> 11) + 0.5) * 0x1.0p-53; }

    // standard normal, Marsaglia polar method
    double norm() {
        if (has_spare_) { has_spare_ = false; return spare_; }
        double u, v, s;
        do {
            u = 2.0 * unif() - 1.0; v = 2.0 * unif() - 1.0; s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double f = std::sqrt(-2.0 * std::log(s) / s);
        spare_ = v * f; has_spare_ = true;
        return u * f;
    }

    // gamma with unit scale, Marsaglia & Tsang (2000); shape < 1 by boosting
    double gamma(double shape) {
        if (shape < 1.0) return gamma(shape + 1.0) * std::pow(unif(), 1.0 / shape);
        const double d = shape - 1.0 / 3.0, c = 1.0 / std::sqrt(9.0 * d);
        while (true) {
            double z, v;
            do { z = norm(); v = 1.0 + c * z; } while (v <= 0.0);
            v = v * v * v;
            if (std::log(unif()) < 0.5 * z * z + d - d * v + d * std::log(v)) return d * v;
        }
    }

    double beta(double a, double b) {
        const double x = gamma(a), y = gamma(b);
        return x / (x + y);
    }

    // N(mu, 1) truncated to [lo, hi]; an empty interval returns its endpoint
    double tnorm(double mu, double lo, double hi) {
        return mu + std_tnorm(lo - mu, hi - mu);
    }

private:
    uint64_t s_[4];
    bool has_spare_;
    double spare_;

    static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

    // Robert (1995): standard normal truncated to [a, b]
    double std_tnorm(double a, double b) {
        if (!(a < b)) return std::isfinite(a) ? a : b;
        if (a >= 0.0) return upper(a, b);
        if (b <= 0.0) return -upper(-b, -a);
        // interval straddles zero
        if (b - a < 2.5066282746310002) {          // narrow: uniform proposal, density max at 0
            while (true) {
                const double z = a + (b - a) * unif();
                if (std::log(unif()) <= -0.5 * z * z) return z;
            }
        }
        while (true) {                              // wide: plain normal proposal
            const double z = norm();
            if (z >= a && z <= b) return z;
        }
    }

    // standard normal truncated to [a, b] with a >= 0
    double upper(double a, double b) {
        const double alpha = 0.5 * (a + std::sqrt(a * a + 4.0));
        const double uniform_cutoff = a + 2.0 / (a + std::sqrt(a * a + 4.0)) *
            std::exp(0.25 * (a * a - a * std::sqrt(a * a + 4.0)) + 0.5);
        if (b <= uniform_cutoff) {                  // narrow: uniform proposal, density max at a
            while (true) {
                const double z = a + (b - a) * unif();
                if (std::log(unif()) <= 0.5 * (a * a - z * z)) return z;
            }
        }
        while (true) {                              // exponential proposal
            const double z = a - std::log(unif()) / alpha;
            if (z > b) continue;
            if (std::log(unif()) <= -0.5 * (z - alpha) * (z - alpha)) return z;
        }
    }
};

} // namespace mmbpfa

#endif
