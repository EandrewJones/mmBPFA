# Loadings Prior and Scalable Sampler: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the IBP with the winning loadings prior from the spike ablations. Then rebuild the sampler so that one sweep costs O(nnz·K + (N+P)·K³), runs entirely in C++, and uses OpenMP. That makes it scale across N, P and sparsity, ready for a later SVI derivation.

**Architecture:**
- **Phase A** (scratchpad spike, throwaway code): finish the ablations, run the Stage-2 stress tests, and pick the algorithm.
- **Phase B** (package, test-first):
  - a sparse data structure, so only observed cells carry latent x;
  - C++ kernels for x, α, ω, the loadings, and an RNG with per-column streams;
  - one C++ driver for the whole chain, behind a new R interface;
  - post-processing per the design;
  - benchmarks and optimization.

**Tech Stack:** R (≥ 4.0), Rcpp, RcppArmadillo (C++14, which R provides by default), OpenMP (optional), testthat edition 3.

**Spec:** `docs/superpowers/specs/2026-10-19-loadings-prior-design.md`

## Global Constraints

- **Statistical tests:** use `expect_z()` for means and `expect_rel()` for variances, both in `tests/testthat/helper-sim.R`. Never use `expect_equal(tolerance =)` on small statistical quantities.
- **Default priors:**
  - slab θ_k ~ InvGamma(1, ½) (the "uniform" preset);
  - ρ_k ~ Beta(1, 1);
  - a_cusp ~ Gamma(1, 1);
  - spike variance `spike_var` is fixed by Task 1, choosing between 0.05², 0.02² and the CUSP default 0.05.
- **Default truncation:** H = min(⌊(P − 1)/2⌋, 50). `K_start` is capped at H, with a message.
- **K for a draw** = the number of on columns (c_k > k) that have at least 3 active items.
- **The sampler output records every prior actually used.** `get_priors()` returns them.
- **The NAMESPACE file is maintained by hand.** C++ helpers that aren't public use dot-prefixed export names.
- **Determinism:** the same `seed` must give identical draws regardless of the number of OpenMP threads.
- **Git:** commit locally on branch `feat/cusp-prior`. Do not merge, tag or push without the author's go-ahead.

## Review Focus

1. **Items with one observed level or no observations** (constant or empty columns). The x-step and thresholds must skip them, not crash or produce NaN. Pinned in Task 4.
2. **Users with no observations** (empty rows). ω_i must be drawn from its N(0, I) prior. Pinned in Task 6.
3. **Continuous items where every value is distinct** (number of levels = n_j). The x-step must stay O(n_j) per column, not O(n_j × levels). Pinned in Task 4.
4. **`K_start > H`, H > P, or P < 7** (so the default H is less than 3). Cap K_start with a message, and error clearly when H < 1. Pinned in Task 8.
5. **Over 99% missing, with some items observed fewer than 3 times.** Those items still get loadings, and their α_j posterior variance must be 1/n_j, not 1/N. Pinned in Task 6.

---

## Phase A: choose the algorithm (spike, scratchpad)

### Task 1: Finish the spike-width ablations and choose the candidate

**Files:** scratchpad `spike/grid3.rds` (E3 with sd 0.02 and 0.01) and `spike/grid4.rds` (C with spike sd 0.05 and 0.02); `spike/summarize.R`.

- [ ] Run `Rscript summarize.R grid2.rds grid3.rds grid4.rds`.
- [ ] Choose by design §7, with the criteria treated as soft. Rank on:
  1. agreement across starts;
  2. median K exactly right;
  3. item-map F1;
  4. correlation error within 0.005 of baseline A (`grid.rds`).
  Tie-break on the share of iterations where K changes. Keep the top two for Task 2.
- [ ] Record the result in design §8 and §9, and fix `spike_var`.

### Task 2: Stage-2 stress tests on diverse data-generating processes

**Files:** create scratchpad `spike/stress_dgp.R` (the generators) and `spike/stress.R` (the runner). Reuse `fit()` from `harness.R`, extended with held-out masking.

**Generators.** Each returns `list(y, L, Z, truth_defined)`. N = 400, with P from 60 to 300.
- binary: balanced, and CTR 2%;
- ordinal: Likert 1–5, and J-shaped;
- nominal: 4 categories, one-hot encoded;
- counts: Poisson, NB(size = 1), zero-inflated (70% zeros), power-law;
- continuous: normal, lognormal, t₃, Beta, censored at the 80th percentile;
- rankings: Plackett–Luce top-10 over 50 items, encoded as rank position with NA when unranked;
- multi-view: rating + click + count on the same items, with one shared and one modality-specific factor;
- missingness: MCAR 90% and 99%, and MNAR exposure (selection probability rising with popularity and with x);
- misspecification: mixture-of-Gaussians ω, t-copula (ν = 3), a near-duplicate item pair, and 5% random-responder rows.

**Metrics.** For each generator × candidate × start ∈ {2, 10, 30}, compute the following, with 10% of observed cells held out:
- agreement across starts;
- K error where the truth is defined;
- item F1;
- copula-correlation error, or Kendall's τ error when the copula isn't Gaussian;
- held-out log predictive probability;
- AUC for binary cells, and NDCG@10 per row for ranking data.

The comparison set is the top two candidates from Task 1, baseline A, and normal-score PCA with parallel analysis for K.

- [ ] Implement the generators and check each one quickly: dimensions, missing rate, and level counts.
- [ ] Run the stress grid in the background on 20 cores. Summarize into design §8, under "Stage-2 results".
- [ ] **Decision:** pick one algorithm. Record failures that come from the structural limits (nominal, rankings, MNAR) as findings; they are not grounds for rejection.

---

## Phase B: implement in the package (branch `feat/cusp-prior`, based on `fix/sampler-correctness`)

### Task 3: Branch setup and baseline reference

**Files:**
- Create: `tests/testthat/fixtures/baseline_ibp_grid.rds` (spike grid A results);
- Create: `tests/testthat/fixtures/README.md` (the source commit and settings).

- [ ] Commit `fix/sampler-correctness` locally: the sampler fixes, the tests, `CLAUDE.md` and the docs.
- [ ] Create `feat/cusp-prior` from that branch.
- [ ] Copy the baseline results file, and record the commit SHA in the README.
- [ ] Commit.

### Task 4: Sparse data structure and the C++ x-step

**Files:**
- Create: `R/prepare_data.R`, `src/data_structure.h`, `src/sample_x_sparse.cpp`;
- Test: `tests/testthat/test-prepare-data.R`, `tests/testthat/test-sample-x.R` (rewritten for the new interface).

**Interfaces:**
- Produces `prepare_data(dat) -> list(n, p, row_idx, col_ptr, level, level_ptr, n_levels, row_ptr, col_idx, cell_of_row)`:
  - a column-compressed (CSC) layout of observed cells, sorted by level within each column (0-based level codes);
  - `level_ptr` holds the start of each level inside a column's segment;
  - a row-compressed (CSR) index for the row-wise updates.
- Produces `.sample_x_sparse(data, x, mu_fun inputs: lambda, omega, alpha, seed, iter) -> numeric(nnz)`. x is stored per observed cell in CSC order.

**Steps:**
- [ ] **Tests first:**
  - `prepare_data` round-trip: rebuilding the dense matrix from CSC gives `dat`, NA included;
  - constant and empty columns produce no levels and no error;
  - a column with all-distinct continuous values has `n_levels == n_j`.
- [ ] **Port the x-step tests to the new interface:**
  - the ordering invariant after every sweep;
  - pooled x is N(0,1) at stationarity under μ = 0 (binary, ordinal with missing, mixed counts);
  - a 5,000-row continuous column runs in under 50 ms per sweep.
- [ ] **Watch them fail, then implement.** Per column, walk the levels in order:
  - lower bound = the running max of the previous level's new draws;
  - upper bound = the min of the next level's current draws, precomputed per level and refreshed as each level is drawn.

  That makes each column O(n_j). Draws use the Task 5 truncated-normal sampler.
- [ ] Run the tests until they pass, then commit.

### Task 5: RNG and truncated-normal kernel

**Files:** create `src/rng.h`, `src/rng_test.cpp`; test `tests/testthat/test-rng.R`.

**Interfaces:**
- Produces a `Rng` class: xoshiro256++ seeded by splitmix64 from (seed, iteration, stream). Methods `unif()`, `norm()`, `gamma(shape)`, `beta(a, b)`, `tnorm(mu, lo, hi)`.
- `tnorm` uses inverse-CDF sampling in log space when the interval carries at least 1e-8 of the mass, and Robert (1995) exponential rejection for one-sided tails.
- Exported for tests only: `.rng_draws(kind, n, seed, params)`.

**Steps:**
- [ ] **Tests first:**
  - moments of `norm`, `gamma(2.5)`, `beta(0.5, 5)` and `tnorm(0, 1, 2)` against analytic values, using `expect_z` / `expect_rel`;
  - extreme tails: `tnorm(0, 8, Inf)` has mean 8.1170, and `tnorm(0, -Inf, -10)` returns finite values below −10;
  - identical output for identical (seed, iteration, stream).
- [ ] Run them to confirm they fail.
- [ ] Implement.
- [ ] Run them to confirm they pass, then commit.

### Task 6: Observed-only α and ω kernels

**Files:** create `src/sample_alpha_omega_sparse.cpp`; test `tests/testthat/test-sample-alpha.R` and `test-sample-omega.R` (rewritten).

**Interfaces:**
- Consumes the Task 4 `data` list and x in CSC order.
- Produces:
  - `.sample_alpha_sparse(data, x, lambda, omega, seed, iter) -> numeric(p)`, with α_j ~ N(mean over observed i of (λ_jᵀω_i − x_ij), 1/n_j);
  - `.sample_omega_sparse(data, x, lambda, alpha, seed, iter) -> matrix(n, H)`, with ω_i ~ N(Σ_i Σ_{j∈obs(i)} λ_j(x_ij + α_j), Σ_i) and Σ_i = (Σ_{j∈obs(i)} λ_jλ_jᵀ + I)⁻¹, drawn via a Cholesky factor per row.

**Steps:**
- [ ] **Tests first.** Analytic-posterior tests, adapted from the current files, with 40% of cells missing:
  - α's variance is `1/n_j` for each column (Review Focus 5);
  - an empty row gives ω_i with mean 0 and covariance I (Review Focus 2).
- [ ] Confirm the tests fail, implement, confirm they pass.
- [ ] Commit.

### Task 7: Loadings sweep for the chosen prior

**Files:** create `src/sample_loadings.cpp`; test `tests/testthat/test-sample-loadings.R`.

**Interfaces:** Produces `.sample_loadings(data, x, omega, alpha, state, prior, seed, iter) -> state`.
- `state` is a list: `lambda` (p×H), `z` (p×H, integer), `rho` (H), `c` (H, integer, 1-based), `v` (H), `a_cusp`, `theta` (H). For E3 it also has `eta`, `xi`, `m`.
- `prior` is a list: `slab_shape`, `slab_rate`, `spike_var`, `a_cusp_shape`, `a_cusp_rate`, `rho_a`, `rho_b`.

**Algorithm:**
- For each column k, keep a running residual over observed cells only. Removing column k adds back λ_jk·ω_ik for the observed cells of j.
- Then, for each item j:
  - prec = Σ_{i∈obs(j)} ω_ik² + 1/θ_k;
  - μ = Σ ω_ik e_ij / prec;
  - log-odds = logit ρ_k + ½(−log θ_k − log prec) + ½·prec·μ².
- Then ρ_k, the column categorical c_k, the v's, a_cusp and θ_k, all as in design §4–5.

**Steps:**
- [ ] **Tests first:**
  - λ given z = 1 against the analytic Gaussian, with K = 1, ρ fixed at 1 and the column on;
  - z frequency against the hand-computed odds;
  - the c_k categorical against hand-computed probabilities for a 3-column fixture with fixed loadings (spike likelihood Π N(0, spike_var) against the multivariate-t slab);
  - θ_k and a_cusp against their Gamma conditionals;
  - **prior invariance:** with no observed cells, the distribution over 20k sweeps of the number of on columns matches simulation from the prior (Kolmogorov–Smirnov p > 0.01 against 20k prior draws). Use a different seed per run.
- [ ] Confirm each test fails, implement, confirm it passes.
- [ ] Commit.

### Task 8: C++ chain driver and the new R interface

**Files:**
- Create: `src/run_chain.cpp`, `R/priors.R`;
- Rewrite: `R/mmBPFA_sampler.R`;
- Delete: `R/sampler_init.R`, `R/remove_inactive.R`, `src/sample_lambda_zeros.cpp`, `src/sample_gamma_k.cpp`, `src/sample_ibp.cpp`, `src/sample_alpha.cpp`, `src/sample_omega.cpp`, `R/sample_x.R`, `R/calculate_dmask.R`, and the IBP tests;
- Test: `tests/testthat/test-sampler-run.R`, `tests/testthat/test-priors.R`.

**Interfaces:**
- Produces `mmBPFA_sampler(dat, priors = "uniform", H = NULL, K_start = 10, warmup = 1000, iter = 5000, thin = 5, chains = 1, n_cores = 1, threads = 1, seed = 101473) -> bpfa.results`.
  - Draws are stored per kept iteration: `lambda` (p×H), `z`, `c`, `omega` (n×H), `alpha`, `theta`, `thresholds` (a list per column), `log_lik`.
  - The `priors` attribute holds the full resolved list.
- Produces `resolve_priors(priors, p) -> list` (the presets: "uniform", and "conservative" with slab rate 0.25) and `get_priors(x)` for both raw and processed results.

**Steps:**
- [ ] **Tests first:**
  - `get_priors(res)` returns every field of the resolved list, including `spike_var`, `H`, and the preset name;
  - raw list entries override the preset;
  - `K_start = 40` with H = 29 gives the message "K_start (40) exceeds H (29); using 29" and runs (Review Focus 4);
  - H < 1 gives a clear error;
  - chains run sequentially and in parallel;
  - the same seed with `threads = 1` and `threads = 4` gives identical draws (this test lands in Task 11);
  - the sampler runs with `iter < 70`.
- [ ] Confirm the tests fail.
- [ ] Implement the driver loop in C++: x, then loadings, then ω, then α, with thresholds and log-likelihood computed on kept iterations.
- [ ] Delete the IBP code and its tests. Update `NAMESPACE`, regenerate `RcppExports` and the Rd files.
- [ ] Confirm the tests pass, then commit.

### Task 9: Post-processing per design §4.4–4.5

**Files:** modify `R/process_results.R` and `R/calculate_loglik.R` (sparse thresholds and cell log-probabilities); test `tests/testthat/test-process-results.R`, `test-likelihood.R`.

**Interfaces:** `process_results(mcmc_object, SE = TRUE, p = c(0.05, 0.95), to_tibble = TRUE)` keeps its signature. Its output adds:
- `"K Distribution"`: a table of K over draws;
- `"Uncounted Correlation Share"`: a numeric value;
- `"Priors"`.

"Factor Loadings", "Copula Loadings" and "Factor Scores" cover counted columns only. The model correlation used by `check_accuracy` uses all H columns. `check_accuracy` takes item levels from `prepare_data()`, since `calculate_d_mask()` is deleted in Task 8.

**Steps:**
- [ ] **Tests first**, using fixtures:
  - K counts only on columns with at least 3 items;
  - the uncounted share is 0 when every column is counted, and positive when a spike column has non-zero loadings;
  - relabelling still works on alternating draws;
  - GMP and log-likelihood still match the hand-computed values on sparse input.
- [ ] Confirm the tests fail, implement, confirm they pass.
- [ ] Commit.

### Task 10: Recovery acceptance test

**Files:** rewrite `tests/testthat/test-recovery.R` (`skip_on_cran`).

- [ ] Port a reduced version of the spike's binary grid: K ∈ {2, 5}, P = 60, starts {2, 30}.
- [ ] Assert the soft criteria as warnings with their sizes, plus hard checks: starts agree within 1, and correlation error is below 0.08.
- [ ] Compare against `fixtures/baseline_ibp_grid.rds`.
- [ ] Run the test and commit.

### Task 11: OpenMP parallelism with deterministic streams

**Files:** modify the `src/*.cpp` kernels from Tasks 4–8 and `src/Makevars*`.

- [ ] **Test first:** `threads = 1` and `threads = 4` give identical draws for the same seed.
- [ ] Parallelize:
  - the x-step over columns;
  - α over columns;
  - ω over rows;
  - the item loop inside each loadings column.

  Each parallel unit gets its own `Rng` stream (seed, iteration, unit). Reductions run in a fixed order.
- [ ] Confirm the test passes, and that the full suite is green.
- [ ] Commit.

### Task 12: Benchmarks and optimization loop

**Files:** create `bench/scaling.R`, which writes `bench/results/*.csv`; record the summary in design §8.

- [ ] **Scaling grid:**
  - N ∈ {200, 2k, 20k};
  - P ∈ {100, 1k, 10k};
  - sparsity ∈ {0, 90%, 99%, 99.9%};
  - K = 5, threads ∈ {1, 8}.

  Measure seconds per sweep, ESS/s for K and the leading copula correlations, and peak memory.
- [ ] Profile the slowest cell with `perf`, or Rprof for the R parts. Optimize the top hotspot, rerun the full test suite, rerun the benchmark, and commit. Repeat until either:
  - the sweep cost is within 2× of the O(nnz·K + (N+P)K³) operation count; or
  - the remaining hotspots are inherent to that count.
- [ ] Rerun the Stage-2 stress tests on the package implementation at larger N and P, and update design §8.

### Task 13: Docs

- [ ] Update `CLAUDE.md` (architecture, sampler layering, test commands), `DESCRIPTION` (the description text and the beta-process wording is the author's call, so leave a TODO for them), and `NEWS.md`: the IBP removed, pointing to the tags; the new priors; and `get_priors()`.
- [ ] Commit.
