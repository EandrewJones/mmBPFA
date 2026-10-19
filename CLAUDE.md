# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`mmBPFA` is an R package (Rcpp/RcppArmadillo) for Bayesian factor analysis of any mix of binary, ordinal, count and continuous items. The margins are left unmodelled: the data are treated as monotone transforms of a latent Gaussian copula, through the extended rank likelihood (Hoff 2007; Murray et al. 2013).

The number of dimensions K is learned with a cumulative shrinkage (CUSP) prior over columns (Legramanti, Durante & Dunson 2020). Within each column, item inclusion is learned separately (Schiavon, Canale & Dunson 2022). The fit therefore says both how many dimensions there are and which items load on each one.

The original IBP sampler was removed. The corrected IBP version is commit 90e383c. The design and its evidence are in `docs/superpowers/specs/2026-10-19-loadings-prior-design.md`, and the implementation plan is in `docs/superpowers/plans/`.

## Commands

There is no CI or vignette directory.

```sh
# Run the test suite from source; NOT_CRAN enables the slower recovery and timing tests
NOT_CRAN=true Rscript -e 'devtools::test()'
# One test file, e.g. tests/testthat/test-sample-loadings.R
NOT_CRAN=true Rscript -e 'devtools::test(filter = "sample-loadings")'

# Regenerate src/RcppExports.cpp and R/RcppExports.R after changing any // [[Rcpp::export]] signature
Rscript -e 'Rcpp::compileAttributes()'
# Regenerate man/*.Rd only (see the NAMESPACE note below)
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'

# Full package check
R CMD build . && R CMD check --no-manual mmBPFA_1.0.tar.gz

# Scaling benchmark (N, P, sparsity, threads); not part of the package build
Rscript bench/scaling.R results.csv
```

**Timing tests need an optimized build.** `devtools::test()` / `pkgload::load_all()` compile at `-O0`, which makes the timing tests fail. Run `pkgbuild::compile_dll(debug = FALSE)` and then `load_all(compile = FALSE)`. Delete `src/*.o src/*.so` first if a debug build is already present.

**Some tests need the installed package.** The socket-cluster test in `test-sampler-run.R` is skipped unless mmBPFA is installed, because PSOCK workers call `library(mmBPFA)`.

## Testing conventions

These are MCMC samplers, so most tests are statistical:
- **Conditional samplers.** Fix every other input at known values, draw many samples, and compare against the analytic conditional posterior.
- **Prior invariance.** With no data, the sampler must reproduce the prior; for example, the CUSP prior on the number of active columns.
- **Recovery.** Simulate from the generative model with known parameters and check recovery of identified quantities. The copula correlation matrix is invariant to rotation, label switching and scale.
- **Determinism.** The same seed must give bit-identical draws for any number of threads.

Use `expect_z()` (a Monte Carlo z-score against the standard error) for means and `expect_rel()` for variances; both are in `helper-sim.R`. Don't use `expect_equal(tolerance =)` for these: in testthat edition 3 it silently switches to an absolute comparison when the expected value is smaller than the tolerance.

## NAMESPACE is hand-maintained

`NAMESPACE` has no roxygen header and uses `exportPattern("^[[:alpha:]]+")`. So every function whose name starts with a letter is exported, including the "internal" ones. C++ test entry points use dot-prefixed export names (`.sample_loadings`, `.run_chain`, ...) so they stay internal.

The `@import` tags in `R/mmBPFA.R` do **not** drive `NAMESPACE`: add dependencies to `NAMESPACE` and `DESCRIPTION` by hand. Don't let `devtools::document()` overwrite `NAMESPACE`.

## Architecture

**Pipeline:** `generate_sim_data()` (optional) → `mmBPFA_sampler()` → `process_results()` → `check_accuracy()` / `check_logliks()`.

**Data layout (`R/prepare_data.R`, `src/data_structure.h`).** Only observed cells are stored:
- **Column layout (CSC):** cells are stored column by column, sorted by value within each column, so the rank-likelihood levels are contiguous.
- **Row index (CSR):** a row-wise index into the same cells, used by the row updates.
- **Missing cells** are treated as missing at random and never imputed, so a sweep costs O(nnz·H + N·H³).
- **Long-format input** (`row`, `col`, `value`) avoids ever building the N×P matrix.
- **Latent x** lives as one value per observed cell. It must start rank-consistent (`initial_x`).

**Sampler.**
- `mmBPFA_sampler()` (R) resolves the priors (`R/priors.R`; `get_priors()` retrieves them) and checks the memory needed for stored draws.
- `run_one_chain()` (R) initializes a chain from `set.seed(seed + chain - 1)` and calls `.run_chain` (C++, `src/run_chain.cpp`).
- `.run_chain` runs the whole chain. Each sweep is `sample_x_sparse` → `sample_loadings` → `sample_omega_sparse`, then the optional scale moves (`moves`, `src/sample_moves.cpp`). On kept draws it also computes K and the log-likelihood, with thresholds taken from the gaps between levels of x.
- **No item intercepts.** Under the rank likelihood an intercept is not identified (shifting x_j and α_j together changes nothing), and with a flat prior its posterior is improper. So the model has none: the latent location is fixed by ΛΩ', and item difficulty is the level thresholds ("Item Thresholds"). Don't reintroduce a non-identified parameter. Any move that touches one breaks invariance (see design doc §9).
- Chains run sequentially, with `mclapply`, or on a PSOCK cluster on Windows.

**Kernels (`src/`).**
- `sample_x_sparse.cpp`: sequential-by-level truncated normals, O(n_j) per column.
- `sample_loadings.cpp`:
  - item inclusion and loadings with a running residual over observed cells; a dense residual matrix is used when the data are fully observed;
  - CUSP column indicators, given the current loadings;
  - the concentration `a_cusp`, slice-sampled with the sticks integrated out. Conditioning on the sticks traps the chain at 0.
- `sample_omega_sparse.cpp`:
  - fully observed rows share one precision matrix and are drawn in one batch;
  - rows that are mostly observed subtract their missing items' outer products from it;
  - sparse rows accumulate their observed items.
- `rng.h`: xoshiro256++ streams keyed by (seed, iteration, kernel, unit), plus Robert's (1995) truncated normal.

**Parallelism (`threads.h`).** OpenMP is used over columns, rows and items, with a fixed-order reduction for the log-likelihood. Never draw from R's RNG or call the R API inside a parallel region. Never share an RNG stream between kernels: use `stream_id(kernel, unit)`.

**Dimension conventions.**
- `lambda` and `z` are P×H, `omega` is N×H, and `c`, `v`, `rho`, `theta` have length H.
- Column k (1-based) is **on** when `c[k] > k`.
- A draw's **K** is the number of on columns with at least 3 active items (an identification requirement).

**Post-processing (`R/process_results.R`).**
- Draws from all chains are pooled. Counted dimensions are matched to a median-K pivot draw by sparsity pattern (`align_dimensions`), then **sign-aligned** to the pivot, so sign flips across chains don't cancel.
- `find_stable_dims` keeps the pivot dimensions present in at least half of the draws.
- Loadings are reported on the sampler scale and on the copula scale (λ_j/√(1 + ‖λ_j‖²), from each draw's full row).
- With several chains, R-hat for K is reported, with a warning above 1.1.
- "Uncounted Correlation Share" measures how much correlation uncounted columns carry.

**Known limits (see the design doc, §7b–8).**
- **Information-poor data mix slowly over K.** With rare positives or very sparse data, chains from extreme starts can take about 15k iterations to agree. Use several chains and check R-hat.
- **Model assumptions:**
  - zeros unrelated to the factors (structural zero-inflation) break the rank likelihood's monotonicity;
  - top-k lists encoded as rank positions are missing not at random;
  - one-hot nominal items create spurious dimensions.
