# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`mmBPFA` is an R package (Rcpp/RcppArmadillo) implementing Multi-Modal Beta-Process Factor Analysis: a Bayesian factor analysis that puts a beta-process / Indian Buffet Process prior on the loadings so the number of latent dimensions K is learned, and treats observed margins (continuous, discrete, or mixed) as manifestations of a latent Gaussian copula via the extended rank likelihood (Hoff 2007; Murray et al. 2013).

## Commands

There is no CI or vignette directory. Development is done from R (an RStudio project is present).

```sh
# Run the test suite from source (compiles src/ via load_all); NOT_CRAN enables the slower recovery test
NOT_CRAN=true Rscript -e 'devtools::test()'
# Run one test file, e.g. tests/testthat/test-sample-x.R
NOT_CRAN=true Rscript -e 'devtools::test(filter = "sample-x")'

# Regenerate src/RcppExports.cpp and R/RcppExports.R after changing any // [[Rcpp::export]] signature
Rscript -e 'Rcpp::compileAttributes()'

# Regenerate man/*.Rd only (see the NAMESPACE note below)
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'

# Full package check
R CMD build . && R CMD check --no-manual mmBPFA_1.0.tar.gz
```

The socket-cluster test in `test-sampler-run.R` is skipped unless mmBPFA is installed, because PSOCK workers call `library(mmBPFA)`.

## Testing conventions

These are MCMC samplers, so most tests are statistical:
- **Conditional samplers:** fix every other input at known values, draw many samples, and compare against the analytic conditional posterior.
- **Integration:** simulate from the generative model with known parameters (`helper-sim.R`) and check recovery of identified quantities. The copula correlation matrix is invariant to rotation, label switching and scale.

Use `expect_z()` (a Monte Carlo z-score against the standard error) for means, and `expect_rel()` for variances. Don't use `expect_equal(tolerance =)` for these: in testthat edition 3 it silently switches to an absolute comparison when the expected value is smaller than the tolerance, which makes small variances and probabilities untestable.

## NAMESPACE is hand-maintained

`NAMESPACE` has no roxygen header and uses `exportPattern("^[[:alpha:]]+")`, so every function whose name starts with a letter is exported, including the "internal" ones. The `@import` tags in `R/mmBPFA.R` do **not** drive `NAMESPACE`. If you add a dependency, add it to `NAMESPACE` and `DESCRIPTION` `Imports` by hand. Don't let `devtools::document()` overwrite `NAMESPACE` unless you mean to switch to roxygen-managed exports.

## Architecture

**Pipeline:** `generate_sim_data()` (optional) → `mmBPFA_sampler()` → `process_results()` → diagnostics (`check_logliks()`, `check_accuracy()`, `process_benchmarks()`).

**Sampler layering (`R/mmBPFA_sampler.R`):**
- `mmBPFA_sampler()` validates input, computes the rank-likelihood mask (`calculate_d_mask()`), initializes state (`initialize_sampler()`), then runs one or more chains (sequential, `mclapply` on Unix, or a PSOCK cluster on Windows). It returns a list of class `bpfa.results`, with run settings stored as attributes (`mode`, `parallel`, `thin`, ...). Downstream code reads those attributes.
- `chain_call()` runs three phases: **pre-sparse** (no IBP prior, z fixed at 1, no new dimensions; helps mixing and keeps K from collapsing to 0), **warmup**, and **main**. Draws are stored every `thin` steps in a tibble of list-columns (`zeros`, `lambda`, `omega`, `alpha`, `gamma_k`, `cutpoints`).
- `sampler_phase()` is one Gibbs sweep, in this order: `sample_x` (R) → `sample_lambda_and_zeros` (C++, which also proposes new dimensions when `infinite`) → `remove_inactive` + `sample_IBP_a/b` (warmup/main only) → `sample_omega` → `sample_alpha` → `sample_gamma_k` → `sample_d` (only if `share_noise`). In the main phase only, it then runs `remove_weakly_active` and `calculate_log_lik`.
- State moves between steps by destructuring named lists with zeallot `c(...) %<-%`. The order of elements in each function's returned `output_list` must match the order of names on the left-hand side at every call site. If you add or reorder outputs, update every destructuring site.

**Dimension conventions:** `dat`/`x` are N×P, `lambda`/`zeros` are P×K, `omega` is N×K, `alpha` has length P, and `gamma_k`/`dc` have length K. K changes between iterations as dimensions are added or removed, so stored draws have ragged K. The `k_*` variables are the copies with weakly active dimensions pruned, and those are what get stored.

**`mode` (`"fixed"`, `"multi"`, `"mixed"`)** controls how margins are handled:
- `"fixed"`: all columns share one set of levels. `d_mask` and `margin_vals` are global.
- Anything else: per-column `margin_vals`/`n_levels` and a nested per-column `d_mask`.

`calculate_d_mask`, `initialize_sampler`, `sample_x`, `process_results`, `check_accuracy`, and `generate_sim_data` all branch on `mode`, so a change to margin handling usually touches all of them.

**C++ (`src/`):** each `.cpp` file holds one conditional sampler, exported with `// [[Rcpp::export(name = "...")]]` and documented with `//'` roxygen comments. `src/RcppExports.cpp` and `R/RcppExports.R` are generated, so don't edit them by hand. `Makevars` builds with C++11 and links OpenMP, LAPACK, and BLAS.

**Post-processing (`R/process_results.R`):**
- Label switching: `align_dimensions` matches each draw's columns to a pivot draw (median K), and `find_stable_dims` keeps the pivot dimensions present in at least half of the draws.
- Loadings are reported on both the sampler scale ("Factor Loadings") and the copula correlation scale ("Copula Loadings", λ_j/√(1 + ‖λ_j‖²) per draw). Λ is never normalized inside the sampler, because that changes the stationary distribution.
- Likelihood and GMP use the ordered-probit probability of each observed cell, with thresholds taken from the gaps between levels of x (`estimate_thresholds`, `cell_log_probs` in `R/calculate_loglik.R`).
- A large commented-out older version of `process_results` sits in the middle of the file; it is dead code.
