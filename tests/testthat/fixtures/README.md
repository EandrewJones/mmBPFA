# Test fixtures

`baseline_ibp_grid.rds`: recovery-grid results for the collapsed two-parameter IBP with the copula-calibrated slab (gamma_k ~ Gamma(1, rate 1/2)), the "candidate A" baseline in the design doc (§7).

- Sampler code: commit 90e383c on `fix/sampler-correctness`. Tag `v1.1-ibp` is pending the author's go-ahead.
- Settings: N = 200; P ∈ {60, 200}; K ∈ {2, 5}; binary and binary+count margins; 5% MCAR. Starts K = 2, 10, 30; truncation 30; 2,000 iterations with 1,000 burn-in.
- Data: data seed 100 + K + P, generated with `sim_data()` from the spike harness.
- Columns:
  - `K_med`: median number of columns with at least 3 active items;
  - `C_err`: mean absolute error of the copula correlation;
  - `F1`: item-map F1;
  - `secs`: runtime in seconds.
