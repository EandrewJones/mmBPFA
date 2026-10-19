# mmBPFA 1.1.0

## Breaking changes

* The IBP loadings prior and its sampler are removed. Factors now use a cumulative shrinkage prior over columns with item-level inclusion within columns (see `?resolve_priors`). The corrected IBP sampler is preserved at commit 90e383c. The original version that produced the published results is commit 761a9ba.
* `mmBPFA_sampler()` takes `priors`, `H`, `K_start`, `warmup`, `iter`, `thin`, `chains`, `n_cores`, `threads`, `seed` and `max_memory_gb`. The `mode`, `pre_sparse`, `share_noise`, `parallel` and `benchmark` arguments are gone: margins are handled per item automatically.
* `process_results()` output: "Slab Variances" replaces "Factor Precisions". New fields are "K Distribution", "R-hat K", "Copula Loadings", "Uncounted Correlation Share" and "Priors".
* `process_benchmarks()` and `calculate_d_mask()` are removed.
* Item intercepts are removed from the model. They are not identified under the rank likelihood, and their flat prior made their posterior improper, so the reported intercepts and cutpoints drifted between chains. Item difficulty is now reported as "Item Thresholds": one per pair of adjacent observed values, with credible intervals for items with at most 10 levels. "Item-level Intercepts" and the stored `alpha` draws are gone.

## New features

* `get_priors()` returns every prior setting used for a fit, including defaults.
* Long-format input (`row`, `col`, `value`) so sparse data never build the N x P matrix.
* The whole chain runs in C++ on observed cells only, with OpenMP threads (`threads`). Results are bit-identical for any thread count.
* R-hat for the number of dimensions across chains, with a warning when chains disagree.
* Optional scale moves (`moves = c("item_scale", "column_scale")`). They leave the posterior unchanged and speed up early mixing slightly.

## Bug fixes

* Many sampler correctness fixes, which are listed in the commit message of 90e383c.
