quiet_run <- function(...) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(mmBPFA_sampler(...)))
  out
}

small_binary <- function(seed = 1) {
  generate_sim_data(n = 60, p = 40, K_true = 2, noise = 0, seed = seed)$data
}

test_that("the sampler runs with fewer than 70 main iterations", {
  res <- quiet_run(small_binary(), mode = "fixed", K_start = 5, pre_sparse = 2,
                   warmup = 2, iter = 10, thin = 1)
  expect_length(res[["Log-Likelihoods"]], 10)
})

test_that("process_results pools chains that were run sequentially", {
  res <- quiet_run(small_binary(), mode = "fixed", K_start = 5, pre_sparse = 2,
                   warmup = 2, iter = 10, thin = 1, chains = 2)
  utils::capture.output(out <- process_results(res, SE = FALSE))
  expect_equal(nrow(out[["log-Likelihoods"]]), 20)
})

test_that("socket-cluster chains (the Windows path) run and return one result per chain", {
  skip_if_not("mmBPFA" %in% rownames(utils::installed.packages()), "workers need the installed package")
  dat <- small_binary()
  margins <- calculate_d_mask(dat, "fixed")
  init <- initialize_sampler(dat, "fixed", margins$n_levels, margins$margin_vals, K_start = 5)
  args <- c(list(pre_sparse = 2, warmup = 2, iter = 4, thin = 1, dat = dat, mode = "fixed",
                 share_noise = FALSE, d_mask = margins$d_mask, n_levels = margins$n_levels,
                 margin_vals = margins$margin_vals, benchmark = FALSE),
            init[c("x", "lambda", "zeros", "omega", "alpha", "gamma_k", "d", "dc", "ibp_a", "ibp_b")])
  utils::capture.output(out <- run_chains_psock(chains = 2, n_cores = 2, chain_args = args))
  expect_length(out, 2)
  expect_length(out[[2]][["Log-Likelihoods"]], 4)
})
