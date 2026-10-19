quiet_run <- function(...) {
  out <- NULL
  utils::capture.output(out <- mmBPFA_sampler(...))
  out
}

small_binary <- function(seed = 1, n = 60, p = 40) {
  generate_sim_data(n = n, p = p, K_true = 2, noise = 0, seed = seed)$data
}

test_that("the sampler stores one draw per kept iteration with the documented shapes", {
  res <- quiet_run(small_binary(), K_start = 5, warmup = 5, iter = 20, thin = 2, seed = 1)
  ch <- res[[1]]
  expect_s3_class(res, "bpfa.results")
  expect_length(ch$K, 10)
  expect_length(ch$log_lik, 10)
  expect_equal(dim(ch$lambda[[1]]), c(40L, 19L))     # H = min(floor(39 / 2), 50)
  expect_equal(dim(ch$omega[[1]]), c(60L, 19L))
  expect_length(ch$alpha[[1]], 40)
  expect_true(all(is.finite(ch$log_lik)))
})

test_that("every prior actually used is retrievable from the results", {
  res <- quiet_run(small_binary(), priors = list(spike_var = 0.04), K_start = 3, warmup = 2, iter = 4, thin = 1)
  pr <- get_priors(res)
  expect_equal(pr$spike_var, 0.04)
  expect_equal(pr$preset, "uniform")
  expect_true(all(c("slab_shape", "slab_rate", "a_cusp_shape", "a_cusp_rate", "rho_a", "rho_b", "H") %in% names(pr)))
})

test_that("K_start above H is capped with a message", {
  expect_message(quiet_run(small_binary(p = 20), K_start = 40, warmup = 1, iter = 2, thin = 1),
                 "K_start \\(40\\) exceeds H \\(9\\); using 9")
})

test_that("the sampler runs with fewer than 70 iterations and with sparse data", {
  y <- small_binary(); y[sample(length(y), round(0.9 * length(y)))] <- NA
  res <- quiet_run(y, K_start = 4, warmup = 3, iter = 10, thin = 1)
  expect_length(res[[1]]$K, 10)
})

test_that("the same seed reproduces the chain; chains differ from each other", {
  y <- small_binary()
  a <- quiet_run(y, K_start = 4, warmup = 3, iter = 6, thin = 1, seed = 9)
  b <- quiet_run(y, K_start = 4, warmup = 3, iter = 6, thin = 1, seed = 9)
  expect_identical(a[[1]]$lambda, b[[1]]$lambda)
  two <- quiet_run(y, K_start = 4, warmup = 3, iter = 6, thin = 1, chains = 2, seed = 9)
  expect_length(two, 2)
  expect_false(identical(two[[1]]$lambda, two[[2]]$lambda))
})

test_that("forked parallel chains match sequential chains for the same seed", {
  skip_on_os("windows")
  y <- small_binary()
  seqr <- quiet_run(y, K_start = 4, warmup = 3, iter = 6, thin = 1, chains = 2, seed = 3)
  par <- quiet_run(y, K_start = 4, warmup = 3, iter = 6, thin = 1, chains = 2, n_cores = 2, seed = 3)
  expect_identical(seqr[[2]]$lambda, par[[2]]$lambda)
})

test_that("socket-cluster chains (the Windows path) run", {
  skip_if_not("mmBPFA" %in% rownames(utils::installed.packages()), "workers need the installed package")
  y <- small_binary()
  out <- run_chains_psock(chains = 2, n_cores = 2,
                          chain_fun = function(chain) mmBPFA:::run_one_chain(prepare_data(y),
                            resolve_priors("uniform", ncol(y)), 4, 3, 6, 1, 3 + chain))
  expect_length(out, 2)
})

test_that("within-chain threads do not change the draws", {
  skip_if_not(.openmp_available(), "built without OpenMP")
  for (miss in c(0, 0.5)) {
    y <- small_binary(n = 120, p = 40); y[matrix(runif(length(y)) < miss, nrow(y))] <- NA
    one <- quiet_run(y, K_start = 4, warmup = 5, iter = 10, thin = 1, threads = 1, seed = 11)
    four <- quiet_run(y, K_start = 4, warmup = 5, iter = 10, thin = 1, threads = 4, seed = 11)
    expect_identical(one[[1]]$lambda, four[[1]]$lambda)
    expect_identical(one[[1]]$omega, four[[1]]$omega)
    expect_identical(one[[1]]$log_lik, four[[1]]$log_lik)
  }
})

test_that("runs whose stored draws would exceed the memory limit stop before sampling", {
  y <- small_binary()
  expect_error(mmBPFA_sampler(y, iter = 1000, thin = 1, max_memory_gb = 1e-3),
               "Stored draws would need about .* GB .*increase `thin`")
  expect_silent(utils::capture.output(mmBPFA_sampler(y, warmup = 1, iter = 2, thin = 1, max_memory_gb = 1e-3)))
})

test_that("each chain reports its sampling time per sweep", {
  res <- quiet_run(small_binary(), K_start = 4, warmup = 5, iter = 10, thin = 1)
  expect_true(is.numeric(res[[1]]$seconds_per_sweep) && res[[1]]$seconds_per_sweep > 0)
})
