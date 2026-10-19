# Parallel tempering. Replica r targets pi_b(theta, x) ∝ p(theta) prod N(x_ij; lambda_j' omega_i, 1)^b 1{x in D},
# i.e. noise variance 1 / b. Because the rank likelihood is scale invariant, x~ = sqrt(b) x and
# lambda~ = sqrt(b) lambda turn it into the untempered model with slab rate and spike variance
# multiplied by b, so every replica runs the ordinary kernels. Adjacent replicas swap states with
# probability min(1, exp((b_r - b_s) (l_s - l_r))), l = -1/2 sum (x - lambda' omega)^2 in the
# original scale. Only the b = 1 replica is kept.

tempering_data <- function(seed, sparse) {
  set.seed(seed); n <- 150; p <- 12
  y <- (outer(rnorm(n), rnorm(p)) + matrix(rnorm(n * p), n, p) > 0) * 1
  if (sparse) y[matrix(runif(n * p) < 0.8, n, p)] <- NA
  y
}

copula_cor_draws <- function(res, p) {
  lapply(res, function(ch) sapply(ch$lambda, function(L) {
    C <- tcrossprod(L) + diag(p); D <- 1 / sqrt(diag(C)); (C * outer(D, D))[upper.tri(C)] }))
}

batch_se <- function(cc) {
  apply(do.call(cbind, lapply(cc, function(m) sapply(split(seq_len(ncol(m)), rep(1:20, each = ncol(m) / 20)),
                                                     function(b) rowMeans(m[, b, drop = FALSE])))), 1, sd) /
    sqrt(20 * length(cc))
}

for (sparse in c(FALSE, TRUE)) {
  test_that(paste("the tempered chain targets the same posterior as the plain chain:",
                  if (sparse) "sparse" else "dense", "data"), {
    skip_on_cran()
    y <- tempering_data(6, sparse)
    fit <- function(temperatures, seed) {
      utils::capture.output(r <- mosaicFA_sampler(y, H = 4, K_start = 2, warmup = 500, iter = 8000, thin = 2,
                                                  chains = 2, seed = seed, temperatures = temperatures))
      cc <- copula_cor_draws(r, 12)
      list(mean = rowMeans(do.call(cbind, cc)), se = batch_se(cc), K = unlist(lapply(r, `[[`, "K")))
    }
    a <- fit(1, 1); b <- fit(c(1, 0.6, 0.35, 0.2), 2)
    expect_z(a$mean, b$mean, sqrt(a$se^2 + b$se^2), max_z = 5)
  })
}

test_that("tempered runs report swap acceptance between adjacent replicas", {
  y <- tempering_data(7, FALSE)
  utils::capture.output(r <- mosaicFA_sampler(y, H = 4, K_start = 2, warmup = 20, iter = 40, thin = 2,
                                              temperatures = c(0.3, 1, 0.6)))
  rate <- r[[1]]$swap_rate
  expect_length(rate, 2)
  expect_true(all(rate >= 0 & rate <= 1))
  expect_equal(attr(r, "temperatures"), c(1, 0.6, 0.3))      # sorted, cold replica first
  expect_length(r[[1]]$K, 20)
})

test_that("temperatures must lie in (0, 1] and include 1", {
  y <- tempering_data(8, FALSE)
  expect_error(mosaicFA_sampler(y, iter = 2, warmup = 0, thin = 1, temperatures = c(0.5, 0.2)), "include 1")
  expect_error(mosaicFA_sampler(y, iter = 2, warmup = 0, thin = 1, temperatures = c(1, 0)), "\\(0, 1\\]")
})

test_that("tempered draws do not depend on the number of threads", {
  skip_if_not(.openmp_available(), "built without OpenMP")
  y <- tempering_data(9, TRUE)
  run <- function(threads) {
    utils::capture.output(r <- mosaicFA_sampler(y, H = 4, K_start = 2, warmup = 5, iter = 10, thin = 1,
                                                temperatures = c(1, 0.5), threads = threads, seed = 4))
    r[[1]]$lambda
  }
  expect_identical(run(1), run(4))
})
