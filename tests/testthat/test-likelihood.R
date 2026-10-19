# Reference implementations of the level thresholds and the ordered-probit
# log-likelihood, written independently in R. The sampler computes both in C++.

ref_thresholds <- function(x, dat) {
  lapply(seq_len(ncol(dat)), function(j) {
    obs <- !is.na(dat[, j]); levels <- sort(unique(dat[obs, j]))
    if (length(levels) < 2) return(numeric(0))
    vapply(seq_len(length(levels) - 1), function(v)
      (max(x[obs & dat[, j] == levels[v], j]) + min(x[obs & dat[, j] == levels[v + 1], j])) / 2, numeric(1))
  })
}

ref_log_lik <- function(lambda, omega, dat, thresholds) {
  mu <- omega %*% t(lambda)
  sum(cell_log_probs(mu, dat, thresholds), na.rm = TRUE)
}

test_that("the reference thresholds sit midway in the gap between adjacent levels", {
  x <- cbind(c(-2, -1, 1, 3), c(-2, 0, 1, 4))
  y <- cbind(c(0, 0, 1, 1), c(0, 1, 1, 2))
  expect_equal(ref_thresholds(x, y), list(0, c(-1, 2.5)))
})

test_that("cell log-probabilities are the ordered-probit probability of the observed level", {
  # mu = (1, -1) for y = (1, 0) with threshold 0: each cell has probability pnorm(1)
  expect_equal(ref_log_lik(matrix(1, 1, 1), matrix(c(1, -1), 2, 1), cbind(c(1, 0)), list(0)),
               2 * log(0.8413447), tolerance = 1e-6)
  # mu = 0, y = (0, 1, 2, NA), thresholds (-1, 1)
  expect_equal(ref_log_lik(matrix(0, 1, 1), matrix(0, 4, 1), cbind(c(0, 1, 2, NA)), list(c(-1, 1))),
               2 * log(0.1586553) + log(0.6826895), tolerance = 1e-6)
})

test_that("the sampler's log-likelihood matches the reference on mixed data with missing cells", {
  set.seed(1)
  n <- 50; p <- 6; K <- 2
  y <- cbind(matrix(rbinom(n * 3, 1, 0.4), n, 3), matrix(rpois(n * 2, 2), n, 2), rnorm(n))
  y[sample(length(y), 40)] <- NA
  d <- prepare_data(y); x <- initial_x(d)
  lam <- matrix(rnorm(p * K, 0, 0.5), p, K); om <- matrix(rnorm(n * K), n, K)
  want <- ref_log_lik(lam, om, y, ref_thresholds(cells_to_matrix(d, x), y))
  expect_equal(.log_likelihood_from_x(d, x, lam, om), want, tolerance = 1e-10)
})

test_that("GMP is the geometric mean predictive probability of the observed data", {
  dat <- cbind(c(1, 0), c(0, 1))
  pred <- cbind(c(1, -1), c(-1, 1))
  expect_equal(gmp(dat, pred, list(c(0, Inf), c(0, Inf))), 0.8413447, tolerance = 1e-6)
})
