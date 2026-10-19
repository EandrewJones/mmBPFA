test_that("thresholds sit midway in the gap between adjacent levels of x", {
  x <- cbind(c(-2, -1, 1, 3), c(-2, 0, 1, 4))
  y <- cbind(c(0, 0, 1, 1), c(0, 1, 1, 2))
  expect_equal(estimate_thresholds(x, y), list(0, c(-1, 2.5)))
})

test_that("thresholds ignore missing cells", {
  x <- cbind(c(-2, 5, 1, 3))
  y <- cbind(c(0, NA, 1, 1))
  expect_equal(estimate_thresholds(x, y), list(-0.5))
})

test_that("log-likelihood is the probit probability of the observed category", {
  # mu = (1, -1) for y = (1, 0) with threshold 0: each cell has probability pnorm(1)
  lambda <- matrix(1, 1, 1); omega <- matrix(c(1, -1), 2, 1)
  ll <- calculate_log_lik(lambda, omega, 0, cbind(c(1, 0)), list(0))
  expect_equal(ll, 2 * log(0.8413447), tolerance = 1e-6)
})

test_that("log-likelihood handles middle categories and skips missing cells", {
  # mu = 0, y = (0, 1, 2, NA), thresholds (-1, 1):
  # pnorm(-1) = 0.1586553, pnorm(1) - pnorm(-1) = 0.6826895, 1 - pnorm(1) = 0.1586553
  lambda <- matrix(0, 1, 1); omega <- matrix(0, 4, 1)
  ll <- calculate_log_lik(lambda, omega, 0, cbind(c(0, 1, 2, NA)), list(c(-1, 1)))
  expect_equal(ll, 2 * log(0.1586553) + log(0.6826895), tolerance = 1e-6)
})

test_that("GMP is the geometric mean predictive probability of the observed data", {
  # two binary items, each cell predicted correctly with probability pnorm(1)
  dat <- cbind(c(1, 0), c(0, 1))
  pred <- cbind(c(1, -1), c(-1, 1))
  expect_equal(gmp(dat, pred, list(c(0, Inf), c(0, Inf))), 0.8413447, tolerance = 1e-6)
})
