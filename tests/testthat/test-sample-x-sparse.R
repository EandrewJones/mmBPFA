# Sparse (observed-cells-only) data layout and the C++ extended-rank-likelihood x-step.

test_that("prepare_data round-trips the observed data, missing cells included", {
  dat <- rbind(c(1, NA, 0.5), c(0, 3, NA), c(1, 2, 0.1), c(NA, 2, 0.9))
  d <- prepare_data(dat)
  expect_equal(d$n, 4L); expect_equal(d$p, 3L); expect_equal(length(d$value), 9L)
  expect_identical(cells_to_matrix(d, d$value), dat)
})

test_that("constant, empty and all-distinct columns get the right number of levels", {
  dat <- cbind(c(1, 1, 1, NA), c(NA, NA, NA, NA), c(0.3, 0.1, 0.2, 0.4))
  d <- prepare_data(dat)
  expect_equal(d$n_levels, c(1L, 0L, 4L))
})

ordered_ok <- function(x, y) {
  all(vapply(seq_len(ncol(y)), function(j) {
    obs <- !is.na(y[, j]); lv <- sort(unique(y[obs, j]))
    if (length(lv) < 2) return(TRUE)
    all(vapply(seq_len(length(lv) - 1), function(v)
      max(x[obs & y[, j] == lv[v], j]) < min(x[obs & y[, j] == lv[v + 1], j]), TRUE))
  }, TRUE))
}

test_that("initial_x is strictly increasing across levels and tied within them", {
  y <- cbind(c(2, 0, 1, 1, NA), c(0.5, 0.1, 0.9, 0.3, 0.7))
  d <- prepare_data(y); x <- cells_to_matrix(d, initial_x(d))
  expect_true(ordered_ok(x, y))
  expect_equal(x[3, 1], x[4, 1])
})

# mean-zero chain: pooled x must be N(0, 1) and the ordering must hold every sweep
run_sparse_chain <- function(y, iters) {
  d <- prepare_data(y)
  x <- initial_x(d)
  pooled <- numeric(0); ok <- ordered_ok(cells_to_matrix(d, x), y)
  for (t in seq_len(iters)) {
    x <- .sample_x_sparse(d, x, matrix(0, d$p, 1), matrix(0, d$n, 1), 1, t)
    ok <- ok && ordered_ok(cells_to_matrix(d, x), y)
    if (t > iters / 2) pooled <- c(pooled, x)
  }
  list(ordered = ok, pooled = pooled)
}

test_that("binary draws keep the rank ordering and are N(0, 1) at stationarity", {
  set.seed(1)
  r <- run_sparse_chain(matrix(rbinom(40 * 30, 1, 0.3), 40, 30), 200)
  expect_true(r$ordered)
  expect_lt(abs(mean(r$pooled)), 0.05)
  expect_rel(var(r$pooled), 1, 0.05)
})

test_that("ordinal draws with missing cells keep the ordering and are N(0, 1)", {
  set.seed(2)
  y <- matrix(sample(0:2, 40 * 30, TRUE), 40, 30); y[sample(length(y), 300)] <- NA
  r <- run_sparse_chain(y, 200)
  expect_true(r$ordered)
  expect_rel(var(r$pooled), 1, 0.05)
})

test_that("count and continuous columns keep the ordering and are N(0, 1)", {
  set.seed(3)
  y <- cbind(matrix(rpois(40 * 10, 3), 40, 10), matrix(rnorm(40 * 10), 40, 10))
  r <- run_sparse_chain(y, 200)
  expect_true(r$ordered)
  expect_rel(var(r$pooled), 1, 0.05)
})

test_that("unconstrained cells are N(lambda_j' omega_i, 1)", {
  # single-level columns carry no ordering, so x - mu must be exactly standard normal
  set.seed(4)
  n <- 500; p <- 4
  y <- matrix(1, n, p); y[sample(n * p, 200)] <- NA
  d <- prepare_data(y)
  lam <- matrix(c(1, -2, 0.5, 0, 0.3, 1, -1, 2), p, 2)
  om <- matrix(rnorm(n * 2), n, 2)
  mu <- om %*% t(lam)
  x <- .sample_x_sparse(d, rep(0, length(d$value)), lam, om, 1, 0)
  r <- (cells_to_matrix(d, x) - mu)[!is.na(y)]
  expect_lt(abs(mean(r)), 4.5 / sqrt(length(r)))
  expect_rel(var(r), 1, 0.08)
})

test_that("the same seed and iteration give identical draws", {
  y <- matrix(rbinom(200, 1, 0.5), 20, 10); d <- prepare_data(y); x0 <- rnorm(length(d$value))
  a <- .sample_x_sparse(d, x0, matrix(0, 10, 1), matrix(0, 20, 1), 7, 3)
  expect_identical(a, .sample_x_sparse(d, x0, matrix(0, 10, 1), matrix(0, 20, 1), 7, 3))
})

test_that("a continuous column of 5000 rows costs O(n) per sweep", {
  y <- matrix(rnorm(5000 * 20), 5000, 20); d <- prepare_data(y)
  x <- rnorm(length(d$value)); lam <- matrix(0, 20, 1); om <- matrix(0, 5000, 1)
  t0 <- proc.time()[["elapsed"]]
  for (t in 1:10) x <- .sample_x_sparse(d, x, lam, om, 1, t)
  expect_lt((proc.time()[["elapsed"]] - t0) / 10, 0.05)
})


test_that("long-format input (row, col, value) gives the same layout as a matrix", {
  set.seed(5)
  y <- matrix(rpois(30 * 8, 2), 30, 8); y[sample(length(y), 150)] <- NA
  obs <- which(!is.na(y), arr.ind = TRUE)
  long <- data.frame(row = obs[, 1], col = obs[, 2], value = y[obs])[sample(nrow(obs)), ]
  expect_identical(prepare_data(long, n = 30, p = 8), prepare_data(y))
})

test_that("long-format input infers n and p, and rejects duplicates and bad indices", {
  long <- data.frame(row = c(1, 3, 2), col = c(2, 1, 2), value = c(1, 0, 1))
  d <- prepare_data(long)
  expect_equal(c(d$n, d$p), c(3L, 2L))
  expect_error(prepare_data(rbind(long, long[1, ])), "duplicate")
  expect_error(prepare_data(data.frame(row = 0, col = 1, value = 1)), "positive")
})
