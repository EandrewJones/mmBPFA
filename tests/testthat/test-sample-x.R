# Under the extended rank likelihood with mean zero, the x's in a column are
# n iid N(0, 1) draws conditioned only on their ranks agreeing with y. The
# rank event does not depend on the order statistics, so the pooled x's are
# marginally N(0, 1), and the ordering must hold after every update.

ordered_ok <- function(x, y) {
  all(vapply(seq_len(ncol(y)), function(j) {
    obs <- !is.na(y[, j])
    lv <- sort(unique(y[obs, j]))
    if (length(lv) < 2) return(TRUE)
    all(vapply(seq_len(length(lv) - 1), function(v)
      max(x[obs & y[, j] == lv[v], j]) < min(x[obs & y[, j] == lv[v + 1], j]), TRUE))
  }, TRUE))
}

run_x_chain <- function(y, mode, iters) {
  n <- nrow(y); p <- ncol(y)
  c(d_mask, n_levels, margin_vals) %<-% calculate_d_mask(y, mode)
  # valid starting point: x = rank-consistent normal scores
  x <- apply(y, 2, function(col) qnorm((rank(col, ties.method = "random", na.last = "keep") - 0.5) / n))
  x[is.na(x)] <- rnorm(sum(is.na(x)))
  pooled <- numeric(0); always_ordered <- TRUE
  for (t in seq_len(iters)) {
    c(x, cuts) %<-% sample_x(y, mode, d_mask, n_levels, margin_vals, x,
                             matrix(0, p, 1), rep(0, p), matrix(0, n, 1))
    always_ordered <- always_ordered && ordered_ok(x, y)
    if (t > iters / 2) pooled <- c(pooled, x)
  }
  list(ordered = always_ordered, pooled = pooled)
}

test_that("binary x draws respect the rank ordering and are N(0, 1) at stationarity", {
  set.seed(1)
  y <- matrix(rbinom(40 * 30, 1, 0.3), 40, 30)
  r <- run_x_chain(y, "fixed", 200)
  expect_true(r$ordered)
  expect_lt(abs(mean(r$pooled)), 0.05)
  expect_rel(var(r$pooled), 1, 0.05)
})

test_that("ordinal x draws respect the rank ordering, with missing data", {
  set.seed(2)
  y <- matrix(sample(0:2, 40 * 30, TRUE), 40, 30)
  y[sample(length(y), 60)] <- NA
  r <- run_x_chain(y, "fixed", 200)
  expect_true(r$ordered)
  expect_rel(var(r$pooled), 1, 0.05)
})

test_that("multi-margin x draws respect the rank ordering", {
  set.seed(3)
  y <- cbind(matrix(rbinom(40 * 10, 1, 0.5), 40, 10),
             matrix(rpois(40 * 10, 3), 40, 10))
  r <- run_x_chain(y, "multi", 200)
  expect_true(r$ordered)
  expect_rel(var(r$pooled), 1, 0.05)
})
