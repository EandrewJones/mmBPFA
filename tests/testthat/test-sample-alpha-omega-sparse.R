# Observed-cells-only conditionals:
#   alpha_j | . ~ N(mean_{i in obs(j)} (lambda_j' omega_i - x_ij), 1 / n_j)
#   omega_i | . ~ N(S_i sum_{j in obs(i)} lambda_j (x_ij + alpha_j), S_i),
#     S_i = (sum_{j in obs(i)} lambda_j lambda_j' + I)^-1

sparse_fixture <- function(n = 60, p = 12, K = 2, miss = 0.4, empty_row = TRUE, seed = 1) {
  set.seed(seed)
  lam <- matrix(rnorm(p * K, 0, 0.7), p, K); om <- matrix(rnorm(n * K, 1, 0.5), n, K)
  al <- rnorm(p)
  xm <- om %*% t(lam) - matrix(al, n, p, byrow = TRUE) + matrix(rnorm(n * p), n, p)
  y <- xm; y[matrix(runif(n * p) < miss, n, p)] <- NA
  if (empty_row) y[n, ] <- NA
  d <- prepare_data(y)
  list(d = d, x = xm[cbind(d$row_idx + 1, d$col_idx + 1)], xm = xm, y = y,
       lam = lam, om = om, al = al)
}

test_that("alpha draws match the observed-only conditional, with variance 1 / n_j", {
  f <- sparse_fixture()
  xm <- cells_to_matrix(f$d, f$x)
  n_j <- colSums(!is.na(f$y))
  post_mean <- colMeans(f$om %*% t(f$lam) - xm, na.rm = TRUE)
  a <- sapply(1:4000, function(t) .sample_alpha_sparse(f$d, f$x, f$lam, f$om, 1, t))
  expect_z(rowMeans(a), post_mean, sqrt(1 / n_j / 4000))
  expect_rel(apply(a, 1, var), 1 / n_j, 0.1)
})

test_that("omega draws match the observed-only conditional for each row", {
  f <- sparse_fixture()
  xm <- cells_to_matrix(f$d, f$x)
  K <- ncol(f$lam)
  draws <- simplify2array(lapply(1:3000, function(t) .sample_omega_sparse(f$d, f$x, f$lam, f$al, 2, t)))
  for (i in c(1, 7, 30)) {
    o <- which(!is.na(f$y[i, ]))
    S <- solve(crossprod(f$lam[o, , drop = FALSE]) + diag(K))
    m <- S %*% crossprod(f$lam[o, , drop = FALSE], xm[i, o] + f$al[o])
    expect_z(rowMeans(draws[i, , ]), as.vector(m), sqrt(diag(S) / 3000))
    expect_rel(diag(cov(t(draws[i, , ]))), diag(S), 0.1)
  }
})

test_that("a row with no observations draws omega from its N(0, I) prior", {
  f <- sparse_fixture()
  n <- nrow(f$y)
  draws <- sapply(1:4000, function(t) .sample_omega_sparse(f$d, f$x, f$lam, f$al, 3, t)[n, ])
  expect_z(rowMeans(draws), c(0, 0), sqrt(1 / 4000))
  expect_rel(apply(draws, 1, var), c(1, 1), 0.1)
})

test_that("an item with no observations gets alpha = 0 (flat prior, nothing to condition on)", {
  f <- sparse_fixture()
  y <- f$y; y[, 3] <- NA
  d <- prepare_data(y); x <- f$xm[cbind(d$row_idx + 1, d$col_idx + 1)]
  a <- .sample_alpha_sparse(d, x, f$lam, f$om, 4, 1)
  expect_true(all(is.finite(a)))
  expect_equal(a[3], 0)
})

test_that("omega draws match the conditional when rows are fully observed or nearly so", {
  for (miss in c(0, 0.03)) {
    f <- sparse_fixture(n = 40, p = 12, miss = miss, empty_row = FALSE, seed = 5)
    xm <- cells_to_matrix(f$d, f$x); K <- ncol(f$lam)
    draws <- simplify2array(lapply(1:3000, function(t) .sample_omega_sparse(f$d, f$x, f$lam, f$al, 6, t)))
    for (i in c(1, 20, 40)) {
      o <- which(!is.na(f$y[i, ]))
      S <- solve(crossprod(f$lam[o, , drop = FALSE]) + diag(K))
      m <- S %*% crossprod(f$lam[o, , drop = FALSE], xm[i, o] + f$al[o])
      expect_z(rowMeans(draws[i, , ]), as.vector(m), sqrt(diag(S) / 3000))
      expect_rel(diag(cov(t(draws[i, , ]))), diag(S), 0.1)
    }
  }
})

test_that("omega costs O(N P H) on dense data and O(nnz H^2) on sparse data", {
  skip_on_cran()
  time_omega <- function(miss) {
    set.seed(1); n <- 2000; p <- 200; H <- 30
    y <- matrix(rbinom(n * p, 1, 0.5), n, p); y[matrix(runif(n * p) < miss, n, p)] <- NA
    d <- prepare_data(y); x <- initial_x(d)
    lam <- matrix(rnorm(p * H, 0, 0.3), p, H); al <- rnorm(p)
    t0 <- proc.time()[["elapsed"]]
    for (r in 1:5) .sample_omega_sparse(d, x, lam, al, 1, r)
    (proc.time()[["elapsed"]] - t0) / 5
  }
  expect_lt(time_omega(0), 0.06)     # was 0.40 s
  expect_lt(time_omega(0.9), 0.03)   # was 0.077 s
})

test_that("omega draws match the conditional for rows with fewer observed items than dimensions", {
  set.seed(7)
  n <- 30; p <- 40; K <- 10
  lam <- matrix(rnorm(p * K, 0, 0.6), p, K); al <- rnorm(p)
  xm <- matrix(rnorm(n * p), n, p)
  y <- matrix(NA_real_, n, p)
  for (i in 1:n) { o <- sample(p, 1 + (i %% 4)); y[i, o] <- xm[i, o] }   # 1 to 4 observed items per row
  d <- prepare_data(y); x <- xm[cbind(d$row_idx + 1, d$col_idx + 1)]
  draws <- simplify2array(lapply(1:4000, function(t) .sample_omega_sparse(d, x, lam, al, 8, t)))
  for (i in c(1, 2, 3, 4)) {
    o <- which(!is.na(y[i, ]))
    S <- solve(crossprod(lam[o, , drop = FALSE]) + diag(K))
    m <- S %*% crossprod(lam[o, , drop = FALSE], xm[i, o] + al[o])
    expect_z(rowMeans(draws[i, , ]), as.vector(m), sqrt(diag(S) / 4000))
    expect_rel(diag(cov(t(draws[i, , ]))), diag(S), 0.1)
    expect_z(cov(t(draws[i, , ]))[1, 2], S[1, 2], sqrt((S[1, 1] * S[2, 2] + S[1, 2]^2) / 4000))
  }
})

test_that("very sparse rows cost O(n_i^2 H), not O(H^3)", {
  skip_on_cran()
  set.seed(1); n <- 20000; p <- 1000; H <- 40
  long <- data.frame(row = rep(1:n, each = 3), col = as.vector(replicate(n, sample(p, 3))), value = rbinom(3 * n, 1, 0.5))
  d <- prepare_data(long, n = n, p = p); x <- initial_x(d)
  lam <- matrix(rnorm(p * H, 0, 0.3), p, H); al <- rnorm(p)
  t0 <- proc.time()[["elapsed"]]
  for (r in 1:3) .sample_omega_sparse(d, x, lam, al, 1, r)
  expect_lt((proc.time()[["elapsed"]] - t0) / 3, 0.05)
})
