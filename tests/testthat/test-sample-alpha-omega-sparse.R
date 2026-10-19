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
