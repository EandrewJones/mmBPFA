# Group moves for mixing (Liu & Wu 1999; Liu & Sabatti 2000). Each move draws a
# scale s from p(s) ∝ pi(g_s(state)) |J_s| / s, which leaves the posterior invariant.
#
# Item scale:   (x_j, lambda_j) -> s (x_j, lambda_j); the ranks of x_j are unchanged and
#   s^2 ~ Gamma((n_j + m_j) / 2, A_j / 2), with
#   A_j = sum_obs (x_ij - lambda_j' omega_i)^2 + sum_k z_jk lambda_jk^2 / theta_k.
# Column scale: (lambda_k, omega_k) -> (s lambda_k, omega_k / s); Lambda Omega' is unchanged
#   and s^2 ~ GIG((m_k - N) / 2, chi = sum omega_ik^2, psi = sum lambda_jk^2 / theta_k).

move_fixture <- function(n = 60, p = 6, H = 3, miss = 0.3, seed = 1) {
  set.seed(seed)
  y <- matrix(rbinom(n * p, 3, 0.4), n, p)
  y[matrix(runif(n * p) < miss, n, p)] <- NA
  d <- prepare_data(y)
  z <- matrix(rbinom(p * H, 1, 0.6), p, H)
  list(d = d, x = initial_x(d), z = z, lambda = z * matrix(rnorm(p * H), p, H),
       omega = matrix(rnorm(n * H), n, H), theta = c(1.5, 0.7, 0.05))
}

test_that("item scale: s^2 is Gamma((n_j + m_j) / 2, A_j / 2)", {
  f <- move_fixture()
  xm <- cells_to_matrix(f$d, f$x)
  mu <- f$omega %*% t(f$lambda)
  n_j <- colSums(!is.na(xm))
  m_j <- rowSums(f$z)
  A <- colSums((xm - mu)^2, na.rm = TRUE) + rowSums(f$z * f$lambda^2 / rep(f$theta, each = 6))
  shape <- (n_j + m_j) / 2; rate <- A / 2
  s2 <- draws(4000, local({ t <- 0; function() { t <<- t + 1
    .item_scale_move(f$d, f$x, f$lambda, f$z, f$omega, f$theta, 7, t)$s^2 } }))
  expect_z(rowMeans(s2), shape / rate, sqrt(shape / rate^2 / 4000))
  expect_rel(apply(s2, 1, var), shape / rate^2, 0.1)
})

test_that("item scale multiplies x_j and lambda_j by the same s and keeps ranks", {
  f <- move_fixture(seed = 2)
  out <- .item_scale_move(f$d, f$x, f$lambda, f$z, f$omega, f$theta, 3, 1)
  s_cell <- out$s[f$d$col_idx + 1]
  expect_equal(out$x, f$x * s_cell)
  expect_equal(out$lambda, f$lambda * out$s)
  expect_true(all(out$s > 0))
})

# moments of GIG(p, chi, psi) by quadrature on log u
gig_moments <- function(p, chi, psi) {
  logf <- function(t) p * t - 0.5 * (psi * exp(t) + chi * exp(-t))
  mode <- optimize(logf, c(-30, 30), maximum = TRUE)$maximum
  f <- function(t, k) exp(k * t + logf(t) - logf(mode))
  m <- sapply(0:2, function(k) integrate(f, mode - 15, mode + 15, k = k, rel.tol = 1e-10)$value)
  c(mean = m[2] / m[1], var = m[3] / m[1] - (m[2] / m[1])^2)
}

test_that("column scale: s^2 is GIG((m_k - N) / 2, sum omega^2, sum lambda^2 / theta)", {
  f <- move_fixture(n = 40)
  m_k <- colSums(f$z)
  mom <- sapply(1:3, function(k) gig_moments((m_k[k] - 40) / 2, sum(f$omega[, k]^2),
                                             sum(f$lambda[, k]^2) / f$theta[k]))
  s2 <- draws(4000, local({ t <- 0; function() { t <<- t + 1
    .column_scale_move(f$lambda, f$z, f$omega, f$theta, 11, t)$s^2 } }))
  expect_z(rowMeans(s2), mom["mean", ], sqrt(mom["var", ] / 4000))
  expect_rel(apply(s2, 1, var), mom["var", ], 0.1)
})

test_that("column scale leaves Lambda Omega' unchanged", {
  f <- move_fixture(seed = 3)
  out <- .column_scale_move(f$lambda, f$z, f$omega, f$theta, 5, 1)
  expect_equal(out$omega %*% t(out$lambda), f$omega %*% t(f$lambda))
  expect_equal(out$lambda, f$lambda * rep(out$s, each = 6))
})

test_that("a column with no active items rescales omega towards its prior", {
  f <- move_fixture(n = 50, seed = 4)
  f$z[, 2] <- 0L; f$lambda[, 2] <- 0
  f$omega[, 2] <- f$omega[, 2] * 5                     # far too spread for N(0, 1)
  sd2 <- sapply(1:500, function(t) sd(.column_scale_move(f$lambda, f$z, f$omega, f$theta, 9, t)$omega[, 2]))
  expect_lt(abs(mean(sd2) - 1), 0.15)
})

# Sparse binary data give the item scales little restoring force, so a move that is
# not invariant (for example one with a wrong Jacobian) drifts the scales and the copula
# correlations away from the plain sampler.
invariance_data <- function(seed) {
  set.seed(seed); n <- 150; p <- 12
  y <- (outer(rnorm(n), rnorm(p)) + matrix(rnorm(n * p), n, p) > 0) * 1
  list(dense = y, sparse = replace(y, matrix(runif(n * p) < 0.8, n, p), NA))
}

for (case in c("dense", "sparse")) {
  test_that(paste("the chain with moves targets the same posterior as without:", case, "data"), {
    skip_on_cran()
    y <- invariance_data(5)[[case]]
    fit <- function(moves, seed) {
      utils::capture.output(r <- mmBPFA_sampler(y, H = 4, K_start = 2, warmup = 500, iter = 8000, thin = 2,
                                                chains = 2, seed = seed, moves = moves))
      cc <- lapply(r, function(ch) sapply(ch$lambda, function(L) {
        C <- tcrossprod(L) + diag(12); D <- 1 / sqrt(diag(C)); (C * outer(D, D))[upper.tri(C)] }))
      list(mean = rowMeans(do.call(cbind, cc)),
           # batch-means standard error over both chains
           se = apply(do.call(cbind, lapply(cc, function(m) sapply(split(seq_len(ncol(m)), rep(1:20, each = ncol(m) / 20)),
                                                                     function(b) rowMeans(m[, b, drop = FALSE])))), 1, sd) / sqrt(40))
    }
    a <- fit(character(), 1); b <- fit(c("item_scale", "column_scale"), 2)
    expect_z(a$mean, b$mean, sqrt(a$se^2 + b$se^2), max_z = 5)
  })
}
