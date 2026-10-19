call_lz <- function(x, lambda, zeros, omega, alpha, gamma_k, a = 1, b = 1,
                    tau = 1, sparse = TRUE, infinite = FALSE) {
  sample_lambda_and_zeros(x, lambda, zeros, omega, alpha, gamma_k,
                          rep(0L, ncol(lambda)), a, b, tau, sparse, infinite)
}

test_that("pre-sparse phase draws lambda from its Gaussian conditional with z fixed at 1", {
  s <- sim_latent(n = 200, p = 6, K = 1)
  gam <- 2
  prec <- sum(s$omega^2) + gam
  resid <- s$x + matrix(s$alpha, 200, 6, byrow = TRUE)
  post_mean <- as.vector(crossprod(s$omega, resid)) / prec

  set.seed(1)
  out <- replicate(4000, simplify = FALSE,
    call_lz(s$x, s$lambda, matrix(1, 6, 1), s$omega, s$alpha, gam, sparse = FALSE))
  z <- sapply(out, function(o) o[[1]])
  l <- sapply(out, function(o) o[[2]])

  expect_true(all(z == 1))
  expect_z(rowMeans(l), post_mean, sqrt(1 / prec / 4000))
  expect_rel(apply(l, 1, var), rep(1 / prec, 6), 0.1)
})

test_that("z_jk is drawn with the two-parameter IBP conditional odds", {
  # P = 2 items, K = 1; item 2 is active so m_{-1,1} = 1 when item 1 is updated
  set.seed(4)
  n <- 50; b <- 3; gam <- 1
  omega <- matrix(rnorm(n), n, 1)
  x <- cbind(0.15 * omega[, 1] + rnorm(n), 0.5 * omega[, 1] + rnorm(n))
  e1 <- x[, 1]                                   # alpha = 0, lambda_11 removed
  prec <- sum(omega^2) + gam
  mu <- sum(omega * e1) / prec
  odds <- sqrt(gam / prec) * exp(0.5 * prec * mu^2) * 1 / (2 + b - 1 - 1)
  want <- odds / (1 + odds)

  set.seed(5)
  z1 <- replicate(6000, call_lz(x, matrix(c(0.1, 0.5), 2, 1), matrix(1, 2, 1),
                                omega, c(0, 0), gam, b = b)[[1]][1, 1])
  expect_z(mean(z1), want, sqrt(want * (1 - want) / 6000))
})

test_that("tiny activation probabilities are not inflated by truncation", {
  # Omega = 0 makes the likelihood ratio 1, so each inactive item joins the
  # single active one with prior odds 1 / (P + b - 2): about one join expected.
  P <- 400000; n <- 4
  zeros <- matrix(0, P, 1); zeros[1] <- 1
  set.seed(6)
  out <- call_lz(matrix(rnorm(n * P), n, P), zeros, zeros, matrix(0, n, 1),
                 rep(0, P), 1, b = 1)
  expect_lte(sum(out[[1]]) - 1, 6)
})

test_that("collapsed likelihood ratio for new features matches hand-derived values", {
  # log v_L = -N/2 log(1 + |l|^2) + |l|^2 / (2 (1 + |l|^2)) * sum(e^2)
  expect_equal(.feature_loglik_ratio(rep(1, 4), 1), -2 * log(2) + 1)
  expect_equal(.feature_loglik_ratio(c(1, 2), c(1, 1)), -log(3) + 5 / 3)
  expect_equal(.feature_loglik_ratio(c(1, 2), numeric(0)), 0)
})

# Run only the IBP moves on a single variable with no observations: the
# likelihood is flat, so the number of features must be distributed Pois(a).
run_prior_chain <- function(iters, a, tau) {
  lambda <- matrix(0, 1, 0); zeros <- matrix(0, 1, 0); k <- integer(iters)
  for (t in seq_len(iters)) {
    K <- ncol(lambda)
    out <- sample_lambda_and_zeros(matrix(0, 0, 1), lambda, zeros, matrix(0, 0, K),
                                   0, rep(1, K), rep(0L, K), a, 1, tau, TRUE, TRUE)
    keep <- colSums(out[[1]]) > 0
    zeros <- out[[1]][, keep, drop = FALSE]
    lambda <- out[[2]][, keep, drop = FALSE]
    k[t] <- ncol(zeros)
    if (k[t] > 100) stop("number of features diverged at iteration ", t)
  }
  k
}

test_that("new-feature moves leave the IBP prior invariant", {
  set.seed(7)
  k <- run_prior_chain(20000, a = 3, tau = 2)
  # draws are autocorrelated, so tolerances allow an effective sample size of a few thousand
  expect_rel(mean(k), 3, 0.03)
  expect_rel(var(k), 3, 0.08)
})

test_that("pure-noise residuals rarely gain features, and new features are singletons", {
  s <- sim_latent(n = 200, p = 60, K = 3)
  set.seed(8)
  added <- replicate(50, {
    out <- call_lz(s$x, s$lambda, matrix(1, 60, 3), s$omega, s$alpha, rep(1, 3),
                   a = 3, tau = 5, infinite = TRUE)
    newz <- out[[1]][, -(1:3), drop = FALSE]
    expect_true(all(colSums(newz) == 1))
    ncol(newz)
  })
  # With no signal the posterior cannot support more singletons than the prior
  # mean P * nu = 60 * 3 / 60 = 3.
  expect_lt(mean(added), 3)
})

test_that("a variable with an idiosyncratic factor gains a new feature", {
  s <- sim_latent(n = 200, p = 60, K = 3)
  set.seed(9)
  s$x[, 1] <- s$x[, 1] + 2 * rnorm(200)          # unmodelled factor on item 1 only
  gained <- replicate(50, {
    out <- call_lz(s$x, s$lambda, matrix(1, 60, 3), s$omega, s$alpha, rep(1, 3),
                   a = 3, tau = 60, infinite = TRUE)
    newz <- out[[1]][, -(1:3), drop = FALSE]
    c(item1 = any(newz[1, ] == 1), others = mean(rowSums(newz[-1, , drop = FALSE]) > 0))
  })
  expect_gt(mean(gained["item1", ]), 0.8)
  expect_lt(mean(gained["others", ]), 0.3)
})
