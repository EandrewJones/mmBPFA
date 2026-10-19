# Shared fixtures: data drawn from the paper's generative model with known
# parameters, x = Omega Lambda' - 1 alpha' + eps, eps ~ N(0, 1).
sim_latent <- function(n = 300, p = 40, K = 2, seed = 1) {
  set.seed(seed)
  lambda <- matrix(rnorm(p * K), p, K)
  omega <- matrix(rnorm(n * K), n, K)
  alpha <- rnorm(p, 0, 0.5)
  mu <- omega %*% t(lambda) - matrix(alpha, n, p, byrow = TRUE)
  list(x = mu + matrix(rnorm(n * p), n, p), mu = mu,
       lambda = lambda, omega = omega, alpha = alpha)
}

# Monte Carlo draws of a sampler that returns a vector; one column per draw.
draws <- function(n_draws, f) vapply(seq_len(n_draws), function(i) f(), f())

# Relative-error expectation. testthat's expect_equal() silently switches to an
# absolute comparison when |expected| < tolerance, which makes small posterior
# variances and probabilities untestable with it.
expect_rel <- function(actual, expected, rel) {
  err <- max(abs(actual - expected) / abs(expected))
  expect(err < rel, sprintf(
    "max relative error %.3g >= %.3g\n  actual: %s\nexpected: %s", err, rel,
    paste(signif(actual, 4), collapse = " "), paste(signif(expected, 4), collapse = " ")))
  invisible(actual)
}

# Monte Carlo z-test: every estimate within max_z standard errors of its target.
expect_z <- function(estimate, expected, se, max_z = 4.5) {
  z <- max(abs(estimate - expected) / se)
  expect(z < max_z, sprintf(
    "max |z| %.2f >= %.2f\n  estimate: %s\n  expected: %s", z, max_z,
    paste(signif(estimate, 4), collapse = " "), paste(signif(expected, 4), collapse = " ")))
  invisible(estimate)
}
