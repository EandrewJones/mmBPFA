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
