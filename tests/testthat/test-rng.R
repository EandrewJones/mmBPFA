# Distribution checks for the C++ RNG used by the sampler kernels. Expected
# values come from closed-form moments computed in base R.

tn_moments <- function(mu, lo, hi) {
  a <- lo - mu; b <- hi - mu
  Z <- pnorm(b) - pnorm(a)
  da <- dnorm(a); db <- dnorm(b)
  ta <- if (is.finite(a)) a * da else 0; tb <- if (is.finite(b)) b * db else 0
  m <- (da - db) / Z
  c(mean = mu + m, var = 1 + (ta - tb) / Z - m^2)
}

test_that("normal draws have mean 0 and variance 1", {
  d <- .rng_draws("norm", 2e5, 1, numeric(0))
  expect_z(mean(d), 0, sqrt(1 / 2e5))
  expect_rel(var(d), 1, 0.02)
})

test_that("uniform draws lie in (0, 1) with mean 1/2", {
  d <- .rng_draws("unif", 2e5, 2, numeric(0))
  expect_true(all(d > 0 & d < 1))
  expect_z(mean(d), 0.5, sqrt(1 / 12 / 2e5))
})

test_that("gamma draws match their moments, including shape below 1", {
  for (shape in c(2.5, 0.3)) {
    d <- .rng_draws("gamma", 2e5, 3, shape)
    expect_z(mean(d), shape, sqrt(shape / 2e5))
    expect_rel(var(d), shape, 0.05)
  }
})

test_that("beta draws match their moments", {
  a <- 0.5; b <- 5
  d <- .rng_draws("beta", 2e5, 4, c(a, b))
  v <- a * b / ((a + b)^2 * (a + b + 1))
  expect_z(mean(d), a / (a + b), sqrt(v / 2e5))
  expect_rel(var(d), v, 0.05)
})

test_that("truncated normal draws match their moments across regimes", {
  cases <- list(c(0, 1, 2), c(0, -Inf, 0.5), c(1, -1, Inf), c(0, 3, 3.5))
  for (cs in cases) {
    d <- .rng_draws("tnorm", 1e5, 5, cs)
    m <- tn_moments(cs[1], cs[2], cs[3])
    expect_true(all(d >= cs[2] & d <= cs[3]))
    expect_z(mean(d), m[["mean"]], sqrt(m[["var"]] / 1e5))
    expect_rel(var(d), m[["var"]], 0.05)
  }
})

test_that("truncated normal is stable far in the tails", {
  up <- .rng_draws("tnorm", 1e5, 6, c(0, 8, Inf))
  expect_true(all(is.finite(up) & up >= 8))
  m <- dnorm(8) / pnorm(8, lower.tail = FALSE)
  expect_z(mean(up), m, sqrt(tn_moments(0, 8, Inf)[["var"]] / 1e5))
  down <- .rng_draws("tnorm", 1e4, 7, c(0, -Inf, -10))
  expect_true(all(is.finite(down) & down <= -10))
  narrow <- .rng_draws("tnorm", 1e4, 8, c(0, 8, 8.01))
  expect_true(all(narrow >= 8 & narrow <= 8.01))
})

test_that("truncated normal with an empty interval returns its endpoint", {
  expect_equal(.rng_draws("tnorm", 3, 9, c(0, 2, 2)), rep(2, 3))
})

test_that("streams are reproducible and distinct", {
  a <- .rng_stream(42, 3, 7, 5)
  expect_identical(a, .rng_stream(42, 3, 7, 5))
  expect_false(identical(a, .rng_stream(42, 3, 8, 5)))
  expect_false(identical(a, .rng_stream(42, 4, 7, 5)))
})
