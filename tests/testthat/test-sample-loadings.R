# Loadings prior: CUSP column shrinkage (Legramanti, Durante & Dunson 2020) with
# item-level inclusion (Schiavon, Canale & Dunson 2022). Every conditional is checked
# against a closed form computed here in base R.

default_prior <- function(spike_var = 0.05^2) {
  list(slab_shape = 1, slab_rate = 0.5, spike_var = spike_var, a_cusp_shape = 1,
       a_cusp_rate = 1, rho_a = 1, rho_b = 1)
}

new_state <- function(p, H, on = H, rho = 0.5, theta = 1) {
  list(lambda = matrix(0, p, H), z = matrix(0L, p, H), rho = rep(rho, H),
       c = as.integer(ifelse(seq_len(H) <= on, H, 1)), v = c(rep(0.3, H - 1), 1),
       a_cusp = 2, theta = rep(theta, H))
}

lz_fixture <- function(n = 150, p = 8, miss = 0.3, seed = 1) {
  set.seed(seed)
  om <- matrix(rnorm(n), n, 1)
  y <- om %*% t(rnorm(p)) + matrix(rnorm(n * p), n, p)
  xm <- y
  y[matrix(runif(n * p) < miss, n, p)] <- NA
  d <- prepare_data(y)
  list(d = d, x = xm[cbind(d$row_idx + 1, d$col_idx + 1)], om = om, y = y)
}

test_that("lambda given z = 1 is Gaussian, with and without missing cells", {
  for (miss in c(0.3, 0)) {
  f <- lz_fixture(miss = miss)
  st <- new_state(8, 1, rho = 1, theta = 2)
  obs <- !is.na(f$y)
  xm <- cells_to_matrix(f$d, f$x)
  prec <- colSums((obs * as.vector(f$om))^2) + 1 / 2
  mu <- colSums(obs * as.vector(f$om) * ifelse(obs, xm, 0)) / prec
  l <- sapply(1:4000, function(t) .sample_loadings(f$d, f$x, f$om, rep(0, 8), st, default_prior(), 1, t)$lambda[, 1])
  expect_z(rowMeans(l), mu, sqrt(1 / prec / 4000))
  expect_rel(apply(l, 1, var), 1 / prec, 0.1)
  }
})

test_that("z_jk is Bernoulli with the collapsed log-odds, with and without missing cells", {
  for (miss in c(0.3, 0)) {
  f <- lz_fixture(n = 40, p = 5, miss = miss, seed = 2)
  st <- new_state(5, 1, rho = 0.3, theta = 1)
  obs <- !is.na(f$y)
  xm <- cells_to_matrix(f$d, f$x)
  prec <- colSums((obs * as.vector(f$om))^2) + 1
  mu <- colSums(obs * as.vector(f$om) * ifelse(obs, xm, 0)) / prec
  lo <- qlogis(0.3) + 0.5 * (log(1) - log(prec)) + 0.5 * prec * mu^2
  want <- plogis(lo)
  z <- sapply(1:6000, function(t) .sample_loadings(f$d, f$x, f$om, rep(0, 5), st, default_prior(), 2, t)$z[, 1])
  expect_z(rowMeans(z), want, sqrt(pmax(want * (1 - want), 1e-4) / 6000))
  }
})

test_that("column indicators follow the CUSP categorical given the loadings", {
  pr <- default_prior()
  lam <- cbind(c(1, -0.8, 0.6, 0), c(0.02, -0.03, 0, 0.01), c(0, 0, 0, 0))
  z <- (lam != 0) * 1L
  v <- c(0.5, 0.4, 1)
  w <- v * c(1, cumprod(1 - v)[-3])
  lmvt <- function(l) { m <- length(l); a <- pr$slab_shape; b <- pr$slab_rate
    lgamma(a + m / 2) - lgamma(a) - m / 2 * log(2 * pi * b) - (a + m / 2) * log1p(sum(l^2) / (2 * b)) }
  want <- t(sapply(1:3, function(h) {
    l <- lam[z[, h] == 1, h]
    lp <- log(w) + ifelse(1:3 <= h, sum(dnorm(l, 0, sqrt(pr$spike_var), log = TRUE)), lmvt(l))
    exp(lp - max(lp)) / sum(exp(lp - max(lp)))
  }))
  draws <- sapply(1:20000, function(t) .sample_column_indicators(lam, z, v, pr, 3, t))
  for (h in 1:3) {
    freq <- tabulate(draws[h, ], 3) / 20000
    expect_z(freq, want[h, ], sqrt(pmax(want[h, ] * (1 - want[h, ]), 1e-4) / 20000))
  }
})

test_that("stick weights and the CUSP concentration follow their conjugate updates", {
  pr <- default_prior()
  cc <- c(4L, 4L, 1L, 2L, 1L)            # H = 5
  H <- 5
  s <- lapply(1:20000, function(t) .sample_sticks(cc, 2, pr, 4, t))
  v1 <- sapply(s, function(r) r$v[1]); v2 <- sapply(s, function(r) r$v[2])
  # v_l ~ Beta(1 + #{c = l}, a + #{c > l}) with a = 2
  for (l in 1:2) {
    a1 <- 1 + sum(cc == l); b1 <- 2 + sum(cc > l)
    vv <- if (l == 1) v1 else v2
    m <- a1 / (a1 + b1); var_ <- a1 * b1 / ((a1 + b1)^2 * (a1 + b1 + 1))
    expect_z(mean(vv), m, sqrt(var_ / 20000))
  }
  expect_true(all(sapply(s, function(r) r$v[H]) == 1))
})

test_that("a_cusp is drawn from p(a | c) with the sticks integrated out", {
  # p(a | c) is proportional to Gamma(a; shape, rate) * prod_{l < H} B(1 + n_l, a + g_l) / B(1, a),
  # n_l = #{c_h = l}, g_l = #{c_h > l}. Conditioning on the sticks instead traps
  # the chain at a = 0 once every column is off.
  pr <- default_prior()
  post <- function(a, cc) {
    H <- length(cc)
    sapply(a, function(ai) exp(dgamma(ai, pr$a_cusp_shape, pr$a_cusp_rate, log = TRUE) +
      sum(sapply(seq_len(H - 1), function(l) lbeta(1 + sum(cc == l), ai + sum(cc > l)) - lbeta(1, ai)))))
  }
  for (cc in list(c(4L, 4L, 1L, 2L, 1L), rep(1L, 8))) {
    Z <- integrate(post, 0, Inf, cc = cc)$value
    m <- integrate(function(a) a * post(a, cc), 0, Inf)$value / Z
    v2 <- integrate(function(a) a^2 * post(a, cc), 0, Inf)$value / Z - m^2
    # slice-sampling chain on a with c fixed; draws are autocorrelated, so the
    # standard error is inflated by 2 (one-dimensional slice sampling mixes fast)
    a <- numeric(20000); cur <- 1
    for (t in 1:20000) { cur <- .sample_a_cusp(cc, cur, pr, 5, t); a[t] <- cur }
    expect_true(all(a > 0))
    expect_z(mean(a), m, 2 * sqrt(v2 / 20000))
  }
})

test_that("slab variances are InvGamma on active loadings; off columns use the spike", {
  pr <- default_prior()
  lam <- cbind(c(1, 1, 0, 2), c(0.5, 0, 0, 0))
  z <- (lam != 0) * 1L
  cc <- c(2L, 1L)                         # column 1 on (c > 1), column 2 off
  th <- sapply(1:20000, function(t) .sample_theta(lam, z, cc, pr, 6, t))
  shape <- 1 + 3 / 2; rate <- 0.5 + (1 + 1 + 4) / 2
  expect_z(mean(1 / th[1, ]), shape / rate, sqrt(shape / rate^2 / 20000))
  expect_true(all(th[2, ] == pr$spike_var))
})

test_that("with no data the sampler leaves the CUSP prior on the number of active columns invariant", {
  H <- 8; p <- 6; n <- 5
  d <- prepare_data(matrix(NA_real_, n, p))
  pr <- default_prior()
  st <- new_state(p, H, on = 3)
  on <- integer(0)
  for (t in 1:30000) {
    st <- .sample_loadings(d, numeric(0), matrix(0, n, H), rep(0, p), st, pr, 7, t)
    if (t %% 10 == 0) on <- c(on, sum(st$c > seq_len(H)))
  }
  set.seed(8)
  prior_on <- replicate(20000, {
    a <- rgamma(1, pr$a_cusp_shape, pr$a_cusp_rate)
    v <- c(rbeta(H - 1, 1, a), 1); w <- v * c(1, cumprod(1 - v)[-H])
    sum(sample.int(H, H, replace = TRUE, prob = w) > seq_len(H))
  })
  expect_rel(mean(on), mean(prior_on), 0.06)
  expect_rel(var(on), var(prior_on), 0.12)
})

test_that("a loadings sweep on fully observed data uses one matrix-vector product per column", {
  skip_on_cran()
  set.seed(1); n <- 2000; p <- 200; H <- 30
  y <- matrix(rbinom(n * p, 1, 0.5), n, p); d <- prepare_data(y); x <- initial_x(d)
  st <- new_state(p, H); st$lambda <- matrix(rnorm(p * H, 0, 0.3), p, H); st$z[] <- 1L
  om <- matrix(rnorm(n * H), n, H)
  t0 <- proc.time()[["elapsed"]]
  for (r in 1:5) .sample_loadings(d, x, om, rep(0, p), st, default_prior(), 1, r)
  # was 0.044 s with per-cell loops; ~0.026 s with the reference BLAS (faster with OpenBLAS/MKL)
  expect_lt((proc.time()[["elapsed"]] - t0) / 5, 0.035)
})
