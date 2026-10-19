# alpha_j | X, Lambda, Omega ~ N(mean_i(lambda_j' omega_i - x_ij), 1/N) under a flat prior.

test_that("sample_alpha draws from the analytic conditional posterior", {
  s <- sim_latent(n = 400, p = 5, K = 2)
  post_mean <- colMeans(s$omega %*% t(s$lambda) - s$x)
  post_var <- 1 / 400

  set.seed(2)
  a <- draws(4000, function() as.vector(sample_alpha(s$x, s$lambda, s$omega)))

  expect_z(rowMeans(a), post_mean, sqrt(post_var / 4000))
  # relative SE of a sample variance is sqrt(2 / 3999) ~ 0.022
  expect_rel(apply(a, 1, var), rep(post_var, 5), 0.1)
})

test_that("sample_alpha recovers item intercepts", {
  s <- sim_latent(n = 2000, p = 5, K = 2)
  set.seed(3)
  a <- as.vector(sample_alpha(s$x, s$lambda, s$omega))
  # posterior sd is 1/sqrt(2000) ~ 0.022
  expect_lt(max(abs(a - s$alpha)), 0.12)
})
