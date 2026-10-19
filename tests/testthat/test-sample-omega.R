# omega_i | X, Lambda, alpha ~ N(Sigma Lambda' (x_i + alpha), Sigma), Sigma = (Lambda' Lambda + I)^-1.

test_that("sample_omega draws from the analytic conditional posterior", {
  set.seed(1)
  n <- 30; p <- 20; K <- 2
  lambda <- matrix(rnorm(p * K, 0, 0.4), p, K)
  alpha <- rnorm(p)
  # ideal points centred away from zero, so the posterior means are too
  omega <- matrix(rnorm(n * K, 1, 0.3), n, K)
  x <- omega %*% t(lambda) - matrix(alpha, n, p, byrow = TRUE) + matrix(rnorm(n * p), n, p)

  Sigma <- solve(crossprod(lambda) + diag(K))
  post_mean <- t(Sigma %*% t(lambda) %*% t(x + matrix(alpha, n, p, byrow = TRUE)))

  set.seed(2)
  d <- replicate(3000, sample_omega(lambda, x, alpha))   # n x K x draws

  expect_z(apply(d, 1:2, mean), post_mean, sqrt(matrix(diag(Sigma), n, K, byrow = TRUE) / 3000))
  # SE of a sample covariance entry: sqrt((s_ii s_jj + s_ij^2) / draws)
  cov_se <- sqrt((outer(diag(Sigma), diag(Sigma)) + Sigma^2) / 3000)
  expect_z(cov(t(d[1, , ])), Sigma, cov_se)
  expect_z(cov(t(d[n, , ])), Sigma, cov_se)
})
