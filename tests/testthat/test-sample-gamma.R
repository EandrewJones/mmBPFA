# gamma_k | Z, Lambda ~ Gamma(shape = c + m_k / 2, rate = d + sum_j z_jk lambda_jk^2 / 2)

test_that("sample_gamma_k draws from the analytic conditional posterior", {
  zeros <- cbind(rep(1, 50), c(rep(1, 10), rep(0, 40)))
  lambda <- cbind(rep(1, 50), c(rep(2, 10), rep(0, 40)))
  # k = 1: shape 1 + 25 = 26, rate 2 + 25 = 27; k = 2: shape 1 + 5 = 6, rate 2 + 20 = 22
  set.seed(1)
  g <- draws(20000, function() as.vector(sample_gamma_k(zeros, lambda, c = 1, d = 2)))
  post_var <- c(26 / 27^2, 6 / 22^2)
  expect_z(rowMeans(g), c(26 / 27, 6 / 22), sqrt(post_var / 20000))
  # relative SE of a Gamma(6) sample variance is ~ sqrt((2 + 1) / 20000) ~ 0.012
  expect_rel(apply(g, 1, var), post_var, 0.06)
})

test_that("the default precision prior is proper", {
  set.seed(1)
  g <- sample_gamma_k(matrix(0, 10, 1), matrix(0, 10, 1))
  expect_gt(g, 0)
})
