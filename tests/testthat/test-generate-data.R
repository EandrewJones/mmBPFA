test_that("simulated latent data have unit-variance noise around the linear predictor", {
  sim <- generate_sim_data(n = 300, p = 100, K_true = 3, noise = 0, seed = 1)
  tv <- sim$true_values
  mu <- tv$omega %*% t(tv$lambda) - matrix(tv$alpha, 300, 100, byrow = TRUE)
  expect_rel(var(as.vector(tv$x - mu)), 1, 0.03)
})

test_that("binary simulated data are the sign of the latent data", {
  sim <- generate_sim_data(n = 100, p = 50, K_true = 3, noise = 0, seed = 2)
  expect_true(all((sim$data > 0) == (sim$true_values$x > 0)))
})

test_that("active loadings have variance 1 / gamma_k and are not normalized", {
  sim <- generate_sim_data(n = 10, p = 4000, K_true = 2, sparsity = 0.5, noise = 0, seed = 3)
  tv <- sim$true_values
  v <- vapply(1:2, function(k) var(tv$lambda[tv$zeros[, k] == 1, k]), numeric(1))
  expect_rel(v, 1 / tv$gamma_k, 0.1)
})

test_that("feature probabilities follow the two-parameter IBP Beta(ab/K, b)", {
  # K = 200, a = 20, b = 5: pi_k ~ Beta(0.5, 5), mean 1/11, sd sqrt(m (1 - m) / 6.5)
  sim <- generate_sim_data(n = 10, p = 20, K_true = 200, a = 20, b = 5, noise = 0, seed = 4)
  m <- 0.5 / 5.5
  expect_z(mean(sim$true_values$pi_k), m, sqrt(m * (1 - m) / 6.5 / 200))
})
