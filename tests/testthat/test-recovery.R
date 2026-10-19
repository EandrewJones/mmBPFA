# End-to-end recovery on data simulated from the paper's model with known
# loadings. The copula correlation matrix C = D^-1/2 (L L' + I) D^-1/2 is
# invariant to rotation, label switching and the scale of the latent data, so
# it is the identified target to compare against.

copula_cor <- function(L) {
  C <- tcrossprod(L) + diag(nrow(L))
  D <- 1 / sqrt(diag(C))
  C * outer(D, D)
}

test_that("the sampler recovers the latent copula correlations of binary data", {
  skip_on_cran()
  set.seed(42)
  n <- 200; p <- 60; K <- 2
  L <- matrix(rnorm(p * K), p, K) * matrix(rbinom(p * K, 1, 0.7), p, K)
  O <- matrix(rnorm(n * K), n, K)
  a <- rnorm(p, 0, 0.5)
  y <- ((O %*% t(L) - matrix(a, n, p, byrow = TRUE) + matrix(rnorm(n * p), n, p)) > 0) * 1

  utils::capture.output(res <- suppressWarnings(mmBPFA_sampler(
    y, mode = "fixed", K_start = 10, pre_sparse = 50, warmup = 300, iter = 300, thin = 3, seed = 1)))

  Cs <- lapply(res[["MCMC Results"]]$lambda, copula_cor)
  C_est <- Reduce(`+`, Cs) / length(Cs)
  C_true <- copula_cor(L)
  off <- upper.tri(C_true)
  # The original sampler, which shrank Lambda every iteration, scores 0.14 here.
  expect_lt(mean(abs(C_est[off] - C_true[off])), 0.08)
  expect_gt(cor(C_est[off], C_true[off]), 0.95)

  # posterior predictive fit beats the base-rate model
  utils::capture.output(out <- process_results(res, SE = FALSE))
  acc <- check_accuracy(y, out, confusion = FALSE)
  base <- exp(mean(log(ifelse(y == 1, colMeans(y)[col(y)], 1 - colMeans(y)[col(y)]))))
  expect_gt(acc$GMP, base + 0.1)
})
