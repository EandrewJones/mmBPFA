# End-to-end recovery on data simulated from the model with known loadings.
# The copula correlation C = D^-1/2 (L L' + I) D^-1/2 is invariant to rotation,
# label switching and the scale of the latent data, so it is the identified
# target. Hard checks: chains from different starts agree on K, K is recovered
# within 1, and C is recovered. Gaps to the IBP baseline are reported (soft).

copula_cor <- function(L) {
  C <- tcrossprod(L) + diag(nrow(L))
  D <- 1 / sqrt(diag(C))
  C * outer(D, D)
}

sim_binary <- function(n, p, K, seed) {   # same generator as the spike's sim_data()
  set.seed(seed)
  Z <- matrix(rbinom(p * K, 1, 0.4), p, K)
  for (k in seq_len(K)) if (sum(Z[, k]) < 3) Z[sample(p, 3), k] <- 1
  L <- Z * matrix(rnorm(p * K), p, K)
  O <- matrix(rnorm(n * K), n, K); a <- rnorm(p, 0, 0.5)
  x <- O %*% t(L) - matrix(a, n, p, byrow = TRUE) + matrix(rnorm(n * p), n, p)
  y <- (x > 0) * 1
  y[sample(n * p, round(0.05 * n * p))] <- NA
  list(y = y, L = L)
}

fit_summary <- function(y, K_start, seed) {
  utils::capture.output(res <- mmBPFA_sampler(y, H = 30, K_start = K_start, warmup = 1000,
                                              iter = 1000, thin = 1, seed = seed))
  ch <- res[[1]]
  Cs <- Reduce(`+`, lapply(ch$lambda, copula_cor)) / length(ch$lambda)
  list(K = stats::median(ch$K), C = Cs)
}

for (K in c(2, 5)) {
  test_that(sprintf("starts at 2 and 30 agree on K and recover C (K = %d, P = 60, binary)", K), {
    skip_on_cran()
    d <- sim_binary(200, 60, K, seed = 100 + K + 60)
    C_true <- copula_cor(d$L); off <- upper.tri(C_true)
    fits <- lapply(c(2, 30), function(ks) fit_summary(d$y, ks, seed = ks))
    Ks <- sapply(fits, `[[`, "K")
    C_err <- sapply(fits, function(f) mean(abs(f$C[off] - C_true[off])))

    expect_lte(diff(range(Ks)), 1)
    expect_true(all(abs(Ks - K) <= 1))
    expect_true(all(C_err < 0.08))

    base <- readRDS(test_path("fixtures", "baseline_ibp_grid.rds"))
    base <- base[base$K == K & base$p == 60 & base$margins == "binary", ]
    gap <- mean(C_err) - mean(base$C_err)
    if (gap > 0.005) warning(sprintf("soft criterion: C error %.3f vs IBP baseline %.3f (K = %d)",
                                     mean(C_err), mean(base$C_err), K))
  })
}
