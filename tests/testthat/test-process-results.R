test_that("dimensions are matched to the pivot draw by their sparsity patterns", {
  pivot <- cbind(c(1, 1, 1, 1, 0, 0), c(0, 0, 1, 1, 1, 1), c(1, 0, 1, 0, 1, 0))
  permuted <- pivot[, c(3, 1, 2)]
  permuted[1, 3] <- 1 - permuted[1, 3]          # one flipped entry
  partial <- pivot[, c(3, 1)]
  perms <- align_dimensions(list(pivot, permuted, partial), pivot = 1)
  expect_equal(perms, list(c(1L, 2L, 3L), c(2L, 3L, 1L), c(2L, NA, 1L)))
})

# A bpfa.results object whose draws store the same two dimensions in
# alternating column order.
fake_results <- function(L, O) {
  n <- nrow(O); p <- nrow(L)
  Z <- (L != 0) * 1
  flip <- function(m, t) if (t %% 2 == 0) m[, 2:1] else m
  draws <- 1:6
  cuts <- list(a = matrix(c(-Inf, 0)[1 + (seq_len(n * p) %% 2)], n, p),
               b = matrix(c(0, Inf)[1 + (seq_len(n * p) %% 2)], n, p))
  store <- tibble::tibble(
    zeros = lapply(draws, function(t) flip(Z, t)),
    lambda = lapply(draws, function(t) flip(L, t)),
    omega = lapply(draws, function(t) flip(O, t)),
    alpha = lapply(draws, function(t) rep(0, p)),
    gamma_k = lapply(draws, function(t) if (t %% 2 == 0) c(2, 1) else c(1, 2)),
    cutpoints = lapply(draws, function(t) cuts)
  )
  res <- list("Number of Dimensions" = rep(2, 6), "Log-Likelihoods" = rep(0, 6),
              "MCMC Results" = store)
  class(res) <- "bpfa.results"
  attr(res, "mode") <- "fixed"; attr(res, "chains") <- 1; attr(res, "parallel") <- FALSE
  res
}

test_that("process_results averages relabelled draws and reports copula-scale loadings", {
  # dimension 1 is denser, so it is reported first
  L <- cbind(c(1, -1, 2, 0.5, 1, 0), c(0, 0, 0, 0, 2, 1))
  O <- cbind(c(1, 2, 3, 4), c(-1, 0, 1, 2))
  utils::capture.output(out <- process_results(fake_results(L, O), SE = FALSE))

  expect_equal(unname(as.matrix(out[["Factor Loadings"]][, -1])), L)
  expect_equal(unname(as.matrix(out[["Factor Scores"]][, -1])), O)
  expect_equal(out[["Factor Precisions"]]$gamma_k, c(1, 2))
  copula <- L / sqrt(1 + rowSums(L^2))
  expect_equal(unname(as.matrix(out[["Copula Loadings"]][, -1])), copula, tolerance = 1e-6)
})

test_that("check_accuracy compares as many dimensions as are both true and stable", {
  n <- 4; p <- 3
  L <- matrix(c(1, -1, 0.5, 0.2, 1, -0.3, 0.7, 0, 1), p, 3)
  O <- matrix(c(1, -1, 0.5, 2, 0, 1, -1, 1, 0.3, -0.2, 1, -1), n, 3)
  est_L <- cbind(L[, 1] + 3, L[, 2:3], 0)            # 4 stable dims; dim 1 is off by 3
  est_O <- cbind(O, 0)
  res <- list(
    "Stable Dimensions" = 1:4,
    "Factor Loadings" = tibble::as_tibble(cbind(row = 1:p, est_L), .name_repair = "minimal"),
    "Factor Scores" = tibble::as_tibble(cbind(row = 1:n, est_O), .name_repair = "minimal"),
    "Factor Precisions" = tibble::tibble(dimension = 1:4, gamma_k = rep(1, 4)),
    "Item-level Intercepts" = tibble::tibble(item = 1:p, alpha = rep(0, p)),
    "Cutpoints" = rep(list(c(0, Inf)), p)
  )
  class(res) <- "mcmc.output.processed"
  attr(res, "Include standard errors") <- FALSE
  attr(res, "Mode") <- "fixed"
  dat <- (O %*% t(L) > 0) * 1
  truth <- list(O %*% t(L), L, (L != 0) * 1, O, rep(0, p), rep(0.5, 3), rep(1, 3))
  acc <- check_accuracy(dat, res, true_vals = truth, confusion = FALSE)
  m <- acc[["Accuracy Metrics"]]
  # squared error 9 on one of three loading columns: sqrt(9 * 3 / 9)
  expect_equal(m$rmse[m$parameters == "lambda"], sqrt(3))
})

test_that("cutpoints are processed per column when every column has the same number of levels", {
  L <- cbind(c(1, -1, 2, 0.5, 1, 0), c(0, 0, 0, 0, 2, 1))
  O <- cbind(c(1, 2, 3, 4), c(-1, 0, 1, 2))
  res <- fake_results(L, O)
  attr(res, "mode") <- "multi"
  utils::capture.output(out <- process_results(res, SE = FALSE))
  expect_length(out[["Cutpoints"]], nrow(L))
  expect_equal(out[["Cutpoints"]][[1]], c(0, Inf))
})
