# Post-processing of sampler output (design §4.4–4.5): counted dimensions,
# relabelling and sign alignment across draws and chains, copula loadings,
# the K distribution, convergence diagnostics and retrievable priors.

test_that("dimensions are matched to the pivot draw by their sparsity patterns", {
  pivot <- cbind(c(1, 1, 1, 1, 0, 0), c(0, 0, 1, 1, 1, 1), c(1, 0, 1, 0, 1, 0))
  permuted <- pivot[, c(3, 1, 2)]
  permuted[1, 3] <- 1 - permuted[1, 3]          # one flipped entry
  partial <- pivot[, c(3, 1)]
  perms <- align_dimensions(list(pivot, permuted, partial), pivot = 1)
  expect_equal(perms, list(c(1L, 2L, 3L), c(2L, 3L, 1L), c(2L, NA, 1L)))
})

# H = 3 columns; columns 1 and 2 are on with >= 3 items, column 3 is off.
# Even-numbered draws store the two counted dimensions in swapped order, and
# chain 2 stores dimension 1 with flipped sign.
fake_results <- function(L, O, n_draws = 6, chains = 1, flip_chain2 = FALSE) {
  p <- nrow(L); n <- nrow(O)
  one_chain <- function(ch) {
    draws <- seq_len(n_draws)
    swap <- function(m, t) if (t %% 2 == 0) m[, c(2, 1, 3), drop = FALSE] else m
    sgn <- if (flip_chain2 && ch == 2) c(-1, 1, 1) else c(1, 1, 1)
    Lf <- cbind(L, 0.01) * rep(sgn, each = p); Of <- cbind(O, rnorm(n)) * rep(sgn, each = n)
    Zf <- (cbind(L, 0) != 0) * 1L
    list(lambda = lapply(draws, function(t) swap(Lf, t)),
         z = lapply(draws, function(t) swap(Zf, t)),
         c = lapply(draws, function(t) c(3L, 3L, 1L)),
         omega = lapply(draws, function(t) swap(Of, t)),
         theta = lapply(draws, function(t) swap(matrix(c(1, 2, 0.05), 1), t)[1, ]),
         a_cusp = rep(1, n_draws), K = rep(2L, n_draws), log_lik = rep(-10, n_draws),
         thresholds = rep(0, p), threshold_draws = matrix(0, p, n_draws))
  }
  res <- lapply(seq_len(chains), one_chain)
  structure(res, class = "bpfa.results", priors = resolve_priors("uniform", p, H = 3),
            n = n, p = p, levels = rep(list(c(0, 1)), p), warmup = 0, iter = n_draws,
            thin = 1, chains = chains)
}

L <- cbind(c(1, -1, 2, 0.5, 1, 0), c(0, 0.7, 0, 0.4, 2, 1))     # dimension 1 is denser
O <- cbind(c(1, 2, 3, 4), c(-1, 0, 1, 2))

test_that("relabelled draws average to the true loadings and scores", {
  set.seed(1)
  out <- process_results(fake_results(L, O), SE = FALSE)
  expect_equal(unname(as.matrix(out[["Factor Loadings"]][, -1])), L)
  expect_equal(unname(as.matrix(out[["Factor Scores"]][, -1])), O)
})

test_that("chains that found a dimension with opposite sign are aligned before averaging", {
  set.seed(2)
  out <- process_results(fake_results(L, O, chains = 2, flip_chain2 = TRUE), SE = FALSE)
  expect_equal(unname(as.matrix(out[["Factor Loadings"]][, -1])), L)
  expect_equal(unname(as.matrix(out[["Factor Scores"]][, -1])), O)
})

test_that("copula loadings use each draw's full loadings row, then keep counted dimensions", {
  set.seed(3)
  out <- process_results(fake_results(L, O), SE = FALSE)
  full <- cbind(L, 0.01)
  copula <- (full / sqrt(1 + rowSums(full^2)))[, 1:2]
  expect_equal(unname(as.matrix(out[["Copula Loadings"]][, -1])), copula, tolerance = 1e-6)
})

test_that("only on columns with at least three items are counted", {
  set.seed(4)
  out <- process_results(fake_results(L, O), SE = FALSE)
  expect_equal(out[["Stable Dimensions"]], 1:2)
  expect_equal(as.vector(out[["K Distribution"]]), 6)        # all six draws have K = 2
  expect_equal(names(out[["K Distribution"]]), "2")
})

test_that("the uncounted share is the correlation carried by columns that are not counted", {
  set.seed(5)
  out <- process_results(fake_results(L, O), SE = FALSE)
  full <- cbind(L, 0.01)
  cc <- function(M) { C <- tcrossprod(M) + diag(nrow(M)); D <- 1 / sqrt(diag(C)); C * outer(D, D) }
  want <- mean(abs(cc(full) - cc(L))[upper.tri(diag(6))])
  expect_equal(out[["Uncounted Correlation Share"]], want, tolerance = 1e-8)
  expect_gt(want, 0)
})

test_that("priors travel with the processed results", {
  out <- process_results(fake_results(L, O), SE = FALSE)
  expect_equal(get_priors(out)$spike_var, 0.05)
})

test_that("chains that disagree on K produce a convergence warning", {
  res <- fake_results(L, O, chains = 2)
  res[[2]]$c <- lapply(res[[2]]$c, function(cc) c(3L, 1L, 1L))   # chain 2: one dimension on
  expect_warning(out <- process_results(res, SE = FALSE), "Chains disagree on the number of dimensions")
  expect_gt(out[["R-hat K"]], 1.1)
})

test_that("check_accuracy compares as many dimensions as are both true and stable", {
  n <- 4; p <- 3
  Lt <- matrix(c(1, -1, 0.5, 0.2, 1, -0.3, 0.7, 0, 1), p, 3)
  Ot <- matrix(c(1, -1, 0.5, 2, 0, 1, -1, 1, 0.3, -0.2, 1, -1), n, 3)
  est_L <- cbind(Lt[, 1] + 3, Lt[, 2:3], 0)           # 4 stable dims; dim 1 is off by 3
  est_O <- cbind(Ot, 0)
  res <- list(
    "Stable Dimensions" = 1:4,
    "Factor Loadings" = tibble::as_tibble(cbind(row = 1:p, est_L), .name_repair = "minimal"),
    "Factor Scores" = tibble::as_tibble(cbind(row = 1:n, est_O), .name_repair = "minimal"),
    "Slab Variances" = tibble::tibble(dimension = 1:4, theta = rep(1, 4)),
    "Cutpoints" = rep(list(c(0, Inf)), p)
  )
  class(res) <- "mcmc.output.processed"
  attr(res, "Include standard errors") <- FALSE
  dat <- (Ot %*% t(Lt) > 0) * 1
  truth <- list(Ot %*% t(Lt), Lt, (Lt != 0) * 1, Ot, rep(0, p), rep(0.5, 3), rep(1, 3))
  acc <- check_accuracy(dat, res, true_vals = truth, confusion = FALSE)
  m <- acc[["Accuracy Metrics"]]
  # squared error 9 on one of three loading columns: sqrt(9 * 3 / 9)
  expect_equal(m$rmse[m$parameters == "lambda"], sqrt(3))
})

test_that("process_results runs on real sampler output, with standard errors", {
  y <- generate_sim_data(n = 60, p = 30, K_true = 2, noise = 0, seed = 1)$data
  utils::capture.output(res <- mmBPFA_sampler(y, K_start = 4, warmup = 50, iter = 40, thin = 2, chains = 2))
  out <- suppressWarnings(process_results(res))
  expect_s3_class(out, "mcmc.output.processed")
  expect_true(all(c("mean", "ci_0.05", "ci_0.95") %in% names(out[["Factor Loadings"]])))
  expect_length(out[["Cutpoints"]], 30)
  acc <- check_accuracy(y, suppressWarnings(process_results(res, SE = FALSE)), confusion = TRUE)
  expect_true(acc$GMP > 0 && acc$GMP < 1)
})

test_that("results with no counted dimensions are processed into empty loadings and scores", {
  res <- fake_results(L, O)
  res[[1]]$c <- lapply(res[[1]]$c, function(cc) c(1L, 1L, 1L))     # every column off
  out <- process_results(res, SE = FALSE)
  expect_equal(out[["Median Number Dimensions"]], 0)
  expect_length(out[["Stable Dimensions"]], 0)
  expect_equal(dim(out[["Factor Loadings"]]), c(6L, 1L))           # the row column only
  expect_equal(dim(out[["Factor Scores"]]), c(4L, 1L))
  out_se <- process_results(res)
  expect_equal(nrow(out_se[["Factor Loadings"]]), 0)
})
