# The model has no item intercepts. Under the rank likelihood an intercept is not
# identified (shifting x_j and alpha_j together changes nothing), so item difficulty
# lives in the level thresholds, which are identified once the latent location is
# fixed by Lambda Omega'. With no factor structure, x_ij ~ N(0, 1 + |lambda_j|^2) with
# lambda near zero, so the threshold of a binary item with P(y = 1) = q is about
# qnorm(1 - q).

test_that("thresholds of independent binary items recover qnorm(1 - q) from any seed", {
  skip_on_cran()
  set.seed(11)
  n <- 600; q <- c(0.2, 0.5, 0.8)
  y <- sapply(rep(q, each = 4), function(qq) rbinom(n, 1, qq))
  utils::capture.output(res <- mosaicFA_sampler(y, K_start = 1, warmup = 1000, iter = 4000, thin = 4,
                                              chains = 2, seed = 3))
  out <- suppressWarnings(process_results(res, SE = FALSE))
  th <- out[["Item Thresholds"]]
  want <- qnorm(1 - colMeans(y))
  expect_equal(th$item, 1:12)
  expect_lt(max(abs(th$mean - want)), 0.15)
  # both chains agree: the thresholds do not drift
  per_chain <- sapply(res, function(ch) ch$thresholds)
  expect_lt(max(abs(per_chain[, 1] - per_chain[, 2])), 0.1)
})

test_that("item thresholds come with credible intervals for items with few levels", {
  set.seed(12)
  y <- cbind(matrix(rbinom(80 * 4, 1, 0.4), 80, 4), matrix(rbinom(80 * 2, 3, 0.5), 80, 2), rnorm(80))
  utils::capture.output(res <- mosaicFA_sampler(y, K_start = 1, warmup = 20, iter = 40, thin = 2, seed = 1))
  out <- suppressWarnings(process_results(res))
  th <- out[["Item Thresholds"]]
  expect_null(out[["Item-level Intercepts"]])
  expect_equal(th$item[1:4], 1:4)                                   # one boundary per binary item
  expect_equal(sum(th$item == 5), length(unique(y[, 5])) - 1)
  expect_true(all(c("lower_value", "upper_value", "mean", "ci_0.05", "ci_0.95") %in% names(th)))
  few <- th$item <= 6
  expect_true(all(th[["ci_0.05"]][few] <= th$mean[few] & th$mean[few] <= th[["ci_0.95"]][few]))
  expect_true(all(is.na(th[["ci_0.05"]][!few])))                    # continuous item: mean only
  expect_null(res[[1]]$alpha)
})
