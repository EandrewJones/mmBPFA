test_that("remove_weakly_active returns zero dimensions when none pass the threshold", {
  # P = 200 items, each dimension active on 1 item: m_k = 1 is not > 0.01 * P = 2
  zeros <- diag(1, 200, 3)
  setTimeLimit(elapsed = 10, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  out <- remove_weakly_active(zeros, zeros, matrix(0, 5, 3), rep(1, 3), rep(20L, 3))
  expect_equal(ncol(out$lambda), 0)
  expect_equal(ncol(out$omega), 0)
  expect_length(out$gamma_k, 0)
})

test_that("remove_weakly_active keeps established dimensions with enough loadings", {
  zeros <- cbind(rep(1, 200), c(1, rep(0, 199)), rep(1, 200))
  dc <- c(20L, 20L, 3L)   # third dimension is too new
  out <- remove_weakly_active(zeros, zeros * 2, matrix(1:15, 5, 3), 1:3, dc)
  expect_equal(out$gamma_k, 1)
  expect_equal(out$omega, matrix(1:5, 5, 1))
})
