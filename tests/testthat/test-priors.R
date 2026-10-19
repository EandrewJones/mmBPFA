test_that("the uniform preset resolves to the copula-calibrated defaults", {
  pr <- resolve_priors("uniform", p = 60)
  expect_equal(pr$preset, "uniform")
  expect_equal(pr[c("slab_shape", "slab_rate", "spike_var", "a_cusp_shape", "a_cusp_rate", "rho_a", "rho_b")],
               list(slab_shape = 1, slab_rate = 0.5, spike_var = 0.05, a_cusp_shape = 1,
                    a_cusp_rate = 1, rho_a = 1, rho_b = 1))
  expect_equal(pr$H, 29)            # min(floor((60 - 1) / 2), 50)
})

test_that("raw values override the preset and H can be set", {
  pr <- resolve_priors(list(preset = "uniform", spike_var = 0.01, slab_rate = 0.25), p = 200, H = 12)
  expect_equal(pr$spike_var, 0.01)
  expect_equal(pr$slab_rate, 0.25)
  expect_equal(pr$H, 12)
  expect_equal(resolve_priors("uniform", p = 500)$H, 50)
})

test_that("unknown presets and fields are rejected", {
  expect_error(resolve_priors("bogus", p = 10), "Unknown prior preset")
  expect_error(resolve_priors(list(slab_rat = 1), p = 10), "Unknown prior field")
})

test_that("H below 1 is an error", {
  expect_error(resolve_priors("uniform", p = 2), "at least one column")
})
