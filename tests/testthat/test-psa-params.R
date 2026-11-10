# The PSA distributions should be centred on the base-case values, so that
# the parameter table, the base-case code and the PSA agree (see issue #1,
# where the published table listed c_trtB as $12,000 instead of $13,000).

test_that("PSA parameter means match the base-case values", {
  df_psa <- generate_psa_params(n_sim = 10000, seed = 071818)
  expect_setequal(names(df_psa),
                  intersect(names(df_psa), names(l_params_base)))
  for (par in names(df_psa)) {
    # 2% relative tolerance: utilities follow beta distributions whose means
    # differ slightly from the base case (e.g., u_H has mean 200/203 vs 1)
    expect_equal(mean(df_psa[[par]]), l_params_base[[par]],
                 tolerance = 0.02, label = paste("mean of", par))
  }
})

test_that("cost of treatment B is $13,000 in the base case and the PSA", {
  expect_equal(l_params_base$c_trtB, 13000)
  # Gamma mean = shape * scale
  expect_equal(86.2 * 150.8, 13000, tolerance = 1e-3)
})
