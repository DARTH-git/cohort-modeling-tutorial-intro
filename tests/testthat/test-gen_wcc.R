# Tests for gen_wcc() in R/Functions.R
# Simpson's rules must integrate polynomials of degree <= 3 exactly on [0, n]:
#   sum_t w_t * t^k == n^(k + 1) / (k + 1)

test_that("Simpson 1/3 weights are assigned by time t, not vector position", {
  expect_equal(gen_wcc(2, "Simpson1/3"), c(1, 4, 1) / 3)
  expect_equal(gen_wcc(4, "Simpson1/3"), c(1, 4, 2, 4, 1) / 3)
  expect_equal(gen_wcc(6, "Simpson1/3"), c(1, 4, 2, 4, 2, 4, 1) / 3)
  # Entry for t = 1 (second position) must be 4/3
  expect_equal(gen_wcc(10)[2], 4/3)
})

test_that("Simpson 1/3 is the default method", {
  expect_identical(gen_wcc(10), gen_wcc(10, "Simpson1/3"))
})

test_that("odd number of cycles uses Simpson 3/8 on the last three cycles", {
  # n = 3: pure Simpson's 3/8 rule
  expect_equal(gen_wcc(3), c(3, 9, 9, 3) / 8)
  # n = 5: Simpson's 1/3 on [0, 2] plus 3/8 on [2, 5]
  expect_equal(gen_wcc(5), c(1/3, 4/3, 1/3 + 3/8, 9/8, 9/8, 3/8))
})

test_that("Simpson weights sum to n and integrate cubics exactly", {
  for (n in c(2:12, 74, 75, 100)) {
    w   <- gen_wcc(n, "Simpson1/3")
    v_t <- 0:n
    expect_length(w, n + 1)
    expect_true(all(w > 0))
    for (k in 0:3) {
      expect_equal(sum(w * v_t^k), n^(k + 1) / (k + 1),
                   label = paste0("n = ", n, ", k = ", k))
    }
  }
})

test_that("reported examples (issue #2, email of 2025-01-10) are fixed", {
  w <- gen_wcc(4, "Simpson1/3")
  expect_equal(sum(w), 4)                 # was 3.33 in v0.2.1
  expect_equal(sum(w * (0:4)^2), 64 / 3)  # was 17.33 in v0.2.1
})

test_that("half-cycle and none methods are unchanged", {
  expect_equal(gen_wcc(4, "half-cycle"), c(0.5, 1, 1, 1, 0.5))
  expect_equal(gen_wcc(1, "half-cycle"), c(0.5, 0.5))
  expect_equal(gen_wcc(4, "none"), rep(1, 5))
  expect_equal(sum(gen_wcc(75, "half-cycle")), 75)
})

test_that("invalid inputs give informative errors", {
  expect_error(gen_wcc(0), "positive integer")
  expect_error(gen_wcc(-3), "positive integer")
  expect_error(gen_wcc(7.5), "positive integer")
  expect_error(gen_wcc(c(4, 6)), "positive integer")
  expect_error(gen_wcc(NA_real_), "positive integer")
  expect_error(gen_wcc("10"), "positive integer")
  expect_error(gen_wcc(1, "Simpson1/3"), "n_cycles` >= 2")
  expect_error(gen_wcc(10, "trapezoid"))
})

test_that("integer-valued doubles are accepted", {
  # e.g., n_cycles <- (100 - 25) / 1
  expect_equal(gen_wcc((100 - 25) / 1), gen_wcc(75L))
})
