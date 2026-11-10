# check_sum_of_transition_array() must accept valid arrays whose cycle sums
# differ from n_states only by floating-point rounding, and reject invalid
# ones. (The exact `==` comparison used up to v0.3.1 made the PSA of the
# state-residence model in cohort-modeling-tutorial-timedep fail.)

test_that("transition array check tolerates floating-point rounding", {
  # Rounding-sized error added explicitly: whether a given set of
  # probabilities sums exactly to n_states depends on the platform
  a_P <- array(rep(c(0.5, 0.3, 0.2), each = 3), dim = c(3, 3, 2))
  a_P[1, 1, ] <- a_P[1, 1, ] + 1e-13
  expect_false(sum(rowSums(a_P[, , 1])) == 3)
  expect_no_error(check_sum_of_transition_array(a_P, n_states = 3, n_cycles = 2,
                                                verbose = FALSE))
  # An invalid array is still rejected
  a_P[1, 1, 2] <- 0.9
  expect_error(check_sum_of_transition_array(a_P, n_states = 3, n_cycles = 2,
                                             verbose = FALSE),
               "not a valid transition array")
})

test_that("transition matrix check rejects rows that do not sum to one", {
  m_P <- matrix(rep(c(0.5, 0.3, 0.2), each = 3), nrow = 3)
  expect_no_error(check_sum_of_transition_array(m_P, n_states = 3, n_cycles = 1,
                                                verbose = FALSE))
  m_P[1, 1] <- 0.9
  expect_error(check_sum_of_transition_array(m_P, n_states = 3, n_cycles = 1,
                                             verbose = FALSE),
               "not a valid transition matrix")
})
