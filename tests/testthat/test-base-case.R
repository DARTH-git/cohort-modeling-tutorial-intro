# Regression test for the base-case cost-effectiveness analysis of the
# time-independent Sick-Sicker model with Simpson's 1/3 (+ 3/8) within-cycle
# correction (v0.3.0). Reference values were produced by
# analysis/cSTM_time_indep.R after the gen_wcc() fix.

test_that("base-case costs and QALYs match reference values", {
  df_ce <- calculate_ce_out(l_params_base)
  expect_equal(df_ce$Strategy,
               c("Standard of care", "Strategy A", "Strategy B", "Strategy AB"))
  expect_equal(df_ce$Cost,
               c(152281.754738, 285889.452912, 260222.079374, 380377.893552),
               tolerance = 1e-9)
  expect_equal(df_ce$Effect,
               c(21.0318903681, 21.8259368701, 22.5049038291, 23.4638309628),
               tolerance = 1e-9)
})

test_that("base-case ICERs and dominance are unchanged qualitatively", {
  df_ce  <- calculate_ce_out(l_params_base)
  df_cea <- suppressWarnings(
    calculate_icers(df_ce$Cost, df_ce$Effect, df_ce$Strategy)
  )
  status <- setNames(df_cea$Status, df_cea$Strategy)
  expect_equal(unname(status[c("Standard of care", "Strategy B",
                               "Strategy AB", "Strategy A")]),
               c("ND", "ND", "ND", "D"))
  icer <- setNames(df_cea$ICER, df_cea$Strategy)
  expect_equal(unname(icer[c("Strategy B", "Strategy AB")]),
               c(73278.5731372, 125302.3404645), tolerance = 1e-9)
})

test_that("PSA wrapper uses the same within-cycle correction as the script", {
  # Total discounted, WCC-corrected life years with all utilities set to 1
  # must equal what the corrected gen_wcc() implies for the SoC cohort trace.
  l_par <- modifyList(l_params_base, list(u_S1 = 1, u_S2 = 1, u_trtA = 1))
  l_m_M <- decision_model(l_par)
  n_cycles <- 75
  v_alive  <- rowSums(l_m_M[[1]][, c("H", "S1", "S2")])
  v_dw     <- 1 / (1.03^(0:n_cycles))
  expect_equal(calculate_ce_out(l_par)$Effect[1],
               sum(v_alive * v_dw * gen_wcc(n_cycles)))
})
