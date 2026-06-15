# Characterization (snapshot) tests for the TAC carryover / borrowing logic
# in future_vpa_R (R/future.r, currently L1008-1094).
#
# Purpose: lock down the *current* numeric behavior before extracting the
# carryover block into apply_TAC_carryover() (refactor-tac-carryover, phase 1).
# These snapshots are the safety net: phase 2 (behavior-preserving extraction)
# must keep them green; phase 4 (intentional bug fixes) will deliberately turn
# them red and the snapshots get reviewed + updated.
#
# Coverage matrix (do_MSE x CV x banking/borrowing/amount), see
# ../../frasyr_memo/frasyr_TAC_carryover_refactor.md phase 1.

context("TAC carryover characterization snapshots")

options(warn = -1)
data(res_vpa_org)
data(res_sr_HSL2)

# --- compact, deterministic base future-data input -------------------------
# nsim small + fixed seed -> reproducible; carryover-relevant settings only.
.carryover_base_input <- function() {
  list(
    res_vpa = res_vpa_org,
    nsim = 5, nyear = 10,
    future_initial_year_name = 2017,
    start_F_year_name = 2018,
    start_biopar_year_name = 2018,
    start_random_rec_year_name = 2018,
    waa_year = 2015:2017, waa = NULL,
    waa_catch_year = 2015:2017, waa_catch = NULL,
    maa_year = 2015:2017, maa = NULL,
    M_year = 2015:2017, M = NULL,
    faa_year = 2015:2017,
    currentF = NULL, futureF = NULL,
    start_ABC_year_name = 2019,
    HCR_beta = 1, HCR_Blimit = -1, HCR_Bban = -1,
    HCR_year_lag = 0, HCR_function_name = "HCR_default",
    res_SR = res_sr_HSL2,
    seed_number = 1,
    resid_type = "lognormal",
    resample_year_range = 0,
    bias_correction = TRUE,
    recruit_intercept = 0,
    Pope = res_vpa_org$input$Pope
  )
}

# build future data with carryover overrides applied
.make_carryover_data <- function(overrides = list()) {
  input <- .carryover_base_input()
  input[names(overrides)] <- overrides
  safe_call(make_future_data, input)
}

# extract the carryover-relevant slices, rounded for portable comparison
.snap_carryover <- function(res, years = 2019:2027) {
  yy <- as.character(years)
  flds <- c("wcatch", "original_ABC", "original_ABC_plus", "reserved_catch")
  list(
    realized      = signif(res$HCR_realized[yy, , flds], 6),
    expect_wcatch = signif(res$HCR_mat[yy, , "expect_wcatch"], 6)
  )
}

test_that("TAC carryover numeric behavior is stable (characterization)", {
  testthat::local_edition(3)

  ## --- non-MSE, no CV --------------------------------------------------
  # (A) banking by rate, denom = original_ABC_plus
  resA <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_rate = 0.1, HCR_TAC_carry_rate = NA,
    HCR_reserve_denom = "original_ABC_plus"))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resA), style = "serialize", tolerance = 1e-6)

  # (B) borrowing by rate, alternating across sims, carry_rate = 1
  resB <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_rate = c(-0.1, 0), HCR_TAC_carry_rate = 1))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resB), style = "serialize", tolerance = 1e-6)

  # (C) banking by amount, carry by amount
  resC <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_amount = 3000, HCR_TAC_carry_amount = 1000,
    HCR_TAC_carry_rate = NA, HCR_TAC_reserve_rate = NA))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resC), style = "serialize", tolerance = 1e-6)

  # (D) borrowing + Blimit guard (SSB < Blimit -> no borrowing)
  resD <- future_vpa(.make_carryover_data(list(
    HCR_Blimit = 26000 * 10,
    HCR_TAC_reserve_rate = c(-0.1, 0), HCR_TAC_carry_rate = 1))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resD), style = "serialize", tolerance = 1e-6)

  ## --- non-MSE + CV limits --------------------------------------------
  # (E) banking by rate + upper/lower CV
  resE <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_rate = 0.1, HCR_TAC_carry_rate = NA,
    HCR_reserve_denom = "original_ABC_plus",
    HCR_TAC_upper_CV = 0.1, HCR_TAC_lower_CV = 0.1))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resE), style = "serialize", tolerance = 1e-6)

  # (F) borrowing by rate + upper/lower CV
  resF <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_rate = c(-0.1, 0), HCR_TAC_carry_rate = 1,
    HCR_TAC_upper_CV = 0.1, HCR_TAC_lower_CV = 0.1))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resF), style = "serialize", tolerance = 1e-6)

  ## --- denom = "original_ABC" (rate_original) path ----------------------
  # (I) banking by rate, denom = original_ABC (hits L1040 + 0.01 floor)
  resI <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_rate = 0.1, HCR_TAC_carry_rate = NA,
    HCR_reserve_denom = "original_ABC"))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resI), style = "serialize", tolerance = 1e-6)

  ## --- amount borrowing -----------------------------------------------
  # (J) borrowing by amount (hits amount borrowing branch + Blimit guard + floor)
  resJ <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_amount = -3000, HCR_TAC_carry_amount = NA,
    HCR_TAC_carry_rate = NA, HCR_TAC_reserve_rate = NA))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resJ), style = "serialize", tolerance = 1e-6)

  ## --- year-to-year banking<->borrowing mix (weakness 1) ---------------
  # (K) alternating banking/borrowing across years, denom = original_ABC.
  #     Locks current (weakness-1, asymmetric A vs P base) behavior; phase 4
  #     will deliberately update this snapshot when fixing the asymmetry.
  #     NOTE: sim-wise sign mixing is out of scope (cannot arise via the args;
  #     vectors recycle over YEARS, identical across sims) and is left for a
  #     phase-4 guard, not snapshotted here.
  resK <- future_vpa(.make_carryover_data(list(
    HCR_TAC_reserve_rate = c(0.1, -0.1), HCR_reserve_denom = "original_ABC",
    HCR_TAC_carry_rate = 1))$data,
    optim_method = "none", multi_init = 1)
  expect_snapshot_value(.snap_carryover(resK), style = "serialize", tolerance = 1e-6)
})

test_that("TAC carryover numeric behavior is stable under MSE (characterization)", {
  testthat::local_edition(3)

  base <- .make_carryover_data(list(
    HCR_TAC_reserve_rate = 0.1, HCR_TAC_carry_rate = NA,
    HCR_reserve_denom = "original_ABC_plus"))

  # (G) MSE banking by rate
  resG <- future_vpa(tmb_data = base$data,
                     optim_method = "none", multi_init = 1, SPRtarget = 0.3,
                     do_MSE = TRUE, MSE_input_data = base, MSE_nsim = 30)
  expect_snapshot_value(.snap_carryover(resG), style = "serialize", tolerance = 1e-6)

  # (H) MSE borrowing + TAC_adjust
  base2 <- .make_carryover_data(list(
    HCR_TAC_reserve_rate = c(-0.1, 0), HCR_TAC_adjust = 1, HCR_TAC_carry_rate = 10))
  resH <- future_vpa(tmb_data = base2$data,
                     optim_method = "none", multi_init = 1, SPRtarget = 0.3,
                     do_MSE = TRUE, MSE_input_data = base2, MSE_nsim = 30)
  expect_snapshot_value(.snap_carryover(resH), style = "serialize", tolerance = 1e-6)
})
