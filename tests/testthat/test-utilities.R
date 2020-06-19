context("Utilities")
test_that("make_kobe_ratio", {

  result_vpa  <- load_data("../../inst/extdata/res_vpa_pma.rda")
  result_msy  <- load_data("../../inst/extdata/res_MSY_pma_pre.rda")

  kobe_ratio <- make_kobe_ratio(result_vpa, result_msy)

  expect_is(kobe_ratio, "data.frame")
  expect_equal(colnames(kobe_ratio), c("year", "Fratio", "Bratio"))
  expect_equal(kobe_ratio$year, as.character(1982:2011))
  expect_is(kobe_ratio$Fratio, "numeric")
  expect_is(kobe_ratio$Bratio, "numeric")
})

test_that("pull single table from table list", {
  # pull_var_from_kobeII_table() is not tested yet
  expect_equal(1, 1)
})

test_that("convert table from kobeIItable to the required format", {
  before <- data.frame(beta = c(1, 0.5, 0),
                       y2019 = seq(12345, 11234, length.out = 3),
                       y2020 = seq(22345, 12345, length.out = 3),
                       y2021 = seq(32345, 23456, length.out = 3))

  after_raw <- format_beta_table(before, divide_by = 1)
  expect_equal(after_raw[, 1], c(1, 0.5, 0))
  expect_equal(after_raw[, "y2019"], c(12345, 11790, 11234))
  expect_equal(after_raw[, "y2020"], c(22345, 17345, 12345))
  expect_equal(after_raw[, "y2021"], c(32345, 27900, 23456))

  after_ton <- format_beta_table(before, divide_by = 1000, round = FALSE)
  expect_equal(after_ton[, 1], c(1, 0.5, 0))
  expect_equal(after_ton[, "y2019"], c(12.345, 11.790, 11.234), tolerance = 5e-4)
  expect_equal(after_ton[, "y2020"], c(22.345, 17.345, 12.345), tolerance = 5e-4)
  expect_equal(after_ton[, "y2021"], c(32.345, 27.900, 23.456), tolerance = 5e-4)

  after_ton_rounded <- format_beta_table(before, divide_by = 1000, round = TRUE)
  expect_equal(after_ton_rounded[, 1], c(1, 0.5, 0))
  expect_equal(after_ton_rounded[, "y2019"], c(12, 12, 11))
  expect_equal(after_ton_rounded[, "y2020"], c(22, 17, 12))
  expect_equal(after_ton_rounded[, "y2021"], c(32, 28, 23))
})

test_that("test for HCR function", {
    res_HCR <- HCR_default(matrix(1:10,2,5),matrix(5,2,5),matrix(1,2,5),matrix(0.8,2,5))
    expect_equal(res_HCR, matrix(c(0,0.2,0.4,0.6,rep(0.8,6)),2,5))
})

