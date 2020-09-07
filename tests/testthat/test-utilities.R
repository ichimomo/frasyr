context("Utilities")

result_vpa  <- load_data("../../inst/extdata/res_vpa_pma.rda")
result_msy  <- load_data("../../inst/extdata/res_MSY_pma_pre.rda")

test_that("make_kobe_ratio", {
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

test_that("calc_future_perSPR accepts list with different length vectors", {
  future_data <- generate_dummy_future_data(result_vpa)

  perspr <- calc_future_perSPR(fout = list(waa       = future_data$data$waa_mat,
                                           maa       = future_data$data$maa_mat,
                                           M         = future_data$data$M_mat,
                                           waa.catch = future_data$data$waa_catch_mat),
                               Fvector = apply_year_colum(result_vpa$faa, 2007:2011),
                               target.year = list(waa       = 2014:2018,
                                                  waa.catch = 2014:2018,
                                                  maa       = 2016:2018,
                                                  M         = 2014:2018))
  expect_is(perspr, "numeric")
})

test_that("test caa.est.mat", {

  expect_catch <- 0.5

  # set usual F => OK
  res <- caa.est.mat(c(1,1,1,1),c(1,1,1,1),c(1,1,1,1),c(0,0,0,0),catch.obs=expect_catch,Pope=TRUE,set_max1=FALSE)
  expect_equal(round(sum(res$caa),3),round(expect_catch,3))

  # set very small F => OK
  res <- caa.est.mat(c(1,1,1,1),c(0.00001,0.00001,0.00001,0.00001),c(1,1,1,1),c(0,0,0,0),
                     catch.obs=expect_catch,Pope=TRUE)
  expect_equal(round(sum(res$caa),3),round(expect_catch,3))

  # naa is very small => warning
  expect_warning(caa.est.mat(c(0.01,0.01,0.01,0.01),c(1,1,1,1),c(1,1,1,1),c(0,0,0,0),
                     catch.obs=expect_catch,Pope=TRUE))

})


test_that("calc.rel.abund"{

  age.test <- c(1:5)
  Fc.test <- rep(1,length(age.test))
  waa.test <- c(1:5)
  maa.test <- c(0,0.5,1,1,1)
  M.test <- rep(0.5,length(age.test))

  calc_rel_abund_popeT_check <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = TRUE)

  calc_rel_abund_popeF_check <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = FALSE)

  #上記計算内容をエクセルで計算したものを読み込み
  calc_rel_abund <- read.csv("./inst/extdata/check_calc_rel_abund.csv",header = T)
  #データ整形
  calc_rel_abund_popeT <- list(calc_rel_abund$rel.abundant,calc_rel_abund$ypr1.popeT,calc_rel_abund$spr)
  names(calc_rel_abund_popeT) <- c("rel.abund","ypr","spr")
  calc_rel_abund_popeF <- list(calc_rel_abund$rel.abundant,calc_rel_abund$ypr1.popeF,calc_rel_abund$spr)
  names(calc_rel_abund_popeF) <- c("rel.abund","ypr","spr")
  # 結果照合
  expect_equal(calc_rel_abund_popeT,calc_rel_abund_popeT_check)
  expect_equal(calc_rel_abund_popeF,calc_rel_abund_popeF_check)

})

test_that("catch_equation"{

  expect_equal(catch_equation(1,1,1,1), 1*(1-exp(-1))*exp(-0.5)*1)
  expect_equal(catch_equation(1,1,1,1,Pope = F), 1*(1-exp(-1-1))*1/(1+1)*1 )

})

test_that("solv.Feq"{

  age.test <- c(1:5)
  faa.test <- rep(1,length(age.test))
  waa.test <- c(1:5)
  naa.test <- rep(3,length(age.test))
  maa.test <- c(0,0.5,1,1,1)
  M.test <- rep(0.5,length(age.test))
  # Baranov eqをもちいてfaa,M,naaをつかってcaaを求める
  caa.test <- faa.test/(faa.test+M.test) *(1-exp(-faa.test- M.test)) *naa.test

  # tolerance=1e-4　でチェック
  expect_equal(faa.test, solv.Feq(cvec = caa.test,nvec = naa.test,mvec = M.test),tolerance=1e-4)

})


test_that("get.SPR" {

})


test_that("convert_faa_perSPR" {

})

