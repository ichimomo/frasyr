library(frasyr)

context("future SRdata")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  SRdata_pma_check <- get.SRdata(res_pma)
  SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
  SRdata0usingPeriodFrom1990To2000 <- get.SRdata(res_pma,years=1990:2000)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRdata0_pma_check <- read.csv(system.file("extdata","future_SRdata0_pma_check.csv",package="frasyr"),row.names=1)
  SRdata0usingPeriodFrom1990To2000_pma_check <- read.csv(system.file("extdata","future_SRdata0usingPeriodFrom1990To2000_pma_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合
  expect_equal(SRdata_pma_check, SRdata_pma)

  expect_equal(SRdata0$year, SRdata0_pma_check$year)
  #expect_equal(SRdata0$SSB, SRdata0_pma_check$SSB)
  #expect_equal(SRdata0$R, SRdata0_pma_check$R)

  expect_equal(SRdata0usingPeriodFrom1990To2000$year, SRdata0usingPeriodFrom1990To2000_pma_check$year)
  expect_equal(SRdata0usingPeriodFrom1990To2000$SSB, SRdata0usingPeriodFrom1990To2000_pma_check$SSB)
  expect_equal(SRdata0usingPeriodFrom1990To2000$R, SRdata0usingPeriodFrom1990To2000_pma_check$R)

})

context("future fitSR")

test_that("oututput value check",{
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRpma_HS_L1_AR0_check <- fit.SR(SRdata_pma,SR="HS",method="L1",AR=0,hessian=FALSE)
  SRpma_HS_L1_AR1_check <- fit.SR(SRdata_pma,SR="HS",method="L1",AR=1,hessian=FALSE)
  SRpma_HS_L2_AR0_check <- fit.SR(SRdata_pma,SR="HS",method="L2",AR=0,hessian=FALSE)
  SRpma_HS_L2_AR1_check <- fit.SR(SRdata_pma,SR="HS",method="L2",AR=1,hessian=FALSE)
  SRpma_BH_L1_AR0_check <- fit.SR(SRdata_pma,SR="BH",method="L1",AR=0,hessian=FALSE)
  SRpma_BH_L1_AR1_check <- fit.SR(SRdata_pma,SR="BH",method="L1",AR=1,hessian=FALSE)
  SRpma_BH_L2_AR0_check <- fit.SR(SRdata_pma,SR="BH",method="L2",AR=0,hessian=FALSE)
  SRpma_BH_L2_AR1_check <- fit.SR(SRdata_pma,SR="BH",method="L2",AR=1,hessian=FALSE)
  SRpma_RI_L1_AR0_check <- fit.SR(SRdata_pma,SR="RI",method="L1",AR=0,hessian=FALSE)
  SRpma_RI_L1_AR1_check <- fit.SR(SRdata_pma,SR="RI",method="L1",AR=1,hessian=FALSE)
  SRpma_RI_L2_AR0_check <- fit.SR(SRdata_pma,SR="RI",method="L2",AR=0,hessian=FALSE)
  SRpma_RI_L2_AR1_check <- fit.SR(SRdata_pma,SR="RI",method="L2",AR=1,hessian=FALSE)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","SRpma_HS_L1_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_HS_L1_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_HS_L2_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_HS_L2_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L1_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L1_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L2_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L2_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L1_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L1_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L2_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L2_AR1.rda",package = "frasyr"))


  #結果の数値を照合
  #HS.L1.AR0
  expect_equal(SRpma_HS_L1_AR0_check, SRpma_HS_L1_AR0)
  #HS.L1.AR1
  expect_equal(SRpma_HS_L1_AR1_check, SRpma_HS_L1_AR1)
  #HS.L2.AR0
  expect_equal(SRpma_HS_L2_AR0_check, SRpma_HS_L2_AR0)
  #HS.L2.AR1
  expect_equal(SRpma_HS_L2_AR1_check, SRpma_HS_L2_AR1)

  #BH.L1.AR0
  expect_equal(SRpma_BH_L1_AR0_check, SRpma_BH_L1_AR0)
  #BH.L1.AR1
  expect_equal(SRpma_BH_L1_AR1_check, SRpma_BH_L1_AR1)
  #BH.L2.AR0
  expect_equal(SRpma_BH_L2_AR0_check, SRpma_BH_L2_AR0)
  #BH.L2.AR1
  expect_equal(SRpma_BH_L2_AR1_check, SRpma_BH_L2_AR1)

  #RI.L1.AR0
  expect_equal(SRpma_RI_L1_AR0_check, SRpma_RI_L1_AR0)
  #RI.L1.AR1
  expect_equal(SRpma_RI_L1_AR1_check, SRpma_RI_L1_AR1)
  #RI.L2.AR0
  expect_equal(SRpma_RI_L2_AR0_check, SRpma_RI_L2_AR0)
  #RI.L2.AR1
  expect_equal(SRpma_RI_L2_AR1_check, SRpma_RI_L2_AR1)

})
