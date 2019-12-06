library(frasyr)

context("future SRdata")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  SRdata_pma_check <- get.SRdata(res_pma)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  #結果の数値を照合
  expect_equal(SRdata_pma_check, SRdata_pma)
})

context("future fitSR")

test_that("oututput value check",{
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRpma_HS_L1_AR0_check <- fit.SR(SRdata_pma,SR="HS",method="L1",AR=0,hessian=FALSE)
  SRpma_HS_L1_AR1_check <- fit.SR(SRdata_pma,SR="HS",method="L1",AR=1,hessian=FALSE,out.AR = TRUE)
  SRpma_HS_L1_AR1_outAR_F_check <- fit.SR(SRdata_pma,SR="HS",method="L1",AR=1,hessian=FALSE,out.AR = FALSE)
  SRpma_HS_L2_AR0_check <- fit.SR(SRdata_pma,SR="HS",method="L2",AR=0,hessian=FALSE)
  SRpma_HS_L2_AR1_check <- fit.SR(SRdata_pma,SR="HS",method="L2",AR=1,hessian=FALSE,out.AR = TRUE)
  SRpma_HS_L2_AR1_outAR_F_check <- fit.SR(SRdata_pma,SR="HS",method="L2",AR=1,hessian=FALSE,out.AR = FALSE)
  SRpma_BH_L1_AR0_check <- fit.SR(SRdata_pma,SR="BH",method="L1",AR=0,hessian=FALSE)
  SRpma_BH_L1_AR1_check <- fit.SR(SRdata_pma,SR="BH",method="L1",AR=1,hessian=FALSE,out.AR = TRUE)
  SRpma_BH_L1_AR1_outAR_F_check <- fit.SR(SRdata_pma,SR="BH",method="L1",AR=1,hessian=FALSE,out.AR = FALSE)
  SRpma_BH_L2_AR0_check <- fit.SR(SRdata_pma,SR="BH",method="L2",AR=0,hessian=FALSE)
  SRpma_BH_L2_AR1_check <- fit.SR(SRdata_pma,SR="BH",method="L2",AR=1,hessian=FALSE,out.AR = TRUE)
  SRpma_BH_L2_AR1_outAR_F_check <- fit.SR(SRdata_pma,SR="BH",method="L2",AR=1,hessian=FALSE,out.AR = FALSE)
  SRpma_RI_L1_AR0_check <- fit.SR(SRdata_pma,SR="RI",method="L1",AR=0,hessian=FALSE)
  SRpma_RI_L1_AR1_check <- fit.SR(SRdata_pma,SR="RI",method="L1",AR=1,hessian=FALSE,out.AR = TRUE)
  SRpma_RI_L1_AR1_outAR_F_check <- fit.SR(SRdata_pma,SR="RI",method="L1",AR=1,hessian=FALSE,out.AR = FALSE)
  SRpma_RI_L2_AR0_check <- fit.SR(SRdata_pma,SR="RI",method="L2",AR=0,hessian=FALSE)
  SRpma_RI_L2_AR1_check <- fit.SR(SRdata_pma,SR="RI",method="L2",AR=1,hessian=FALSE,out.AR = TRUE)
  SRpma_RI_L2_AR1_outAR_F_check <- fit.SR(SRdata_pma,SR="RI",method="L2",AR=1,hessian=FALSE,out.AR = FALSE)

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


  #結果の数値を照合($par)----
  #HS.L1.AR0
  expect_equal(SRpma_HS_L1_AR0_check$par, SRpma_HS_L1_AR0$par)
  #HS.L1.AR1
  expect_equal(SRpma_HS_L1_AR1_outAR_F_check$par, SRpma_HS_L1_AR1$par)
  #HS.L2.AR0
  expect_equal(SRpma_HS_L2_AR0_check$par, SRpma_HS_L2_AR0$par)
  #HS.L2.AR1
  expect_equal(SRpma_HS_L2_AR1_outAR_F_check$par, SRpma_HS_L2_AR1$par)

  #BH.L1.AR0
  expect_equal(SRpma_BH_L1_AR0_check$par, SRpma_BH_L1_AR0$par)
  #BH.L1.AR1
  expect_equal(SRpma_BH_L1_AR1_outAR_F_check$par, SRpma_BH_L1_AR1$par)
  #BH.L2.AR0
  expect_equal(SRpma_BH_L2_AR0_check$par, SRpma_BH_L2_AR0$par)
  #BH.L2.AR1
  expect_equal(SRpma_BH_L2_AR1_outAR_F_check$par, SRpma_BH_L2_AR1$par)

  #RI.L1.AR0
  expect_equal(SRpma_RI_L1_AR0_check$par, SRpma_RI_L1_AR0$par)
  #RI.L1.AR1
  expect_equal(SRpma_RI_L1_AR1_outAR_F_check$par, SRpma_RI_L1_AR1$par)
  #RI.L2.AR0
  expect_equal(SRpma_RI_L2_AR0_check$par, SRpma_RI_L2_AR0$par)
  #RI.L2.AR1
  expect_equal(SRpma_RI_L2_AR1_outAR_F_check$par, SRpma_RI_L2_AR1$par)

  #結果の数値を照合($resid)----
  #HS.L1.AR0
  expect_equal(SRpma_HS_L1_AR0_check$resid, SRpma_HS_L1_AR0$resid)
  #HS.L1.AR1
  expect_equal(SRpma_HS_L1_AR1_outAR_F_check$resid, SRpma_HS_L1_AR1$resid)
  #HS.L2.AR0
  expect_equal(SRpma_HS_L2_AR0_check$resid, SRpma_HS_L2_AR0$resid)
  #HS.L2.AR1
  expect_equal(SRpma_HS_L2_AR1_outAR_F_check$resid, SRpma_HS_L2_AR1$resid)

  #BH.L1.AR0
  expect_equal(SRpma_BH_L1_AR0_check$resid, SRpma_BH_L1_AR0$resid)
  #BH.L1.AR1
  expect_equal(SRpma_BH_L1_AR1_outAR_F_check$resid, SRpma_BH_L1_AR1$resid)
  #BH.L2.AR0
  expect_equal(SRpma_BH_L2_AR0_check$resid, SRpma_BH_L2_AR0$resid)
  #BH.L2.AR1
  expect_equal(SRpma_BH_L2_AR1_outAR_F_check$resid, SRpma_BH_L2_AR1$resid)

  #RI.L1.AR0
  expect_equal(SRpma_RI_L1_AR0_check$resid, SRpma_RI_L1_AR0$resid)
  #RI.L1.AR1
  expect_equal(SRpma_RI_L1_AR1_outAR_F_check$resid, SRpma_RI_L1_AR1$resid)
  #RI.L2.AR0
  expect_equal(SRpma_RI_L2_AR0_check$resid, SRpma_RI_L2_AR0$resid)
  #RI.L2.AR1
  expect_equal(SRpma_RI_L2_AR1_outAR_F_check$resid, SRpma_RI_L2_AR1$resid)


  #結果の数値を照合($resid2)----
  #HS.L1.AR0
  expect_equal(SRpma_HS_L1_AR0_check$resid2, SRpma_HS_L1_AR0$resid2)
  #HS.L1.AR1
  expect_equal(SRpma_HS_L1_AR1_outAR_F_check$resid2, SRpma_HS_L1_AR1$resid2)
  #HS.L2.AR0
  expect_equal(SRpma_HS_L2_AR0_check$resid2, SRpma_HS_L2_AR0$resid2)
  #HS.L2.AR1
  expect_equal(SRpma_HS_L2_AR1_outAR_F_check$resid2, SRpma_HS_L2_AR1$resid2)

  #BH.L1.AR0
  expect_equal(SRpma_BH_L1_AR0_check$resid2, SRpma_BH_L1_AR0$resid2)
  #BH.L1.AR1
  expect_equal(SRpma_BH_L1_AR1_outAR_F_check$resid2, SRpma_BH_L1_AR1$resid2)
  #BH.L2.AR0
  expect_equal(SRpma_BH_L2_AR0_check$resid2, SRpma_BH_L2_AR0$resid2)
  #BH.L2.AR1
  expect_equal(SRpma_BH_L2_AR1_outAR_F_check$resid2, SRpma_BH_L2_AR1$resid2)

  #RI.L1.AR0
  expect_equal(SRpma_RI_L1_AR0_check$resid2, SRpma_RI_L1_AR0$resid2)
  #RI.L1.AR1
  expect_equal(SRpma_RI_L1_AR1_outAR_F_check$resid2, SRpma_RI_L1_AR1$resid2)
  #RI.L2.AR0
  expect_equal(SRpma_RI_L2_AR0_check$resid2, SRpma_RI_L2_AR0$resid2)
  #RI.L2.AR1
  expect_equal(SRpma_RI_L2_AR1_outAR_F_check$resid2, SRpma_RI_L2_AR1$resid2)


  #結果の数値を照合($opt)----
  #HS.L1.AR0
  expect_equal(SRpma_HS_L1_AR0_check$opt, SRpma_HS_L1_AR0$opt)
  #HS.L1.AR1
  expect_equal(SRpma_HS_L1_AR1_outAR_F_check$opt, SRpma_HS_L1_AR1$opt)
  #HS.L2.AR0
  expect_equal(SRpma_HS_L2_AR0_check$opt, SRpma_HS_L2_AR0$opt)
  #HS.L2.AR1
  expect_equal(SRpma_HS_L2_AR1_outAR_F_check$opt, SRpma_HS_L2_AR1$opt)

  #BH.L1.AR0
  expect_equal(SRpma_BH_L1_AR0_check$opt, SRpma_BH_L1_AR0$opt)
  #BH.L1.AR1
  expect_equal(SRpma_BH_L1_AR1_outAR_F_check$opt, SRpma_BH_L1_AR1$opt)
  #BH.L2.AR0
  expect_equal(SRpma_BH_L2_AR0_check$opt, SRpma_BH_L2_AR0$opt)
  #BH.L2.AR1
  expect_equal(SRpma_BH_L2_AR1_outAR_F_check$opt, SRpma_BH_L2_AR1$opt)

  #RI.L1.AR0
  expect_equal(SRpma_RI_L1_AR0_check$opt, SRpma_RI_L1_AR0$opt)
  #RI.L1.AR1
  expect_equal(SRpma_RI_L1_AR1_outAR_F_check$opt, SRpma_RI_L1_AR1$opt)
  #RI.L2.AR0
  expect_equal(SRpma_RI_L2_AR0_check$opt, SRpma_RI_L2_AR0$opt)
  #RI.L2.AR1
  expect_equal(SRpma_RI_L2_AR1_outAR_F_check$opt, SRpma_RI_L2_AR1$opt)


})

