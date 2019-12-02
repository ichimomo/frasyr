library(frasyr)

context("future ref.F")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  res_ref_f_pma_check <- ref.F(res_pma,Fcurrent=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                     waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                     d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                     iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )

  res.ref.f_2times <- ref.F(res_pma,Fcurrent=res_pma$Fc.at.age*2,
                            waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                     waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                     d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,
                     pSPR = c(seq(10,90,by=10),res_ref_f_pma_check$currentSPR$perSPR*100),
                     iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )
  times2_check <- as.numeric(table(unlist(res_ref_f_pma_check$summary/res.ref.f_2times$summary[,1:16])))
  expect_equal(times2_check[1],2)
  expect_equal(times2_check[2],31)
  expect_equal(times2_check[3],15)
  expect_equal(round(res.ref.f_2times$summary[3,17],2),0.5)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","res_ref_f_pma.rda",package = "frasyr"))

  #結果の数値を照合
  expect_equal(res_ref_f_pma_check, res_ref_f_pma)

  })


context("future future.vpa HS")

test_that("oututput value check (iteration # of future sim is fixed as 2) ",{#このテストはfutrue.rにおいてHS.recAR関数内１行目にset.seed(0)を設定した結果を参照する。デフォルト設定では一致しない。
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_HS_L1_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_HS_L1_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_HS_L2_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_HS_L2_AR1.rda",package = "frasyr"))

  fres_pma_HS_L1_AR0_check <- future.vpa(res0=res_pma,
                        multi=1, # res.pma$Fc.at.ageに掛ける乗数
                        nyear=50, # 将来予測の年数
                        start.year=2012, # 将来予測の開始年
                        N=2, # 確率的計算の繰り返し回数
                        ABC.year=2013, # ABCを計算する年
                        waa.year=2009:2011, # 生物パラメータの参照年
                        maa.year=2009:2011,
                        M.year=2009:2011,
                        is.plot=TRUE, # 結果をプロットするかどうか
                        seed=1,
                        silent=TRUE,
                        recfunc=HS.recAR, # 再生産関係の関数
                        # recfuncに対する引数
                        rec.arg=list(a=SRpma_HS_L1_AR0$pars$a,b=SRpma_HS_L1_AR0$pars$b,
                                     rho=SRpma_HS_L1_AR0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                     sd=SRpma_HS_L1_AR0$pars$sd,resid=SRpma_HS_L1_AR0$resid)
                        )

  fres_pma_HS_L2_AR0_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=HS.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_HS_L2_AR0$pars$a,b=SRpma_HS_L2_AR0$pars$b,
                                                      rho=SRpma_HS_L2_AR0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_HS_L2_AR0$pars$sd,resid=SRpma_HS_L2_AR0$resid)
  )

  fres_pma_HS_L1_AR1_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=HS.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_HS_L1_AR1$pars$a,b=SRpma_HS_L1_AR1$pars$b,
                                                      rho=SRpma_HS_L1_AR1$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_HS_L1_AR1$pars$sd,resid=SRpma_HS_L1_AR1$resid)
  )

  fres_pma_HS_L2_AR1_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=HS.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_HS_L2_AR1$pars$a,b=SRpma_HS_L2_AR1$pars$b,
                                                      rho=SRpma_HS_L2_AR1$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_HS_L2_AR1$pars$sd,resid=SRpma_HS_L2_AR1$resid)
  )

  #上記引数での計算結果を読み込み
  # HS L1
  load(system.file("extdata","fres_pma_HS_L1_AR0.rda",package = "frasyr"))
  # HS L2
  load(system.file("extdata","fres_pma_HS_L2_AR0.rda",package = "frasyr"))

  # HS.AR L1
  load(system.file("extdata","fres_pma_HS_L1_AR1.rda",package = "frasyr"))
  # HS.AR L2
  load(system.file("extdata","fres_pma_HS_L2_AR1.rda",package = "frasyr"))

  #結果の数値を照合
  # HS L1
  #expect_equal(fres_pma_HS_L1_AR0, fres_pma_HS_L1_AR0_check)
  # HS L2
  #expect_equal(fres_pma_HS_L2_AR0, fres_pma_HS_L2_AR0_check)

  # HS.AR L1
  #expect_equal(fres_pma_HS_L1_AR1, fres_pma_HS_L1_AR1_check)
  # HS.AR L2
  #expect_equal(fres_pma_HS_L2_AR1, fres_pma_HS_L2_AR1_check)

})

context("future future.vpa BH")

test_that("oututput value check (iteration # of future sim is fixed as 2) ",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L1_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L1_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L2_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_BH_L2_AR1.rda",package = "frasyr"))

  fres_pma_BH_L1_AR0_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=BH.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_BH_L1_AR0$pars$a,b=SRpma_BH_L1_AR0$pars$b,
                                                      rho=SRpma_BH_L1_AR0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_BH_L1_AR0$pars$sd,resid=SRpma_BH_L1_AR0$resid)
  )

  fres_pma_BH_L2_AR0_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=BH.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_BH_L2_AR0$pars$a,b=SRpma_BH_L2_AR0$pars$b,
                                                      rho=SRpma_BH_L2_AR0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_BH_L2_AR0$pars$sd,resid=SRpma_BH_L2_AR0$resid)
  )

  fres_pma_BH_L1_AR1_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=BH.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_BH_L1_AR1$pars$a,b=SRpma_BH_L1_AR1$pars$b,
                                                      rho=SRpma_BH_L1_AR1$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_BH_L1_AR1$pars$sd,resid=SRpma_BH_L1_AR1$resid)
  )

  fres_pma_BH_L2_AR1_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=BH.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_BH_L2_AR1$pars$a,b=SRpma_BH_L2_AR1$pars$b,
                                                      rho=SRpma_BH_L2_AR1$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_BH_L2_AR1$pars$sd,resid=SRpma_BH_L2_AR1$resid)
  )

  #上記引数での計算結果を読み込み
  # BH L1
  load(system.file("extdata","fres_pma_BH_L1_AR0.rda",package = "frasyr"))
  # BH L2
  load(system.file("extdata","fres_pma_BH_L2_AR0.rda",package = "frasyr"))

  # BH AR L1
  load(system.file("extdata","fres_pma_BH_L1_AR1.rda",package = "frasyr"))
  # BH AR L2
  load(system.file("extdata","fres_pma_BH_L2_AR1.rda",package = "frasyr"))

  #結果の数値を照合
  # BH L1
  #expect_equal(fres_pma_BH_L1_AR0, fres_pma_BH_L1_AR0_check)
  # BH L2
  #expect_equal(fres_pma_BH_L2_AR0, fres_pma_BH_L2_AR0_check)

  # BH AR L1
  #expect_equal(fres_pma_BH_L1_AR1, fres_pma_BH_L1_AR1_check)
  # BH L2
  #expect_equal(fres_pma_BH_L2_AR1, fres_pma_BH_L2_AR1_check)

})

context("future future.vpa RI")

test_that("oututput value check (iteration # of future sim is fixed as 2) ",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L1_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L1_AR1.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L2_AR0.rda",package = "frasyr"))
  load(system.file("extdata","SRpma_RI_L2_AR1.rda",package = "frasyr"))

  fres_pma_RI_L1_AR0_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=RI.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_RI_L1_AR0$pars$a,b=SRpma_RI_L1_AR0$pars$b,
                                                      rho=SRpma_RI_L1_AR0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_RI_L1_AR0$pars$sd,resid=SRpma_RI_L1_AR0$resid)
  )

  fres_pma_RI_L2_AR0_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=RI.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_RI_L2_AR0$pars$a,b=SRpma_RI_L2_AR0$pars$b,
                                                      rho=SRpma_RI_L2_AR0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_RI_L2_AR0$pars$sd,resid=SRpma_RI_L2_AR0$resid)
  )

  fres_pma_RI_L1_AR1_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=RI.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_RI_L1_AR1$pars$a,b=SRpma_RI_L1_AR1$pars$b,
                                                      rho=SRpma_RI_L1_AR1$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_RI_L1_AR1$pars$sd,resid=SRpma_RI_L1_AR1$resid)
  )

  fres_pma_RI_L2_AR1_check <- future.vpa(res0=res_pma,
                                         multi=1, # res.pma$Fc.at.ageに掛ける乗数
                                         nyear=50, # 将来予測の年数
                                         start.year=2012, # 将来予測の開始年
                                         N=2, # 確率的計算の繰り返し回数
                                         ABC.year=2013, # ABCを計算する年
                                         waa.year=2009:2011, # 生物パラメータの参照年
                                         maa.year=2009:2011,
                                         M.year=2009:2011,
                                         is.plot=TRUE, # 結果をプロットするかどうか
                                         seed=1,
                                         silent=TRUE,
                                         recfunc=RI.recAR, # 再生産関係の関数
                                         # recfuncに対する引数
                                         rec.arg=list(a=SRpma_RI_L2_AR1$pars$a,b=SRpma_RI_L2_AR1$pars$b,
                                                      rho=SRpma_RI_L2_AR1$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                      sd=SRpma_RI_L2_AR1$pars$sd,resid=SRpma_RI_L2_AR1$resid)
  )

  #上記引数での計算結果を読み込み
  # RI L1
  load(system.file("extdata","fres_pma_RI_L1_AR0.rda",package = "frasyr"))
  # RI L2
  load(system.file("extdata","fres_pma_RI_L2_AR0.rda",package = "frasyr"))

  # RI AR L1
  load(system.file("extdata","fres_pma_RI_L1_AR1.rda",package = "frasyr"))
  # RI AR L2
  load(system.file("extdata","fres_pma_RI_L2_AR1.rda",package = "frasyr"))

  #結果の数値を照合
  # RI L1
  #expect_equal(fres_pma_RI_L1_AR0, fres_pma_RI_L1_AR0_check)
  # RI L2
  #expect_equal(fres_pma_RI_L2_AR0, fres_pma_RI_L2_AR0_check)

  # RI AR L1
  #expect_equal(fres_pma_RI_L1_AR1, fres_pma_RI_L1_AR1_check)
  # RI L2
  #expect_equal(fres_pma_RI_L2_AR1, fres_pma_RI_L2_AR1_check)

})

context("future future.vpa (option of futureF)")

test_that("oututput value check (iteration # of future sim is fixed as 2) ",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  SRdata <- get.SRdata(res.pma)

  HS.par0 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
  HS.par1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,hessian=FALSE)

  currentF.test <- 1:4/10
  futureF.test <- 5:8/10
  fres.HS.check <- future.vpa(res.pma,
                              multi=2,
                              currentF=currentF.test,
                              futureF=5:8/10,
                              nyear=50,
                              start.year=2012,
                              N=2, ABC.year=2013,
                        waa.year=2009:2011,
                        maa.year=2009:2011,
                        M.year=2009:2011,
                        is.plot=TRUE,
                        seed=1,
                        silent=TRUE,
                        recfunc=HS.recAR,
                        rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                     rho=HS.par0$pars$rho,
                                     sd=HS.par0$pars$sd,resid=HS.par0$resid)
                        )


  # faaが想定通りに入っていればOKくらいのテストです
  for(i in 1:4){
      expect_true(mean(fres.HS.check$faa[i,dimnames(fres.HS.check$faa)[[2]]==2012,])==currentF.test[i])
      expect_true(mean(fres.HS.check$faa[i,dimnames(fres.HS.check$faa)[[2]]>2012,])==futureF.test[i]*2)
  }
})


