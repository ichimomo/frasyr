library(frasyr)

context("future ref.F")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  res_ref_f_pma_check <- ref.F(res_vpa_pma,Fcurrent=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                               waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                               d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                               iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )

  res.ref.f_2times <- ref.F(res_vpa_pma,Fcurrent=res_vpa_pma$Fc.at.age*2,
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

  #照合内容
  testcontents <-c("sel","max.age","min.age","rps.q","spr.q","Fcurrent","Fmed","Flow","Fhigh","Fmax","F0.1","Fmean","rps.data","FpSPR","summary","ypr.spr","waa","waa.catch","maa","spr0")

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("res_ref_f_pma$",testcontents[i]))),eval(parse(text=paste("res_ref_f_pma_check$",testcontents[i]))))
  }

})

context("future SRdata")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
  SRdata0usingPeriodFrom1990To2000 <- get.SRdata(res_vpa_pma,years=1990:2000)

  #上記引数での計算結果を読み込み
  SRdata0_pma_check <- read.csv(system.file("extdata","future_SRdata0_pma_check.csv",package="frasyr"),row.names=1)
  SRdata0usingPeriodFrom1990To2000_pma_check <- read.csv(system.file("extdata","future_SRdata0usingPeriodFrom1990To2000_pma_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合

  expect_equal(SRdata0$year, SRdata0_pma_check$year)
  #expect_equal(SRdata0$SSB, SRdata0_pma_check$SSB)
  #expect_equal(SRdata0$R, SRdata0_pma_check$R)

  expect_equal(SRdata0usingPeriodFrom1990To2000$year, SRdata0usingPeriodFrom1990To2000_pma_check$year)
  expect_equal(SRdata0usingPeriodFrom1990To2000$SSB, SRdata0usingPeriodFrom1990To2000_pma_check$SSB)
  expect_equal(SRdata0usingPeriodFrom1990To2000$R, SRdata0usingPeriodFrom1990To2000_pma_check$R)

})

context("future future.vpa")

test_that("oututput value check (iteration for future sim is fixed as 2) ",{

  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
  SR.list <- list()

  for (i in 1:nrow(SRmodel.list)) {
    SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                           AR = SRmodel.list$AR.type[i], out.AR=SRmodel.list$out.AR[i], hessian = FALSE)
  }

  for (i in 1:nrow(SRmodel.list)) {
    infile.name <- sprintf("SRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    resfitSR <- load(system.file("extdata",infile.name,package = "frasyr"))

    fittedSR <- eval(parse(text=resfitSR))
    fres_pma_recarg_list <-list(a=fittedSR$pars$a,b=fittedSR$pars$b,
                                rho=fittedSR$pars$rho, # ここではrho=0なので指定しなくてもOK
                                sd=fittedSR$pars$sd,resid=fittedSR$resid)

    selectedrecSR <- sprintf("%s.recAR",fittedSR$input$SR[1])

    fres_pma_check <-future.vpa(res0=res_vpa_pma,
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
                                recfunc=eval(parse(text=selectedrecSR)), # 再生産関係の関数
                                # recfuncに対する引数
                                rec.arg=fres_pma_recarg_list)

    assign(sprintf("fres_pma_%s_%s_AR%d_outAR%d_check",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),fres_pma_check)

  }

  #照合内容
  testcontents <-c("faa","naa","wcaa","M","maa","vbiom","recruit","eaa","alpha","thisyear.ssb","waa","waa.catch","currentF","vssb","vwcaa","naa_all","years","fyear.year","ABC","waa.year","maa.year","multi","multi.year")

  # HS L1 AR0 ----
  load(system.file("extdata","fres_pma_HS_L1_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_HS_L1_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_HS_L1_AR0_outAR0_check$",testcontents[i]))), label=testcontents[i])
  }

  # HS L1 AR1 outAR F(0) ----
  load(system.file("extdata","fres_pma_HS_L1_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
      expect_equal(eval(parse(text=paste("fres_pma_HS_L1_AR1_outAR0$",testcontents[i]))),
                   eval(parse(text=paste("fres_pma_HS_L1_AR1_outAR0_check$",testcontents[i]))),
                   label=testcontents[i])
  }

  # HS L1 AR1 outAR T(1) ----
  load(system.file("extdata","fres_pma_HS_L1_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_HS_L1_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("fres_pma_HS_L1_AR1_outAR1_check$",testcontents[i]))))
  }

  # HS L2 AR0 ----
  load(system.file("extdata","fres_pma_HS_L2_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_HS_L2_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_HS_L2_AR0_outAR0_check$",testcontents[i]))))
  }
  # HS L2 AR1 outAR F(0) ----
  load(system.file("extdata","fres_pma_HS_L2_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_HS_L2_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_HS_L2_AR1_outAR0_check$",testcontents[i]))))
  }
  # HS L2 AR1 outAR T(1) ----
  load(system.file("extdata","fres_pma_HS_L2_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_HS_L2_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("fres_pma_HS_L2_AR1_outAR1_check$",testcontents[i]))))
  }


  # BH L1 AR0 ----
  load(system.file("extdata","fres_pma_BH_L1_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_BH_L1_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_BH_L1_AR0_outAR0_check$",testcontents[i]))))
  }

  # BH L1 AR1 outAR F(0) ----
  load(system.file("extdata","fres_pma_BH_L1_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_BH_L1_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_BH_L1_AR1_outAR0_check$",testcontents[i]))))
  }

  # BH L1 AR1 outAR T(1) ----
  load(system.file("extdata","fres_pma_BH_L1_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_BH_L1_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("fres_pma_BH_L1_AR1_outAR1_check$",testcontents[i]))))
  }

  # BH L2 AR0 ----
  load(system.file("extdata","fres_pma_BH_L2_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_BH_L2_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_BH_L2_AR0_outAR0_check$",testcontents[i]))))
  }
  # BH L2 AR1 outAR F(0) ----
  load(system.file("extdata","fres_pma_BH_L2_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_BH_L2_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_BH_L2_AR1_outAR0_check$",testcontents[i]))))
  }
  # BH L2 AR1 outAR T(1) ----
  load(system.file("extdata","fres_pma_BH_L2_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_BH_L2_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("fres_pma_BH_L2_AR1_outAR1_check$",testcontents[i]))))
  }




  # RI L1 AR0 ----
  load(system.file("extdata","fres_pma_RI_L1_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_RI_L1_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_RI_L1_AR0_outAR0_check$",testcontents[i]))))
  }

  # RI L1 AR1 outAR F(0) ----
  load(system.file("extdata","fres_pma_RI_L1_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_RI_L1_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_RI_L1_AR1_outAR0_check$",testcontents[i]))))
  }

  # RI L1 AR1 outAR T(1) ----
  load(system.file("extdata","fres_pma_RI_L1_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_RI_L1_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("fres_pma_RI_L1_AR1_outAR1_check$",testcontents[i]))))
  }

  # RI L2 AR0 ----
  load(system.file("extdata","fres_pma_RI_L2_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_RI_L2_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_RI_L2_AR0_outAR0_check$",testcontents[i]))))
  }
  # RI L2 AR1 outAR F(0) ----
  load(system.file("extdata","fres_pma_RI_L2_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_RI_L2_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("fres_pma_RI_L2_AR1_outAR0_check$",testcontents[i]))))
  }
  # RI L2 AR1 outAR T(1) ----
  load(system.file("extdata","fres_pma_RI_L2_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("fres_pma_RI_L2_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("fres_pma_RI_L2_AR1_outAR1_check$",testcontents[i]))))
  }



     })

context("future est MSY")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
  SR.list <- list()

  for (i in 1:nrow(SRmodel.list)) {
    SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                           AR = SRmodel.list$AR.type[i], out.AR =SRmodel.list$out.AR[i], hessian = FALSE)
  }

  SRmodel.list$AICc <- sapply(SR.list, function(x) x$AICc)
  SRmodel.list$delta.AIC <- SRmodel.list$AICc - min(SRmodel.list$AICc)
  SR.list <- SR.list[order(SRmodel.list$AICc)]  # AICの小さい順に並べたもの
  (SRmodel.list <- SRmodel.list[order(SRmodel.list$AICc), ]) # 結果

  SRmodel.base <- SR.list[[1]] # AIC最小モデルを今後使っていく

  selectedSR <- sprintf("%s.recAR",SRmodel.base$input$SR[1])
 # future fcurrent ----
  res_future_Fcurrent_pma <- future.vpa(res_vpa_pma,
                                        multi=1,
                                        nyear=50, # 将来予測の年数
                                        start.year=2012, # 将来予測の開始年
                                        N=100, # 確率的計算の繰り返し回数
                                        ABC.year=2013, # ABCを計算する年
                                        waa.year=2009:2011, # 生物パラメータの参照年
                                        maa.year=2009:2011,
                                        M.year=2009:2011,
                                        is.plot=FALSE, # 結果をプロットするかどうか
                                        seed=1,
                                        silent=TRUE,
                                        recfunc=eval(parse(text=selectedSR))
                                        , # 再生産関係の関数
                                        # recfuncに対する引数
                                        rec.arg=list(a=SRmodel.base$pars$a,b=SRmodel.base$pars$b,
                                                     rho=SRmodel.base$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                     sd=SRmodel.base$pars$sd,resid=SRmodel.base$resid))

  # est MSY ----
  nyear <- round(Generation.Time(res_vpa,
                                 maa.year=2009:2011,
                                 M.year=2009:2011)*20)
  res_MSY_pma_check <- est.MSY(res_vpa_pma, # VPAの計算結果
                         res_future_Fcurrent_pma$input, # 将来予測で使用した引数
                         seed=res_future_Fcurrent_pma$input$seed,
                         N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                         calc.yieldcurve=TRUE,
                         PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                         onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                         B0percent=c(0.2,0.3,0.4),
                         nyear=nyear,
                         Bempirical=c(round(tail(colSums(res_vpa_pma$ssb),n=1)),
                                      round(max(colSums(res_vpa_pma$ssb))),
                                      24000, # 現行Blimit
                                      SRmodel.base$pars$b) # HSの折れ点
  )

  # 上記設定の結果を読み込み ----
  load(system.file("extdata","res_MSY_pma.rda",package = "frasyr"))
  load(system.file("extdata","res_MSY_pma_pre.rda",package = "frasyr"))

  #照合内容
  #  testcontents <-c("summary","summaryAR","summary_tb","F.msy","all.stat","all.statAR","all.stat_tb","trace","ssb.ar.mean","SPR.msy")
  testcontents <-c("summary","summary_tb","F.msy","all.stat","trace")

  #読み込んだ結果と照合 future2.1.r + utility.r(future-vpa ver.)との比較
  for(i in 1:length(testcontents)){
      expect_equal(eval(parse(text=paste("res_MSY_pma$",testcontents[i]))),
                   eval(parse(text=paste("res_MSY_pma_check$",testcontents[i]))),
                   label=testcontents[i])
  }
  #読み込んだ結果と照合 future.r + utility.r(2019/12/18以前ver.)との比較
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("res_MSY_pma_pre$",testcontents[i]))),
                 eval(parse(text=paste("res_MSY_pma_check$",testcontents[i]))),
                 label=testcontents[i])
  }

})


context("future future.vpa (option of futureF)")

test_that("oututput value check (iteration for future sim is fixed as 2) ",{
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


