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

context("future SRdata")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  SRdata_pma_check <- get.SRdata(res_pma)
  SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
  SRdata0usingPeriodFrom1990To2000 <- get.SRdata(res_pma,years=1990:2000)

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

test_that("oututput value check (iteration for future sim is fixed as 2) ",{#デフォルト設定では将来予測でHS.recAR関数の乱数生成により一致しない。

  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

for(i in 1:3){ # SR function type
  if(i==1) SRtype = "HS"
  if(i==2) SRtype = "BH"
  if(i==3) SRtype = "RI"
  for(j in 1:2){ # Estimation method
    if(j==1) Est="L1"
    if(j==2) Est="L2"
    for(k in 1:2){ # AR option
      if(k==1) ARopt="AR0"
      if(k==2) ARopt="AR1"
      if(k==3) ARopt="AR1_outAR_F"

      infile.name <- sprintf("SRpma_%s_%s_%s.rda",SRtype,Est,ARopt)
      resfitSR <- load(system.file("extdata",infile.name,package = "frasyr"))
      fittedSR <- eval(parse(text=resfitSR)
      )

      fres_pma_recarg_list <-list(a=fittedSR$pars$a,b=fittedSR$pars$b,
                                  rho=fittedSR$pars$rho, # ここではrho=0なので指定しなくてもOK
                                  sd=fittedSR$pars$sd,resid=fittedSR$resid)

      fres_pma_check <-future.vpa(res0=res_pma,
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
                                  rec.arg=fres_pma_recarg_list)

      checkfile.name <- sprintf("fres_pma_%s_%s_%s.rda",SRtype,Est,ARopt)
      resfuturevpa <- load(system.file("extdata",checkfile.name,package = "frasyr"))
      fres_pma <- eval(parse(text=resfuturevpa))

      expect_equal(fres_pma$faa, fres_pma_check$faa)
      #expect_equal(fres_pma$naa, fres_pma_check$naa)
      #expect_equal(fres_pma$biom, fres_pma_check$biom)
      #expect_equal(fres_pma$baa, fres_pma_check$baa)
      #expect_equal(fres_pma$ssb, fres_pma_check$ssb)
      #expect_equal(fres_pma$wcaa, fres_pma_check$wcaa)
      #expect_equal(fres_pma$caa, fres_pma_check$caa)
      expect_equal(fres_pma$M, fres_pma_check$M)
      #expect_equal(fres_pma$rps, fres_pma_check$rps)
      #expect_equal(fres_pma$recruit, fres_pma_check$recurit)
      expect_equal(fres_pma$maa, fres_pma_check$maa)
      #expect_equal(fres_pma$vbiom, fres_pma_check$vbiom)
      #expect_equal(fres_pma$eaa, fres_pma_check$eaa)
      expect_equal(fres_pma$alpha, fres_pma_check$alpha)
      #expect_equal(fres_pma$thisyear.ssb, fres_pma_check$thisyear.ssb)
      expect_equal(fres_pma$waa, fres_pma_check$waa)
      expect_equal(fres_pma$waa.catch, fres_pma_check$waa.catch)
      expect_equal(fres_pma$currentF, fres_pma_check$currentF)
      expect_equal(fres_pma$futureF, fres_pma_check$futureF)
      #expect_equal(fres_pma$vssb, fres_pma_check$vssb)
      #expect_equal(fres_pma$vwcaa, fres_pma_check$vwcaa)
      #expect_equal(fres_pma$naa_all, fres_pma_check$naa_all)
      expect_equal(fres_pma$years, fres_pma_check$years)
      expect_equal(fres_pma$fyear.year, fres_pma_check$fyear.year)
      #expect_equal(fres_pma$ABC, fres_pma_check$ABC)
      #expect_equal(fres_pma$recfunc, fres_pma_check$recfunc)
      #expect_equal(fres_pma$rec.arg, fres_pma_check$rec.arg)
      expect_equal(fres_pma$waa.year, fres_pma_check$waa.year)
      expect_equal(fres_pma$maa.year, fres_pma_check$maa.year)
      expect_equal(fres_pma$multi, fres_pma_check$multi)
      expect_equal(fres_pma$multi.year, fres_pma_check$multi.year)
      expect_equal(fres_pma$Frec, fres_pma_check$Frec)
      expect_equal(fres_pma$rec.new, fres_pma_check$rec.new)
      expect_equal(fres_pma$pre.catch, fres_pma_check$pre.catch)
      #expect_equal(fres_pma$input, fres_pma_check$input)
    }
  }
}

})

context("future est MSY")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), L.type = c("L1", "L2"))
SR.list <- list()
for (i in 1:nrow(SRmodel.list)) {
  SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                         AR = SRmodel.list$AR.type[i], hessian = FALSE)
}

SRmodel.list$AICc <- sapply(SR.list, function(x) x$AICc)
SRmodel.list$delta.AIC <- SRmodel.list$AICc - min(SRmodel.list$AICc)
SR.list <- SR.list[order(SRmodel.list$AICc)]  # AICの小さい順に並べたもの
(SRmodel.list <- SRmodel.list[order(SRmodel.list$AICc), ]) # 結果

SRmodel.base <- SR.list[[1]] # AIC最小モデルを今後使っていく

selectedSR <- sprintf("%s.recAR",SRmodel.base$input$SR[1])

res_future_Fcurrent_pma <- future.vpa(res_pma,
                                  multi=1,
                                  nyear=50, # 将来予測の年数
                                  start.year=2012, # 将来予測の開始年
                                  N=100, # 確率的計算の繰り返し回数
                                  ABC.year=2013, # ABCを計算する年
                                  waa.year=2009:2011, # 生物パラメータの参照年
                                  maa.year=2009:2011,
                                  M.year=2009:2011,
                                  is.plot=TRUE, # 結果をプロットするかどうか
                                  seed=1,
                                  silent=TRUE,
                                  recfunc=eval(parse(text=selectedSR))
, # 再生産関係の関数
                                  # recfuncに対する引数
                                  rec.arg=list(a=SRmodel.base$pars$a,b=SRmodel.base$pars$b,
                                               rho=SRmodel.base$pars$rho, # ここではrho=0なので指定しなくてもOK
                                               sd=SRmodel.base$pars$sd,resid=SRmodel.base$resid))

res_MSY_pma_check <- est.MSY(res_pma, # VPAの計算結果
                   res_future_Fcurrent_pma$input, # 将来予測で使用した引数
                   resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                   N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                   calc.yieldcurve=TRUE,
                   PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                   onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                   B0percent=c(0.2,0.3,0.4),
                   Bempirical=c(round(tail(colSums(res_vpa$ssb),n=1)),
                                round(max(colSums(res_vpa$ssb))),
                                24000, # 現行Blimit
                                SRmodel.base$pars$b) # HSの折れ点
)

load(system.file("extdata","res_MSY_pma.rda",package = "frasyr"))

#expect_equal(res_MSY_pma,res_MSY_pma_check)

})


context("future HCR")

test_that("oututput value check",{
  load(system.file("extdata","res_MSY_pma.rda",package = "frasyr"))
  refs_all_pma <- res_MSY_pma$summary

  refs_all_pma$RP.definition[refs_all_pma$RP_name=="B0-20%" & refs_all_pma$AR==FALSE] <- "Btarget1"  # たとえばBtargetの代替値としてB020%も候補に残しておきたい場合
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="PGY_0.95_lower" & refs_all_pma$AR==FALSE] <- "Btarget2"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-19431" & refs_all_pma$AR==FALSE] <- "Bcurrent"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-63967" & refs_all_pma$AR==FALSE] <- "Bmax"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-24000" & refs_all_pma$AR==FALSE] <- "Blimit1"
  refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-51882" & refs_all_pma$AR==FALSE] <- "B_HS"

  refs_all_pma %>% dplyr::select(RP_name,RP.definition)

  refs_base_pma_check <- refs_all_pma %>%
    dplyr::filter(!is.na(RP.definition)) %>% # RP.definitionがNAでないものを抽出
    arrange(desc(SSB)) %>% # SSBを大きい順に並び替え
    dplyr::select(RP.definition,RP_name,SSB,SSB2SSB0,Catch,Catch.CV,U,Fref2Fcurrent)

  load(system.file("extdata","refs_base_pma.rda",package = "frasyr"))

  expect_equal(refs_base_pma,refs_base_pma_check)


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


