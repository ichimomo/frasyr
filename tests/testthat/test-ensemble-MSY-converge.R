library(frasyr)

context("check future est MSY ensemble")

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
                                        nyear=58, # 将来予測の年数
                                        start.year=2012, # 将来予測の開始年
                                        N=1000, # 確率的計算の繰り返し回数
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
  #このテストでは future.rのest.MSY関数内でfarg.tmp$det.run <- TRUE としたうえで実行
  #ただし、farg.tmp$det.run <- FALSEでもres_MSY_pma_check$summaryとres_MSY_pma_check$F.msyはN=999 (future2.1.rではN=1000)でテストクリアするはず
  res_MSY_pma_check <- est.MSY(res_vpa_pma, # VPAの計算結果
                               res_future_Fcurrent_pma$input, # 将来予測で使用した引数
                               seed=res_future_Fcurrent_pma$input$seed,
                               N=999, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
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

  for(n in 2:4){
    sample.size <- 10^n
    filename <- sprintf("res_MSY_pma_mean_samplesize%d.rda",sample.size)
    load(system.file("extdata",filename,package = "frasyr"))
    assign(sprintf("res_MSY_pma_mean_samplesize%d",sample.size),res_MSY_pma_mean)
    filename <- sprintf("res_MSY_pma_sd_samplesize%d.rda",sample.size)
    load(system.file("extdata",filename,package = "frasyr"))
    assign(sprintf("res_MSY_pma_sd_samplesize%d",sample.size),res_MSY_pma_sd)
  }

  # testで収束をチェックするのはF.msyのみ
  # stochastic sim iteration = nsim(100,1000,10000)につき、独立に100回シミュレーションし、grand.meanとそのsdをgenerate_res_MSY_for_ensemble.Rで求めたので、grand.meanに対して3σにテストデータの結果が入るかをチェック

  check <- as.data.frame(eval(parse(text=paste("res_MSY_pma_check$F.msy"))))
  grand.mean <- as.data.frame(eval(parse(text=paste("res_MSY_pma_mean_samplesize100$F.msy"))))
  grand.mean.sd <- as.data.frame(eval(parse(text=paste("res_MSY_pma_sd_samplesize100$F.msy"))))
  for(i in nrow(check)){
    expect_equal(check[i,1], grand.mean[i,1], tolerance=3*grand.mean.sd[i,1])
  }

  check <- as.data.frame(eval(parse(text=paste("res_MSY_pma_check$F.msy"))))
  grand.mean <- as.data.frame(eval(parse(text=paste("res_MSY_pma_mean_samplesize1000$F.msy"))))
  grand.mean.sd <- as.data.frame(eval(parse(text=paste("res_MSY_pma_sd_samplesize1000$F.msy"))))
  for(i in nrow(check)){
    expect_equal(check[i,1], grand.mean[i,1], tolerance=3*grand.mean.sd[i,1])
  }

  check <- as.data.frame(eval(parse(text=paste("res_MSY_pma_check$F.msy"))))
  grand.mean <- as.data.frame(eval(parse(text=paste("res_MSY_pma_mean_samplesize10000$F.msy"))))
  grand.mean.sd <- as.data.frame(eval(parse(text=paste("res_MSY_pma_sd_samplesize10000$F.msy"))))
  for(i in nrow(check)){
    expect_equal(check[i,1], grand.mean[i,1], tolerance=3*grand.mean.sd[i,1])
  }

  })
