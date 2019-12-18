library(frasyr)

context("stock-recruitment SRdata")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  SRdata_pma_check <- get.SRdata(res_vpa_pma)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  #結果の数値を照合
  expect_equal(SRdata_pma_check$year, SRdata_pma$year)
  expect_equal(SRdata_pma_check$SSB, SRdata_pma$SSB)
  expect_equal(SRdata_pma_check$R, SRdata_pma$R)
})

context("stock-recruitment fitSR")

test_that("oututput value check",{
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
  SR.list <- list()
  for (i in 1:nrow(SRmodel.list)){
    SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                           AR = SRmodel.list$AR.type[i], out.AR =SRmodel.list$out.AR[i], hessian = FALSE)
  }

  # SRタイプ、L1L2、自己相関タイプごとに異なるオブジェクトへ格納
  for (i in 1:nrow(SRmodel.list)) {
    assign(sprintf("SRpma_%s_%s_AR%d_outAR%d_check",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),SR.list[[i]])
  }

  # HS L1 AR0
  load(system.file("extdata","SRpma_HS_L1_AR0_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_HS_L1_AR0_outAR0_check$opt$par)){
    expect_equal(SRpma_HS_L1_AR0_outAR0$opt$par[i],SRpma_HS_L1_AR0_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_HS_L1_AR0_outAR0$opt$value,SRpma_HS_L1_AR0_outAR0_check$opt$value)
  for(i in 1:length(SRpma_HS_L1_AR0_outAR0_check$opt$counts)){
    expect_equal(SRpma_HS_L1_AR0_outAR0$opt$counts[i],SRpma_HS_L1_AR0_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_HS_L1_AR0_outAR0$opt$convergence,SRpma_HS_L1_AR0_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_HS_L1_AR0_outAR0_check$resid)){
    expect_equal(SRpma_HS_L1_AR0_outAR0$resid,SRpma_HS_L1_AR0_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_HS_L1_AR0_outAR0_check$resid2)){
    expect_equal(SRpma_HS_L1_AR0_outAR0$resid2,SRpma_HS_L1_AR0_outAR0_check$resid2)
  }
  expect_equal(SRpma_HS_L1_AR0_outAR0$pars,SRpma_HS_L1_AR0_outAR0_check$pars)
  expect_equal(SRpma_HS_L1_AR0_outAR0$loglik,SRpma_HS_L1_AR0_outAR0_check$loglik)
  expect_equal(SRpma_HS_L1_AR0_outAR0$pred,SRpma_HS_L1_AR0_outAR0_check$pred)
  expect_equal(SRpma_HS_L1_AR0_outAR0$k,SRpma_HS_L1_AR0_outAR0_check$k)
  expect_equal(SRpma_HS_L1_AR0_outAR0$AIC,SRpma_HS_L1_AR0_outAR0_check$AIC)
  expect_equal(SRpma_HS_L1_AR0_outAR0$AICc,SRpma_HS_L1_AR0_outAR0_check$AICc)
  expect_equal(SRpma_HS_L1_AR0_outAR0$BIC,SRpma_HS_L1_AR0_outAR0_check$BIC)


  # HS L1 AR1 outAR False
  load(system.file("extdata","SRpma_HS_L1_AR1_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_HS_L1_AR1_outAR0_check$opt$par)){
    expect_equal(SRpma_HS_L1_AR1_outAR0$opt$par[i],SRpma_HS_L1_AR1_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_HS_L1_AR1_outAR0$opt$value,SRpma_HS_L1_AR1_outAR0_check$opt$value)
  for(i in 1:length(SRpma_HS_L1_AR1_outAR0_check$opt$counts)){
    expect_equal(SRpma_HS_L1_AR1_outAR0$opt$counts[i],SRpma_HS_L1_AR1_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_HS_L1_AR1_outAR0$opt$convergence,SRpma_HS_L1_AR1_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_HS_L1_AR1_outAR0_check$resid)){
    expect_equal(SRpma_HS_L1_AR1_outAR0$resid,SRpma_HS_L1_AR1_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_HS_L1_AR1_outAR0_check$resid2)){
    expect_equal(SRpma_HS_L1_AR1_outAR0$resid2,SRpma_HS_L1_AR1_outAR0_check$resid2)
  }
  expect_equal(SRpma_HS_L1_AR1_outAR0$pars,SRpma_HS_L1_AR1_outAR0_check$pars)
  expect_equal(SRpma_HS_L1_AR1_outAR0$loglik,SRpma_HS_L1_AR1_outAR0_check$loglik)
  expect_equal(SRpma_HS_L1_AR1_outAR0$pred,SRpma_HS_L1_AR1_outAR0_check$pred)
  expect_equal(SRpma_HS_L1_AR1_outAR0$k,SRpma_HS_L1_AR1_outAR0_check$k)
  expect_equal(SRpma_HS_L1_AR1_outAR0$AIC,SRpma_HS_L1_AR1_outAR0_check$AIC)
  expect_equal(SRpma_HS_L1_AR1_outAR0$AICc,SRpma_HS_L1_AR1_outAR0_check$AICc)
  expect_equal(SRpma_HS_L1_AR1_outAR0$BIC,SRpma_HS_L1_AR1_outAR0_check$BIC)

  # HS L1 AR1 outAR True
  load(system.file("extdata","SRpma_HS_L1_AR1_outAR1.rda",package = "frasyr"))
  for(i in 1:length(SRpma_HS_L1_AR1_outAR1_check$opt$par)){
    expect_equal(SRpma_HS_L1_AR1_outAR1$opt$par[i],SRpma_HS_L1_AR1_outAR1_check$opt$par[i])
  }
  expect_equal(SRpma_HS_L1_AR1_outAR1$opt$value,SRpma_HS_L1_AR1_outAR1_check$opt$value)
  for(i in 1:length(SRpma_HS_L1_AR1_outAR1_check$opt$counts)){
    expect_equal(SRpma_HS_L1_AR1_outAR1$opt$counts[i],SRpma_HS_L1_AR1_outAR1_check$opt$counts[i])
  }
  expect_equal(SRpma_HS_L1_AR1_outAR1$opt$convergence,SRpma_HS_L1_AR1_outAR1_check$opt$convergence)
  for(i in 1:length(SRpma_HS_L1_AR1_outAR1_check$resid)){
    expect_equal(SRpma_HS_L1_AR1_outAR1$resid,SRpma_HS_L1_AR1_outAR1_check$resid)
  }
  for(i in 1:length(SRpma_HS_L1_AR1_outAR1_check$resid2)){
    expect_equal(SRpma_HS_L1_AR1_outAR1$resid2,SRpma_HS_L1_AR1_outAR1_check$resid2)
  }
  expect_equal(SRpma_HS_L1_AR1_outAR1$pars,SRpma_HS_L1_AR1_outAR1_check$pars)
  expect_equal(SRpma_HS_L1_AR1_outAR1$loglik,SRpma_HS_L1_AR1_outAR1_check$loglik)
  expect_equal(SRpma_HS_L1_AR1_outAR1$pred,SRpma_HS_L1_AR1_outAR1_check$pred)
  expect_equal(SRpma_HS_L1_AR1_outAR1$k,SRpma_HS_L1_AR1_outAR1_check$k)
  expect_equal(SRpma_HS_L1_AR1_outAR1$AIC,SRpma_HS_L1_AR1_outAR1_check$AIC)
  expect_equal(SRpma_HS_L1_AR1_outAR1$AICc,SRpma_HS_L1_AR1_outAR1_check$AICc)
  expect_equal(SRpma_HS_L1_AR1_outAR1$BIC,SRpma_HS_L1_AR1_outAR1_check$BIC)


  # HS L2 AR0
  load(system.file("extdata","SRpma_HS_L2_AR0_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_HS_L2_AR0_outAR0_check$opt$par)){
    expect_equal(SRpma_HS_L2_AR0_outAR0$opt$par[i],SRpma_HS_L2_AR0_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_HS_L2_AR0_outAR0$opt$value,SRpma_HS_L2_AR0_outAR0_check$opt$value)
  for(i in 1:length(SRpma_HS_L2_AR0_outAR0_check$opt$counts)){
    expect_equal(SRpma_HS_L2_AR0_outAR0$opt$counts[i],SRpma_HS_L2_AR0_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_HS_L2_AR0_outAR0$opt$convergence,SRpma_HS_L2_AR0_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_HS_L2_AR0_outAR0_check$resid)){
    expect_equal(SRpma_HS_L2_AR0_outAR0$resid,SRpma_HS_L2_AR0_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_HS_L2_AR0_outAR0_check$resid2)){
    expect_equal(SRpma_HS_L2_AR0_outAR0$resid2,SRpma_HS_L2_AR0_outAR0_check$resid2)
  }
  expect_equal(SRpma_HS_L2_AR0_outAR0$pars,SRpma_HS_L2_AR0_outAR0_check$pars)
  expect_equal(SRpma_HS_L2_AR0_outAR0$loglik,SRpma_HS_L2_AR0_outAR0_check$loglik)
  expect_equal(SRpma_HS_L2_AR0_outAR0$pred,SRpma_HS_L2_AR0_outAR0_check$pred)
  expect_equal(SRpma_HS_L2_AR0_outAR0$k,SRpma_HS_L2_AR0_outAR0_check$k)
  expect_equal(SRpma_HS_L2_AR0_outAR0$AIC,SRpma_HS_L2_AR0_outAR0_check$AIC)
  expect_equal(SRpma_HS_L2_AR0_outAR0$AICc,SRpma_HS_L2_AR0_outAR0_check$AICc)
  expect_equal(SRpma_HS_L2_AR0_outAR0$BIC,SRpma_HS_L2_AR0_outAR0_check$BIC)


  # HS L2 AR1 outAR False
  load(system.file("extdata","SRpma_HS_L2_AR1_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_HS_L2_AR1_outAR0_check$opt$par)){
    expect_equal(SRpma_HS_L2_AR1_outAR0$opt$par[i],SRpma_HS_L2_AR1_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_HS_L2_AR1_outAR0$opt$value,SRpma_HS_L2_AR1_outAR0_check$opt$value)
  for(i in 1:length(SRpma_HS_L2_AR1_outAR0_check$opt$counts)){
    expect_equal(SRpma_HS_L2_AR1_outAR0$opt$counts[i],SRpma_HS_L2_AR1_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_HS_L2_AR1_outAR0$opt$convergence,SRpma_HS_L2_AR1_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_HS_L2_AR1_outAR0_check$resid)){
    expect_equal(SRpma_HS_L2_AR1_outAR0$resid,SRpma_HS_L2_AR1_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_HS_L2_AR1_outAR0_check$resid2)){
    expect_equal(SRpma_HS_L2_AR1_outAR0$resid2,SRpma_HS_L2_AR1_outAR0_check$resid2)
  }
  expect_equal(SRpma_HS_L2_AR1_outAR0$pars,SRpma_HS_L2_AR1_outAR0_check$pars)
  expect_equal(SRpma_HS_L2_AR1_outAR0$loglik,SRpma_HS_L2_AR1_outAR0_check$loglik)
  expect_equal(SRpma_HS_L2_AR1_outAR0$pred,SRpma_HS_L2_AR1_outAR0_check$pred)
  expect_equal(SRpma_HS_L2_AR1_outAR0$k,SRpma_HS_L2_AR1_outAR0_check$k)
  expect_equal(SRpma_HS_L2_AR1_outAR0$AIC,SRpma_HS_L2_AR1_outAR0_check$AIC)
  expect_equal(SRpma_HS_L2_AR1_outAR0$AICc,SRpma_HS_L2_AR1_outAR0_check$AICc)
  expect_equal(SRpma_HS_L2_AR1_outAR0$BIC,SRpma_HS_L2_AR1_outAR0_check$BIC)

  # HS L2 AR1 outAR True
  load(system.file("extdata","SRpma_HS_L2_AR1_outAR1.rda",package = "frasyr"))
  for(i in 1:length(SRpma_HS_L2_AR1_outAR1_check$opt$par)){
    expect_equal(SRpma_HS_L2_AR1_outAR1$opt$par[i],SRpma_HS_L2_AR1_outAR1_check$opt$par[i])
  }
  expect_equal(SRpma_HS_L2_AR1_outAR1$opt$value,SRpma_HS_L2_AR1_outAR1_check$opt$value)
  for(i in 1:length(SRpma_HS_L2_AR1_outAR1_check$opt$counts)){
    expect_equal(SRpma_HS_L2_AR1_outAR1$opt$counts[i],SRpma_HS_L2_AR1_outAR1_check$opt$counts[i])
  }
  expect_equal(SRpma_HS_L2_AR1_outAR1$opt$convergence,SRpma_HS_L2_AR1_outAR1_check$opt$convergence)
  for(i in 1:length(SRpma_HS_L2_AR1_outAR1_check$resid)){
    expect_equal(SRpma_HS_L2_AR1_outAR1$resid,SRpma_HS_L2_AR1_outAR1_check$resid)
  }
  for(i in 1:length(SRpma_HS_L2_AR1_outAR1_check$resid2)){
    expect_equal(SRpma_HS_L2_AR1_outAR1$resid2,SRpma_HS_L2_AR1_outAR1_check$resid2)
  }
  expect_equal(SRpma_HS_L2_AR1_outAR1$pars,SRpma_HS_L2_AR1_outAR1_check$pars)
  expect_equal(SRpma_HS_L2_AR1_outAR1$loglik,SRpma_HS_L2_AR1_outAR1_check$loglik)
  expect_equal(SRpma_HS_L2_AR1_outAR1$pred,SRpma_HS_L2_AR1_outAR1_check$pred)
  expect_equal(SRpma_HS_L2_AR1_outAR1$k,SRpma_HS_L2_AR1_outAR1_check$k)
  expect_equal(SRpma_HS_L2_AR1_outAR1$AIC,SRpma_HS_L2_AR1_outAR1_check$AIC)
  expect_equal(SRpma_HS_L2_AR1_outAR1$AICc,SRpma_HS_L2_AR1_outAR1_check$AICc)
  expect_equal(SRpma_HS_L2_AR1_outAR1$BIC,SRpma_HS_L2_AR1_outAR1_check$BIC)


  # BH L1 AR0
  load(system.file("extdata","SRpma_BH_L1_AR0_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_BH_L1_AR0_outAR0_check$opt$par)){
    expect_equal(SRpma_BH_L1_AR0_outAR0$opt$par[i],SRpma_BH_L1_AR0_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_BH_L1_AR0_outAR0$opt$value,SRpma_BH_L1_AR0_outAR0_check$opt$value)
  for(i in 1:length(SRpma_BH_L1_AR0_outAR0_check$opt$counts)){
    expect_equal(SRpma_BH_L1_AR0_outAR0$opt$counts[i],SRpma_BH_L1_AR0_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_BH_L1_AR0_outAR0$opt$convergence,SRpma_BH_L1_AR0_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_BH_L1_AR0_outAR0_check$resid)){
    expect_equal(SRpma_BH_L1_AR0_outAR0$resid,SRpma_BH_L1_AR0_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_BH_L1_AR0_outAR0_check$resid2)){
    expect_equal(SRpma_BH_L1_AR0_outAR0$resid2,SRpma_BH_L1_AR0_outAR0_check$resid2)
  }
  expect_equal(SRpma_BH_L1_AR0_outAR0$pars,SRpma_BH_L1_AR0_outAR0_check$pars)
  expect_equal(SRpma_BH_L1_AR0_outAR0$loglik,SRpma_BH_L1_AR0_outAR0_check$loglik)
  expect_equal(SRpma_BH_L1_AR0_outAR0$pred,SRpma_BH_L1_AR0_outAR0_check$pred)
  expect_equal(SRpma_BH_L1_AR0_outAR0$k,SRpma_BH_L1_AR0_outAR0_check$k)
  expect_equal(SRpma_BH_L1_AR0_outAR0$AIC,SRpma_BH_L1_AR0_outAR0_check$AIC)
  expect_equal(SRpma_BH_L1_AR0_outAR0$AICc,SRpma_BH_L1_AR0_outAR0_check$AICc)
  expect_equal(SRpma_BH_L1_AR0_outAR0$BIC,SRpma_BH_L1_AR0_outAR0_check$BIC)


  # BH L1 AR1 outAR False
  load(system.file("extdata","SRpma_BH_L1_AR1_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_BH_L1_AR1_outAR0_check$opt$par)){
    expect_equal(SRpma_BH_L1_AR1_outAR0$opt$par[i],SRpma_BH_L1_AR1_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_BH_L1_AR1_outAR0$opt$value,SRpma_BH_L1_AR1_outAR0_check$opt$value)
  for(i in 1:length(SRpma_BH_L1_AR1_outAR0_check$opt$counts)){
    expect_equal(SRpma_BH_L1_AR1_outAR0$opt$counts[i],SRpma_BH_L1_AR1_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_BH_L1_AR1_outAR0$opt$convergence,SRpma_BH_L1_AR1_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_BH_L1_AR1_outAR0_check$resid)){
    expect_equal(SRpma_BH_L1_AR1_outAR0$resid,SRpma_BH_L1_AR1_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_BH_L1_AR1_outAR0_check$resid2)){
    expect_equal(SRpma_BH_L1_AR1_outAR0$resid2,SRpma_BH_L1_AR1_outAR0_check$resid2)
  }
  expect_equal(SRpma_BH_L1_AR1_outAR0$pars,SRpma_BH_L1_AR1_outAR0_check$pars)
  expect_equal(SRpma_BH_L1_AR1_outAR0$loglik,SRpma_BH_L1_AR1_outAR0_check$loglik)
  expect_equal(SRpma_BH_L1_AR1_outAR0$pred,SRpma_BH_L1_AR1_outAR0_check$pred)
  expect_equal(SRpma_BH_L1_AR1_outAR0$k,SRpma_BH_L1_AR1_outAR0_check$k)
  expect_equal(SRpma_BH_L1_AR1_outAR0$AIC,SRpma_BH_L1_AR1_outAR0_check$AIC)
  expect_equal(SRpma_BH_L1_AR1_outAR0$AICc,SRpma_BH_L1_AR1_outAR0_check$AICc)
  expect_equal(SRpma_BH_L1_AR1_outAR0$BIC,SRpma_BH_L1_AR1_outAR0_check$BIC)

  # BH L1 AR1 outAR True
  load(system.file("extdata","SRpma_BH_L1_AR1_outAR1.rda",package = "frasyr"))
  for(i in 1:length(SRpma_BH_L1_AR1_outAR1_check$opt$par)){
    expect_equal(SRpma_BH_L1_AR1_outAR1$opt$par[i],SRpma_BH_L1_AR1_outAR1_check$opt$par[i])
  }
  expect_equal(SRpma_BH_L1_AR1_outAR1$opt$value,SRpma_BH_L1_AR1_outAR1_check$opt$value)
  for(i in 1:length(SRpma_BH_L1_AR1_outAR1_check$opt$counts)){
    expect_equal(SRpma_BH_L1_AR1_outAR1$opt$counts[i],SRpma_BH_L1_AR1_outAR1_check$opt$counts[i])
  }
  expect_equal(SRpma_BH_L1_AR1_outAR1$opt$convergence,SRpma_BH_L1_AR1_outAR1_check$opt$convergence)
  for(i in 1:length(SRpma_BH_L1_AR1_outAR1_check$resid)){
    expect_equal(SRpma_BH_L1_AR1_outAR1$resid,SRpma_BH_L1_AR1_outAR1_check$resid)
  }
  for(i in 1:length(SRpma_BH_L1_AR1_outAR1_check$resid2)){
    expect_equal(SRpma_BH_L1_AR1_outAR1$resid2,SRpma_BH_L1_AR1_outAR1_check$resid2)
  }
  expect_equal(SRpma_BH_L1_AR1_outAR1$pars,SRpma_BH_L1_AR1_outAR1_check$pars)
  expect_equal(SRpma_BH_L1_AR1_outAR1$loglik,SRpma_BH_L1_AR1_outAR1_check$loglik)
  expect_equal(SRpma_BH_L1_AR1_outAR1$pred,SRpma_BH_L1_AR1_outAR1_check$pred)
  expect_equal(SRpma_BH_L1_AR1_outAR1$k,SRpma_BH_L1_AR1_outAR1_check$k)
  expect_equal(SRpma_BH_L1_AR1_outAR1$AIC,SRpma_BH_L1_AR1_outAR1_check$AIC)
  expect_equal(SRpma_BH_L1_AR1_outAR1$AICc,SRpma_BH_L1_AR1_outAR1_check$AICc)
  expect_equal(SRpma_BH_L1_AR1_outAR1$BIC,SRpma_BH_L1_AR1_outAR1_check$BIC)


  # BH L2 AR0
  load(system.file("extdata","SRpma_BH_L2_AR0_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_BH_L2_AR0_outAR0_check$opt$par)){
    expect_equal(SRpma_BH_L2_AR0_outAR0$opt$par[i],SRpma_BH_L2_AR0_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_BH_L2_AR0_outAR0$opt$value,SRpma_BH_L2_AR0_outAR0_check$opt$value)
  for(i in 1:length(SRpma_BH_L2_AR0_outAR0_check$opt$counts)){
    expect_equal(SRpma_BH_L2_AR0_outAR0$opt$counts[i],SRpma_BH_L2_AR0_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_BH_L2_AR0_outAR0$opt$convergence,SRpma_BH_L2_AR0_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_BH_L2_AR0_outAR0_check$resid)){
    expect_equal(SRpma_BH_L2_AR0_outAR0$resid,SRpma_BH_L2_AR0_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_BH_L2_AR0_outAR0_check$resid2)){
    expect_equal(SRpma_BH_L2_AR0_outAR0$resid2,SRpma_BH_L2_AR0_outAR0_check$resid2)
  }
  expect_equal(SRpma_BH_L2_AR0_outAR0$pars,SRpma_BH_L2_AR0_outAR0_check$pars)
  expect_equal(SRpma_BH_L2_AR0_outAR0$loglik,SRpma_BH_L2_AR0_outAR0_check$loglik)
  expect_equal(SRpma_BH_L2_AR0_outAR0$pred,SRpma_BH_L2_AR0_outAR0_check$pred)
  expect_equal(SRpma_BH_L2_AR0_outAR0$k,SRpma_BH_L2_AR0_outAR0_check$k)
  expect_equal(SRpma_BH_L2_AR0_outAR0$AIC,SRpma_BH_L2_AR0_outAR0_check$AIC)
  expect_equal(SRpma_BH_L2_AR0_outAR0$AICc,SRpma_BH_L2_AR0_outAR0_check$AICc)
  expect_equal(SRpma_BH_L2_AR0_outAR0$BIC,SRpma_BH_L2_AR0_outAR0_check$BIC)


  # BH L2 AR1 outAR False
  load(system.file("extdata","SRpma_BH_L2_AR1_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_BH_L2_AR1_outAR0_check$opt$par)){
    expect_equal(SRpma_BH_L2_AR1_outAR0$opt$par[i],SRpma_BH_L2_AR1_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_BH_L2_AR1_outAR0$opt$value,SRpma_BH_L2_AR1_outAR0_check$opt$value)
  for(i in 1:length(SRpma_BH_L2_AR1_outAR0_check$opt$counts)){
    expect_equal(SRpma_BH_L2_AR1_outAR0$opt$counts[i],SRpma_BH_L2_AR1_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_BH_L2_AR1_outAR0$opt$convergence,SRpma_BH_L2_AR1_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_BH_L2_AR1_outAR0_check$resid)){
    expect_equal(SRpma_BH_L2_AR1_outAR0$resid,SRpma_BH_L2_AR1_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_BH_L2_AR1_outAR0_check$resid2)){
    expect_equal(SRpma_BH_L2_AR1_outAR0$resid2,SRpma_BH_L2_AR1_outAR0_check$resid2)
  }
  expect_equal(SRpma_BH_L2_AR1_outAR0$pars,SRpma_BH_L2_AR1_outAR0_check$pars)
  expect_equal(SRpma_BH_L2_AR1_outAR0$loglik,SRpma_BH_L2_AR1_outAR0_check$loglik)
  expect_equal(SRpma_BH_L2_AR1_outAR0$pred,SRpma_BH_L2_AR1_outAR0_check$pred)
  expect_equal(SRpma_BH_L2_AR1_outAR0$k,SRpma_BH_L2_AR1_outAR0_check$k)
  expect_equal(SRpma_BH_L2_AR1_outAR0$AIC,SRpma_BH_L2_AR1_outAR0_check$AIC)
  expect_equal(SRpma_BH_L2_AR1_outAR0$AICc,SRpma_BH_L2_AR1_outAR0_check$AICc)
  expect_equal(SRpma_BH_L2_AR1_outAR0$BIC,SRpma_BH_L2_AR1_outAR0_check$BIC)

  # BH L2 AR1 outAR True
  load(system.file("extdata","SRpma_BH_L2_AR1_outAR1.rda",package = "frasyr"))
  for(i in 1:length(SRpma_BH_L2_AR1_outAR1_check$opt$par)){
    expect_equal(SRpma_BH_L2_AR1_outAR1$opt$par[i],SRpma_BH_L2_AR1_outAR1_check$opt$par[i])
  }
  expect_equal(SRpma_BH_L2_AR1_outAR1$opt$value,SRpma_BH_L2_AR1_outAR1_check$opt$value)
  for(i in 1:length(SRpma_BH_L2_AR1_outAR1_check$opt$counts)){
    expect_equal(SRpma_BH_L2_AR1_outAR1$opt$counts[i],SRpma_BH_L2_AR1_outAR1_check$opt$counts[i])
  }
  expect_equal(SRpma_BH_L2_AR1_outAR1$opt$convergence,SRpma_BH_L2_AR1_outAR1_check$opt$convergence)
  for(i in 1:length(SRpma_BH_L2_AR1_outAR1_check$resid)){
    expect_equal(SRpma_BH_L2_AR1_outAR1$resid,SRpma_BH_L2_AR1_outAR1_check$resid)
  }
  for(i in 1:length(SRpma_BH_L2_AR1_outAR1_check$resid2)){
    expect_equal(SRpma_BH_L2_AR1_outAR1$resid2,SRpma_BH_L2_AR1_outAR1_check$resid2)
  }
  expect_equal(SRpma_BH_L2_AR1_outAR1$pars,SRpma_BH_L2_AR1_outAR1_check$pars)
  expect_equal(SRpma_BH_L2_AR1_outAR1$loglik,SRpma_BH_L2_AR1_outAR1_check$loglik)
  expect_equal(SRpma_BH_L2_AR1_outAR1$pred,SRpma_BH_L2_AR1_outAR1_check$pred)
  expect_equal(SRpma_BH_L2_AR1_outAR1$k,SRpma_BH_L2_AR1_outAR1_check$k)
  expect_equal(SRpma_BH_L2_AR1_outAR1$AIC,SRpma_BH_L2_AR1_outAR1_check$AIC)
  expect_equal(SRpma_BH_L2_AR1_outAR1$AICc,SRpma_BH_L2_AR1_outAR1_check$AICc)
  expect_equal(SRpma_BH_L2_AR1_outAR1$BIC,SRpma_BH_L2_AR1_outAR1_check$BIC)


  # RI L1 AR0
  load(system.file("extdata","SRpma_RI_L1_AR0_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_RI_L1_AR0_outAR0_check$opt$par)){
    expect_equal(SRpma_RI_L1_AR0_outAR0$opt$par[i],SRpma_RI_L1_AR0_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_RI_L1_AR0_outAR0$opt$value,SRpma_RI_L1_AR0_outAR0_check$opt$value)
  for(i in 1:length(SRpma_RI_L1_AR0_outAR0_check$opt$counts)){
    expect_equal(SRpma_RI_L1_AR0_outAR0$opt$counts[i],SRpma_RI_L1_AR0_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_RI_L1_AR0_outAR0$opt$convergence,SRpma_RI_L1_AR0_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_RI_L1_AR0_outAR0_check$resid)){
    expect_equal(SRpma_RI_L1_AR0_outAR0$resid,SRpma_RI_L1_AR0_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_RI_L1_AR0_outAR0_check$resid2)){
    expect_equal(SRpma_RI_L1_AR0_outAR0$resid2,SRpma_RI_L1_AR0_outAR0_check$resid2)
  }
  expect_equal(SRpma_RI_L1_AR0_outAR0$pars,SRpma_RI_L1_AR0_outAR0_check$pars)
  expect_equal(SRpma_RI_L1_AR0_outAR0$loglik,SRpma_RI_L1_AR0_outAR0_check$loglik)
  expect_equal(SRpma_RI_L1_AR0_outAR0$pred,SRpma_RI_L1_AR0_outAR0_check$pred)
  expect_equal(SRpma_RI_L1_AR0_outAR0$k,SRpma_RI_L1_AR0_outAR0_check$k)
  expect_equal(SRpma_RI_L1_AR0_outAR0$AIC,SRpma_RI_L1_AR0_outAR0_check$AIC)
  expect_equal(SRpma_RI_L1_AR0_outAR0$AICc,SRpma_RI_L1_AR0_outAR0_check$AICc)
  expect_equal(SRpma_RI_L1_AR0_outAR0$BIC,SRpma_RI_L1_AR0_outAR0_check$BIC)


  # RI L1 AR1 outAR False
  load(system.file("extdata","SRpma_RI_L1_AR1_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_RI_L1_AR1_outAR0_check$opt$par)){
    expect_equal(SRpma_RI_L1_AR1_outAR0$opt$par[i],SRpma_RI_L1_AR1_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_RI_L1_AR1_outAR0$opt$value,SRpma_RI_L1_AR1_outAR0_check$opt$value)
  for(i in 1:length(SRpma_RI_L1_AR1_outAR0_check$opt$counts)){
    expect_equal(SRpma_RI_L1_AR1_outAR0$opt$counts[i],SRpma_RI_L1_AR1_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_RI_L1_AR1_outAR0$opt$convergence,SRpma_RI_L1_AR1_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_RI_L1_AR1_outAR0_check$resid)){
    expect_equal(SRpma_RI_L1_AR1_outAR0$resid,SRpma_RI_L1_AR1_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_RI_L1_AR1_outAR0_check$resid2)){
    expect_equal(SRpma_RI_L1_AR1_outAR0$resid2,SRpma_RI_L1_AR1_outAR0_check$resid2)
  }
  expect_equal(SRpma_RI_L1_AR1_outAR0$pars,SRpma_RI_L1_AR1_outAR0_check$pars)
  expect_equal(SRpma_RI_L1_AR1_outAR0$loglik,SRpma_RI_L1_AR1_outAR0_check$loglik)
  expect_equal(SRpma_RI_L1_AR1_outAR0$pred,SRpma_RI_L1_AR1_outAR0_check$pred)
  expect_equal(SRpma_RI_L1_AR1_outAR0$k,SRpma_RI_L1_AR1_outAR0_check$k)
  expect_equal(SRpma_RI_L1_AR1_outAR0$AIC,SRpma_RI_L1_AR1_outAR0_check$AIC)
  expect_equal(SRpma_RI_L1_AR1_outAR0$AICc,SRpma_RI_L1_AR1_outAR0_check$AICc)
  expect_equal(SRpma_RI_L1_AR1_outAR0$BIC,SRpma_RI_L1_AR1_outAR0_check$BIC)

  # RI L1 AR1 outAR True
  load(system.file("extdata","SRpma_RI_L1_AR1_outAR1.rda",package = "frasyr"))
  for(i in 1:length(SRpma_RI_L1_AR1_outAR1_check$opt$par)){
    expect_equal(SRpma_RI_L1_AR1_outAR1$opt$par[i],SRpma_RI_L1_AR1_outAR1_check$opt$par[i])
  }
  expect_equal(SRpma_RI_L1_AR1_outAR1$opt$value,SRpma_RI_L1_AR1_outAR1_check$opt$value)
  for(i in 1:length(SRpma_RI_L1_AR1_outAR1_check$opt$counts)){
    expect_equal(SRpma_RI_L1_AR1_outAR1$opt$counts[i],SRpma_RI_L1_AR1_outAR1_check$opt$counts[i])
  }
  expect_equal(SRpma_RI_L1_AR1_outAR1$opt$convergence,SRpma_RI_L1_AR1_outAR1_check$opt$convergence)
  for(i in 1:length(SRpma_RI_L1_AR1_outAR1_check$resid)){
    expect_equal(SRpma_RI_L1_AR1_outAR1$resid,SRpma_RI_L1_AR1_outAR1_check$resid)
  }
  for(i in 1:length(SRpma_RI_L1_AR1_outAR1_check$resid2)){
    expect_equal(SRpma_RI_L1_AR1_outAR1$resid2,SRpma_RI_L1_AR1_outAR1_check$resid2)
  }
  expect_equal(SRpma_RI_L1_AR1_outAR1$pars,SRpma_RI_L1_AR1_outAR1_check$pars)
  expect_equal(SRpma_RI_L1_AR1_outAR1$loglik,SRpma_RI_L1_AR1_outAR1_check$loglik)
  expect_equal(SRpma_RI_L1_AR1_outAR1$pred,SRpma_RI_L1_AR1_outAR1_check$pred)
  expect_equal(SRpma_RI_L1_AR1_outAR1$k,SRpma_RI_L1_AR1_outAR1_check$k)
  expect_equal(SRpma_RI_L1_AR1_outAR1$AIC,SRpma_RI_L1_AR1_outAR1_check$AIC)
  expect_equal(SRpma_RI_L1_AR1_outAR1$AICc,SRpma_RI_L1_AR1_outAR1_check$AICc)
  expect_equal(SRpma_RI_L1_AR1_outAR1$BIC,SRpma_RI_L1_AR1_outAR1_check$BIC)


  # RI L2 AR0
  load(system.file("extdata","SRpma_RI_L2_AR0_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_RI_L2_AR0_outAR0_check$opt$par)){
    expect_equal(SRpma_RI_L2_AR0_outAR0$opt$par[i],SRpma_RI_L2_AR0_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_RI_L2_AR0_outAR0$opt$value,SRpma_RI_L2_AR0_outAR0_check$opt$value)
  for(i in 1:length(SRpma_RI_L2_AR0_outAR0_check$opt$counts)){
    expect_equal(SRpma_RI_L2_AR0_outAR0$opt$counts[i],SRpma_RI_L2_AR0_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_RI_L2_AR0_outAR0$opt$convergence,SRpma_RI_L2_AR0_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_RI_L2_AR0_outAR0_check$resid)){
    expect_equal(SRpma_RI_L2_AR0_outAR0$resid,SRpma_RI_L2_AR0_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_RI_L2_AR0_outAR0_check$resid2)){
    expect_equal(SRpma_RI_L2_AR0_outAR0$resid2,SRpma_RI_L2_AR0_outAR0_check$resid2)
  }
  expect_equal(SRpma_RI_L2_AR0_outAR0$pars,SRpma_RI_L2_AR0_outAR0_check$pars)
  expect_equal(SRpma_RI_L2_AR0_outAR0$loglik,SRpma_RI_L2_AR0_outAR0_check$loglik)
  expect_equal(SRpma_RI_L2_AR0_outAR0$pred,SRpma_RI_L2_AR0_outAR0_check$pred)
  expect_equal(SRpma_RI_L2_AR0_outAR0$k,SRpma_RI_L2_AR0_outAR0_check$k)
  expect_equal(SRpma_RI_L2_AR0_outAR0$AIC,SRpma_RI_L2_AR0_outAR0_check$AIC)
  expect_equal(SRpma_RI_L2_AR0_outAR0$AICc,SRpma_RI_L2_AR0_outAR0_check$AICc)
  expect_equal(SRpma_RI_L2_AR0_outAR0$BIC,SRpma_RI_L2_AR0_outAR0_check$BIC)


  # RI L2 AR1 outAR False
  load(system.file("extdata","SRpma_RI_L2_AR1_outAR0.rda",package = "frasyr"))
  for(i in 1:length(SRpma_RI_L2_AR1_outAR0_check$opt$par)){
    expect_equal(SRpma_RI_L2_AR1_outAR0$opt$par[i],SRpma_RI_L2_AR1_outAR0_check$opt$par[i])
  }
  expect_equal(SRpma_RI_L2_AR1_outAR0$opt$value,SRpma_RI_L2_AR1_outAR0_check$opt$value)
  for(i in 1:length(SRpma_RI_L2_AR1_outAR0_check$opt$counts)){
    expect_equal(SRpma_RI_L2_AR1_outAR0$opt$counts[i],SRpma_RI_L2_AR1_outAR0_check$opt$counts[i])
  }
  expect_equal(SRpma_RI_L2_AR1_outAR0$opt$convergence,SRpma_RI_L2_AR1_outAR0_check$opt$convergence)
  for(i in 1:length(SRpma_RI_L2_AR1_outAR0_check$resid)){
    expect_equal(SRpma_RI_L2_AR1_outAR0$resid,SRpma_RI_L2_AR1_outAR0_check$resid)
  }
  for(i in 1:length(SRpma_RI_L2_AR1_outAR0_check$resid2)){
    expect_equal(SRpma_RI_L2_AR1_outAR0$resid2,SRpma_RI_L2_AR1_outAR0_check$resid2)
  }
  expect_equal(SRpma_RI_L2_AR1_outAR0$pars,SRpma_RI_L2_AR1_outAR0_check$pars)
  expect_equal(SRpma_RI_L2_AR1_outAR0$loglik,SRpma_RI_L2_AR1_outAR0_check$loglik)
  expect_equal(SRpma_RI_L2_AR1_outAR0$pred,SRpma_RI_L2_AR1_outAR0_check$pred)
  expect_equal(SRpma_RI_L2_AR1_outAR0$k,SRpma_RI_L2_AR1_outAR0_check$k)
  expect_equal(SRpma_RI_L2_AR1_outAR0$AIC,SRpma_RI_L2_AR1_outAR0_check$AIC)
  expect_equal(SRpma_RI_L2_AR1_outAR0$AICc,SRpma_RI_L2_AR1_outAR0_check$AICc)
  expect_equal(SRpma_RI_L2_AR1_outAR0$BIC,SRpma_RI_L2_AR1_outAR0_check$BIC)

  # RI L2 AR1 outAR True
  load(system.file("extdata","SRpma_RI_L2_AR1_outAR1.rda",package = "frasyr"))
  for(i in 1:length(SRpma_RI_L2_AR1_outAR1_check$opt$par)){
    expect_equal(SRpma_RI_L2_AR1_outAR1$opt$par[i],SRpma_RI_L2_AR1_outAR1_check$opt$par[i])
  }
  expect_equal(SRpma_RI_L2_AR1_outAR1$opt$value,SRpma_RI_L2_AR1_outAR1_check$opt$value)
  for(i in 1:length(SRpma_RI_L2_AR1_outAR1_check$opt$counts)){
    expect_equal(SRpma_RI_L2_AR1_outAR1$opt$counts[i],SRpma_RI_L2_AR1_outAR1_check$opt$counts[i])
  }
  expect_equal(SRpma_RI_L2_AR1_outAR1$opt$convergence,SRpma_RI_L2_AR1_outAR1_check$opt$convergence)
  for(i in 1:length(SRpma_RI_L2_AR1_outAR1_check$resid)){
    expect_equal(SRpma_RI_L2_AR1_outAR1$resid,SRpma_RI_L2_AR1_outAR1_check$resid)
  }
  for(i in 1:length(SRpma_RI_L2_AR1_outAR1_check$resid2)){
    expect_equal(SRpma_RI_L2_AR1_outAR1$resid2,SRpma_RI_L2_AR1_outAR1_check$resid2)
  }
  expect_equal(SRpma_RI_L2_AR1_outAR1$pars,SRpma_RI_L2_AR1_outAR1_check$pars)
  expect_equal(SRpma_RI_L2_AR1_outAR1$loglik,SRpma_RI_L2_AR1_outAR1_check$loglik)
  expect_equal(SRpma_RI_L2_AR1_outAR1$pred,SRpma_RI_L2_AR1_outAR1_check$pred)
  expect_equal(SRpma_RI_L2_AR1_outAR1$k,SRpma_RI_L2_AR1_outAR1_check$k)
  expect_equal(SRpma_RI_L2_AR1_outAR1$AIC,SRpma_RI_L2_AR1_outAR1_check$AIC)
  expect_equal(SRpma_RI_L2_AR1_outAR1$AICc,SRpma_RI_L2_AR1_outAR1_check$AICc)
  expect_equal(SRpma_RI_L2_AR1_outAR1$BIC,SRpma_RI_L2_AR1_outAR1_check$BIC)

})

