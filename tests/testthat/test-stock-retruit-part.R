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

  #照合内容
  testcontents <-c("opt$par","opt$value","$opt$counts","opt$convergence","resid","resid2","pars","loglik","pred","k","AIC","AICc","BIC")


  # HS L1 AR0 ----
  load(system.file("extdata","SRpma_HS_L1_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_HS_L1_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_HS_L1_AR0_outAR0_check$",testcontents[i]))))
  }

  # HS L1 AR1 outAR False ----
  load(system.file("extdata","SRpma_HS_L1_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_HS_L1_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_HS_L1_AR1_outAR0_check$",testcontents[i]))))
  }

  # HS L1 AR1 outAR True ----
  load(system.file("extdata","SRpma_HS_L1_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_HS_L1_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("SRpma_HS_L1_AR1_outAR1_check$",testcontents[i]))))
  }

  # HS L2 AR0 ----
  load(system.file("extdata","SRpma_HS_L2_AR0_outAR0.rda",package = "frasyr"))
  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_HS_L2_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_HS_L2_AR0_outAR0_check$",testcontents[i]))))
  }

  # HS L2 AR1 outAR False ----
  load(system.file("extdata","SRpma_HS_L2_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_HS_L2_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_HS_L2_AR1_outAR0_check$",testcontents[i]))))
  }

  # HS L2 AR1 outAR True ----
  load(system.file("extdata","SRpma_HS_L2_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_HS_L2_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("SRpma_HS_L2_AR1_outAR1_check$",testcontents[i]))))
  }

  # BH L1 AR0 ----
  load(system.file("extdata","SRpma_BH_L1_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_BH_L1_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_BH_L1_AR0_outAR0_check$",testcontents[i]))))
  }

  # BH L1 AR1 outAR False ----
  load(system.file("extdata","SRpma_BH_L1_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_BH_L1_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_BH_L1_AR1_outAR0_check$",testcontents[i]))))
  }

  # BH L1 AR1 outAR True ----
  load(system.file("extdata","SRpma_BH_L1_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_BH_L1_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("SRpma_BH_L1_AR1_outAR1_check$",testcontents[i]))))
  }

  # BH L2 AR0 ----
  load(system.file("extdata","SRpma_BH_L2_AR0_outAR0.rda",package = "frasyr"))
  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_BH_L2_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_BH_L2_AR0_outAR0_check$",testcontents[i]))))
  }

  # BH L2 AR1 outAR False ----
  load(system.file("extdata","SRpma_BH_L2_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_BH_L2_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_BH_L2_AR1_outAR0_check$",testcontents[i]))))
  }

  # BH L2 AR1 outAR True ----
  load(system.file("extdata","SRpma_BH_L2_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_BH_L2_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("SRpma_BH_L2_AR1_outAR1_check$",testcontents[i]))))
  }

  # RI L1 AR0 ----
  load(system.file("extdata","SRpma_RI_L1_AR0_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_RI_L1_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_RI_L1_AR0_outAR0_check$",testcontents[i]))))
  }

  # RI L1 AR1 outAR False ----
  load(system.file("extdata","SRpma_RI_L1_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_RI_L1_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_RI_L1_AR1_outAR0_check$",testcontents[i]))))
  }

  # RI L1 AR1 outAR True ----
  load(system.file("extdata","SRpma_RI_L1_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_RI_L1_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("SRpma_RI_L1_AR1_outAR1_check$",testcontents[i]))))
  }

  # RI L2 AR0 ----
  load(system.file("extdata","SRpma_RI_L2_AR0_outAR0.rda",package = "frasyr"))
  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_RI_L2_AR0_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_RI_L2_AR0_outAR0_check$",testcontents[i]))))
  }

  # RI L2 AR1 outAR False ----
  load(system.file("extdata","SRpma_RI_L2_AR1_outAR0.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_RI_L2_AR1_outAR0$",testcontents[i]))),eval(parse(text=paste("SRpma_RI_L2_AR1_outAR0_check$",testcontents[i]))))
  }

  # RI L2 AR1 outAR True ----
  load(system.file("extdata","SRpma_RI_L2_AR1_outAR1.rda",package = "frasyr"))

  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("SRpma_RI_L2_AR1_outAR1$",testcontents[i]))),eval(parse(text=paste("SRpma_RI_L2_AR1_outAR1_check$",testcontents[i]))))
  }


})

