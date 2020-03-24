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
    # L1のsdの出力変更のための暫定コード
    if(SRmodel.list$L.type[i]=="L1"){
        # pars$sdは残差のRMSEに一致
        expect_equal(SR.list[[i]]$pars$sd, sqrt(mean((SR.list[[i]]$resid)^2)),tolerance=0.001)
        # sd.predとして出力されるものは過去のsdと一致        
        SR.list[[i]]$pars$sd <- SR.list[[i]]$sd.pred # sdはsd.predに置き換え
    }
  }

  # SRタイプ、L1L2、自己相関タイプごとに異なるオブジェクトへ格納
  for (i in 1:nrow(SRmodel.list)) {
    assign(sprintf("SRpma_%s_%s_AR%d_outAR%d_check",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),SR.list[[i]])
  }

  #照合内容
  testcontents <-c("opt$par","opt$value",#"opt$counts", # 最適化したときの繰り返し回数は環境によって異なることがあるためコメントアウト
                   "opt$convergence","resid","resid2","pars","loglik","pred","k","AIC","AICc","BIC")


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

test_that("tentative test for sd of L1 and L2",{  
      SRdata <- list(year=2000:2019,SSB=c(rep(1:5,2),3),R=1*exp(c(rep(-1,4),-2,rep(1,4),2,0)))
      res_L2 <- fit.SR(SRdata,method="L2",AR=0) 
      #res_L2$pars
      res_L1 <- fit.SR(SRdata,method="L1",AR=0) 
      #res_L1$pars
      #plot_SRdata(SRdata)
      #points(res_L1$pred,type="l")
      #points(res_L2$pred,type="l",col=2)

      # パラメータa,bはすべて1になるはず
      for_test <- as.numeric(c(res_L1$pars[c("a","b")],res_L2$pars[c("a","b")])) %>% round(4)
      expect_equal(for_test,rep(1,4))

      # L1のresidのRMSEがres_L1$pars$sdと一致するはず
      expect_equal(sqrt(mean((res_L1$resid)^2)),res_L1$pars$sd)
      # L2のpars$sdとL1のpars$sdも一致するはず      
      expect_equal(res_L2$pars$sd,res_L1$pars$sd)      
})

