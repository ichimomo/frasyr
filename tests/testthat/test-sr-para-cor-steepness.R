library(frasyr)

context("check the correlation between parameters and the steepness")

test_that("output value check",{

  # vpaとget.SRdataの結果を読み込み
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))

  res.SRfit <- list()

  # fit.SR計算
  for (i in 1:nrow(SRmodel.list)){
    res.SRfit[[i]] <- fit.SR(SRdata = SRdata_pma,SR=SRmodel.list$SR.rel[i],method=SRmodel.list$L.type[i],AR=SRmodel.list$AR.type[i],out.AR=SRmodel.list$out.AR[i])
    }

  year <- as.character(max(res_vpa_pma$input$rec.year))

  #各SR推定の設定でテスト

  # テスト項目corSR
  testcontents.corSR <-c("hessian","cov","cor")
  # テスト項目steepness
  testcontents.h <-c("SPR0","SB0","R0","h")
  SPR0_by_getSPR <- get.SPR(res_vpa_pma)

  # test the returens of corSR and calc_steepness for each SR based on the SRmodel.list array
  i <- 3 #skip first 3 (AR=0 outAR=T)
  # HS L1 AR1 outAR True ----
  i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # BH L1 AR1 outAR True ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # RI L1 AR1 outAR True ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # HS L1 AR0 ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # BH L1 AR0 ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # RI L1 AR0 ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # HS L1 AR1 outAR False ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # BH L1 AR1 outAR False ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # RI L1 AR1 outAR False ----
    i<-i+1
    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み(corSR)
    corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",corSR_file,package = "frasyr"))
    #変数名変更
    corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合corSR
    for(j in 1:length(testcontents.corSR)){
      expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
    }

    #上記結果の読み込み(steepness)
    steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
    load(system.file("extdata",steepness_file,package = "frasyr"))
    #変数名変更
    steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #読み込んだ結果と照合h
    for(j in 1:length(testcontents.h)){
      expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
    }
    # get.SPRからのSPR0と照合h
    expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  i <- i+3 #skip first 3 for L.type=L2 (AR=0 outAR=T)
  # HS L2 AR1 outAR True ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # BH L2 AR1 outAR True ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # RI L2 AR1 outAR True ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # HS L2 AR0 ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # BH L2 AR0 ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # RI L2 AR0 ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # HS L2 AR1 outAR False ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])

  # BH L2 AR1 outAR False ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  # RI L2 AR1 outAR False ----
  i<-i+1
  #再生産関係のパラメータ間相関
  corSR_pma_check <- corSR(res.SRfit[[i]])
  #スティープネス
  steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み(corSR)
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  #変数名変更
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合corSR
  for(j in 1:length(testcontents.corSR)){
    expect_equal(corSR_pma$testcontents.corSR[j],corSR_pma_check$testcontents.corSR[j])
  }

  #上記結果の読み込み(steepness)
  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  #変数名変更
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  #読み込んだ結果と照合h
  for(j in 1:length(testcontents.h)){
    expect_equal(steepness_pma$testcontents[j],steepness_pma_check$testcontents[j])
  }
  # get.SPRからのSPR0と照合h
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])


  })
