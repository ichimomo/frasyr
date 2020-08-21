library(frasyr)

context("stock-recruitment, check the correlation between parameters and the steepness")

test_that("output value check",{

  # vpaとget.SRdataの結果を読み込み
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
  res.SRfit <- list()
  for (i in 1:nrow(SRmodel.list)){

    SRfile <- sprintf("SRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])

    load(system.file("extdata",SRfile,package = "frasyr"))

    res.SRfit[[i]] <- eval(parse(text=paste(sprintf("SRpma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

    #再生産関係のパラメータ間相関
    corSR_pma_check <- corSR(res.SRfit[[i]])
    year <- as.character(max(res_vpa_pma$input$rec.year))
    #スティープネス
    steepness_pma_check <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

    #上記結果の読み込み
  corSR_file <- sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",corSR_file,package = "frasyr"))
  corSR_pma <- eval(parse(text=paste(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  # テスト項目corSR
  testcontens <-c("hessian","cov","cor")
  #読み込んだ結果と照合
    for(i in 1:length(testcontents)){
      expect_equal(corSR_pma$testcontents[i],corSR_pma_check$testcontents[i])
    }

  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  steepness_pma <- eval(parse(text=paste(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))
  }

  SPR0_by_getSPR <- get.SPR(res_vpa_pma)

  # テスト項目steepness
  testcontens <-c("SPR0","SB0","R0","h")
  #読み込んだ結果と照合
  for(i in 1:length(testcontents)){
    expect_equal(steepness_pma$testcontents[i],steepness_pma_check$testcontents[i])
  }
  # get.SPRからのSPR0と照合
  expect_equal(steepness_pma_check$SPR0, SPR0_by_getSPR$ysdata[1,4])
})
