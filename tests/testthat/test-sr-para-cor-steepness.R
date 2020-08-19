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

  corSR_pma <- corSR(res.SRfit[[i]])
  year = max(res_vpa_pma$input$rec.year)
  steepness_pma <- calc_steepness(SR=res.SRfit[[i]]$input$SR,rec_pars=res.SRfit[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  #上記結果の読み込み
  corSR_file <- sprintf("corSRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])

  load(system.file("extdata",corSR_file,package = "frasyr"))
  corSR_pma_check <- eval(parse(text=paste(sprintf("corSRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  steepness_file <- sprintf("steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  load(system.file("extdata",steepness_file,package = "frasyr"))
  steepness_pma_check <- eval(parse(text=paste(sprintf("steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]))))

  }


})
