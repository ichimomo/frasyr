library(frasyr)

context("stock-recruitment SRdata")

test_that("oututput value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
  SRdata_pma_check <- get.SRdata(res_pma)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  #結果の数値を照合
  expect_equal(SRdata_pma_check, SRdata_pma)
})

context("stock-recruitment fitSR")

test_that("oututput value check",{
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1, 1), out.AR=c(TRUE,TRUE,FALSE), L.type = c("L1", "L2"))
  SR.list <- list()
  for (i in 1:nrow(SRmodel.list)) {

    SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                           AR = SRmodel.list$AR.type[i], out.AR =SRmodel.list$out.AR[i], hessian = FALSE)
  }

  # SRタイプ、L1L2、自己相関タイプごとに異なるオブジェクトへ格納
  for (i in 1:nrow(SRmodel.list)) {
    assign(sprintf("SRpma_%s_%s_AR%d_outAR%d_check",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),SR.list[[i]])
  }

  #上記引数での計算結果を読み込み
  for(i in 1:nrow(SRmodel.list)){
    checkfile.name <- sprintf("SRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i],SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])

    resfitSR <- load(system.file("extdata",checkfile.name,package = "frasyr"))
    fitSR <- eval(parse(text=resfitSR))

    resfitSR_check =paste("SRpma_",SRmodel.list$SR.rel[i],"_",SRmodel.list$L.type[i],"_AR", SRmodel.list$AR.type[i],"_outAR",as.numeric(SRmodel.list$out.AR[i]), sep="")
    fitSR_check <- eval(parse(text=resfitSR_check))

    expect_equal(fitSR,fitSR_check)
   }


})

