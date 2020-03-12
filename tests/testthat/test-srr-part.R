library(frasyr)

context("SRR model selection by AICc")

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

  SRmodel.select <- SR.list[[1]] # AIC最小モデルを今後使っていく

  # pmaではHS L2 AR0 (outAR=Tのほうが選ばれているが、Fでも同じはず
  # 正規性チェック----
  srr_check1_pma_check <- shapiro.test(SRmodel.select$resid)
  srr_check2_pma_check <- ks.test(SRmodel.select$resid,y="pnorm")

  #上記条件での結果の読み込み
  load(system.file("extdata","srr_check1_pma.rda",package = "frasyr"))
  load(system.file("extdata","srr_check2_pma.rda",package = "frasyr"))

  #結果の照合
  expect_equal(srr_check1_pma,srr_check1_pma_check)
  expect_equal(srr_check2_pma,srr_check2_pma_check)

  # residual trend & acf ----
  ac_res_pma_check <- acf(SRmodel.select$resid2,plot=FALSE)
  #上記条件での結果の読み込み
  load(system.file("extdata","ac_res_pma.rda",package = "frasyr"))
  #結果の照合
  expect_equal(ac_res_pma,ac_res_pma_check)
})

context("SRR bootstrap")

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

  SRmodel.select <- SR.list[[1]] # AIC最小モデルを今後使っていく

  # boot strap ----
  boot_res_pma_check <- boot.SR(SRmodel.select)

  boot_res_pma_median_pars_a <-median(sapply(1:boot_res_pma_check$input$n, function(i) boot_res_pma_check[[i]]$pars$a))
  boot_res_pma_median_pars_b <-median(sapply(1:boot_res_pma_check$input$n, function(i) boot_res_pma_check[[i]]$pars$b))
  boot_res_pma_median_pars_sd <-median(sapply(1:boot_res_pma_check$input$n, function(i) boot_res_pma_check[[i]]$pars$sd))
  boot_res_pma_median_pars_rho <-median(sapply(1:boot_res_pma_check$input$n, function(i) boot_res_pma_check[[i]]$pars$rho))

  # boot_res_pma_median_pars_check <- cbind(boot_res_pma_median_pars_a,boot_res_pma_median_pars_b,boot_res_pma_median_pars_sd,boot_res_pma_median_pars_rho)
  boot_res_pma_median_pars_check <- c(boot_res_pma_median_pars_a,boot_res_pma_median_pars_b,boot_res_pma_median_pars_sd)

  #上記条件での結果の読み込み
  load(system.file("extdata","boot_res_pma_median_pars.rda",package = "frasyr"))
  
  #結果の照合
  expect_equal(boot_res_pma_median_pars_check/as.numeric(boot_res_pma_median_pars[,-4]),rep(1,3),tolerance=0.01,
               scale=0.01)
})

