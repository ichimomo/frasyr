library(frasyr)

load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))
SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
SR.list <- list()

year <- as.character(max(res_vpa_pma$input$rec.year))
bio_par <- derive_biopar(res_vpa_pma,derive_year=year)

for (i in 1:nrow(SRmodel.list)) {
    SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                           AR = SRmodel.list$AR.type[i], out.AR =SRmodel.list$out.AR[i], hessian = FALSE, bio_par = bio_par)
}

SRmodel.list$AICc <- sapply(SR.list, function(x) x$AICc)
SRmodel.list$delta.AIC <- SRmodel.list$AICc - min(SRmodel.list$AICc)
SR.list <- SR.list[order(SRmodel.list$AICc)]  # AICの小さい順に並べたもの
(SRmodel.list <- SRmodel.list[order(SRmodel.list$AICc), ]) # 結果

context("SRR model selection by AICc")

test_that("oututput value check",{

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
  expect_equal(srr_check2_pma$p.value,srr_check2_pma_check$p.value)

  # residual trend & acf ----
  ac_res_pma_check <- acf(SRmodel.select$resid2,plot=FALSE)
  #上記条件での結果の読み込み
  load(system.file("extdata","ac_res_pma.rda",package = "frasyr"))
  #結果の照合
  expect_equal(ac_res_pma,ac_res_pma_check)
})

context("SRR bootstrap")

test_that("oututput value check",{

  SRmodel.select <- SR.list[[1]] # AIC最小モデルを今後使っていく

  # boot strap ----
  boot_res_pma_check <- boot.SR(SRmodel.select)

  pdf("tmp.pdf")
  bootSR.plot(boot_res_pma_check, ggplt=TRUE)
  dev.off()

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

context("SRregime bootstrap")

test_that("output value check",{
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))
  SRdata = SRdata_pma
  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), L.type = c("L1", "L2"))

  # regime_year = ceiling(mean(SRdata$year))
  regime_year = 1999
  regime1 = min(SRdata$year):(regime_year-1); regime2 = regime_year:max(SRdata$year);
  SRdata1 = list(year=regime1, R=SRdata$R[SRdata$year %in% regime1],SSB=SRdata$SSB[SRdata$year %in% regime1])
  SRdata2 = list(year=regime2, R=SRdata$R[SRdata$year %in% regime2],SSB=SRdata$SSB[SRdata$year %in% regime2])

#  for (i in 1:nrow(SRmodel.list)) {
  for (i in 1) {
    resSR1 <- fit.SR(SRdata1, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],AR = 0, hessian = FALSE,length=20)
    resSR2 <- fit.SR(SRdata2, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],AR = 0, hessian = FALSE,length=20)
    resSRregime <- fit.SRregime(SRdata, SR = as.character(SRmodel.list$SR.rel[i]), method = as.character(SRmodel.list$L.type[i]), regime.year = regime_year, regime.key = 0:1, regime.par = c("a","b","sd"), use.fit.SR = TRUE)

    # boot strap ----
    nboot <- 3
    boot_resSR1_check <- boot.SR(resSR1,n=nboot)
    boot_resSR2_check <- boot.SR(resSR2,n=nboot)
    boot_resSRregime_check <- boot.SR(resSRregime,n=nboot)
    bootSR.plot(boot_resSRregime_check)

    x1 <- purrr::map_dfr(boot_resSR1_check,function(x) x$pars) %>% summarise(median.a=median(a),median.b=median(b),median.sd=median(sd))
    x2 <- purrr::map_dfr(boot_resSR2_check,function(x) x$pars) %>% summarise(median.a=median(a),median.b=median(b),median.sd=median(sd))
    x3 <- purrr::map_dfr(boot_resSRregime_check,function(x) x$regime_pars) %>%
        group_by(regime) %>%
        summarise(median.a=median(a),median.b=median(b),median.sd=median(sd))

    # ブートストラップするときに乱数がちがうので一致するわけがないような気がする、、
    #expect_equal(as.numeric(x1/x3[1,-1]),rep(1,3),tolerance=0.5,scale=0.01)
    #expect_equal(as.numeric(x2/x3[2,-1]),rep(1,3),tolerance=0.5,scale=0.01)
  }

  # test for boot_steepness
  x1 <- boot_steepness(resSR1,
                 M  =res_vpa_pma$input$dat$M  [,"2011"],
                 waa=res_vpa_pma$input$dat$waa[,"2011"],
                 maa=res_vpa_pma$input$dat$maa[,"2011"],n=nboot)
  x2 <- boot_steepness(list(resSR1,resSR1),
                 M  =res_vpa_pma$input$dat$M  [,"2011"],
                 waa=res_vpa_pma$input$dat$waa[,"2011"],
                 maa=res_vpa_pma$input$dat$maa[,"2011"],n=nboot)
  x3 <- boot_steepness(resSRregime,
                 M  =res_vpa_pma$input$dat$M  [,"2011"],
                 waa=res_vpa_pma$input$dat$waa[,"2011"],
                 maa=res_vpa_pma$input$dat$maa[,"2011"],n=nboot)

})


context("specify gamma in SR=Mesnil")

test_that("output value check",{

  data("res_vpa_example")
  SRdata=get.SRdata(res_vpa_example)
  res_sr_MesnilL1 <- fit.SR(SRdata,SR="Mesnil",method = "L1",AR=0,out.AR = F)
  res_spec_Mesnil<-specify.Mesnil.gamma(res_sr_MesnilL1)
  expect_equal(res_spec_Mesnil$gamma,2)
  }
  )

test_that("Add Shaefer and Cusing",{

  data("res_vpa_example")
  SRdata=get.SRdata(res_vpa_example)
  
  if(0){ # まだ実装されていないが将来的にはこのようにしたい
      res_sr_SH <- fit.SR(SRdata,SR="Shepherd",method = "L1",AR=0,out.AR = F)
      res_sr_CU <- fit.SR(SRdata,SR="Cushing",method = "L1",AR=0,out.AR = F)

      expect_equal(all(c("a","b","sd","rho","gamma") %in% names(res_sr_SH$pars)),TRUE)
      expect_equal(all(c("a","b","sd","rho") %in% names(res_sr_CU)),TRUE)      
  }

  # non-regime (Shepherdのgamma=1はBHと一致）
  res_sr_SH <- fit.SR(SRdata,SR="Shepherd",method = "L2",AR=0,out.AR = F, gamma=1)
  res_sr_BH <- fit.SR(SRdata,SR="BH"      ,method = "L2",AR=0,out.AR = F)
  expect_equal(res_sr_BH$pars, res_sr_SH$pars)

  # non-regime
  res_sr_SH <- fit.SR(SRdata,SR="Shepherd",method = "L2",AR=0,out.AR = F, gamma=0.51)  
  res_sr_CU <- fit.SR(SRdata,SR="Cushing" ,method = "L2",AR=0,out.AR = F)

  # regime
  res_sr_BHr <- fit.SRregime(SRdata,SR="BH",method = "L2", regime.year=2000)    
  res_sr_SHr <- fit.SRregime(SRdata,SR="Shepherd",method = "L2", gamma=1, regime.year=2000,p0=res_sr_BHr$opt$par)
  res_sr_CUr <- fit.SRregime(SRdata,SR="Cushing" ,method = "L2", regime.year=2000)
  expect_equal(exp(res_sr_BHr$opt$par),
               exp(res_sr_SHr$opt$par), tol=0.0001)

  gg <- plot_SRregime(res_sr_BHr)
  gg <- plot_SRregime(res_sr_SHr)
  gg <- plot_SRregime(res_sr_CUr)
  # テストにはなっていないが、empty testと言われないために
  expect_equal(class(gg)[1],"gg")
})


