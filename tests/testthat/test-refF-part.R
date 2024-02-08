library(frasyr)

context("check ref.F") # ----

test_that("ref.F (level 2)",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  res_ref_f_pma_check <- ref.F(res_vpa_pma,Fcurrent=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                               waa.year=NULL,maa.year=NULL,rps.year = as.numeric(colnames(res_vpa_pma$naa)),
                               max.age = Inf,min.age = 0,
                               d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                               iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )

  res.ref.f_2times <- ref.F(res_vpa_pma,Fcurrent=res_vpa_pma$Fc.at.age*2,
                            waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                            waa.year=NULL,maa.year=NULL,rps.year = as.numeric(colnames(res_vpa_pma$naa)),
                            max.age = Inf,min.age = 0,
                            d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,
                            pSPR = c(seq(10,90,by=10),res_ref_f_pma_check$currentSPR$perSPR*100),
                            iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )


  times2_check <- as.numeric(table(unlist(res_ref_f_pma_check$summary/res.ref.f_2times$summary[,1:16])))

  expect_equal(times2_check[1],2)
  expect_equal(times2_check[2],31)
  expect_equal(times2_check[3],15)
  expect_equal(round(res.ref.f_2times$summary[3,17],2),0.5)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","res_ref_f_pma.rda",package = "frasyr"))

  #照合内容
  testcontents <-c("sel","max.age","min.age","rps.q","spr.q","Fcurrent","Fmed","Flow","Fhigh","Fmax","F0.1","Fmean","rps.data","FpSPR","summary", #"ypr.spr",
                   "waa","waa.catch","maa","spr0")

  for(i in 1:length(testcontents)){
      tmp1 <- eval(parse(text=paste("res_ref_f_pma$",testcontents[i])))
      tmp2 <- eval(parse(text=paste("res_ref_f_pma_check$",testcontents[i])))
      expect_equivalent(tmp1,tmp2,tolerance=1e-4,label=testcontents[i])
      # %SPRを計算するところで、初期値が変わると1e-4以下の誤差で値がずれるので1e-4をtoleranceに入れる
      # toleranceのつづりが間違っていても誰も教えてくれない（無言でtoleranceを無視する）ため注意
  }

  # vpa結果を与えない場合
  res_ref_independent <- ref.F(res=NULL,
                               Fcurrent=res_vpa_pma$Fc.at.age,
                               waa=res_vpa_pma$input$dat$waa$"2011",
                               maa=res_vpa_pma$input$dat$maa$"2011",
                               M  =res_vpa_pma$input$dat$M$"2011",
                               waa.catch=res_vpa_pma$input$dat$waa$"2011",
                               rps.vector=res_vpa_pma$naa[1,]/colSums(res_vpa_pma$ssb),
                               max.age = Inf,
                               min.age = 0,
                               d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                               iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )
  testthat::expect_equal(res_ref_independent$summary,res_ref_f_pma_check$summary, tol=0.001)

  # vpa結果を与えない場合 (biological parameterを単純なベクトルで与える)
  xx <- as.numeric(unlist(res_vpa_pma$Fc.at.age))
  res_ref_independent <- ref.F(res=NULL,
                               Fcurrent=xx,
                               waa=res_vpa_pma$input$dat$waa$"2011",
                               maa=res_vpa_pma$input$dat$maa$"2011",
                               M  =res_vpa_pma$input$dat$M$"2011",
                               waa.catch=res_vpa_pma$input$dat$waa$"2011",
                               rps.vector=res_vpa_pma$naa[1,]/colSums(res_vpa_pma$ssb),
                               max.age = Inf,
                               min.age = 0,
                               d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                               iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )
  

  # 同じ機能を持つcalc_Fratioとの整合性をチェック=> ref.Fとcalc_Fratioは同じ機能を提供
  Fratio_test <- purrr::map_dfr((1:4) * 10, function(x)
    calc_Fratio(res_vpa_pma$Fc.at.age,
                rev(res_vpa_pma$input$dat$waa)[,1],
                rev(res_vpa_pma$input$dat$maa)[,1],
                rev(res_vpa_pma$input$dat$M  )[,1],
                SPRtarget=x,
                plus_group=TRUE,
                waa.catch=NULL,
                Pope=res_vpa_pma$input$Pope,
                return_SPR=TRUE))
  for_test_tmp <- 1/res_ref_f_pma_check$summary[str_c("FpSPR.",1:4 * 10,".SPR")][3,] %>%
    unlist() %>% as.numeric()
  expect_equal(for_test_tmp,Fratio_test$Fratio,tol=0.0001)
  expect_equal(1:4 * 10,Fratio_test$SPR_est,tol=0.0001)

#  MSY_perSPR1 <- calc_future_perSPR(res_future_0.8HCR,res_vpa=res_vpa_example,Fvector=res_MSY$F.msy)
#  MSY_perSPR2 <- calc_future_perSPR(res_future_0.8HCR,res_vpa=res_vpa_example,Fvector=res_MSY$F.msy,target.year=2040:2045)
#  MSY_perSPR3 <- calc_future_perSPR(res_future_0.8HCR,res_vpa=res_vpa_example,Fvector=res_MSY$F.msy,target.col=30)

  MSY_perSPR1 <- calc_future_perSPR(res_future_HSL1,res_vpa=res_vpa_example,Fvector=res_MSY_HSL1$F.msy)
  MSY_perSPR2 <- calc_future_perSPR(res_future_HSL1,res_vpa=res_vpa_example,Fvector=res_MSY_HSL1$F.msy,target.year=2030:2035)
  MSY_perSPR3 <- calc_future_perSPR(res_future_HSL1,res_vpa=res_vpa_example,Fvector=res_MSY_HSL1$F.msy,target.col=30)
  
  #  expect_equal(MSY_perSPR1,0.2307774,tol=1e-4)
  expect_equal(MSY_perSPR1,0.251,tol=1e-3)  
  expect_equal(MSY_perSPR1,MSY_perSPR2,tol=1e-4)
  expect_equal(MSY_perSPR1,MSY_perSPR3,tol=1e-4)

  # calc_future_perSPR with VPA
  faa <- res_vpa_pma$Fc.at.age %>% as.numeric()
  tmp1 <- calc_future_perSPR(fout=NULL,res_vpa=res_vpa_pma,Fvector=faa,
                             plus_group=res_vpa_example$input$plus.group,target.year=2011:2010)
  expect_equal(res_ref_f_pma_check$currentSPR$perSPR,tmp1)

  # プラスグループを考慮しない場合
  res_vpa_pma_noplus <- res_vpa_pma
  res_vpa_pma_noplus$input$plus.group <- FALSE
  res_ref_f_pma_check_noplus <- ref.F(res_vpa_pma_noplus,Fcurrent=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                                      waa.year=NULL,maa.year=NULL,rps.year = as.numeric(colnames(res_vpa_pma$naa)),
                                      max.age = Inf,min.age = 0,
                                      d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                                      iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )
  tmp2 <- calc_future_perSPR(fout=NULL,res_vpa=res_vpa_pma_noplus,Fvector=faa,
                             target.year=2011:2010)
  expect_equal(res_ref_f_pma_check_noplus$currentSPR$perSPR,tmp2)

  # check convert_Fvector
  tmp1 <- convert_Fvector(res_vpa = res_vpa_pma,        faa_vector=faa, faa_bio_year=2011,
                          saa_vector = rep(1,length(faa)), saa_bio_year=2011)
  tmp2 <- convert_Fvector(res_vpa = res_vpa_pma_noplus, faa_vector=faa, faa_bio_year=2011,
                          saa_vector = rep(1,length(faa)), saa_bio_year=2011)
  tmp3 <- convert_Fvector(res_vpa = res_vpa_pma,        faa_vector=tmp1$Fvector, faa_bio_year=2011,
                          saa_vector = faa, saa_bio_year=2011)
  expect_equal(mean(faa/tmp3$Fvector),1,tol=0.0001)
  expect_equal(tmp1$faa_perSPR, res_ref_f_pma_check$currentSPR$perSPR)
  expect_equal(tmp2$faa_perSPR, res_ref_f_pma_check_noplus$currentSPR$perSPR)  
  expect_equal(mean(faa/tmp3$Fvector),1,tol=0.0001)
  
  # just for run
#  MSY_Fratio <- calc_perspr(res_future_0.8HCR,res_vpa=res_vpa_example,Fvector=res_MSY$F.msy,SPRtarget=0.3)

})

context("check get.SRdata") # ----

test_that("get.SRdata (level 2)",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
  SRdata0usingPeriodFrom1990To2000 <- get.SRdata(res_vpa_pma,years=1990:2000)

  #上記引数での計算結果を読み込み
  SRdata0_pma_check <- read.csv(system.file("extdata","future_SRdata0_pma_check.csv",package="frasyr"),row.names=1)
  SRdata0usingPeriodFrom1990To2000_pma_check <- read.csv(system.file("extdata","future_SRdata0usingPeriodFrom1990To2000_pma_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合

  expect_equal(SRdata0$year, SRdata0_pma_check$year)
  #expect_equal(SRdata0$SSB, SRdata0_pma_check$SSB)
  #expect_equal(SRdata0$R, SRdata0_pma_check$R)

  expect_equal(SRdata0usingPeriodFrom1990To2000$year, SRdata0usingPeriodFrom1990To2000_pma_check$year)
  expect_equal(SRdata0usingPeriodFrom1990To2000$SSB, SRdata0usingPeriodFrom1990To2000_pma_check$SSB)
  expect_equal(SRdata0usingPeriodFrom1990To2000$R, SRdata0usingPeriodFrom1990To2000_pma_check$R)

})
