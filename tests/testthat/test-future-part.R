library(frasyr)

# ref.F test ----
context("future ref.F") 

test_that("output value check",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  res_ref_f_pma_check <- ref.F(res_vpa_pma,Fcurrent=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                               waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                               d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                               iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )

  res.ref.f_2times <- ref.F(res_vpa_pma,Fcurrent=res_vpa_pma$Fc.at.age*2,
                            waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                            waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
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
  testcontents <-c("sel","max.age","min.age","rps.q","spr.q","Fcurrent","Fmed","Flow","Fhigh","Fmax","F0.1","Fmean","rps.data","FpSPR","summary","ypr.spr","waa","waa.catch","maa","spr0")

  for(i in 1:length(testcontents)){
      tmp1 <- eval(parse(text=paste("res_ref_f_pma$",testcontents[i])))
      tmp2 <- eval(parse(text=paste("res_ref_f_pma_check$",testcontents[i])))
      expect_equivalent(tmp1,tmp2,tolerance=1e-4,label=testcontents[i])
      # %SPRを計算するところで、初期値が変わると1e-4以下の誤差で値がずれるので1e-4をtoleranceに入れる
      # toleranceのつづりが間違っていても誰も教えてくれない（無言でtoleranceを無視する）ため注意
  }
})

# check SR data ----
context("future SRdata")

test_that("oututput value check",{
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



