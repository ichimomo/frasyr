library(frasyr)

context("future ref.F")

test_that("inputput value check",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  res.ref.f <- ref.F(res.pma,sel=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                     waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                     d = 0.001,Fspr.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                     iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )
  #上記引数での計算結果を読み込み
  res_ref_f_sel_pma_check <- read.csv(system.file("extdata","future_ref_F_sel_check.csv",package="frasyr"),row.names=1)
  res_ref_f_max_age_pma_check <- read.csv(system.file("extdata","future_ref_F_max_age_check.csv",package="frasyr"),row.names=1)
  res_ref_f_min_age_pma_check <- read.csv(system.file("extdata","future_ref_F_min_age_check.csv",package="frasyr"),row.names=1)
  res_ref_f_rps_q_pma_check <- read.csv(system.file("extdata","future_ref_F_rps-q_check.csv",package="frasyr"),row.names=1)
  res_ref_f_spr_q_pma_check <- read.csv(system.file("extdata","future_ref_F_spr-q_check.csv",package="frasyr"),row.names=1)
  res_ref_f_Fcurrent_pma_check <- read.csv(system.file("extdata","future_ref_F_Fcurrent_check.csv",package="frasyr"),row.names=1)
  res_ref_f_Fmed_pma_check <- read.csv(system.file("extdata","future_ref_F_Fmed_check.csv",package="frasyr"),row.names=1)
  res_ref_f_Flow_pma_check <- read.csv(system.file("extdata","future_ref_F_Flow_check.csv",package="frasyr"),row.names=1)
  res_ref_f_Fhigh_pma_check <- read.csv(system.file("extdata","future_ref_F_Fhigh_check.csv",package="frasyr"),row.names=1)
  res_ref_f_Fmax_pma_check <- read.csv(system.file("extdata","future_ref_F_Fmax_check.csv",package="frasyr"),row.names=1)
  res_ref_f_F01_pma_check <- read.csv(system.file("extdata","future_ref_F_F01_check.csv",package="frasyr"),row.names=1)
  res_ref_f_Fmean_pma_check <- read.csv(system.file("extdata","future_ref_F_Fmean_check.csv",package="frasyr"),row.names=1)
  res_ref_f_rpsdata_pma_check <- read.csv(system.file("extdata","future_ref_F_rps-data_check.csv",package="frasyr"),row.names=1)
  res_ref_f_FpSPR_pma_check <- read.csv(system.file("extdata","future_ref_F_FpSPR_check.csv",package="frasyr"),row.names=1)
  res_ref_f_summary_pma_check <- read.csv(system.file("extdata","future_ref_F_summary_check.csv",package="frasyr"),row.names=1)
  res_ref_f_ypr_spr_pma_check <- read.csv(system.file("extdata","future_ref_F_ypr-spr_check.csv",package="frasyr"),row.names=1)
  res_ref_f_waa_pma_check <- read.csv(system.file("extdata","future_ref_F_waa_check.csv",package="frasyr"),row.names=1)
  res_ref_f_waa_catch_pma_check <- read.csv(system.file("extdata","future_ref_F_waa-catch_check.csv",package="frasyr"),row.names=1)
  res_ref_f_maa_pma_check <- read.csv(system.file("extdata","future_ref_F_maa_check.csv",package="frasyr"),row.names=1)
  res_ref_f_spr0_pma_check <- read.csv(system.file("extdata","future_ref_F_spr0_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合
  for(i in 1:nrow(res_ref_f_sel_pma_check)){
    expect_equal(as.numeric(res.ref.f$sel[i]), res_ref_f_sel_pma_check[i,])
  }
  expect_equal(res.ref.f$max.age,res_ref_f_max_age_pma_check[1,])
  expect_equal(res.ref.f$min.age,res_ref_f_min_age_pma_check[1,])
  for(i in 1:nrow(res_ref_f_rps_q_pma_check)){
    expect_equal(as.numeric(res.ref.f$rps.q[i]), res_ref_f_rps_q_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_spr_q_pma_check)){
    expect_equal(as.numeric(res.ref.f$spr.q[i]), res_ref_f_spr_q_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_Fcurrent_pma_check)){
    expect_equal(as.numeric(res.ref.f$Fcurrent[i]), res_ref_f_Fcurrent_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_Fmed_pma_check)){
    expect_equal(as.numeric(res.ref.f$Fmed[i]), res_ref_f_Fmed_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_Flow_pma_check)){
    expect_equal(as.numeric(res.ref.f$Flow[i]), res_ref_f_Flow_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_Fhigh_pma_check)){
    expect_equal(as.numeric(res.ref.f$Fhigh[i]), res_ref_f_Fhigh_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_Fmax_pma_check)){
    expect_equal(as.numeric(res.ref.f$Fmax[i]), res_ref_f_Fmax_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_F01_pma_check)){
    expect_equal(as.numeric(res.ref.f$F0.1[i]), res_ref_f_F01_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_Fmean_pma_check)){
    expect_equal(as.numeric(res.ref.f$Fmean[i]), res_ref_f_Fmean_pma_check[i,])
  }
  expect_equal(res.ref.f$rps.data,res_ref_f_rpsdata_pma_check)
  for(i in 1:nrow(res_ref_f_FpSPR_pma_check)){
    for(j in 1:ncol(res_ref_f_FpSPR_pma_check)){
      expect_equal(res.ref.f$FpSPR[i,j],res_ref_f_FpSPR_pma_check[i,j])
    }
  }
  expect_equal(res.ref.f$summary,res_ref_f_summary_pma_check)

  expect_equal(res.ref.f$ypr.spr,res_ref_f_ypr_spr_pma_check)
  for(i in 1:nrow(res_ref_f_waa_pma_check)){
    expect_equal(as.numeric(res.ref.f$waa[i]), res_ref_f_waa_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_waa_catch_pma_check)){
    expect_equal(as.numeric(res.ref.f$waa.catch[i]), res_ref_f_waa_catch_pma_check[i,])
  }
  for(i in 1:nrow(res_ref_f_maa_pma_check)){
    expect_equal(as.numeric(res.ref.f$maa[i]), res_ref_f_maa_pma_check[i,])
  }
  expect_equal(res.ref.f$spr0,res_ref_f_spr0_pma_check[1,])
})

