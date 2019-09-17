library(frasyr)

context("future ref.F")

test_that("oututput value check",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  res.ref.f <- ref.F(res.pma,sel=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                     waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                     d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
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
  res_ref_f_FpSPR_pma_check <- read.csv(system.file("extdata","future_ref_F_FpSPR_check.csv",package="frasyr"),row.names=1) # summaryとFpSPRの内容は全く同じなので、本来ならテストは必要ない
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

context("future SRdata")

test_that("oututput value check",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  SRdata <- get.SRdata(res.pma)
  SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
  SRdata0usingPeriodFrom1990To2000 <- get.SRdata(res.pma,years=1990:2000)

  #上記引数での計算結果を読み込み
  SRdata_pma_check <- read.csv(system.file("extdata","future_SRdata_pma_check.csv",package="frasyr"),row.names=1)
  SRdata0_pma_check <- read.csv(system.file("extdata","future_SRdata0_pma_check.csv",package="frasyr"),row.names=1)
  SRdata0usingPeriodFrom1990To2000_pma_check <- read.csv(system.file("extdata","future_SRdata0usingPeriodFrom1990To2000_pma_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合
  expect_equal(SRdata$year, SRdata_pma_check$year)
  expect_equal(SRdata$SSB, SRdata_pma_check$SSB)
  expect_equal(SRdata$R, SRdata_pma_check$R)

  expect_equal(SRdata0$year, SRdata0_pma_check$year)
  #expect_equal(SRdata0$SSB, SRdata0_pma_check$SSB)
  #expect_equal(SRdata0$R, SRdata0_pma_check$R)

  expect_equal(SRdata0usingPeriodFrom1990To2000$year, SRdata0usingPeriodFrom1990To2000_pma_check$year)
  expect_equal(SRdata0usingPeriodFrom1990To2000$SSB, SRdata0usingPeriodFrom1990To2000_pma_check$SSB)
  expect_equal(SRdata0usingPeriodFrom1990To2000$R, SRdata0usingPeriodFrom1990To2000_pma_check$R)

})

context("future fitSR")

test_that("oututput value check",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  SRdata <- get.SRdata(res.pma)

  HS.par0 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
  HS.par1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,hessian=FALSE)
  BH.par0 <- fit.SR(SRdata,SR="BH",method="L2",AR=0,hessian=FALSE)
  BH.par1 <- fit.SR(SRdata,SR="BH",method="L2",AR=1,hessian=FALSE)
  RI.par0 <- fit.SR(SRdata,SR="RI",method="L2",AR=0,hessian=FALSE)
  RI.par1 <- fit.SR(SRdata,SR="RI",method="L2",AR=1,hessian=FALSE)

  #上記引数での計算結果を読み込み
  HSpar0_opt_par_pma_check <- read.csv(system.file("extdata","future_HSpar0_opt_par_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_opt_value_pma_check <- read.csv(system.file("extdata","future_HSpar0_opt_value_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_opt_counts_pma_check <- read.csv(system.file("extdata","future_HSpar0_opt_counts_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_opt_convergence_pma_check <- read.csv(system.file("extdata","future_HSpar0_opt_convergence_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_resid_pma_check <- read.csv(system.file("extdata","future_HSpar0_resid_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_resid2_pma_check <- read.csv(system.file("extdata","future_HSpar0_resid2_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_pars_pma_check <- read.csv(system.file("extdata","future_HSpar0_pars_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_loglik_pma_check <- read.csv(system.file("extdata","future_HSpar0_loglik_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_pred_pma_check <- read.csv(system.file("extdata","future_HSpar0_pred_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_k_pma_check <- read.csv(system.file("extdata","future_HSpar0_k_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_AIC_pma_check <- read.csv(system.file("extdata","future_HSpar0_AIC_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_AICc_pma_check <- read.csv(system.file("extdata","future_HSpar0_AICc_pma_check.csv",package="frasyr"),row.names=1)
  HSpar0_BIC_pma_check <- read.csv(system.file("extdata","future_HSpar0_BIC_pma_check.csv",package="frasyr"),row.names=1)

  HSpar1_opt_par_pma_check <- read.csv(system.file("extdata","future_HSpar1_opt_par_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_opt_value_pma_check <- read.csv(system.file("extdata","future_HSpar1_opt_value_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_opt_counts_pma_check <- read.csv(system.file("extdata","future_HSpar1_opt_counts_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_opt_convergence_pma_check <- read.csv(system.file("extdata","future_HSpar1_opt_convergence_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_resid_pma_check <- read.csv(system.file("extdata","future_HSpar1_resid_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_resid2_pma_check <- read.csv(system.file("extdata","future_HSpar1_resid2_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_pars_pma_check <- read.csv(system.file("extdata","future_HSpar1_pars_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_loglik_pma_check <- read.csv(system.file("extdata","future_HSpar1_loglik_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_pred_pma_check <- read.csv(system.file("extdata","future_HSpar1_pred_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_k_pma_check <- read.csv(system.file("extdata","future_HSpar1_k_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_AIC_pma_check <- read.csv(system.file("extdata","future_HSpar1_AIC_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_AICc_pma_check <- read.csv(system.file("extdata","future_HSpar1_AICc_pma_check.csv",package="frasyr"),row.names=1)
  HSpar1_BIC_pma_check <- read.csv(system.file("extdata","future_HSpar1_BIC_pma_check.csv",package="frasyr"),row.names=1)

  BHpar0_opt_par_pma_check <- read.csv(system.file("extdata","future_BHpar0_opt_par_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_opt_value_pma_check <- read.csv(system.file("extdata","future_BHpar0_opt_value_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_opt_counts_pma_check <- read.csv(system.file("extdata","future_BHpar0_opt_counts_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_opt_convergence_pma_check <- read.csv(system.file("extdata","future_BHpar0_opt_convergence_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_resid_pma_check <- read.csv(system.file("extdata","future_BHpar0_resid_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_resid2_pma_check <- read.csv(system.file("extdata","future_BHpar0_resid2_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_pars_pma_check <- read.csv(system.file("extdata","future_BHpar0_pars_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_loglik_pma_check <- read.csv(system.file("extdata","future_BHpar0_loglik_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_pred_pma_check <- read.csv(system.file("extdata","future_BHpar0_pred_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_k_pma_check <- read.csv(system.file("extdata","future_BHpar0_k_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_AIC_pma_check <- read.csv(system.file("extdata","future_BHpar0_AIC_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_AICc_pma_check <- read.csv(system.file("extdata","future_BHpar0_AICc_pma_check.csv",package="frasyr"),row.names=1)
  BHpar0_BIC_pma_check <- read.csv(system.file("extdata","future_BHpar0_BIC_pma_check.csv",package="frasyr"),row.names=1)

  BHpar1_opt_par_pma_check <- read.csv(system.file("extdata","future_BHpar1_opt_par_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_opt_value_pma_check <- read.csv(system.file("extdata","future_BHpar1_opt_value_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_opt_counts_pma_check <- read.csv(system.file("extdata","future_BHpar1_opt_counts_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_opt_convergence_pma_check <- read.csv(system.file("extdata","future_BHpar1_opt_convergence_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_resid_pma_check <- read.csv(system.file("extdata","future_BHpar1_resid_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_resid2_pma_check <- read.csv(system.file("extdata","future_BHpar1_resid2_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_pars_pma_check <- read.csv(system.file("extdata","future_BHpar1_pars_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_loglik_pma_check <- read.csv(system.file("extdata","future_BHpar1_loglik_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_pred_pma_check <- read.csv(system.file("extdata","future_BHpar1_pred_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_k_pma_check <- read.csv(system.file("extdata","future_BHpar1_k_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_AIC_pma_check <- read.csv(system.file("extdata","future_BHpar1_AIC_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_AICc_pma_check <- read.csv(system.file("extdata","future_BHpar1_AICc_pma_check.csv",package="frasyr"),row.names=1)
  BHpar1_BIC_pma_check <- read.csv(system.file("extdata","future_BHpar1_BIC_pma_check.csv",package="frasyr"),row.names=1)

  RIpar0_opt_par_pma_check <- read.csv(system.file("extdata","future_RIpar0_opt_par_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_opt_value_pma_check <- read.csv(system.file("extdata","future_RIpar0_opt_value_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_opt_counts_pma_check <- read.csv(system.file("extdata","future_RIpar0_opt_counts_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_opt_convergence_pma_check <- read.csv(system.file("extdata","future_RIpar0_opt_convergence_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_resid_pma_check <- read.csv(system.file("extdata","future_RIpar0_resid_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_resid2_pma_check <- read.csv(system.file("extdata","future_RIpar0_resid2_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_pars_pma_check <- read.csv(system.file("extdata","future_RIpar0_pars_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_loglik_pma_check <- read.csv(system.file("extdata","future_RIpar0_loglik_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_pred_pma_check <- read.csv(system.file("extdata","future_RIpar0_pred_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_k_pma_check <- read.csv(system.file("extdata","future_RIpar0_k_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_AIC_pma_check <- read.csv(system.file("extdata","future_RIpar0_AIC_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_AICc_pma_check <- read.csv(system.file("extdata","future_RIpar0_AICc_pma_check.csv",package="frasyr"),row.names=1)
  RIpar0_BIC_pma_check <- read.csv(system.file("extdata","future_RIpar0_BIC_pma_check.csv",package="frasyr"),row.names=1)

  RIpar1_opt_par_pma_check <- read.csv(system.file("extdata","future_RIpar1_opt_par_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_opt_value_pma_check <- read.csv(system.file("extdata","future_RIpar1_opt_value_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_opt_counts_pma_check <- read.csv(system.file("extdata","future_RIpar1_opt_counts_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_opt_convergence_pma_check <- read.csv(system.file("extdata","future_RIpar1_opt_convergence_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_resid_pma_check <- read.csv(system.file("extdata","future_RIpar1_resid_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_resid2_pma_check <- read.csv(system.file("extdata","future_RIpar1_resid2_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_pars_pma_check <- read.csv(system.file("extdata","future_RIpar1_pars_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_loglik_pma_check <- read.csv(system.file("extdata","future_RIpar1_loglik_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_pred_pma_check <- read.csv(system.file("extdata","future_RIpar1_pred_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_k_pma_check <- read.csv(system.file("extdata","future_RIpar1_k_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_AIC_pma_check <- read.csv(system.file("extdata","future_RIpar1_AIC_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_AICc_pma_check <- read.csv(system.file("extdata","future_RIpar1_AICc_pma_check.csv",package="frasyr"),row.names=1)
  RIpar1_BIC_pma_check <- read.csv(system.file("extdata","future_RIpar1_BIC_pma_check.csv",package="frasyr"),row.names=1)


  #結果の数値を照合
  #HS.par0

  for(i in 1:nrow(HSpar0_opt_par_pma_check)){
    expect_equal(HS.par0$opt$par[i],HSpar0_opt_par_pma_check[i,])
  }
  expect_equal(HS.par0$opt$value,as.numeric(HSpar0_opt_value_pma_check))
  for(i in 1:nrow(HSpar0_opt_counts_pma_check)){
    expect_equal(as.numeric(HS.par0$opt$counts[i]),HSpar0_opt_counts_pma_check[i,])
  }
  expect_equal(HS.par0$opt$convergence,as.numeric(HSpar0_opt_convergence_pma_check))
  for(i in 1:nrow(HSpar0_resid_pma_check)){
    expect_equal(HS.par0$resid[i],HSpar0_resid_pma_check[i,])
  }
  for(i in 1:nrow(HSpar0_resid2_pma_check)){
    expect_equal(HS.par0$resid2[i],HSpar0_resid2_pma_check[i,])
  }
  expect_equal(HS.par0$pars,HSpar0_pars_pma_check)
  expect_equal(HS.par0$loglik,as.numeric(HSpar0_loglik_pma_check))
  expect_equal(HS.par0$pred,HSpar0_pred_pma_check)
  expect_equal(HS.par0$k,as.numeric(HSpar0_k_pma_check))
  expect_equal(HS.par0$AIC,as.numeric(HSpar0_AIC_pma_check))
  expect_equal(HS.par0$AICc,as.numeric(HSpar0_AICc_pma_check))
  expect_equal(HS.par0$BIC,as.numeric(HSpar0_BIC_pma_check))

  #HS.par1

  for(i in 1:nrow(HSpar1_opt_par_pma_check)){
    expect_equal(HS.par1$opt$par[i],HSpar1_opt_par_pma_check[i,])
  }
  expect_equal(HS.par1$opt$value,as.numeric(HSpar1_opt_value_pma_check))
  for(i in 1:nrow(HSpar1_opt_counts_pma_check)){
    expect_equal(as.numeric(HS.par1$opt$counts[i]),HSpar1_opt_counts_pma_check[i,])
  }
  expect_equal(HS.par1$opt$convergence,as.numeric(HSpar1_opt_convergence_pma_check))
  for(i in 1:nrow(HSpar1_resid_pma_check)){
    expect_equal(HS.par1$resid[i],HSpar1_resid_pma_check[i,])
  }
  for(i in 1:nrow(HSpar1_resid2_pma_check)){
    expect_equal(HS.par1$resid2[i],HSpar1_resid2_pma_check[i,])
  }
  expect_equal(HS.par1$pars,HSpar1_pars_pma_check)
  expect_equal(HS.par1$loglik,as.numeric(HSpar1_loglik_pma_check))
  expect_equal(HS.par1$pred,HSpar1_pred_pma_check)
  expect_equal(HS.par1$k,as.numeric(HSpar1_k_pma_check))
  expect_equal(HS.par1$AIC,as.numeric(HSpar1_AIC_pma_check))
  expect_equal(HS.par1$AICc,as.numeric(HSpar1_AICc_pma_check))

  #BHpar0

  for(i in 1:nrow(BHpar0_opt_par_pma_check)){
    expect_equal(BH.par0$opt$par[i],BHpar0_opt_par_pma_check[i,])
  }
  expect_equal(BH.par0$opt$value,as.numeric(BHpar0_opt_value_pma_check))
  for(i in 1:nrow(BHpar0_opt_counts_pma_check)){
    expect_equal(as.numeric(BH.par0$opt$counts[i]),BHpar0_opt_counts_pma_check[i,])
  }
  expect_equal(BH.par0$opt$convergence,as.numeric(BHpar0_opt_convergence_pma_check))
  for(i in 1:nrow(BHpar0_resid_pma_check)){
    expect_equal(BH.par0$resid[i],BHpar0_resid_pma_check[i,])
  }
  for(i in 1:nrow(BHpar0_resid2_pma_check)){
    expect_equal(BH.par0$resid2[i],BHpar0_resid2_pma_check[i,])
  }
  expect_equal(BH.par0$pars,BHpar0_pars_pma_check)
  expect_equal(BH.par0$loglik,as.numeric(BHpar0_loglik_pma_check))
  expect_equal(BH.par0$pred,BHpar0_pred_pma_check)
  expect_equal(BH.par0$k,as.numeric(BHpar0_k_pma_check))
  expect_equal(BH.par0$AIC,as.numeric(BHpar0_AIC_pma_check))
  expect_equal(BH.par0$AICc,as.numeric(BHpar0_AICc_pma_check))
  expect_equal(BH.par0$BIC,as.numeric(BHpar0_BIC_pma_check))

  #BHpar1

  for(i in 1:nrow(BHpar1_opt_par_pma_check)){
    expect_equal(BH.par1$opt$par[i],BHpar1_opt_par_pma_check[i,])
  }
  expect_equal(BH.par1$opt$value,as.numeric(BHpar1_opt_value_pma_check))
  for(i in 1:nrow(BHpar1_opt_counts_pma_check)){
      #    expect_equal(as.numeric(BH.par1$opt$counts[i]),BHpar1_opt_counts_pma_check[i,])
      # 環境が変わるとこのカウント数が変わるので、とりあえずコメントアウト
  }
  expect_equal(BH.par1$opt$convergence,as.numeric(BHpar1_opt_convergence_pma_check))
  for(i in 1:nrow(BHpar1_resid_pma_check)){
    expect_equal(BH.par1$resid[i],BHpar1_resid_pma_check[i,])
  }
  for(i in 1:nrow(BHpar1_resid2_pma_check)){
    expect_equal(BH.par1$resid2[i],BHpar1_resid2_pma_check[i,])
  }
  expect_equal(BH.par1$pars,BHpar1_pars_pma_check)
  expect_equal(BH.par1$loglik,as.numeric(BHpar1_loglik_pma_check))
  expect_equal(BH.par1$pred,BHpar1_pred_pma_check)
  expect_equal(BH.par1$k,as.numeric(BHpar1_k_pma_check))
  expect_equal(BH.par1$AIC,as.numeric(BHpar1_AIC_pma_check))
  expect_equal(BH.par1$AICc,as.numeric(BHpar1_AICc_pma_check))

  #RIpar0


  for(i in 1:nrow(RIpar0_opt_par_pma_check)){
    expect_equal(RI.par0$opt$par[i],RIpar0_opt_par_pma_check[i,])
  }
  expect_equal(RI.par0$opt$value,as.numeric(RIpar0_opt_value_pma_check))
  for(i in 1:nrow(RIpar0_opt_counts_pma_check)){
    expect_equal(as.numeric(RI.par0$opt$counts[i]),RIpar0_opt_counts_pma_check[i,])
  }
  expect_equal(RI.par0$opt$convergence,as.numeric(RIpar0_opt_convergence_pma_check))
  for(i in 1:nrow(RIpar0_resid_pma_check)){
    expect_equal(RI.par0$resid[i],RIpar0_resid_pma_check[i,])
  }
  for(i in 1:nrow(RIpar0_resid2_pma_check)){
    expect_equal(RI.par0$resid2[i],RIpar0_resid2_pma_check[i,])
  }
  expect_equal(RI.par0$pars,RIpar0_pars_pma_check)
  expect_equal(RI.par0$loglik,as.numeric(RIpar0_loglik_pma_check))
  expect_equal(RI.par0$pred,RIpar0_pred_pma_check)
  expect_equal(RI.par0$k,as.numeric(RIpar0_k_pma_check))
  expect_equal(RI.par0$AIC,as.numeric(RIpar0_AIC_pma_check))
  expect_equal(RI.par0$AICc,as.numeric(RIpar0_AICc_pma_check))
  expect_equal(RI.par0$BIC,as.numeric(RIpar0_BIC_pma_check))

  #RIpar1

  for(i in 1:nrow(RIpar1_opt_par_pma_check)){
    expect_equal(RI.par1$opt$par[i],RIpar1_opt_par_pma_check[i,])
  }
  expect_equal(RI.par1$opt$value,as.numeric(RIpar1_opt_value_pma_check))
  for(i in 1:nrow(RIpar1_opt_counts_pma_check)){
    expect_equal(as.numeric(RI.par1$opt$counts[i]),RIpar1_opt_counts_pma_check[i,])
  }
  expect_equal(RI.par1$opt$convergence,as.numeric(RIpar1_opt_convergence_pma_check))
  for(i in 1:nrow(RIpar1_resid_pma_check)){
    expect_equal(RI.par1$resid[i],RIpar1_resid_pma_check[i,])
  }
  for(i in 1:nrow(RIpar1_resid2_pma_check)){
    expect_equal(RI.par1$resid2[i],RIpar1_resid2_pma_check[i,])
  }
  expect_equal(RI.par1$pars,RIpar1_pars_pma_check)
  expect_equal(RI.par1$loglik,as.numeric(RIpar1_loglik_pma_check))
  expect_equal(RI.par1$pred,RIpar1_pred_pma_check)
  expect_equal(RI.par1$k,as.numeric(RIpar1_k_pma_check))
  expect_equal(RI.par1$AIC,as.numeric(RIpar1_AIC_pma_check))
  expect_equal(RI.par1$AICc,as.numeric(RIpar1_AICc_pma_check))

})

context("future future.vpa HS")

test_that("oututput value check (iteration for future sim is fixed as 2) ",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  SRdata <- get.SRdata(res.pma)

  HS.par0 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
  HS.par1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,hessian=FALSE)

  fres.HS <- future.vpa(res.pma,
                        multi=1, # res.pma$Fc.at.ageに掛ける乗数
                        nyear=50, # 将来予測の年数
                        start.year=2012, # 将来予測の開始年
                        N=2, # 確率的計算の繰り返し回数
                        ABC.year=2013, # ABCを計算する年
                        waa.year=2009:2011, # 生物パラメータの参照年
                        maa.year=2009:2011,
                        M.year=2009:2011,
                        is.plot=TRUE, # 結果をプロットするかどうか
                        seed=1,
                        silent=TRUE,
                        recfunc=HS.recAR, # 再生産関係の関数
                        # recfuncに対する引数
                        rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                     rho=HS.par0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                     sd=HS.par0$pars$sd,resid=HS.par0$resid)
                        )


  #上記引数での計算結果を読み込み
    # HS
  fresHS_faa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_faa1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_faa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_faa2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_faa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_faa3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_naa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_naa1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_naa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_naa2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_naa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_naa3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_biom1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_biom1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_biom2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_biom2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_biom3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_biom3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_baa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_baa1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_baa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_baa2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_baa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_baa3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_ssb1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_ssb1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_ssb2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_ssb2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_ssb3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_ssb3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_wcaa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_wcaa1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_wcaa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_wcaa2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_wcaa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_wcaa3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_caa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_caa1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_caa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_caa2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_caa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_caa3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_M1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_M1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_M2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_M2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_M3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_M3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_rps_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_rps_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_maa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_maa1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_maa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_maa2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_maa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_maa3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_vbiom_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_vbiom_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_recruit_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_recruit_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_eaa_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_eaa_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_alpha_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_alpha_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_thisyear_ssb_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_thisyear-ssb_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_waa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_waa1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_waa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_waa2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_waa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_waa3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_waacatch1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_waa-catch1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_waacatch2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_waa-catch2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_waacatch3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_waa-catch3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_currentF_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_currentF_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_vssb_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_vssb_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_vwcaa_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_vwcaa_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_naa_all1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_naa_all1_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_naa_all2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_naa_all2_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_naa_all3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_naa_all3_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_years_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_years_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_fyear_year_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_fyear_year_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_ABC_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_ABC_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_recarg_a_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_rec-arg_a_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_recarg_b_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_rec-arg_b_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_recarg_rho_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_rec-arg_rho_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_recarg_sd_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_rec-arg_sd_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_recarg_resid_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_rec-arg_resid_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_recarg_sd2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_rec-arg_sd2_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_waayear_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_waa_year_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_maayear_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_maa_year_pma_check.csv",package="frasyr"),row.names=1)

  fresHS_multi_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_multi_pma_check.csv",package="frasyr"),row.names=1)
  fresHS_multiyear_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_HS_multi_year_pma_check.csv",package="frasyr"),row.names=1)

    # HS.AR

  #結果の数値を照合
    # HS
  for( i in 1:nrow(fresHS_faa1_pma_check)){
    for( j in 1:ncol(fresHS_faa1_pma_check)){
      expect_equal(fres.HS$faa[i,j,1], fresHS_faa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_faa2_pma_check)){
    for( j in 1:ncol(fresHS_faa2_pma_check)){
      expect_equal(fres.HS$faa[i,j,2], fresHS_faa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_faa3_pma_check)){
    for( j in 1:ncol(fresHS_faa3_pma_check)){
      expect_equal(fres.HS$faa[i,j,3], fresHS_faa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_naa1_pma_check)){
    for( j in 1:ncol(fresHS_naa1_pma_check)){
      expect_equal(fres.HS$naa[i,j,1], fresHS_naa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_naa2_pma_check)){
    for( j in 1:ncol(fresHS_naa2_pma_check)){
      expect_equal(fres.HS$naa[i,j,2], fresHS_naa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_naa3_pma_check)){
    for( j in 1:ncol(fresHS_naa3_pma_check)){
      expect_equal(fres.HS$naa[i,j,3], fresHS_naa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_biom1_pma_check)){
    for( j in 1:ncol(fresHS_biom1_pma_check)){
      expect_equal(fres.HS$biom[i,j,1], fresHS_biom1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_biom2_pma_check)){
    for( j in 1:ncol(fresHS_biom2_pma_check)){
      expect_equal(fres.HS$biom[i,j,2], fresHS_biom2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_biom3_pma_check)){
    for( j in 1:ncol(fresHS_biom3_pma_check)){
      expect_equal(fres.HS$biom[i,j,3], fresHS_biom3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_baa1_pma_check)){
    for( j in 1:ncol(fresHS_baa1_pma_check)){
      expect_equal(fres.HS$baa[i,j,1], fresHS_baa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_baa2_pma_check)){
    for( j in 1:ncol(fresHS_baa2_pma_check)){
      expect_equal(fres.HS$baa[i,j,2], fresHS_baa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_baa3_pma_check)){
    for( j in 1:ncol(fresHS_baa3_pma_check)){
      expect_equal(fres.HS$baa[i,j,3], fresHS_baa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_ssb1_pma_check)){
    for( j in 1:ncol(fresHS_ssb1_pma_check)){
      expect_equal(fres.HS$ssb[i,j,1], fresHS_ssb1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_ssb2_pma_check)){
    for( j in 1:ncol(fresHS_ssb2_pma_check)){
      expect_equal(fres.HS$ssb[i,j,2], fresHS_ssb2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_ssb3_pma_check)){
    for( j in 1:ncol(fresHS_ssb3_pma_check)){
      expect_equal(fres.HS$ssb[i,j,3], fresHS_ssb3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_wcaa1_pma_check)){
    for( j in 1:ncol(fresHS_wcaa1_pma_check)){
      expect_equal(fres.HS$wcaa[i,j,1], fresHS_wcaa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_wcaa2_pma_check)){
    for( j in 1:ncol(fresHS_wcaa2_pma_check)){
      expect_equal(fres.HS$wcaa[i,j,2], fresHS_wcaa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_wcaa3_pma_check)){
    for( j in 1:ncol(fresHS_wcaa3_pma_check)){
      expect_equal(fres.HS$wcaa[i,j,3], fresHS_wcaa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_wcaa1_pma_check)){
    for( j in 1:ncol(fresHS_wcaa1_pma_check)){
      expect_equal(fres.HS$wcaa[i,j,1], fresHS_wcaa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_wcaa2_pma_check)){
    for( j in 1:ncol(fresHS_wcaa2_pma_check)){
      expect_equal(fres.HS$wcaa[i,j,2], fresHS_wcaa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_wcaa3_pma_check)){
    for( j in 1:ncol(fresHS_wcaa3_pma_check)){
      expect_equal(fres.HS$wcaa[i,j,3], fresHS_wcaa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_M1_pma_check)){
    for( j in 1:ncol(fresHS_M1_pma_check)){
      expect_equal(fres.HS$M[i,j,1], fresHS_M1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_M2_pma_check)){
    for( j in 1:ncol(fresHS_M2_pma_check)){
      expect_equal(fres.HS$M[i,j,2], fresHS_M2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_M3_pma_check)){
    for( j in 1:ncol(fresHS_M3_pma_check)){
      expect_equal(fres.HS$M[i,j,3], fresHS_M3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_rps_pma_check)){
    for( j in 1:ncol(fresHS_rps_pma_check)){
      expect_equal(fres.HS$rps[i,j], fresHS_rps_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_maa1_pma_check)){
    for( j in 1:ncol(fresHS_maa1_pma_check)){
      expect_equal(fres.HS$maa[i,j,1], fresHS_maa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_maa2_pma_check)){
    for( j in 1:ncol(fresHS_maa2_pma_check)){
      expect_equal(fres.HS$maa[i,j,2], fresHS_maa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_maa3_pma_check)){
    for( j in 1:ncol(fresHS_maa3_pma_check)){
      expect_equal(fres.HS$maa[i,j,3], fresHS_maa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_vbiom_pma_check)){
    for( j in 1:ncol(fresHS_vbiom_pma_check)){
      expect_equal(fres.HS$vbiom[i,j], fresHS_vbiom_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_recruit_pma_check)){
    for( j in 1:ncol(fresHS_recruit_pma_check)){
      expect_equal(fres.HS$recruit[i,j], fresHS_recruit_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_eaa_pma_check)){
    for( j in 1:ncol(fresHS_eaa_pma_check)){
      expect_equal(fres.HS$eaa[i,j], fresHS_eaa_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_alpha_pma_check)){
    for( j in 1:ncol(fresHS_alpha_pma_check)){
      expect_equal(fres.HS$alpha[i,j], fresHS_alpha_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_thisyear_ssb_pma_check)){
    for( j in 1:ncol(fresHS_thisyear_ssb_pma_check)){
      expect_equal(fres.HS$thisyear.ssb[i,j], fresHS_thisyear_ssb_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_waa1_pma_check)){
    for( j in 1:ncol(fresHS_waa1_pma_check)){
      expect_equal(fres.HS$waa[i,j,1], fresHS_waa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_waa2_pma_check)){
    for( j in 1:ncol(fresHS_waa2_pma_check)){
      expect_equal(fres.HS$waa[i,j,2], fresHS_waa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_waa3_pma_check)){
    for( j in 1:ncol(fresHS_waa3_pma_check)){
      expect_equal(fres.HS$waa[i,j,3], fresHS_waa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_waacatch1_pma_check)){
    for( j in 1:ncol(fresHS_waacatch1_pma_check)){
      expect_equal(fres.HS$waa.catch[i,j,1], fresHS_waacatch1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_waacatch2_pma_check)){
    for( j in 1:ncol(fresHS_waacatch2_pma_check)){
      expect_equal(fres.HS$waa.catch[i,j,2], fresHS_waacatch2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_waacatch3_pma_check)){
    for( j in 1:ncol(fresHS_waacatch3_pma_check)){
      expect_equal(fres.HS$waa.catch[i,j,3], fresHS_waacatch3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_currentF_pma_check)){
    expect_equal(as.numeric(fres.HS$currentF[i]), fresHS_currentF_pma_check[i,])
  }

  for( i in 1:nrow(fresHS_vssb_pma_check)){
    for( j in 1:ncol(fresHS_vssb_pma_check)){
      expect_equal(fres.HS$vssb[i,j], fresHS_vssb_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_vwcaa_pma_check)){
    for( j in 1:ncol(fresHS_vwcaa_pma_check)){
      expect_equal(fres.HS$vwcaa[i,j], fresHS_vwcaa_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_naa_all1_pma_check)){
    for( j in 1:ncol(fresHS_naa_all1_pma_check)){
      expect_equal(fres.HS$naa_all[i,j,1], fresHS_naa_all1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_naa_all2_pma_check)){
    for( j in 1:ncol(fresHS_naa_all2_pma_check)){
      expect_equal(fres.HS$naa_all[i,j,2], fresHS_naa_all2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresHS_naa_all3_pma_check)){
    for( j in 1:ncol(fresHS_naa_all3_pma_check)){
      expect_equal(fres.HS$naa_all[i,j,3], fresHS_naa_all3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresHS_years_pma_check)){
    expect_equal(as.numeric(fres.HS$years[i]), fresHS_years_pma_check[i,])
  }

  for( i in 1:nrow(fresHS_fyear_year_pma_check)){
    expect_equal(as.numeric(fres.HS$fyear.year[i]), fresHS_fyear_year_pma_check[i,])
  }

  for( i in 1:nrow(fresHS_ABC_pma_check)){
    expect_equal(as.numeric(fres.HS$ABC[i]), fresHS_ABC_pma_check[i,])
  }

  expect_equal(fres.HS$rec.arg$a, fresHS_recarg_a_pma_check[1,])
  expect_equal(fres.HS$rec.arg$b, fresHS_recarg_b_pma_check[1,])
  expect_equal(fres.HS$rec.arg$rho, fresHS_recarg_rho_pma_check[1,])
  for( i in 1:nrow(fresHS_recarg_sd_pma_check)){
    expect_equal(fres.HS$rec.arg$sd[i], as.numeric(fresHS_recarg_sd_pma_check[i,]))
  }
  for( i in 1:nrow(fresHS_recarg_resid_pma_check)){
    expect_equal(as.numeric(fres.HS$rec.arg$resid[i]), fresHS_recarg_resid_pma_check[i,])
  }
  for( i in 1:nrow(fresHS_recarg_sd2_pma_check)){
    expect_equal(fres.HS$rec.arg$sd2[i], as.numeric(fresHS_recarg_sd2_pma_check[i,]))
  }
  for( i in 1:nrow(fresHS_waayear_pma_check)){
    expect_equal(fres.HS$waa.year[i], as.numeric(fresHS_waayear_pma_check[i,]))
  }
  for( i in 1:nrow(fresHS_maayear_pma_check)){
    expect_equal(fres.HS$maa.year[i], as.numeric(fresHS_maayear_pma_check[i,]))
  }

  expect_equal(fres.HS$multi, fresHS_multi_pma_check[1,])
  expect_equal(fres.HS$multi.year, fresHS_multiyear_pma_check[1,])

    # HS.AR

})

context("future future.vpa BH")

test_that("oututput value check (iteration for future sim is fixed as 2) ",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  SRdata <- get.SRdata(res.pma)

  BH.par0 <- fit.SR(SRdata,SR="BH",method="L2",AR=0,hessian=FALSE)
  BH.par1 <- fit.SR(SRdata,SR="BH",method="L2",AR=1,hessian=FALSE)

  fres.BH <- future.vpa(res.pma,
                        multi=1,
                        nyear=50, # 将来予測の年数
                        start.year=2012, # 将来予測の開始年
                        N=2, # 確率的計算の繰り返し回数
                        ABC.year=2013, # ABCを計算する年
                        waa.year=2009:2011, # 生物パラメータの参照年
                        maa.year=2009:2011,
                        M.year=2009:2011,
                        is.plot=TRUE, # 結果をプロットするかどうか
                        seed=1,
                        silent=TRUE,
                        recfunc=BH.recAR, # 再生産関係の関数
                        # recfuncに対する引数
                        rec.arg=list(a=BH.par0$pars$a,b=BH.par0$pars$b,rho=BH.par0$rho,
                                     sd=BH.par0$pars$sd,resid=BH.par0$resid))
  #上記引数での計算結果を読み込み
  # BH
  fresBH_faa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_faa1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_faa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_faa2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_faa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_faa3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_naa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_naa1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_naa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_naa2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_naa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_naa3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_biom1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_biom1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_biom2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_biom2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_biom3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_biom3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_baa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_baa1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_baa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_baa2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_baa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_baa3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_ssb1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_ssb1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_ssb2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_ssb2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_ssb3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_ssb3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_wcaa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_wcaa1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_wcaa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_wcaa2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_wcaa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_wcaa3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_caa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_caa1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_caa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_caa2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_caa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_caa3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_M1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_M1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_M2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_M2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_M3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_M3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_rps_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_rps_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_maa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_maa1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_maa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_maa2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_maa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_maa3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_vbiom_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_vbiom_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_recruit_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_recruit_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_eaa_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_eaa_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_alpha_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_alpha_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_thisyear_ssb_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_thisyear-ssb_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_waa1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_waa1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_waa2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_waa2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_waa3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_waa3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_waacatch1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_waa-catch1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_waacatch2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_waa-catch2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_waacatch3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_waa-catch3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_currentF_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_currentF_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_vssb_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_vssb_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_vwcaa_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_vwcaa_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_naa_all1_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_naa_all1_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_naa_all2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_naa_all2_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_naa_all3_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_naa_all3_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_years_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_years_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_fyear_year_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_fyear_year_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_ABC_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_ABC_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_recarg_a_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_rec-arg_a_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_recarg_b_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_rec-arg_b_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_recarg_rho_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_rec-arg_rho_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_recarg_sd_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_rec-arg_sd_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_recarg_resid_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_rec-arg_resid_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_recarg_sd2_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_rec-arg_sd2_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_waayear_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_waa_year_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_maayear_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_maa_year_pma_check.csv",package="frasyr"),row.names=1)

  fresBH_multi_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_multi_pma_check.csv",package="frasyr"),row.names=1)
  fresBH_multiyear_pma_check <- read.csv(system.file("extdata", "future_vpa_fres_BH_multi_year_pma_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合
  # BH

  for( i in 1:nrow(fresBH_faa1_pma_check)){
    for( j in 1:ncol(fresBH_faa1_pma_check)){
      expect_equal(fres.BH$faa[i,j,1], fresBH_faa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_faa2_pma_check)){
    for( j in 1:ncol(fresBH_faa2_pma_check)){
      expect_equal(fres.BH$faa[i,j,2], fresBH_faa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_faa3_pma_check)){
    for( j in 1:ncol(fresBH_faa3_pma_check)){
      expect_equal(fres.BH$faa[i,j,3], fresBH_faa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_naa1_pma_check)){
    for( j in 1:ncol(fresBH_naa1_pma_check)){
      expect_equal(fres.BH$naa[i,j,1], fresBH_naa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_naa2_pma_check)){
    for( j in 1:ncol(fresBH_naa2_pma_check)){
      expect_equal(fres.BH$naa[i,j,2], fresBH_naa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_naa3_pma_check)){
    for( j in 1:ncol(fresBH_naa3_pma_check)){
      expect_equal(fres.BH$naa[i,j,3], fresBH_naa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_biom1_pma_check)){
    for( j in 1:ncol(fresBH_biom1_pma_check)){
      expect_equal(fres.BH$biom[i,j,1], fresBH_biom1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_biom2_pma_check)){
    for( j in 1:ncol(fresBH_biom2_pma_check)){
      expect_equal(fres.BH$biom[i,j,2], fresBH_biom2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_biom3_pma_check)){
    for( j in 1:ncol(fresBH_biom3_pma_check)){
      expect_equal(fres.BH$biom[i,j,3], fresBH_biom3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_baa1_pma_check)){
    for( j in 1:ncol(fresBH_baa1_pma_check)){
      expect_equal(fres.BH$baa[i,j,1], fresBH_baa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_baa2_pma_check)){
    for( j in 1:ncol(fresBH_baa2_pma_check)){
      expect_equal(fres.BH$baa[i,j,2], fresBH_baa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_baa3_pma_check)){
    for( j in 1:ncol(fresBH_baa3_pma_check)){
      expect_equal(fres.BH$baa[i,j,3], fresBH_baa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_ssb1_pma_check)){
    for( j in 1:ncol(fresBH_ssb1_pma_check)){
      expect_equal(fres.BH$ssb[i,j,1], fresBH_ssb1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_ssb2_pma_check)){
    for( j in 1:ncol(fresBH_ssb2_pma_check)){
      expect_equal(fres.BH$ssb[i,j,2], fresBH_ssb2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_ssb3_pma_check)){
    for( j in 1:ncol(fresBH_ssb3_pma_check)){
      expect_equal(fres.BH$ssb[i,j,3], fresBH_ssb3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_wcaa1_pma_check)){
    for( j in 1:ncol(fresBH_wcaa1_pma_check)){
      expect_equal(fres.BH$wcaa[i,j,1], fresBH_wcaa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_wcaa2_pma_check)){
    for( j in 1:ncol(fresBH_wcaa2_pma_check)){
      expect_equal(fres.BH$wcaa[i,j,2], fresBH_wcaa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_wcaa3_pma_check)){
    for( j in 1:ncol(fresBH_wcaa3_pma_check)){
      expect_equal(fres.BH$wcaa[i,j,3], fresBH_wcaa3_pma_check[i,j])
    }
  }


  for( i in 1:nrow(fresBH_wcaa1_pma_check)){
    for( j in 1:ncol(fresBH_wcaa1_pma_check)){
      expect_equal(fres.BH$wcaa[i,j,1], fresBH_wcaa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_wcaa2_pma_check)){
    for( j in 1:ncol(fresBH_wcaa2_pma_check)){
      expect_equal(fres.BH$wcaa[i,j,2], fresBH_wcaa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_wcaa3_pma_check)){
    for( j in 1:ncol(fresBH_wcaa3_pma_check)){
      expect_equal(fres.BH$wcaa[i,j,3], fresBH_wcaa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_M1_pma_check)){
    for( j in 1:ncol(fresBH_M1_pma_check)){
      expect_equal(fres.BH$M[i,j,1], fresBH_M1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_M2_pma_check)){
    for( j in 1:ncol(fresBH_M2_pma_check)){
      expect_equal(fres.BH$M[i,j,2], fresBH_M2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_M3_pma_check)){
    for( j in 1:ncol(fresBH_M3_pma_check)){
      expect_equal(fres.BH$M[i,j,3], fresBH_M3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_rps_pma_check)){
    for( j in 1:ncol(fresBH_rps_pma_check)){
      expect_equal(fres.BH$rps[i,j], fresBH_rps_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_maa1_pma_check)){
    for( j in 1:ncol(fresBH_maa1_pma_check)){
      expect_equal(fres.BH$maa[i,j,1], fresBH_maa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_maa2_pma_check)){
    for( j in 1:ncol(fresBH_maa2_pma_check)){
      expect_equal(fres.BH$maa[i,j,2], fresBH_maa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_maa3_pma_check)){
    for( j in 1:ncol(fresBH_maa3_pma_check)){
      expect_equal(fres.BH$maa[i,j,3], fresBH_maa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_vbiom_pma_check)){
    for( j in 1:ncol(fresBH_vbiom_pma_check)){
      expect_equal(fres.BH$vbiom[i,j], fresBH_vbiom_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_recruit_pma_check)){
    for( j in 1:ncol(fresBH_recruit_pma_check)){
      expect_equal(fres.BH$recruit[i,j], fresBH_recruit_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_eaa_pma_check)){
    for( j in 1:ncol(fresBH_eaa_pma_check)){
      expect_equal(fres.BH$eaa[i,j], fresBH_eaa_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_alpha_pma_check)){
    for( j in 1:ncol(fresBH_alpha_pma_check)){
      expect_equal(fres.BH$alpha[i,j], fresBH_alpha_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_thisyear_ssb_pma_check)){
    for( j in 1:ncol(fresBH_thisyear_ssb_pma_check)){
      expect_equal(fres.BH$thisyear.ssb[i,j], fresBH_thisyear_ssb_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_waa1_pma_check)){
    for( j in 1:ncol(fresBH_waa1_pma_check)){
      expect_equal(fres.BH$waa[i,j,1], fresBH_waa1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_waa2_pma_check)){
    for( j in 1:ncol(fresBH_waa2_pma_check)){
      expect_equal(fres.BH$waa[i,j,2], fresBH_waa2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_waa3_pma_check)){
    for( j in 1:ncol(fresBH_waa3_pma_check)){
      expect_equal(fres.BH$waa[i,j,3], fresBH_waa3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_waacatch1_pma_check)){
    for( j in 1:ncol(fresBH_waacatch1_pma_check)){
      expect_equal(fres.BH$waa.catch[i,j,1], fresBH_waacatch1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_waacatch2_pma_check)){
    for( j in 1:ncol(fresBH_waacatch2_pma_check)){
      expect_equal(fres.BH$waa.catch[i,j,2], fresBH_waacatch2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_waacatch3_pma_check)){
    for( j in 1:ncol(fresBH_waacatch3_pma_check)){
      expect_equal(fres.BH$waa.catch[i,j,3], fresBH_waacatch3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_currentF_pma_check)){
    expect_equal(as.numeric(fres.BH$currentF[i]), fresBH_currentF_pma_check[i,])
  }

  for( i in 1:nrow(fresBH_vssb_pma_check)){
    for( j in 1:ncol(fresBH_vssb_pma_check)){
      expect_equal(fres.BH$vssb[i,j], fresBH_vssb_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_vwcaa_pma_check)){
    for( j in 1:ncol(fresBH_vwcaa_pma_check)){
      expect_equal(fres.BH$vwcaa[i,j], fresBH_vwcaa_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_naa_all1_pma_check)){
    for( j in 1:ncol(fresBH_naa_all1_pma_check)){
      expect_equal(fres.BH$naa_all[i,j,1], fresBH_naa_all1_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_naa_all2_pma_check)){
    for( j in 1:ncol(fresBH_naa_all2_pma_check)){
      expect_equal(fres.BH$naa_all[i,j,2], fresBH_naa_all2_pma_check[i,j])
    }
  }
  for( i in 1:nrow(fresBH_naa_all3_pma_check)){
    for( j in 1:ncol(fresBH_naa_all3_pma_check)){
      expect_equal(fres.BH$naa_all[i,j,3], fresBH_naa_all3_pma_check[i,j])
    }
  }

  for( i in 1:nrow(fresBH_years_pma_check)){
    expect_equal(as.numeric(fres.BH$years[i]), fresBH_years_pma_check[i,])
  }

  for( i in 1:nrow(fresBH_fyear_year_pma_check)){
    expect_equal(as.numeric(fres.BH$fyear.year[i]), fresBH_fyear_year_pma_check[i,])
  }

  for( i in 1:nrow(fresBH_ABC_pma_check)){
    expect_equal(as.numeric(fres.BH$ABC[i]), fresBH_ABC_pma_check[i,])
  }

  expect_equal(fres.BH$rec.arg$a, fresBH_recarg_a_pma_check[1,])
  expect_equal(fres.BH$rec.arg$b, fresBH_recarg_b_pma_check[1,])
  expect_equal(fres.BH$rec.arg$rho, fresBH_recarg_rho_pma_check[1,])
  for( i in 1:nrow(fresBH_recarg_sd_pma_check)){
    expect_equal(fres.BH$rec.arg$sd[i], as.numeric(fresBH_recarg_sd_pma_check[i,]))
  }
  for( i in 1:nrow(fresBH_recarg_resid_pma_check)){
    expect_equal(as.numeric(fres.BH$rec.arg$resid[i]), fresBH_recarg_resid_pma_check[i,])
  }
  for( i in 1:nrow(fresBH_recarg_sd2_pma_check)){
    expect_equal(fres.BH$rec.arg$sd2[i], as.numeric(fresBH_recarg_sd2_pma_check[i,]))
  }
  for( i in 1:nrow(fresBH_waayear_pma_check)){
    expect_equal(fres.BH$waa.year[i], as.numeric(fresBH_waayear_pma_check[i,]))
  }
  for( i in 1:nrow(fresBH_maayear_pma_check)){
    expect_equal(fres.BH$maa.year[i], as.numeric(fresBH_maayear_pma_check[i,]))
  }

  expect_equal(fres.BH$multi, fresBH_multi_pma_check[1,])
  expect_equal(fres.BH$multi.year, fresBH_multiyear_pma_check[1,])

})

context("future future.vpa (option of futureF)")

test_that("oututput value check (iteration for future sim is fixed as 2) ",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  SRdata <- get.SRdata(res.pma)

  HS.par0 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
  HS.par1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,hessian=FALSE)

  currentF.test <- 1:4/10
  futureF.test <- 5:8/10
  fres.HS.check <- future.vpa(res.pma,
                              multi=2,
                              currentF=currentF.test,
                              futureF=5:8/10,
                              nyear=50, 
                              start.year=2012, 
                              N=2, ABC.year=2013, 
                        waa.year=2009:2011, 
                        maa.year=2009:2011,
                        M.year=2009:2011,
                        is.plot=TRUE, 
                        seed=1,
                        silent=TRUE,
                        recfunc=HS.recAR, 
                        rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                     rho=HS.par0$pars$rho, 
                                     sd=HS.par0$pars$sd,resid=HS.par0$resid)
                        )


  # faaが想定通りに入っていればOKくらいのテストです
  for(i in 1:4){
      expect_true(mean(fres.HS.check$faa[i,dimnames(fres.HS.check$faa)[[2]]==2012,])==currentF.test[i])
      expect_true(mean(fres.HS.check$faa[i,dimnames(fres.HS.check$faa)[[2]]>2012,])==futureF.test[i]*2)      
  }
})


