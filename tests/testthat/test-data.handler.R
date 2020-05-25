context("generate vpa results with various data type (with dummy data)") 

test_that("future_vpa function (with dummy vpa data) (level 2-3?)",{
  # read various data ----
  # data with caa=maa=waa=1, M=0
  data_base <- readr::read_csv(system.file("extdata","all_dummy_data_base.csv",package="frasyr")) 
  # data with caa=maa=waa=1, M=0 but plus group have changed 
  data_pgc <- readr::read_csv(system.file("extdata","all_dummy_data_plus_group_change.csv",package="frasyr")) 
  # data with caa=maa=waa=1, M=0 but first recruit age is 1
  data_rec <- readr::read_csv(system.file("extdata","all_dummy_data_rec.csv",package="frasyr")) 

  # create various vpa data ----
  vpadat_base0 <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
                               waa=to_vpa_data(data_base, label_name="waa"),
                               maa=to_vpa_data(data_base, label_name="maa"),
                               M  = 0,
                               index = to_vpa_data(data_base, label_name="abund"),
                               maa.tune = NULL,
                               waa.catch = NULL,
                               catch.prop = NULL)
  # waa.catch is given 
  vpadat_base1 <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
                               waa=to_vpa_data(data_base, label_name="waa"),
                               maa=to_vpa_data(data_base, label_name="maa"),
                               M  = 0,
                               index = to_vpa_data(data_base, label_name="abund"),
                               maa.tune = NULL,
                               waa.catch = to_vpa_data(data_base, label_name="waa")*2,
                               catch.prop = NULL)
  
  vpadat_pgc0 <- data.handler(caa=to_vpa_data(data_pgc, label_name="caa"),
                              waa=to_vpa_data(data_pgc, label_name="waa"),
                              maa=to_vpa_data(data_pgc, label_name="maa"),
                              M  = 0,
                              index = to_vpa_data(data_pgc, label_name="abund"),
                              maa.tune = NULL,
                              waa.catch = NULL,
                              catch.prop = NULL)
  
  vpadat_rec0 <- data.handler(caa=to_vpa_data(data_rec, label_name="caa"),
                              waa=to_vpa_data(data_rec, label_name="waa"),
                              maa=to_vpa_data(data_rec, label_name="maa"),
                              M  = 0,
                              index = to_vpa_data(data_rec, label_name="abund"),
                              maa.tune = NULL,
                              waa.catch = NULL,
                              catch.prop = NULL)

  # vpa (no tuning) ----

  # 普通のVPAの場合、0-3歳の尾数は4,3,2,2になる。それをtrue.numberとして格納しておく
  true_number <- c(4,3,2,2)
  # また、indexは左右に対照なデータになるので、sd=mean(log(CPUE))になる
  true_sd <- vpadat_base0$index %>% log() %>% unlist() %>% mean(na.rm=TRUE)
  
  res_vpa_base0_nontune <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_base0_nontune$naa)), 
               true_number)
  
  res_vpa_base1_nontune <- vpa(vpadat_base1, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_base1_nontune$naa)), 
               true_number)

  # プラスグループが変わる場合はtrue_numberには一致しない；どうテストすべきか？
  res_vpa_pgc0_nontune <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(unlist(res_vpa_pgc0_nontune$naa["2017"])), 
               c(3,2,2,NA), tol=0.0001)
  
  res_vpa_rec0_nontune <- vpa(vpadat_rec0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_rec0_nontune$naa)), 
               true_number)
  
  # catch計算用のwaaを２倍にしているbase1データでは漁獲量が倍になる
  expect_equal(res_vpa_base0_nontune$wcaa*2,
               res_vpa_base1_nontune$wcaa)
  
  # Part1: dataset is "vpadat_base0" ----
  # 1-1: 二段階法によるtuning ----
 
  sel.f1<-res_vpa_base0_nontune$saa$"2017"
  round(sel.f1,3)
  
  #二段階法：est.method=最小二乗法による推定
  res_vpa_base0_tune1l <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ls", b.est=FALSE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(rowMeans(res_vpa_base0_tune1l$naa)),true_number,tol=0.0001)
  expect_equal(as.numeric(res_vpa_base0_tune1l$sigma),true_sd, tol=0.0001)

  #二段階法：est.method=最尤法による推定
  res_vpa_base0_tune1m <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ml", b.est=FALSE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(rowMeans(res_vpa_base0_tune1m$naa)),true_number, tol=0.0001)
  expect_equal(as.numeric(res_vpa_base0_tune1m$sigma),rep(true_sd,2), tol=0.0001)

  #二段階法：est.method=最尤法による推定(2つの指数のsdが異なる場合）
  vpadat_base0_index_change <- vpadat_base0
  vpadat_base0_index_change$index[1,is.na(vpadat_base0_index_change$index[1,])] <- exp(mean(log(c(1,2))))
  res_vpa_base0_index_change_tune1m <- vpa(vpadat_base0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ml", b.est=FALSE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(round(as.numeric(rowMeans(res_vpa_base0_index_change_tune1m$naa)),2),
               true_number+c(0.04,0.02,0.01,0.01),
               tol=0.01)
  expect_equal(round(as.numeric(res_vpa_base0_index_change_tune1m$sigma),2),
               c(0.25,0.35))  
  
  #二段階法：est.method=最小二乗法による推定＋指標値の非線形性bの推定
  # 現状のvpaでのb,sd,naa推定値
  b_est_tmp <- c(-0.31,0.31)
  sd_est_tmp <- 0.34
  naa_est_tmp<-c(3.70,2.83,1.91,1.91)
  
  res_vpa_base0_tune1l_b <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ls", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune1l_b$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune1l_b$b,2)),b_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune1l_b$sigma,2)),sd_est_tmp)
  
  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定
  res_vpa_base0_tune1m_b <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune1m_b$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune1m_b$b,2)),b_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune1m_b$sigma,2)),rep(sd_est_tmp,2))
  
  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定＋一部のbは固定
  res_vpa_base0_tune1m_b_fix <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA),b.fix=c(NA,-0.1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune1m_b_fix$naa),2)),
               c(5.20,3.73,2.39,2.39))
  expect_equal(as.numeric(round(res_vpa_base0_tune1m_b_fix$b,2)),c(0.14,-0.10))
  expect_equal(as.numeric(round(res_vpa_base0_tune1m_b_fix$sigma,2)),rep(sd_est_tmp,2))

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定(2つの指数のsdが異なる場合）
  res_vpa_base0_index_change_tune1m_b <- vpa(vpadat_base0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_index_change_tune1m_b$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune1m_b$b,2)),c(-0.28,0.31))
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune1m_b$sigma,2)),c(0.25,0.34))
  
  
  #1-2: 選択率更新法によるtuning ----
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最小二乗法
  res_vpa_base0_tune2l <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ls", b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune2l$sigma,2)),true_sd,tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l$saa),2)),c(0.42,0.58,1.00,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法
  res_vpa_base0_tune2m <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune2m$sigma,2)),rep(true_sd,2), tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m$saa),2)),c(0.42,0.58,1.00,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最小二乗法．b推定する
  res_vpa_base0_tune2l_b <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ls", b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l_b$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune2l_b$b,2)),b_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune2l_b$sigma,2)),sd_est_tmp)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l_b$saa),2)),c(0.44,0.59,1.00,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する
  res_vpa_base0_tune2m_b <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m_b$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune2m_b$b,2)),b_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune2m_b$sigma,2)),rep(sd_est_tmp,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m_b$saa),2)),c(0.44,0.59,1.00,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する(2つの指数のsdが異なる場合)
  res_vpa_base0_index_change_tune2m_b <- vpa(vpadat_base0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_index_change_tune2m_b$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune2m_b$b,2)),c(-0.28,0.31))
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune2m_b$sigma,2)),c(0.25,sd_est_tmp))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_index_change_tune2m_b$saa),2)),c(0.44,0.59,1.00,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する．1部のｂは固定
  res_vpa_base0_tune2m_b_fix <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),b.fix=c(NA,-0.1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m_b_fix$naa),2)),c(5.48,3.79,2.40,2.40))
  expect_equal(as.numeric(round(res_vpa_base0_tune2m_b_fix$b,2)),c(0.14,-0.10))
  expect_equal(as.numeric(round(res_vpa_base0_tune2m_b_fix$sigma,2)),rep(sd_est_tmp,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m_b_fix$saa),2)),c(0.40,0.58,1.00,1.00))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最小二乗法．
  res_vpa_base0_tune3l <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ls" ,b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3l$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune3l$sigma,2)),c(0.42))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3l$saa),2)),c(0.14,0.19,0.33,0.33))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法
  res_vpa_base0_tune3m <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3m$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune3m$sigma,2)),c(0.46,0.38))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3m$saa),2)),c(0.14,0.19,0.33,0.33))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最小二乗法．ｂを推定する
  res_vpa_base0_tune3lb <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ls" ,b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3lb$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune3lb$b,2)),b_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune3lb$sigma,2)),sd_est_tmp)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3lb$saa),2)),c(0.14,0.19,0.33,0.33))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する
  res_vpa_base0_tune3mb <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3mb$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune3mb$b,2)),b_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_tune3mb$sigma,2)),rep(sd_est_tmp,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3mb$saa),2)),c(0.14,0.19,0.33,0.33))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する(2つの指数のsdが異なる場合)
  res_vpa_base0_index_change_tune3mb <- vpa(vpadat_base0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_index_change_tune3mb$naa),2)),naa_est_tmp)
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune3mb$b,2)),c(-0.28,0.31))
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune3mb$sigma,2)),c(0.25,sd_est_tmp))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_index_change_tune3mb$saa),2)),c(0.14,0.19,0.33,0.33))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する．一部のbを固定する
  res_vpa_base0_tune3mb_fix <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=TRUE,abund=c("SSB","SSB"),b.fix=c(NA,-0.1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3mb_fix $naa),2)),c(5.48,3.79,2.40,2.40))
  expect_equal(as.numeric(round(res_vpa_base0_tune3mb_fix $b,2)),c(0.14,-0.10))
  expect_equal(as.numeric(round(res_vpa_base0_tune3mb_fix $sigma,2)),rep(sd_est_tmp,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3mb_fix $saa),2)),c(0.13,0.19,0.34,0.34))
  
  
  #1-3: 全F法によるtuning ----
  
  #  全F推定法．最小二乗法．
  res_vpa_base0_tune4l <- vpa(vpadat_base0,  last.catch.zero = FALSE,tf.year=2015:2016,
                              Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ls" ,b.est=FALSE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4l$naa),2)),c(4.07,3.00,2.00,2.00))
  expect_equal(as.numeric(round(res_vpa_base0_tune4l$sigma,2)),c(0.35))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4l$saa),2)),c(0.14,0.20,0.33,0.33))
  
  #  全F推定法．最尤法
  res_vpa_base0_tune4m <- vpa(vpadat_base0, last.catch.zero = FALSE, 
                              Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=FALSE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m$naa),2)),c(4.07,3.00,2.00,2.00))
  expect_equal(as.numeric(round(res_vpa_base0_tune4m$sigma,2)),rep(0.35,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m$saa),2)),c(0.14,0.20,0.33,0.33))
  
  #  全F推定法．最小二乗法．b推定あり
  res_vpa_base0_tune4l_b <- vpa(vpadat_base0, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ls" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4l_b$naa),2)),c(4.32,3.26,1.91,1.91))
  expect_equal(as.numeric(round(res_vpa_base0_tune4l_b$b,2)),c(0.56,-0.56))
  expect_equal(as.numeric(round(res_vpa_base0_tune4l_b$sigma,2)),c(0.33))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4l_b$saa),2)),c(0.13,0.19,0.34,0.34))
  
  #  全F推定法．最尤法．b推定あり
  res_vpa_base0_tune4m_b <- vpa(vpadat_base0, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m_b$naa),2)),c(4.32,3.26,1.91,1.91))
  expect_equal(as.numeric(round(res_vpa_base0_tune4m_b$b,2)),c(0.56,-0.56))
  expect_equal(as.numeric(round(res_vpa_base0_tune4m_b$sigma,2)),rep(0.33,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m_b$saa),2)),c(0.13,0.19,0.34,0.34))
  
  #  全F推定法．最尤法．b推定あり(2つの指数のsdが異なる場合)
  res_vpa_base0_index_change_tune4m_b <- vpa(vpadat_base0_index_change, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_index_change_tune4m_b$naa),2)),c(4.33,3.26,1.91,1.91))
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune4m_b$b,2)),c(0.56,-0.56))
  expect_equal(as.numeric(round(res_vpa_base0_index_change_tune4m_b$sigma,2)),c(0.24,0.33))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_index_change_tune4m_b$saa),2)),c(0.13,0.19,0.34,0.34))
  
  # 　全F推定法．最尤法．b推定あり， Pope=FALSE
  res_vpa_base0_tune4m_b_baranov <- vpa(vpadat_base0, last.catch.zero = FALSE, 
                                Pope = FALSE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m_b_baranov$naa),2)),c(4.32,3.26,1.91,1.91))
  expect_equal(as.numeric(round(res_vpa_base0_tune4m_b_baranov$b,2)),c(0.56,-0.56))
  expect_equal(as.numeric(round(res_vpa_base0_tune4m_b_baranov$sigma,2)),rep(0.33,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m_b_baranov$saa),2)),c(0.13,0.19,0.34,0.34))
  
  # 　全F推定法．最尤法．b推定あり，alpha=0.3(最高齢と最高齢―1のFの比：Fa=alpha*Fa-1)
  res_vpa_base0_tune4m_b_alpha <- vpa(vpadat_base0, last.catch.zero = FALSE, 
                                        Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"),alpha=0.3)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m_b_alpha$naa),2)),c(4.09,3.03,1.98,2.03))
  expect_equal(as.numeric(round(res_vpa_base0_tune4m_b_alpha$b,2)),c(4.68,-4.68))
  expect_equal(as.numeric(round(res_vpa_base0_tune4m_b_alpha$sigma,2)),rep(0.31,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m_b_alpha$saa),2)),c(0.18,0.25,0.44,0.13))

  
  #1-4: Ridge VPA ----
  
  # set lambda + 全F推定法．最尤法．b推定あり
  res_vpa_base0_tune5m_b <- vpa(vpadat_base0, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"), lambda=0.02)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune5m_b$naa),2)),c(13.41,3.09,1.96,1.96))
  expect_equal(as.numeric(round(res_vpa_base0_tune5m_b$b,2)),c(1.38,-1.38))
  expect_equal(as.numeric(round(res_vpa_base0_tune5m_b$sigma,2)),rep(0.33,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune5m_b$saa),2)),c(0.13,0.19,0.34,0.34))
  
  # TMB true + set lambda + 全F推定法．最尤法．b推定あり
  library(TMB)
  use_rvpa_tmb()
  res_vpa_base0_tune6m_b <- vpa(vpadat_base0, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"), lambda=0.01, TMB=TRUE)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune6m_b$naa),2)),c(15.68,3.13,1.95,1.95))
  expect_equal(as.numeric(round(res_vpa_base0_tune6m_b$b,2)),c(0.95,-0.95))
  expect_equal(as.numeric(round(res_vpa_base0_tune6m_b$sigma,2)),rep(0.33,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune6m_b$saa),2)),c(0.13,0.19,0.34,0.34))
  
  
  # Part2: dataset is "vpadat_pgc0" ----
  # 2-1: 二段階法によるtuning ---- 
  
  b_est_tmp2 <- c(-0.17,0.17)
 
  naa_est_tmp2<-c(3.55,2.48,2.14,NA)
  
  sel.f1<-res_vpa_pgc0_nontune$saa$"2017"
  round(sel.f1,3)
  
  #二段階法：est.method=最小二乗法による推定
  res_vpa_pgc0_tune1l <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ls", b.est=FALSE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1l$naa),2)),naa_est_tmp2)
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1l$sigma,2)),0.4)
  
  #二段階法：est.method=最尤法による推定
  res_vpa_pgc0_tune1m <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=FALSE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1m$naa),2)),naa_est_tmp2)
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1m$sigma,2)),c(0.41,0.40))
  
  #二段階法：est.method=最尤法による推定(2つの指数のsdが異なる場合）
  vpa_pgc0_index_change <- vpadat_pgc0
  vpa_pgc0_index_change$index[1,is.na(vpa_pgc0_index_change$index[1,])] <- exp(mean(log(c(1,2))))
  res_vpa_pgc0_index_change_tune1m <- vpa(vpa_pgc0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                                           Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=FALSE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(round(as.numeric(rowMeans(res_vpa_pgc0_index_change_tune1m$naa)),2),
               c(3.58,2.49,2.16,NA))
  expect_equal(round(as.numeric(res_vpa_pgc0_index_change_tune1m$sigma),2),
               c(0.33,0.40))  
  
  #二段階法：est.method=最小二乗法による推定＋指標値の非線形性bの推定
  res_vpa_pgc0_tune1l_b <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ls", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1l_b$naa),2)),c(3.22,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1l_b$b,2)),b_est_tmp2)
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1l_b$sigma,2)),0.34)
  
  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定
  res_vpa_pgc0_tune1m_b <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1m_b$naa),2)),c(3.22,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1m_b$b,2)),b_est_tmp2)
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1m_b$sigma,2)),rep(0.34,2))
  
  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定＋一部のbは固定
  res_vpa_pgc0_tune1m_b_fix <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                    Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA),b.fix=c(NA,-0.1))
  expect_equal(as.numeric(round(rowMeans( res_vpa_pgc0_tune1m_b_fix$naa),2)),
               c(5.53,3.54,3.20,NA))
  expect_equal(as.numeric(round( res_vpa_pgc0_tune1m_b_fix$b,2)),c(0.09,-0.10))
  expect_equal(as.numeric(round( res_vpa_pgc0_tune1m_b_fix$sigma,2)),rep(0.34,2))
  
  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定(2つの指数のsdが異なる場合）
  res_vpa_pgc0_index_change_tune1m_b <- vpa(vpa_pgc0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                                             Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),min.age=c(NA,NA),max.age=c(NA,NA))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_index_change_tune1m_b$naa),2)),c(3.22,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_index_change_tune1m_b$b,2)),c(-0.11,0.17))
  expect_equal(as.numeric(round(res_vpa_pgc0_index_change_tune1m_b$sigma,2)),c(0.25,0.34))
  
  #2-2: 選択率更新法によるtuning ----
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最小二乗法
  res_vpa_pgc0_tune2l <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ls", b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2l$naa),2)),c(3.56,2.48,2.14,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2l$sigma,2)),0.4)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2l$saa),2)),c(0.39,0.62,0.76,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法
  res_vpa_pgc0_tune2m <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2m$naa),2)),c(3.56,2.48,2.14,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2m$sigma,2)),c(0.41,0.40))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2m$saa),2)),c(0.39,0.62,0.76,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最小二乗法．b推定する
  res_vpa_pgc0_tune2l_b <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ls", b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2l_b$naa),2)),c(7.14,4.28,3.95,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2l_b$b,2)),c(0.07,-0.07))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2l_b$sigma,2)),0.34)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2l_b$saa),2)),c(0.33,0.51,0.66,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する
  res_vpa_pgc0_tune2m_b <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                Pope = TRUE, p.init = 1, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2m_b$naa),2)),c(3.22,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2m_b$b,2)),c(-0.17,0.17))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2m_b$sigma,2)),rep(0.34,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2m_b$saa),2)),c(0.43,0.68,0.82,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する(2つの指数のsdが異なる場合)
  res_vpa_pgc0_index_change_tune2m_b <- vpa(vpa_pgc0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                                             Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_index_change_tune2m_b$naa),2)),c(6.92,4.18,3.84,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_index_change_tune2m_b$b,2)),c(0.05,-0.07))
  expect_equal(as.numeric(round(res_vpa_pgc0_index_change_tune2m_b$sigma,2)),c(0.25,0.34))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_index_change_tune2m_b$saa),2)),c(0.33,0.52,0.66,1.00))
  
  #選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する．1部のｂは固定
  res_vpa_pgc0_tune2m_b_fix <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                    Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("SSB","SSB"),b.fix=c(NA,-0.1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2m_b_fix$naa),2)),c(5.66,3.54,3.20,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2m_b_fix$b,2)),c(0.09,-0.10))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune2m_b_fix$sigma,2)),rep(0.34,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune2m_b_fix$saa),2)),c(0.34,0.53,0.68,1.00))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最小二乗法．
  res_vpa_pgc0_tune3l <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ls" ,b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3l$naa),2)),c(3.22,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune3l$sigma,2)),c(0.5))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3l$saa),2)),c(0.14,0.22,0.27,0.36))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法
  res_vpa_pgc0_tune3m <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=FALSE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3m$naa),2)),c(3.22,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune3m$sigma,2)),c(0.55,0.46))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3m$saa),2)),c(0.14,0.22,0.27,0.36))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最小二乗法．ｂを推定する
  res_vpa_pgc0_tune3lb <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ls" ,b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3lb$naa),2)),c(7.14,4.28,3.95,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune3lb$b,2)),c(0.07,-0.07))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune3lb$sigma,2)),0.34)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3lb$saa),2)),c(0.12,0.20,0.24,0.44))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する
  res_vpa_pgc0_tune3mb <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3mb$naa),2)),c(3.70,2.83,1.91,1.91))
  expect_equal(as.numeric(round(res_vpa_base0_tune3mb$b,2)),c(-0.31,0.31))
  expect_equal(as.numeric(round(res_vpa_base0_tune3mb$sigma,2)),rep(0.34,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3mb$saa),2)),c(0.14,0.19,0.33,0.33))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する(2つの指数のsdが異なる場合)
  res_vpa_pgc0_index_change_tune3mb <- vpa(vpa_pgc0_index_change, tf.year=2015:2016, last.catch.zero = FALSE, 
                                            Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=TRUE,abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans( res_vpa_pgc0_index_change_tune3mb$naa),2)),c(6.92,4.18,3.84,NA))
  expect_equal(as.numeric(round( res_vpa_pgc0_index_change_tune3mb$b,2)),c(0.05,-0.07))
  expect_equal(as.numeric(round( res_vpa_pgc0_index_change_tune3mb$sigma,2)),c(0.25,0.34))
  expect_equal(as.numeric(round(rowMeans( res_vpa_pgc0_index_change_tune3mb$saa),2)),c(0.12,0.20,0.24,0.44))
  
  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する．一部のbを固定する
  res_vpa_pgc0_tune3mb_fix <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                                   Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=TRUE,abund=c("SSB","SSB"),b.fix=c(NA,-0.1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3mb_fix $naa),2)),c(5.66,3.54,3.20,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune3mb_fix $b,2)),c(0.09,-0.10))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune3mb_fix $sigma,2)),rep(0.34,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune3mb_fix $saa),2)),c(0.13,0.20,0.25,0.42))
  
  #2-3: 全F法によるtuning ----
  
  #  全F推定法．最小二乗法．
  res_vpa_pgc0_tune4l <- vpa(vpadat_pgc0,  last.catch.zero = FALSE, 
                              Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ls" ,b.est=FALSE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4l$naa),2)),c(3.59,2.48,2.14,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4l$sigma,2)),c(0.4))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4l$saa),2)),c(0.19,0.32,0.36,NA))
  
  #  全F推定法．最尤法
  res_vpa_pgc0_tune4m <- vpa(vpadat_pgc0, last.catch.zero = FALSE, 
                              Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=FALSE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4m$naa),2)),c(3.59,2.48,2.14,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4m$sigma,2)),c(0.41,0.40))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4m$saa),2)),c(0.19,0.32,0.36,NA))
  
  #  全F推定法．最小二乗法．b推定あり
  res_vpa_pgc0_tune4l_b <- vpa(vpadat_pgc0, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ls" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4l_b$naa),2)),c(3.41,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4l_b$b,2)),c(-0.17,0.17))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4l_b$sigma,2)),c(0.34))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4l_b$saa),2)),c(0.19,0.32,0.37,NA))
  
  #  全F推定法．最尤法．b推定あり
  res_vpa_pgc0_tune4m_b <- vpa(vpadat_pgc0, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans( res_vpa_pgc0_tune4m_b$naa),2)),c(5.39,4.28,3.95,NA))
  expect_equal(as.numeric(round( res_vpa_pgc0_tune4m_b$b,2)),c(0.07,-0.07))
  expect_equal(as.numeric(round( res_vpa_pgc0_tune4m_b$sigma,2)),rep(0.34,2))
  expect_equal(as.numeric(round(rowMeans( res_vpa_pgc0_tune4m_b$saa),2)),c(0.21,0.31,0.35,NA))
  
  #  全F推定法．最尤法．b推定あり(2つの指数のsdが異なる場合)
  res_vpa_pgc0_index_change_tune4m_b <- vpa(vpa_pgc0_index_change, last.catch.zero = FALSE, 
                                             Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_index_change_tune4m_b$naa),2)),c(3.41,2.30,1.97,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_index_change_tune4m_b$b,2)),c(-0.11,0.17))
  expect_equal(as.numeric(round(res_vpa_pgc0_index_change_tune4m_b$sigma,2)),c(0.25,0.34))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_index_change_tune4m_b$saa),2)),c(0.19,0.32,0.37,NA))
  
  # 　全F推定法．最尤法．b推定あり，alpha=0.3(最高齢と最高齢―1のFの比：Fa=alpha*Fa-1)
  res_vpa_pgc0_tune4m_b_alpha <- vpa(vpadat_pgc0, last.catch.zero = FALSE, 
                                      Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(2,2,2,2),abund=c("SSB","SSB"),alpha=0.3)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4m_b_alpha$naa),2)),c(3.23,2.31,1.99,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4m_b_alpha$b,2)),c(-0.18,0.18))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4m_b_alpha$sigma,2)),rep(0.34,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4m_b_alpha$saa),2)),c(0.26,0.43,0.26,NA))
  
  #2-4: Ridge VPA ----
  
  # set lambda + 全F推定法．最尤法．b推定あり
  res_vpa_pgc0_tune5m_b <- vpa(vpadat_pgc0, last.catch.zero = FALSE, 
                                Pope = TRUE,  tune=TRUE, term.F="all",sel.def="mean", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB"), lambda=0.02)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune5m_b$naa),2)),c(15.52,4.29,3.95,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune5m_b$b,2)),c(0.07,-0.07))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune5m_b$sigma,2)),rep(0.34,2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune5m_b$saa),2)),c(0.18,0.32,0.37,NA))
  
  # TMB true is not set to work for vpadat_pgc0
  
  save(res_vpa_base0_nontune,
       res_vpa_base1_nontune,
       res_vpa_pgc0_nontune,
       res_vpa_rec0_nontune,
       file="res_vpa_files.rda")
})
