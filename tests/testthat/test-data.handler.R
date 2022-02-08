context("generate vpa results with various data type (with dummy data)")

test_that("vpa function (with dummy data) (level 2-3?)",{
  # read various data ----
  # data with caa=maa=waa=1, M=0
  data_base <- readr::read_csv(system.file("extdata","all_dummy_data_base.csv",package="frasyr"))
  # data with caa=maa=waa=1, M=0 but plus group have changed
  data_pgc <- readr::read_csv(system.file("extdata","all_dummy_data_plus_group_change.csv",package="frasyr"))
  # data with caa=maa=waa=1, M=0 but first recruit age is 1
  data_rec <- readr::read_csv(system.file("extdata","all_dummy_data_rec.csv",package="frasyr"))
  # data from "https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw"with some modification
  data_estb <-readr::read_csv(system.file("extdata","all_dummy_data_estb.csv",package="frasyr"))
  # plus group change of the above
  data_pgc_estb<-readr::read_csv(system.file("extdata","all_dummy_data_pgc_estb.csv",package="frasyr"))
  # data with caa=maa=waa=1, M=0 but first recruit age is 2
  data_rec2 <- readr::read_csv(system.file("extdata","all_dummy_data_rec2.csv",package="frasyr"))
  # data with caa=maa=waa=0 for the last year
  data_lst0 <- readr::read_csv(system.file("extdata","all_dummy_data_lst0.csv",package="frasyr"))

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

  vpadat_estb <- data.handler(caa=to_vpa_data(data_estb, label_name="caa"),
                               waa=to_vpa_data(data_estb, label_name="waa"),
                               maa=to_vpa_data(data_estb, label_name="maa"),
                               M  = 0.4,
                               index = to_vpa_data(data_estb, label_name="abund"),
                               maa.tune = NULL,
                               waa.catch = NULL,
                               catch.prop = NULL)

  vpadat_pgc0_estb <- data.handler(caa=to_vpa_data(data_pgc_estb, label_name="caa"),
                              waa=to_vpa_data(data_pgc_estb, label_name="waa"),
                              maa=to_vpa_data(data_pgc_estb, label_name="maa"),
                              M  = 0.4,
                              index = to_vpa_data(data_pgc_estb, label_name="abund"),
                              maa.tune = NULL,
                              waa.catch = NULL,
                              catch.prop = NULL)

  vpadat_rec2 <- data.handler(caa=to_vpa_data(data_rec2, label_name="caa"),
                              waa=to_vpa_data(data_rec2, label_name="waa"),
                              maa=to_vpa_data(data_rec2, label_name="maa"),
                              M  = 0,
                              index = to_vpa_data(data_rec2, label_name="abund"),
                              maa.tune = NULL,
                              waa.catch = NULL,
                              catch.prop = NULL)

  vpadat_lst0 <- data.handler(caa=to_vpa_data(data_lst0, label_name="caa"),
                              waa=to_vpa_data(data_lst0, label_name="waa"),
                              maa=to_vpa_data(data_lst0, label_name="maa"),
                              M  = 0.4,
                              index = to_vpa_data(data_lst0, label_name="abund"),
                              maa.tune = NULL,
                              waa.catch = NULL,
                              catch.prop = NULL)

  # 放流データを含むVPAデータ
  data_release <- to_vpa_data(data_base, label_name="caa")[1,]
  release.number <- 2
  data_release[] <- release.number

  vpadat_base0_release <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
                               waa=to_vpa_data(data_base, label_name="waa"),
                               maa=to_vpa_data(data_base, label_name="maa"),
                               M  = 0,
                               index = to_vpa_data(data_base, label_name="abund"),
                               maa.tune = NULL,
                               waa.catch = NULL,
                               release.dat = data_release,
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

  res_vpa_base0_nontune_release <- vpa(vpadat_base0_release, tf.year=2015:2016, last.catch.zero = FALSE,
                                       Pope = TRUE, p.init = 0.5)
  expect_equal(as.numeric(rowMeans(res_vpa_base0_nontune_release$naa)),
               true_number)

  res_vpa_base1_nontune <- vpa(vpadat_base1, tf.year=2015:2016, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5)
  expect_equal(as.numeric(rowMeans(res_vpa_base1_nontune$naa)),
               true_number)

  #last.catch.zero=TRUEのときでも，最終年（Ｃａｔｃｈが0の年）のSSBが他の値があれば計算されることの確認
  res_vpa_lst0 <- vpa(vpadat_lst0, tf.year=2007:2008, last.catch.zero = TRUE,
                               Pope = TRUE, p.init = 0.5, plus.group=TRUE)
  expect_equal(as.numeric(round(rowMeans(res_vpa_lst0$ssb),2)),c(NA,16.40,37.92,30.63))

  # プラスグループが変わる場合はtrue_numberには一致しない；どうテストすべきか？
  res_vpa_pgc0_nontune <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5)
  expect_equal(as.numeric(unlist(res_vpa_pgc0_nontune$naa["2017"])),
               c(3,2,2,NA), tol=0.0001)

  res_vpa_rec0_nontune <- vpa(vpadat_rec0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5)
  expect_equal(as.numeric(rowMeans(res_vpa_rec0_nontune$naa)),
               true_number)

  res_vpa_rec2_nontune <- vpa(vpadat_rec2, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5)

  data_SR <- get.SRdata(res_vpa_rec2_nontune)
  expect_equal(nrow(data_SR), ncol(res_vpa_rec2_nontune$naa)-2)

  # catch計算用のwaaを２倍にしているbase1データでは漁獲量が倍になる
  expect_equal(res_vpa_base0_nontune$wcaa*2,
               res_vpa_base1_nontune$wcaa)

  #b推定するデータを用いたとき
  res_vpa_estb_nontune <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5,fc.year=1998:2000)

  res_vpa_pgc_estb_nontune <- vpa(vpadat_pgc0_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5,fc.year=1998:2000)


  # Part1: dataset is "vpadat_base0" for b.est=FALSE, "vpadat_estb" for b.est=TRUE----
  # 1-1: 二段階法によるtuning ----

  sel.f1<-res_vpa_base0_nontune$saa$"2017"
  round(sel.f1,3)

  #二段階法：est.method=最小二乗法による推定
  res_vpa_base0_tune1l <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ls", b.est=FALSE,abund=c("B","B"))
  expect_equal(as.numeric(rowMeans(res_vpa_base0_tune1l$naa)),true_number,tol=0.0001)
  expect_equal(as.numeric(res_vpa_base0_tune1l$sigma),true_sd, tol=0.0001)

  #二段階法：est.method=最尤法による推定
  res_vpa_base0_tune1m <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ml", b.est=FALSE,abund=c("B","B"))
  expect_equal(as.numeric(rowMeans(res_vpa_base0_tune1m$naa)),true_number, tol=0.0001)
  expect_equal(as.numeric(res_vpa_base0_tune1m$sigma),rep(true_sd,2), tol=0.0001)

  #二段階法：est.method=最尤法による推定(2つの指数のsdが異なる場合）
  vpadat_base0_index_change <- vpadat_base0
  vpadat_base0_index_change$index[1,is.na(vpadat_base0_index_change$index[1,])] <- exp(mean(log(c(1,2))))
  res_vpa_base0_index_change_tune1m <- vpa(vpadat_base0_index_change, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1,est.method="ml", b.est=FALSE,abund=c("B","B"))
  expect_equal(round(as.numeric(rowMeans(res_vpa_base0_index_change_tune1m$naa)),2),
               true_number+c(0.02,0.01,0.01,0.01),
               tol=0.01)
  expect_equal(round(as.numeric(res_vpa_base0_index_change_tune1m$sigma),2),
               c(0.25,0.35))

  #二段階法：est.method=最小二乗法による推定＋指標値の非線形性bの推定
  sel.f2<-res_vpa_estb_nontune$saa$"2000"
  round(sel.f2,3)

  res_vpa_estb_tune1l_b <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f2,est.method="ls", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune1l_b$naa),2)),c(709.95,346.35,146.32,63.81))
  expect_equal(as.numeric(round(res_vpa_estb_tune1l_b$b,2)),c(1.01,0.53,0.73,0.51,0.67,1.16))
  expect_equal(as.numeric(round(res_vpa_estb_tune1l_b$sigma,2)),0.18)

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定
  res_vpa_estb_tune1m_b <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f2,est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune1m_b$naa),2)),c(715.80,348.57,147.19,64.23))
  expect_equal(as.numeric(round(res_vpa_estb_tune1m_b$b,2)),c(1.06,0.57,0.76,0.54,0.70,1.21))
  expect_equal(as.numeric(round(res_vpa_estb_tune1m_b$sigma,2)),c(0.18,0.16,0.11,0.13,0.17,0.28))

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定+指標２と３のシグマは同じ
  res_vpa_estb_tune1m_b_sigma <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f2,est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_estb_tune1m_b_sigma$sigma[[2]],res_vpa_estb_tune1m_b_sigma$sigma[[3]])


  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定＋一部のbは固定
  res_vpa_estb_tune1m_b_fix <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f2,est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),b.fix=c(NA,0.7,NA,NA,NA,1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune1m_b_fix$naa),2)),c(718.23,349.50,147.55,64.41))
  expect_equal(as.numeric(round(res_vpa_estb_tune1m_b_fix$b,2)),c(1.09,0.70,0.77,0.55,0.72,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune1m_b_fix$sigma,2)),c(0.18,0.16,0.11,0.13,0.17,0.29))

  #二段階法: last.catch.zero=TRUEで，rec.new=NULLのとき，最新年の0歳はNAとなる
  res_vpa_base0_lst0_ni1 <- vpa(vpadat_lst0, tf.year=2007:2008, last.catch.zero = TRUE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f2, est.method="ml", b.est=FALSE, abund=c("N","N","N","N","N","N"),min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),rec.new=NULL)
  expect_equal(as.numeric(res_vpa_base0_lst0_ni1$naa[1,10]),NA_real_)

  #二段階法: last.catch.zero=TRUEで，rec.new=1000のとき，最新年の0歳は1000となる
  res_vpa_base0_lst0_ni2 <- vpa(vpadat_lst0, tf.year=2007:2009, last.catch.zero = TRUE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f2, est.method="ml", b.est=FALSE, abund=c("N","N","N","N","N","N"),min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),rec.new=1000)
  expect_equal(as.numeric(res_vpa_base0_lst0_ni2$naa[1,10]),1000)


  #1-2: 選択率更新法によるtuning ----

  #選択率更新法：選択率の計算方法は，一番高いFを１とする．最小二乗法
  res_vpa_base0_tune2l <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ls", b.est=FALSE,abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune2l$sigma,2)),true_sd,tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l$saa),2)),c(0.42,0.58,1.00,1.00))

  #選択率更新法：選択率の計算方法は，最高齢を１とする．最小二乗法
  res_vpa_base0_tune2l_mxa <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="maxage",sel.update=TRUE, est.method="ls", b.est=FALSE,abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l_mxa $naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune2l_mxa $sigma,2)),true_sd,tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2l_mxa $saa),2)),c(0.42,0.58,1.00,1.00))


  #選択率更新法：選択率の計算方法は，一番高いFを１とする．最尤法
  res_vpa_base0_tune2m <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=FALSE,abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune2m$sigma,2)),rep(true_sd,2), tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune2m$saa),2)),c(0.42,0.58,1.00,1.00))

  #選択率更新法：選択率の計算方法は，一番高いFを１とする．最小二乗法．b推定する
  res_vpa_estb_tune2l_b <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ls", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b$naa),2)),c(705.34,343.43,145.06,63.21))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b$b,2)),c(0.97,0.51,0.71,0.50,0.64,1.11))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b$sigma,2)),0.18)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b$saa),2)),c(0.46,0.65,1.00,1.00))

  #選択率更新法：選択率の計算方法は，最高齢のFを１とする．最小二乗法．b推定する
  res_vpa_estb_tune2l_b_mxa <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="maxage",sel.update=TRUE, est.method="ls", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b_mxa$naa),2)),c(705.34,343.43,145.06,63.21))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b_mxa$b,2)),c(0.97,0.51,0.71,0.50,0.64,1.11))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b_mxa$sigma,2)),0.18)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b_mxa$saa),2)),c(0.46,0.65,1.00,1.00))

  #選択率更新法：選択率の計算方法は，最高齢のFを１とする．最小二乗法．b推定する.aveS=FALSE
  res_vpa_estb_tune2l_b_mxa_aveS <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                   Pope = TRUE,  p.init = 0.5, tune=TRUE, term.F="max",sel.def="maxage",sel.update=TRUE, est.method="ls", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),ave_S=FALSE)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b_mxa_aveS$naa),2)),c(705.26,343.41,145.03,63.19))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b_mxa_aveS$b,2)),c(0.97,0.51,0.71,0.50,0.64,1.11))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b_mxa_aveS$sigma,2)),0.18)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b_mxa_aveS$saa),2)),c(0.46,0.65,1.00,1.00))


  #選択率更新法：選択率の計算方法は，一番高いFを１とする．最尤法．b推定する
  res_vpa_estb_tune2m_b <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b$naa),2)),c(711.09,345.43,145.83,63.58))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b$b,2)),c(1.02,0.54,0.73,0.52,0.67,1.16))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b$sigma,2)),c(0.18,0.16,0.11,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b$saa),2)),c(0.46,0.65,1.00,1.00))

  #選択率更新法：選択率の計算方法は，最高齢のFを１とする．最尤法．b推定する
  res_vpa_estb_tune2m_b_mxa <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="maxage",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b$naa),2)),c(711.09,345.43,145.83,63.58))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b$b,2)),c(1.02,0.54,0.73,0.52,0.67,1.16))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b$sigma,2)),c(0.18,0.16,0.11,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b$saa),2)),c(0.46,0.65,1.00,1.00))


  #選択率更新法：選択率の計算方法は，一番高いFを１とする．最尤法．b推定する,alphaを与える
  res_vpa_estb_tune2m_b_alpha <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),alpha=0.3)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_alpha $naa),2)),c(935.10,501.99,257.51,291.44))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_alpha $b,2)),c(1.54,0.77,0.82,0.54,1.02,1.84))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_alpha $sigma,2)),c(0.17,0.17,0.10,0.15,0.16,0.24))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_alpha $saa),2)),c(0.72,0.86,0.98,0.30))
  expect_equal(as.numeric(res_vpa_estb_tune2m_b_alpha $faa[3,]*0.3),as.numeric(res_vpa_estb_tune2m_b_alpha $faa[4,]))

  #選択率更新法：選択率の計算方法は，一番高いFを１とする．最尤法．b推定する.指標の2番目と3番目のシグマは同じ
  res_vpa_estb_tune2m_b_sigma <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_estb_tune2m_b_sigma$sigma[[2]],res_vpa_estb_tune2m_b_sigma$sigma[[3]])


  #選択率更新法：選択率の計算方法は，一番高いFを１とする．最尤法．b推定する．1部のｂは固定
  res_vpa_estb_tune2m_b_fix <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),b.fix=c(NA,0.7,NA,NA,NA,1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_fix$naa),2)),c(714.73,346.69,146.31,63.81))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_fix$b,2)),c(1.05,0.70,0.75,0.54,0.69,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_fix$sigma,2)),c(0.18,0.17,0.11,0.13,0.17,0.29))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_fix$saa),2)),c(0.46,0.65,1.00,1.00))

  #選択率更新法：平均値Fで割る（sel.def="mean")．最小二乗法．
  res_vpa_base0_tune3l <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 1, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ls" ,b.est=FALSE,abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3l$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune3l$sigma,2)),true_sd,tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3l$saa),2)),c(0.14,0.19,0.33,0.33))

  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法
  res_vpa_base0_tune3m <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 1, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=FALSE,abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3m$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune3m$sigma,2)),rep(true_sd,2), tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune3m$saa),2)),c(0.14,0.19,0.33,0.33))

  #選択率更新法：平均値Fで割る（sel.def="mean")．最小二乗法．ｂを推定する
  res_vpa_estb_tune3lb <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 1, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ls" ,b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune3lb$naa),2)),c(705.26,343.39,145.04,63.20))
  expect_equal(as.numeric(round(res_vpa_estb_tune3lb$b,2)),c(0.97,0.51,0.71,0.50,0.64,1.11))
  expect_equal(as.numeric(round(res_vpa_estb_tune3lb$sigma,2)),0.18)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune3lb$saa),2)),c(0.15,0.21,0.32,0.32))

  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する
  res_vpa_estb_tune3mb <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                              Pope = TRUE, p.init =1, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml" ,b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune3mb$naa),2)),c(711.01,345.39,145.80,63.56))
  expect_equal(as.numeric(round(res_vpa_estb_tune3mb$b,2)),c(1.02,0.54,0.73,0.52,0.67,1.16))
  expect_equal(as.numeric(round(res_vpa_estb_tune3mb$sigma,2)),c(0.18,0.16,0.11,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune3mb$saa),2)),c(0.15,0.21,0.32,0.32))

  #選択率更新法：平均値Fで割る（sel.def="mean")．最尤法．ｂを推定する．一部のbを固定する
  res_vpa_estb_tune3mb_fix <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 1, tune=TRUE, term.F="max",sel.def="mean",sel.update=TRUE, est.method="ml",b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),b.fix=c(NA,0.7,NA,NA,NA,1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune3mb_fix $naa),2)),c(714.67,346.66,146.28,63.80))
  expect_equal(as.numeric(round(res_vpa_estb_tune3mb_fix $b,2)),c(1.05,0.70,0.75,0.54,0.69,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune3mb_fix $sigma,2)),c(0.18,0.17,0.11,0.13,0.17,0.29))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune3mb_fix $saa),2)),c(0.15,0.21,0.32,0.32))

  #選択率更新法：Baranov. 選択率の計算方法は，一番高いFを１とする．最小二乗法．b推定する
  res_vpa_estb_tune2l_b_baranov <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = FALSE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ls", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b_baranov$naa),2)),c(687.09,333.39,140.30,61.13))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b_baranov$b,2)),c(0.98,0.52,0.71,0.50,0.65,1.12))
  expect_equal(as.numeric(round(res_vpa_estb_tune2l_b_baranov$sigma,2)),0.18)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2l_b_baranov$saa),2)),c(0.46,0.65,1.00,1.00))

  #選択率更新法: last.catch.zero=TRUEで，rec.new=NULLのとき，最新年の0歳はNAとなる
  res_vpa_base0_lst0 <- vpa(vpadat_lst0, tf.year=2007:2008, last.catch.zero = TRUE,
                            Pope = TRUE, p.init = 1, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=FALSE, abund=c("N","N","N","N","N","N"),min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),rec.new=NULL)
  expect_equal(as.numeric(res_vpa_base0_lst0$naa[1,10]),NA_real_)

  #選択率更新法: last.catch.zero=TRUEで，rec.new=1000のとき，最新年の0歳は1000となる
  res_vpa_base0_lst0_2 <- vpa(vpadat_lst0, tf.year=2007:2009, last.catch.zero = TRUE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=FALSE, abund=c("N","N","N","N","N","N"),min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),rec.new=1000)
  expect_equal(as.numeric(res_vpa_base0_lst0_2$naa[1,10]),1000)


  #1-3: 全F法によるtuning ----

  #  全F推定法．最小二乗法．
  res_vpa_base0_tune4l <-vpa(vpadat_base0, last.catch.zero = FALSE,
                            Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=FALSE,p.init=c(0.2,0.3,0.6,0.6),abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4l$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune4l$sigma,2)),true_sd,tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4l$saa),2)),c(0.42,0.58,1.00,1.00))

  #  全F推定法．最尤法
  res_vpa_base0_tune4m <- vpa(vpadat_base0, last.catch.zero = FALSE,
                              Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=FALSE,p.init=c(0.2,0.3,0.6,0.6),abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_tune4m$sigma,2)),rep(true_sd,2),tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_tune4m$saa),2)),c(0.42,0.58,1.00,1.00))

  #  全F推定法．最小二乗法．b推定あり
  res_vpa_estb_tune4l_b <- vpa(vpadat_estb,last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4l_b$naa),2)),c(727.24,355.51,153.09,67.08))
  expect_equal(as.numeric(round(res_vpa_estb_tune4l_b$b,2)),c(1.16,0.62,0.75,0.53,0.75,1.35))
  expect_equal(as.numeric(round(res_vpa_estb_tune4l_b$sigma,2)),c(0.18))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4l_b$saa),2)),c(0.49,0.68,0.99,0.99))

  #  全F推定法．最尤法．b推定あり
  res_vpa_estb_tune4m_b <- vpa(vpadat_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4m_b$naa),2)),c(694.83,335.61,142.37,61.91))
  expect_equal(as.numeric(round(res_vpa_estb_tune4m_b$b,2)),c(0.88,0.45,0.64,0.45,0.57,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune4m_b$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.29))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4m_b$saa),2)),c(0.46,0.66,1.00,1.00))

  #  全F推定法．最尤法．b推定あり,指標2と3のシグマは同じとする
  res_vpa_estb_tune4m_b_sigma <- vpa(vpadat_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000, sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_estb_tune4m_b_sigma$sigma[[2]],res_vpa_estb_tune4m_b_sigma$sigma[[3]])

  # 　全F推定法．最小二乗法．b推定あり， Pope=FALSE
  res_vpa_estb_tune4l_b_baranov <- vpa(vpadat_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = FALSE,  tune=TRUE, term.F="all",est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4l_b_baranov$naa),2)),c(710.08,346.11,148.54,65.10))
  expect_equal(as.numeric(round(res_vpa_estb_tune4l_b_baranov$b,2)),c(1.19,0.64,0.76,0.54,0.77,1.38))
  expect_equal(as.numeric(round(res_vpa_estb_tune4l_b_baranov$sigma,2)),0.18)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4l_b_baranov$saa),2)),c(0.49,0.68,1.00,1.00))

  # 　全F推定法．最尤法．b推定あり， Pope=FALSE, 指標２と３のシグマは同じ
  res_vpa_estb_tune4m_b_baranov_sigma <- vpa(vpadat_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                       Pope = FALSE,  tune=TRUE, term.F="all",est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000,sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_estb_tune4m_b_baranov_sigma$sigma[[2]],res_vpa_estb_tune4m_b_baranov_sigma$sigma[[3]])

  # 　全F推定法．最尤法．b推定あり，alpha=0.8(最高齢と最高齢―1のFの比：Fa=alpha*Fa-1)
  res_vpa_estb_tune4m_b_alpha <- vpa(vpadat_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                        Pope = TRUE,  tune=TRUE, term.F="all",est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),alpha=0.8,fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4m_b_alpha$naa),2)),c(712.60,346.43,149.16,72.54))
  expect_equal(as.numeric(round(res_vpa_estb_tune4m_b_alpha$b,2)),c(0.98,0.51,0.69,0.48,0.63,1.11))
  expect_equal(as.numeric(round(res_vpa_estb_tune4m_b_alpha$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune4m_b_alpha$saa),2)),c(0.48,0.68,1.00,0.80))

  #   全F推定法: last.catch.zero=TRUEで，rec.new=NULLのとき，最新年の0歳はNAとなる
  res_vpa_base0_lst0_zen1 <- vpa(vpadat_lst0, tf.year=2007:2009, last.catch.zero = TRUE,
                                 Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="all", est.method="ml", b.est=FALSE, abund=c("N","N","N","N","N","N"),min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),rec.new=NULL)
  expect_equal(as.numeric(res_vpa_base0_lst0_zen1$naa[1,10]),NA_real_)
  expect_equal(as.numeric(round(res_vpa_base0_lst0_zen1$naa[,10],2)),c(NA_real_,538.18,25.87,150.63))

  #   全F推定法: last.catch.zero=TRUEで，rec.new=NULLのとき，plus.group=FALSE
  res_vpa_base0_lst0_zen11 <- vpa(vpadat_lst0, tf.year=2007:2009, last.catch.zero = TRUE,
                                 Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="all", est.method="ml", b.est=FALSE, abund=c("N","N","N","N","N","N"),min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),rec.new=NULL,plus.group=FALSE)

  expect_equal(as.numeric(round(res_vpa_base0_lst0_zen11$naa[,10],2)),c(NA_real_,685.22,86.13,118.49))


  #   全F推定法: last.catch.zero=TRUEで，rec.new=1000のとき，最新年の0歳は1000となる
  res_vpa_base0_lst0_zen2 <- vpa(vpadat_lst0, tf.year=2007:2009, last.catch.zero = TRUE,
                                 Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="all", est.method="ml", b.est=FALSE, abund=c("N","N","N","N","N","N"),min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),rec.new=1000)
  expect_equal(as.numeric(res_vpa_base0_lst0_zen2$naa[1,10]),1000)


  #1-4: Ridge VPA ----

  # set lambda + 全F推定法．最尤法．b推定あり,penalty="p", eta=NULL
  res_vpa_estb_tune5m_b <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="p")
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_b$naa),2)),c(695.60,336.02,142.58,62.01))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_b$b,2)),c(0.88,0.46,0.65,0.45,0.57,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_b$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_b$saa),2)),c(0.46,0.66,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_b$logLik,3)),c(23.116))
  
  # set lambda + 全F推定法．最尤法．b推定あり,penalty="f", eta=NULL(add:2021/11)
  res_vpa_estb_tune5m_bf <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="f",tf.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_bf$naa),2)),c(694.79,335.59,142.36,61.91))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_bf$b,2)),c(0.88,0.45,0.64,0.45,0.57,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_bf$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.29))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_bf$saa),2)),c(0.46,0.66,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_bf$logLik,3)),c(23.127))
  
  # set lambda + 全F推定法．最尤法．b推定あり,penalty="s", eta=NULL(add:2021/11)
  res_vpa_estb_tune5m_bs <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="s",tf.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_bs$naa),2)),c(694.86,335.62,142.38,61.91))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_bs$b,2)),c(0.88,0.45,0.64,0.45,0.57,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_bs$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.29))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_bs$saa),2)),c(0.46,0.66,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_bs$logLik,3)),c(23.128))
  

  # set lambda + 全F推定法．最尤法．b推定あり,penalty="p", eta=0.3, eta.age=0
  res_vpa_estb_tune5m_b_e <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_b_e$naa),2)),c(695.37,335.90,142.52,61.98))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_b_e$b,2)),c(0.88,0.45,0.65,0.45,0.57,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_b_e$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune5m_b_e$saa),2)),c(0.46,0.66,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune5m_b_e$logLik,3)),c(23.12))

  #set lambda + 選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する (eta=NULL) だけど　p_by_age=FALSEのまま（つまりターミナルＦのものだけペナルテイーがかかる）
  res_vpa_estb_tune2m_b_r <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),lambda=0.02)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_r$naa),2)),c(711.10,345.44,145.83,63.58))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_r$b,2)),c(1.02,0.54,0.73,0.52,0.67,1.16))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_r$sigma,2)),c(0.18,0.16,0.11,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_r$saa),2)),c(0.46,0.65,1.00,1.00))
  expect_equal(as.numeric(round( res_vpa_estb_tune2m_b_r$logLik,3)),c(22.698))

  #set lambda + 選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する (eta=NULL) p_by_age=TRUE, penalty_age=3なら一個前の方法と結果は同じになることの確かめ
  res_vpa_estb_tune2m_b_r_p <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                 Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),lambda=0.02, p_by_age=TRUE, penalty_age=3)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_r$naa),2)),c(711.10,345.44,145.83,63.58))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_r$b,2)),c(1.02,0.54,0.73,0.52,0.67,1.16))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_r$sigma,2)),c(0.18,0.16,0.11,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune2m_b_r$saa),2)),c(0.46,0.65,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune2m_b_r_p $logLik,3)),c(22.698))

  #set lambda + 選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する (eta=NULL) p_by_age=TRUE, penalty_age=0:2なら尤度は変化することの確かめ
  res_vpa_estb_tune2m_b_r_p2 <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                   Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),lambda=0.02, p_by_age=TRUE, penalty_age=0:2)
  expect_equal(as.numeric(round( res_vpa_estb_tune2m_b_r_p2$logLik,3)),c(22.695))


  #set lambda + 選択率更新法：選択率の計算方法は，最高齢を１とする．最尤法．b推定する (eta=0) p_by_age=TRUE, eta.age=0, no_eta_age=3
  res_vpa_estb_tune2m_b_r_p_e <- vpa(vpadat_estb, tf.year=1997:1999, last.catch.zero = FALSE,
                                   Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=TRUE, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),lambda=0.02, p_by_age=TRUE, eta=0.1, eta.age=0, no_eta_age=1:3)

  expect_equal(as.numeric(round(  res_vpa_estb_tune2m_b_r_p_e$logLik,3)),c(22.693))


  # set lambda + 全F推定法．最尤法．b推定あり+指標2と３のシグマは同じ
  res_vpa_estb_tune5m_b_sigma <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_estb_tune5m_b_sigma$sigma[[2]],res_vpa_estb_tune5m_b_sigma$sigma[[3]])

  # TMB true + set lambda + 全F推定法．最尤法．b推定あり----
  library(TMB)
  use_rvpa_tmb()
  res_vpa_estb_tune6m_b <- vpa(vpadat_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = TRUE,  tune=TRUE, term.F="all",est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, TMB=TRUE,fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b$naa),2)),c(695.60,336.02,142.58,62.01))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b$b,2)),c(0.88,0.46,0.65,0.45,0.57,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b$saa),2)),c(0.46,0.66,1.00,1.00))

  # set lambda + 全F推定法．最尤法．b推定あり,penalty="p", eta=0.3, eta.age=0, TMB=TRUE
  res_vpa_estb_tune6m_b_e <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                 Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=TRUE)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_e$naa),2)),c(695.37,335.90,142.52,61.98))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_e$b,2)),c(0.88,0.45,0.65,0.45,0.57,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_e$sigma,2)),c(0.18,0.16,0.10,0.13,0.17,0.28))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_e$saa),2)),c(0.46,0.66,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_e$logLik,3)),c(23.12))

  # TMB true + set lambda + 全F推定法．最尤法．b推定あり,指標２と３のシグマは同じ
  res_vpa_estb_tune6m_b_sigma <- vpa(vpadat_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all",est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, TMB=TRUE,fc.year=1998:2000,sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_estb_tune6m_b_sigma$sigma[[2]],res_vpa_estb_tune6m_b_sigma$sigma[[3]])

  # 最小二乗法，set lambda + 全F推定法．b推定あり,penalty="p", eta=0.3, eta.age=0, TMB=TRUE (追加:2021/09/16)
  res_vpa_estb_tune6m_b_maiwashi1 <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                 Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=TRUE)
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi1$naa),2)),c(728.36,356.15,153.39,67.22))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi1$b,2)),c(1.17,0.63,0.75,0.54,0.76,1.36))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi1$sigma,2)),c(0.18))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi1$saa),2)),c(0.49,0.68,0.99,0.99))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi1$logLik,3)),c(18.41))

  # b.fix, 最小二乗法，set lambda + 全F推定法,penalty="p", eta=0.3, eta.age=0, TMB=TRUE (追加:2021/09/16)
  res_vpa_estb_tune6m_b_maiwashi2 <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                         Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=TRUE, b.fix=c(1,NA,1,NA,1,1))
  expect_equal(as.numeric(round(rowMeans( res_vpa_estb_tune6m_b_maiwashi2$naa),2)),c(723.72,352.64,149.99,65.58))
  expect_equal(as.numeric(round( res_vpa_estb_tune6m_b_maiwashi2$b,2)),c(1.00,0.61,1.00,0.56,1.00,1.00))
  expect_equal(as.numeric(round( res_vpa_estb_tune6m_b_maiwashi2$sigma,2)),c(0.19))
  expect_equal(as.numeric(round(rowMeans( res_vpa_estb_tune6m_b_maiwashi2$saa),2)),c(0.47,0.66,1.00,1.00))
  expect_equal(as.numeric(round( res_vpa_estb_tune6m_b_maiwashi2$logLik,3)),c(15.336))

  # abund=B, b.fix, 最小二乗法，set lambda + 全F推定法,penalty="p", eta=0.3, eta.age=0, TMB=TRUE (追加:2021/09/16)
  res_vpa_estb_tune6m_b_maiwashi3 <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                         Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("B","B","B","B","B","B"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=TRUE, b.fix=c(1,NA,1,NA,1,1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi3$naa),2)),c(718.52,346.02,145.28,63.31))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi3$b,2)),c(1.00,0.51,1.00,0.56,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi3$sigma,2)),c(0.22))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi3$saa),2)),c(0.45,0.64,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi3$logLik,3)),c(6.252))

  res_vpa_estb_tune6m_b_maiwashi3notmb <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                         Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("B","B","B","B","B","B"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=FALSE, b.fix=c(1,NA,1,NA,1,1))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi3notmb$naa),2)),c(718.52,346.02,145.28,63.31))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi3notmb$b,2)),c(1.00,0.51,1.00,0.56,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi3notmb$sigma,2)),c(0.22))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi3notmb$saa),2)),c(0.45,0.64,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi3notmb$logLik,3)),c(6.252))

  # abund=SSB, b.fix, 最小二乗法，set lambda + 全F推定法,penalty="p", eta=0.3, eta.age=0, TMB=TRUE (追加:2021/09/16)
  res_vpa_estb_tune6m_b_maiwashi4 <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                         Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB","N","N","SSB","SSB"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=TRUE, b.fix=c(1,NA,1,NA,1,1),
                                         sdreport=FALSE)

  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi4$naa),2)),c(718.62,345.71,144.19,62.79))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi4$b,2)),c(1.00,0.43,1.00,0.57,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi4$sigma,2)),c(0.24))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi4$saa),2)),c(0.44,0.63,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi4$logLik,3)),c(1.050))

  #上記ケースでtmb=FALSEのもの
  res_vpa_estb_tune6m_b_maiwashi4notmb <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                         Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE,
                                         p.init=c(0.2,0.3,0.6,0.6),abund=c("SSB","SSB","N","N","SSB","SSB"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=FALSE, b.fix=c(1,NA,1,NA,1,1),no.est=FALSE)

  # tmb = TRUE/FALSEで結果が一致するか
  expect_equal(as.numeric(res_vpa_estb_tune6m_b_maiwashi4$term.f),as.numeric(res_vpa_estb_tune6m_b_maiwashi4notmb$term.f))
  expect_equal(as.numeric(res_vpa_estb_tune6m_b_maiwashi4$logLik),as.numeric(res_vpa_estb_tune6m_b_maiwashi4notmb$logLik))

  # abund=Bs, b.fix, 最小二乗法，set lambda + 全F推定法,penalty="p", eta=0.3, eta.age=0, TMB=TRUE (追加:2021/09/16)
  res_vpa_estb_tune6m_b_maiwashi5 <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                         Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6),abund=c("Bs","Bs","N","N","Bs","Bs"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=TRUE, b.fix=c(1,NA,1,NA,1,1),
                                         sdreport=FALSE)

  expect_equal(as.numeric(round(rowMeans( res_vpa_estb_tune6m_b_maiwashi5$naa),2)),c(721.88, 347.66, 145.54, 63.44))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi5$b,2)),c(1.00,0.49,1.00,0.58,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi5$sigma,2)),c(0.23))
  expect_equal(as.numeric(round(rowMeans(res_vpa_estb_tune6m_b_maiwashi5$saa),2)),c(0.45,0.64,1.00,1.00))
  expect_equal(as.numeric(round(res_vpa_estb_tune6m_b_maiwashi5$logLik,3)),c(3.603))

  #上記ケースでtmb=FALSEのもの
  res_vpa_estb_tune6m_b_maiwashi5notmb <- vpa(vpadat_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                         Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE, p.init=c(0.2,0.3,0.6),abund=c("Bs","Bs","N","N","Bs","Bs"), lambda=0.02, fc.year=1998:2000,penalty="p", eta=0.3, eta.age=0, TMB=FALSE, b.fix=c(1,NA,1,NA,1,1))

  # tmb = TRUE/FALSEで結果が一致するか
  expect_equal(as.numeric(res_vpa_estb_tune6m_b_maiwashi5$term.f),as.numeric(res_vpa_estb_tune6m_b_maiwashi5notmb$term.f))
  expect_equal(as.numeric(res_vpa_estb_tune6m_b_maiwashi5$logLik),as.numeric(res_vpa_estb_tune6m_b_maiwashi5notmb$logLik))

  #1-5: test abund.extractor function----
  naa_base0<-res_vpa_base0_nontune$naa
  faa_base0<-res_vpa_base0_nontune$faa
  omega_base0<-matrix(0.5,nrow=4,ncol=23)
  colnames(omega_base0)<-c(1995:2017)
  rownames(omega_base0)<-c(0:3)
  naa_base2<-res_vpa_rec2_nontune$naa
  faa_base2<-res_vpa_rec2_nontune$faa

  true_abun_N<-rep(11,length(naa_base0))
  true_abun_Nm<-rep(8.742,length(naa_base0))
  true_abun_B<-rep(0.011,length(naa_base0))
  true_abun_Bm<-rep(0.009,length(naa_base0))
  true_abun_SSB<-rep(0.007,length(naa_base0))
  true_abun_Bs<-rep(0.0074,length(naa_base0))
  true_abun_Bo<-rep(0.0025,length(naa_base0))
  true_abun_Ns<-rep(7.415,length(naa_base0))
  true_abun_SSBm<-rep(0.0053,length(naa_base0))
  true_abun_SSBmsj<-rep(0.0070,length(naa_base0))
  true_abun_N1sj<-c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
  true_abun_N0sj<-c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)

  # when abund="N" :abundance
  abund_vpadat_base0_N <- abund.extractor(dat=vpadat_base0,abund="N", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(abund_vpadat_base0_N),true_abun_N)

  # when abund="Nm":abundance at the middle of the year
  abund_vpadat_base0_Nm <- abund.extractor(dat=vpadat_base0,abund="Nm", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_Nm,3)),true_abun_Nm,tol=0.001)

  # when abund="B":biomass
  abund_vpadat_base0_B <- abund.extractor(dat=vpadat_base0,abund="B", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_B,3)),true_abun_B,tol=0.001)

  # when abund="Bm":biomass at the middle of the year
  abund_vpadat_base0_Bm <- abund.extractor(dat=vpadat_base0,abund="Bm", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_Bm,3)),true_abun_Bm,tol=0.001)

  # when abund="SSB"
  abund_vpadat_base0_SSB <- abund.extractor(dat=vpadat_base0,abund="SSB", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_SSB,3)),true_abun_SSB,tol=0.001)

  # when abund="Bs":biomass multiplied by normal selectivity
  abund_vpadat_base0_Bs <- abund.extractor(dat=vpadat_base0,abund="Bs", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_Bs,4)),true_abun_Bs,tol=0.0001)

  # when abund="Bo":adjust selectivity with par omega (=ratio of catch between different fishery) since the tuning index is only from some certain fishery. Use that omega tuned selectivity and multiply with biomass
  abund_vpadat_base0_Bo <- abund.extractor(dat=vpadat_base0,abund="Bo", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3, omega=omega_base0)
  expect_equal(as.numeric(round(abund_vpadat_base0_Bo,4)),true_abun_Bo,tol=0.0001)

  # when abund="Ns":number multipled by normal selectivity
  abund_vpadat_base0_Ns<- abund.extractor(dat=vpadat_base0,abund="Ns", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_Ns,4)),true_abun_Ns,tol=0.0001)
  
  # when abund="SSBm":SSB at the middle of the year
  abund_vpadat_base0_SSBm<- abund.extractor(dat=vpadat_base0,abund="SSBm", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_SSBm,4)),true_abun_SSBm,tol=0.0001)
  
  # when abund="SSBmsj"
  abund_vpadat_base0_SSBmsj<- abund.extractor(dat=vpadat_base0,abund="SSBmsj", naa=naa_base0, faa=faa_base0, min.age=0, max.age=3)
  expect_equal(as.numeric(round(abund_vpadat_base0_SSBmsj,4)),true_abun_SSBmsj,tol=0.0001)

  # when abund="N1sj": This is a special case for Sukesoudara-Japan Sea stock. Enter fishery from 2years old, but the abundance index is for 0 and 1 year olds. So here estimates the number of age 1
  abund_vpadat_rec2_N1sj<- abund.extractor(dat=vpadat_rec2,abund="N1sj", naa=naa_base2, faa=faa_base2, min.age=2, max.age=4)
  expect_equal(as.numeric(abund_vpadat_rec2_N1sj),true_abun_N1sj)

  # when abund="N0sj": This is a special case for Sukesoudara-Japan Sea stock. Enter fishery from 2years old, but the abundance index is for 0 and 1 year olds. So here estimates the number of age 0
  abund_vpadat_rec2_N0sj<- abund.extractor(dat=vpadat_rec2,abund="N0sj", naa=naa_base2, faa=faa_base2, min.age=2, max.age=4)
  expect_equal(as.numeric(abund_vpadat_rec2_N0sj),true_abun_N0sj)


  #1-6: madara法によるtuning ----

  #  madara法．最小二乗法．
  res_vpa_base0_madara1 <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=FALSE, est.method="ls", b.est=FALSE,abund=c("B","B"),madara=TRUE)

  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_madara1$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_madara1$sigma,2)),true_sd,tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_madara1$saa),2)),c(0.42,0.58,1.00,1.00))

  #  madara法．最尤法．
  res_vpa_base0_madara2 <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=FALSE, est.method="ml", b.est=FALSE,abund=c("B","B"),madara=TRUE)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_madara2$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_madara2$sigma,2)),rep(true_sd,2), tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_madara2$saa),2)),c(0.42,0.58,1.00,1.00))

  #  madara法．最小二乗法実数．
  res_vpa_base0_madara3 <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE,
                               Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.def="max",sel.update=FALSE, est.method="ls_nolog", b.est=FALSE,abund=c("B","B"),madara=TRUE)

  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_madara3$naa),2)),true_number,tol=0.0001)
  expect_equal(as.numeric(round(res_vpa_base0_madara3$sigma,2)),0.5,tol=0.01)
  expect_equal(as.numeric(round(rowMeans(res_vpa_base0_madara3$saa),2)),c(0.42,0.58,1.00,1.00))

  # Part2: dataset is "vpadat_pgc0" for b.est=F, and "vpadat_pgc0_estb" for b.est=T (plus group changes in some years)----
  #現状では，Pope=FALSE（Baranovの方程式）の場合には対応していないのでPope=FALSEのテストは省略

  # 2-1: 二段階法によるtuning ----

  sel.f1<-res_vpa_pgc0_nontune$saa$"2017"
  round(sel.f1,3)

  sel.f3<-res_vpa_pgc_estb_nontune$saa$"2000"
  round(sel.f3,3)

  #二段階法：est.method=最小二乗法による推定 (use.equ="new")
  res_vpa_pgc0_tune1l <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ls", b.est=FALSE,abund=c("B","B"),use.equ="new")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1l$naa),2)),c(3.42,2.41,2.07,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1l$sigma,2)),0.36)

  #二段階法：est.method=最小二乗法による推定 (use.equ="old")
  res_vpa_pgc0_tune1l_o <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE,
                             Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ls", b.est=FALSE,abund=c("B","B"),use.equ="old")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1l_o$naa),2)),c(3.42,2.41,2.07,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1l_o$sigma,2)),0.36)

  #二段階法：est.method=最尤法による推定 (use.equ="new")
  res_vpa_pgc0_tune1m <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE,
                              Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=FALSE,abund=c("B","B"),use.equ="new")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1m$naa),2)),c(3.41,2.40,2.07,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1m$sigma,2)),c(0.39,0.34))

  #二段階法：est.method=最尤法による推定 (use.equ="old")
  res_vpa_pgc0_tune1m_o <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE,
                             Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=FALSE,abund=c("B","B"),use.equ="old")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune1m_o$naa),2)),c(3.41,2.40,2.07,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune1m_o$sigma,2)),c(0.39,0.34))

  #二段階法：est.method=最尤法による推定(2つの指数のsdが異なる場合） (use.equ="new")
  vpa_pgc0_index_change <- vpadat_pgc0
  vpa_pgc0_index_change$index[1,is.na(vpa_pgc0_index_change$index[1,])] <- exp(mean(log(c(1,2))))
  res_vpa_pgc0_index_change_tune1m <- vpa(vpa_pgc0_index_change, tf.year=2015:2016, last.catch.zero = FALSE,
                                           Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=FALSE,abund=c("B","B"), use.equ="new")
  expect_equal(round(as.numeric(rowMeans(res_vpa_pgc0_index_change_tune1m$naa)),2),c(3.45,2.42,2.09,NA))
  expect_equal(round(as.numeric(res_vpa_pgc0_index_change_tune1m$sigma),2),c(0.29,0.35))

  #二段階法：est.method=最尤法による推定(2つの指数のsdが異なる場合） (use.equ="old")
  vpa_pgc0_index_change <- vpadat_pgc0
  vpa_pgc0_index_change$index[1,is.na(vpa_pgc0_index_change$index[1,])] <- exp(mean(log(c(1,2))))
  res_vpa_pgc0_index_change_tune1m_o <- vpa(vpa_pgc0_index_change, tf.year=2015:2016, last.catch.zero = FALSE,
                                          Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f1, est.method="ml", b.est=FALSE,abund=c("B","B"), use.equ="old")
  expect_equal(round(as.numeric(rowMeans(res_vpa_pgc0_index_change_tune1m_o$naa)),2),c(3.45,2.42,2.09,NA))
  expect_equal(round(as.numeric(res_vpa_pgc0_index_change_tune1m_o$sigma),2),c(0.29,0.35))

  #二段階法：est.method=最小二乗法による推定＋指標値の非線形性bの推定(use.equ="new")
  res_vpa_pgc0_estb_tune1l_b <- vpa(vpadat_pgc0_estb,last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3,est.method="ls", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000, use.equ="new")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune1l_b$naa),2)),c(656.24,313.74,160.64,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1l_b$b,2)),c(0.70,0.34,0.55,0.38,0.46,0.78))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1l_b$sigma,2)),0.2)

  #二段階法：est.method=最小二乗法による推定＋指標値の非線形性bの推定(use.equ="old")
  res_vpa_pgc0_estb_tune1l_b_o <- vpa(vpadat_pgc0_estb,last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3,est.method="ls", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000, use.equ="old")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune1l_b_o$naa),2)),c(653.18,311.59,156.09,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1l_b_o$b,2)),c(0.68,0.33,0.53,0.37,0.45,0.76))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1l_b_o$sigma,2)),0.2)

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定(use.equ="new")
  res_vpa_pgc0_estb_tune1m_b <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000, use.equ="new")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune1m_b$naa),2)),c(670.37,318.86,163.79,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b$b,2)),c(0.79,0.39,0.61,0.44,0.51,0.86))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b$sigma,2)),c(0.21,0.18,0.12,0.14,0.18,0.32))

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定(use.equ="old")
  res_vpa_pgc0_estb_tune1m_b_o <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000, use.equ="old")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune1m_b_o$naa),2)),c(667.07,316.62,159.18,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b_o$b,2)),c(0.76,0.37,0.59,0.42,0.50,0.84))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b_o$sigma,2)),c(0.21,0.18,0.12,0.14,0.18,0.32))

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定+指標２と３のシグマは同じ(use.equ="new")
  res_vpa_pgc0_estb_tune1m_b_sigma <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,sigma.constraint=c(1,2,2,3,4,5),use.equ="new")
  expect_equal( res_vpa_pgc0_estb_tune1m_b_sigma$sigma[[2]], res_vpa_pgc0_estb_tune1m_b_sigma$sigma[[3]])

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定+指標２と３のシグマは同じ(use.equ="old")
  res_vpa_pgc0_estb_tune1m_b_sigma_o <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                          Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),fc.year=1998:2000,sigma.constraint=c(1,2,2,3,4,5),use.equ="old")
  expect_equal( res_vpa_pgc0_estb_tune1m_b_sigma_o$sigma[[2]], res_vpa_pgc0_estb_tune1m_b_sigma_o$sigma[[3]])

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定＋一部のbは固定(use.equ="new")
  res_vpa_pgc0_estb_tune1m_b_fix <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),b.fix=c(NA,0.7,NA,NA,NA,1),fc.year=1998:2000,use.equ="new")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune1m_b_fix$naa),2)),c(685.94,324.50,167.26,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b_fix$b,2)),c(0.88,0.70,0.68,0.50,0.57,1.00))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b_fix$sigma,2)),c(0.21,0.20,0.13,0.13,0.18,0.33))

  #二段階法：est.method=最尤法による推定＋指標値の非線形性bの推定＋一部のbは固定(use.equ="old")
  res_vpa_pgc0_estb_tune1m_b_fix_o <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                        Pope = TRUE, p.init = 0.5, tune=TRUE, term.F="max",sel.f=sel.f3, est.method="ml", b.est=TRUE,abund=c("N","N","N","N","N","N"),b.fix=c(NA,0.7,NA,NA,NA,1),fc.year=1998:2000,use.equ="old")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune1m_b_fix_o$naa),2)),c(683.77,322.68,162.90,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b_fix_o$b,2)),c(0.86,0.70,0.66,0.48,0.56,1.00))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune1m_b_fix_o$sigma,2)),c(0.21,0.20,0.13,0.14,0.18,0.33))

  #2-2: 選択率更新法によるtuning ----
  #現状では，＋グループが途中で変わる場合の計算には対応していないのでテストは省略

  #2-3: 全F法によるtuning ----

  #  全F推定法．最小二乗法．
  res_vpa_pgc0_tune4l <- vpa(vpadat_pgc0,  last.catch.zero = FALSE,
                             Pope = TRUE,  tune=TRUE, term.F="all",est.method="ls" ,b.est=FALSE,p.init=0.5,abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4l$naa),2)),c(3.43,2.42,2.09,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4l$sigma,2)),0.36)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4l$saa),2)),c(0.54,0.86,1.00,NA))

  #  全F推定法．最尤法
  res_vpa_pgc0_tune4m <- vpa(vpadat_pgc0, last.catch.zero = FALSE,
                              Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=FALSE,p.init=c(0.2,0.3,0.6,0.6),abund=c("B","B"))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4m$naa),2)),c(3.43,2.43,2.10,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_tune4m$sigma,2)),c(0.40,0.33))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_tune4m$saa),2)),c(0.54,0.86,1.00,NA))

  #  全F推定法．最小二乗法．b推定あり
  res_vpa_pgc0_estb_tune4l_b <- vpa(vpadat_pgc0_estb,last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all", est.method="ls" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune4l_b$naa),2)),c(629.77,296.55,150.07,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune4l_b$b,2)),c(0.55,0.26,0.45,0.31,0.35,0.61))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune4l_b$sigma,2)),c(0.2))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune4l_b$saa),2)),c(0.55,0.89,1.00,NA))

  #  全F推定法．最尤法．b推定あり
  res_vpa_pgc0_estb_tune4m_b <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans( res_vpa_pgc0_estb_tune4m_b$naa),2)),c(633.41,296.91,150.29,NA))
  expect_equal(as.numeric(round( res_vpa_pgc0_estb_tune4m_b$b,2)),c(0.57,0.27,0.46,0.31,0.36,0.62))
  expect_equal(as.numeric(round( res_vpa_pgc0_estb_tune4m_b$sigma,2)),c(0.20,0.18,0.11,0.14,0.18,0.31))
  expect_equal(as.numeric(round(rowMeans( res_vpa_pgc0_estb_tune4m_b$saa),2)),c(0.55,0.89,1.00,NA))

  #  全F推定法．最尤法．b推定あり,指標２と３のシグマは同じ
  res_vpa_pgc0_estb_tune4m_b_sigma <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000,sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_pgc0_estb_tune4m_b_sigma$sigma[[2]],res_vpa_pgc0_estb_tune4m_b_sigma$sigma[[3]])

  # 　全F推定法．最尤法．b推定あり，alpha=0.8(最高齢と最高齢―1のFの比：Fa=alpha*Fa-1)
  res_vpa_pgc0_estb_tune4m_b_alpha <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                     Pope = TRUE,  tune=TRUE, term.F="all",est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),alpha=0.8,fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune4m_b_alpha$naa),2)),c(633.25,298.36,162.34,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune4m_b_alpha$b,2)),c(0.56,0.27,0.46,0.31,0.35,0.61))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune4m_b_alpha$sigma,2)),c(0.19,0.18,0.11,0.14,0.18,0.31))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune4m_b_alpha$saa),2)),c(0.55,0.89,0.88,NA))

  #2-4: Ridge VPA ----
  # set lambda + 全F推定法．最尤法．b推定あり
  res_vpa_pgc0_estb_tune5m_b <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                               Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_b$naa),2)),c(638.85,299.71,152.02,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_b$b,2)),c(0.60,0.29,0.48,0.33,0.38,0.65))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_b$sigma,2)),c(0.20,0.18,0.11,0.14,0.18,0.31))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_b$saa),2)),c(0.56,0.89,1.00,NA))
  
  # set lambda + 全F推定法．最尤法．b推定あり_penalty=p
  res_vpa_pgc0_estb_tune5m_bp <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="p")
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_bp$naa),2)),c(638.85,299.71,152.02,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_bp$b,2)),c(0.60,0.29,0.48,0.33,0.38,0.65))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_bp$sigma,2)),c(0.20,0.18,0.11,0.14,0.18,0.31))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_bp$saa),2)),c(0.56,0.89,1.00,NA))
  
  # set lambda + 全F推定法．最尤法．b推定あり_penalty=f
  res_vpa_pgc0_estb_tune5m_bf <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="f",tf.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_bf$naa),2)),c(636.09,298.29,151.14,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_bf$b,2)),c(0.58,0.28,0.47,0.32,0.37,0.63))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_bf$sigma,2)),c(0.20,0.18,0.11,0.14,0.18,0.31))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_bf$saa),2)),c(0.56,0.89,1.00,NA))
  
  # set lambda + 全F推定法．最尤法．b推定あり_penalty=s
  res_vpa_pgc0_estb_tune5m_bs <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                     Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,penalty="s",tf.year=1998:2000)
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_bs$naa),2)),c(633.63,297.02,150.36,NA))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_bs$b,2)),c(0.57,0.27,0.46,0.31,0.36,0.62))
  expect_equal(as.numeric(round(res_vpa_pgc0_estb_tune5m_bs$sigma,2)),c(0.20,0.18,0.11,0.14,0.18,0.31))
  expect_equal(as.numeric(round(rowMeans(res_vpa_pgc0_estb_tune5m_bs$saa),2)),c(0.55,0.89,1.00,NA))

  # set lambda + 全F推定法．最尤法．b推定あり, 指標２と３のシグマは同じ
  res_vpa_pgc0_estb_tune5m_b_sigma <- vpa(vpadat_pgc0_estb, last.catch.zero = FALSE,  min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                                    Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE, p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"), lambda=0.02, fc.year=1998:2000,sigma.constraint=c(1,2,2,3,4,5))
  expect_equal(res_vpa_pgc0_estb_tune5m_b_sigma$sigma[[2]],res_vpa_pgc0_estb_tune5m_b_sigma$sigma[[3]])

  #現状では，TMBの計算は＋グループが途中で変わるケースには対応していないのでテストは省略する．

  save(res_vpa_base0_nontune,
       res_vpa_base1_nontune,
       res_vpa_pgc0_nontune,
       res_vpa_rec0_nontune,
       res_vpa_base0_nontune_release,
       file="res_vpa_files.rda")
})

