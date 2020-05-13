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
  # この結果が4,3,2,2になるのはなんとなくそんな感じする             
  res_vpa_base0_nontune <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_base0_nontune$naa)), 
               c(4,3,2,2))
  
  res_vpa_base1_nontune <- vpa(vpadat_base1, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_base1_nontune$naa)), 
               c(4,3,2,2))
  
  res_vpa_pgc0_nontune <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(unlist(res_vpa_pgc0_nontune$naa["2017"])), 
               c(3,2,2,NA), tol=0.0001)
  
  res_vpa_rec0_nontune <- vpa(vpadat_rec0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_rec0_nontune$naa)), 
               c(4,3,2,2))
  
  # catch計算用のwaaを２倍にしているbase1データでは漁獲量が倍になる
  expect_equal(res_vpa_base0_nontune$wcaa*2,
               res_vpa_base1_nontune$wcaa)
  
  # vpa (tuning) ----
  # この記述は正しいか不明。計算しただけ。
  res_vpa_base0_tune <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                            Pope = TRUE, p.init = 0.5, tune=TRUE, sel.update=TRUE)

  save(res_vpa_base0_nontune,
       res_vpa_base1_nontune,
       res_vpa_pgc0_nontune,
       res_vpa_rec0_nontune,
       res_vpa_base0_tune,file="res_vpa_files.rda")
})
