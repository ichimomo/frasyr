data_base <- readr::read_csv(system.file("extdata","all_dummy_data_base.csv",package="frasyr"))

test_that("read_vpa with release data",{

  # read.vpaによる放流データ利用
  # vpa with relase data
  res_SRdat_R1 <- read.vpa(system.file("extdata","res_vpa_dummy_release_fish.csv",package="frasyr")) %>%
      get.SRdata(weight.year = 1981:2018)
  expect_equal(names(res_SRdat_R1),
               c("year","SSB","R","release","weight"))
  
  # normal vpa
  res_SRdat_R0 <- read.vpa(system.file("extdata","res_vpa_dummy.csv",package="frasyr")) %>%
      get.SRdata()
  expect_equal(names(res_SRdat_R0),
               c("year","SSB","R"))

  # normal vpa with weight data
  res_SRdat_R0_8118 <- read.vpa(system.file("extdata","res_vpa_dummy.csv",package="frasyr")) %>%
      get.SRdata(weight.year = 1981:2018)  
  expect_equal(names(res_SRdat_R0_8118),
               c("year","SSB","R","weight"))
  
  # data.handlerによる放流データ利用
  data_release <- to_vpa_data(data_base, label_name="caa")[1,]
  release.number <- 2
  data_release[] <- release.number
  
  res_SRdat_D1 <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
               waa=to_vpa_data(data_base, label_name="waa"),
               maa=to_vpa_data(data_base, label_name="maa"),
               M  = 0.4,
               index = to_vpa_data(data_base, label_name="abund"),
               maa.tune = NULL,
               waa.catch = NULL,
               catch.prop = NULL,
               release.dat = data_release) %>%
    vpa(tf.year=2015:2016, last.catch.zero = FALSE, 
        Pope = TRUE, p.init = 0.5) %>%
      get.SRdata(weight.year=1990:2100)
  expect_equal(names(res_SRdat_D1),
               c("year","SSB","R","release","weight"))
  
  res_SRdat_D0 <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
                             waa=to_vpa_data(data_base, label_name="waa"),
                             maa=to_vpa_data(data_base, label_name="maa"),
                             M  = 0.4,
                             index = to_vpa_data(data_base, label_name="abund"),
                             maa.tune = NULL,
                             waa.catch = NULL,
                             catch.prop = NULL, release.dat = NULL) %>%
    vpa(tf.year=2015:2016, last.catch.zero = FALSE, 
        Pope = TRUE, p.init = 0.5) %>%
      get.SRdata(weight.year=1990:2100)
  expect_equal(names(res_SRdat_D0),
               c("year","SSB","R","weight"))  

  expect_equal(mean(res_SRdat_D0$R-res_SRdat_D1$R), release.number)
  expect_equal(mean(res_SRdat_R0$R-res_SRdat_R1$R), 990.491, tol=0.0001) 

  # get.SRdataの段階でRが差し引かれているためこのあと、fit.SRは変更しなくてOK
  SR0.par <- fit.SR(SRdata=res_SRdat_D0, AR=0)
  SR1.par <- fit.SR(SRdata=res_SRdat_D1, AR=0)

  expect_equal(SR0.par$pars$b*SR0.par$pars$a,mean(res_SRdat_D0$R),tol=0.001)
  expect_equal(SR1.par$pars$b*SR1.par$pars$a,mean(res_SRdat_D1$R),tol=0.001)



  
})
