test_that("read_vpa with release data",{

  # read.vpaによる放流データ利用
  # vpa with relase data (only with alived relased fish; old version)
  res_vpa_R1 <- read.vpa(system.file("extdata","res_vpa_dummy_release_fish.csv",package="frasyr"))
  expect_equal(res_vpa_R1$input$dat$release.alive %>% as.numeric(), rep(1000, length(1973:2019)))
  SRdata_R1 <- res_vpa_R1 %>% get.SRdata(weight.year = 1981:2018)
  expect_equal(names(SRdata_R1),
               c("year","SSB","R","weight","release_alive"))

  # vpa with relase data (with alive rate and  relased fish)
  res_vpa_R2 <- read.vpa(system.file("extdata","res_vpa_dummy_release_fish_new.csv",package="frasyr"))
  expect_equal(res_vpa_R2$input$dat$release.all %>% as.numeric(), rep(1000, length(1973:2019)))
  expect_equal(res_vpa_R2$input$dat$release.alive %>% as.numeric(), rep(100, length(1973:2019)))
  expect_equal(res_vpa_R2$input$dat$release.ratealive %>% as.numeric(), rep(0.1, length(1973:2019)))
  
  SRdata_R2 <- res_vpa_R2 %>% get.SRdata()
  expect_silent(plot_SRdata(SRdata_R2))  
  expect_equal(names(SRdata_R2),
               c("year","SSB","R","weight","release_alive","release_all","release_ratealive"))  
  
  # normal vpa
  SRdata_R0 <- read.vpa(system.file("extdata","res_vpa_dummy.csv",package="frasyr")) %>%
      get.SRdata()
  expect_equal(names(SRdata_R0), c("year","SSB","R","weight"))

  # normal vpa with weight data
  SRdata_R0_weight <- read.vpa(system.file("extdata","res_vpa_dummy.csv",package="frasyr")) %>%
      get.SRdata(weight.year = 1981:2018)  
  expect_equal(names(SRdata_R0_weight),
               c("year","SSB","R","weight"))
  
  # data.handlerによる放流データ利用
  load("res_vpa_files.rda")
  load("data_future_test.rda")  
  SRdata_from_vpa_release <- res_vpa_base0_nontune_release %>%
      get.SRdata(weight.year=1990:2100)
  res_SR_release <- fit.SR(SRdata=SRdata_from_vpa_release, AR=0)

  SRdata_normal<- res_vpa_base0_nontune %>%
      get.SRdata(weight.year=1990:2100)
  res_SR_normal <- fit.SR(SRdata=SRdata_normal, AR=0)  

  expect_equal(names(SRdata_from_vpa_release),
               c("year","SSB","R","weight","release_alive"))
  expect_equal(res_SR_release$pars$a*res_SR_release$pars$b,2,tol=0.001)

  # 過去・将来に放流がある場合の将来予測
  res_future_release <- redo_future(data_future_test,
                                    list(res_SR=res_SR_release, recruit_intercept=2))
  res_future <- redo_future(data_future_test,list(res_SR=res_SR_normal))
  plot.SR(res_SR_release, recruit_intercept=2)  
  expect_equal(select(res_future$summary,-intercept),select(res_future_release$summary,-intercept))

  # 過去・将来に放流があり、SD>0、自己相関>0の場合
  res_vpa_base0_nontune_release2 <- res_vpa_base0_nontune_release
  nyear <- ncol(res_vpa_base0_nontune_release2$naa)
  res_vpa_base0_nontune_release2$naa[1,nyear] <- (res_vpa_base0_nontune_release2$naa[1,nyear] - res_vpa_base0_nontune_release2$input$dat$release[nyear]) * exp(log(2)) + res_vpa_base0_nontune_release2$input$dat$release[nyear]
  res_SR_release2 <- res_SR_release
  rho_assumption <- 0.3  
  res_SR_release2$pars$rho <- rho_assumption
  res_SR_release2$pars$sd <- 0.1
  data_future_release2 <- redo_future(data_future_test,
                                     list(res_SR=res_SR_release2, recruit_intercept=2,
                                          res_vpa=res_vpa_base0_nontune_release2),
                                     only_data=T)
  res_future_release2 <- test_sd0_future(data_future_release2)  
  future_deviance1 <- res_future_release2$res2$SR_mat[as.character(2017:2026),1,"deviance"]
  expected_deviance <- future_deviance1["2017"]*(rho_assumption^(0:9))
  expect_equal(as.numeric(exp(future_deviance1)[1]), 2, tol=0.0001)
  expect_equal(as.numeric(future_deviance1/expected_deviance), rep(1,10), tol=0.0001)

  future_deviance2 <- res_future_release2$res1$summary %>% filter(year%in%2017:2026) %>%
      select(deviance) %>% unlist %>% as.numeric
  expected_deviance <- future_deviance2[1]*(rho_assumption^(0:9))    
  expect_equal(future_deviance2-expected_deviance, rep(0,10), tol=0.01)
  expect_equal(future_deviance2-as.numeric(future_deviance1), rep(0,10), tol=0.01)

  # より詳細な放流シナリオの設定: それぞれ過去年のみを設定
  res_SR_R2 <- fit.SR(SRdata_R2)
  expect_silent(plot_SR(res_SR_R2))
  
  setting_release <- list(
      number             = tibble(year=2015:2017), # or value= # or value= & year=
      rate               = tibble(year=2015:2017)) # or value=
  data_future_test$input$setting_release <- setting_release
  data_future_R2 <- redo_future(data_future_test,
                                list(res_SR = res_SR_R2, 
                                     res_vpa= res_vpa_R2),
                                only_data=T)
  data_future_R2$data$SR_mat[,,"intercept"] %>% mean() %>%
      expect_equal(mean(res_SR_R2$input$SRdata$release_alive))
  res_future_R2 <- future_vpa(data_future_R2$data)
  expect_equal(mean(res_future_R2$summary$intercept),100)

  res_SR_R2 <- fit.SR(SRdata_R2)
  expect_silent(plot_SR(res_SR_R2))

  # 値で指定
  setting_release <- list(
      number             = tibble(value=200), # or value= # or value= & year=
      rate               = tibble(value=c(0.2,0.5,0.8)))
  
  data_future_test$input$setting_release <- setting_release
  data_future_R3 <- redo_future(data_future_test,
                                list(res_SR = res_SR_R2, 
                                     res_vpa= res_vpa_R2),
                                only_data=T)
  res_intercept <- data_future_R3$data$SR_mat[,,"intercept"]
  res_intercept[as.character(data_future_test$input$start_random_rec_year_name:tail(dimnames(res_intercept)[[1]],n=1)),] %>%
      as.numeric() %>% unique() %>% sort() %>%
      expect_equal(setting_release$number$value * setting_release$rate$value)

  # 年と値で指定
  setting_release <- list(
      number             = tibble(year=2018:(2018+199), value=c(0,100,rep(200,198))), 
      rate               = tibble(value=c(0.2,0.5,0.8)))
  data_future_test$input$setting_release <- setting_release
  data_future_R3 <- redo_future(data_future_test,
                                list(res_SR = res_SR_R2, 
                                     res_vpa= res_vpa_R2),
                                only_data=T)
  res_intercept <- data_future_R3$data$SR_mat[,,"intercept"]
#  res_intercept[as.character(data_future_test$input$start_random_rec_year_name:tail(dimnames(res_intercept)[[1]],n=1)),]
  expect_equal(as.numeric(res_intercept["2018",]), rep(0,10))
  expect_equal(unique(as.numeric(res_intercept["2019",])) %>% sort(),
               c(0.2,0.5,0.8) * 100)
  expect_equal(unique(as.numeric(res_intercept[as.character(2020:2030),])) %>% sort(),
               c(0.2,0.5,0.8) * 200)

  # VPA結果が１年更新される場合
  res_vpa_R3 <- create_dummy_vpa(res_vpa_R2)
  res_vpa_R3$input$dat$release.all["2020"][] <- 3000
  res_vpa_R3$input$dat$release.alive["2020"][] <- 300

  setting_release <- list(
      number             = tibble(year=2020), # or value= # or value= & year=
      rate               = tibble(year=2015:2017),
      data_source="VPA") # or value=
  data_future_test$input$setting_release <- setting_release  
  
  data_future_R2 <- redo_future(data_future_test,
                                list(res_SR = res_SR_R2, 
                                     res_vpa= res_vpa_R3,
                                     start_random_rec_year_name=2021),
                                only_data=T)
  data_future_R2$data$SR_mat[as.character(2020:2037),,"intercept"] %>% mean() %>%
      expect_equal(300)
  
  res_future_R2 <- future_vpa(data_future_R2$data)
  res_future_R2$summary %>% dplyr::filter(year>2019) %>% select(intercept) %>% unlist() %>% mean() %>% expect_equal(300)

  res_SR_R2 <- fit.SR(SRdata_R2)
  expect_silent(plot_SR(res_SR_R2))  
  
})
