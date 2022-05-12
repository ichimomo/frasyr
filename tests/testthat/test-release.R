test_that("read_vpa with release data",{

  # read.vpaによる放流データ利用
  # vpa with relase data (only with alived relased fish; old version)
  res_vpa_R1 <- read.vpa(system.file("extdata","res_vpa_dummy_release_fish.csv",package="frasyr"))
  expect_equal(res_vpa_R1$input$dat$release.alive %>% as.numeric(), rep(1000, length(1973:2019)))
  res_SRdat_R1 <- res_vpa_R1 %>% get.SRdata(weight.year = 1981:2018)
  expect_equal(names(res_SRdat_R1),
               c("year","SSB","R","weight","release_alive"))

  # vpa with relase data (with alive rate and  relased fish)
  res_vpa_R2 <- read.vpa(system.file("extdata","res_vpa_dummy_release_fish_new.csv",package="frasyr"))
  
  expect_equal(res_vpa_R2$input$dat$release.all %>% as.numeric(), rep(10000, length(1973:2019)))
  expect_equal(res_vpa_R2$input$dat$release.alive %>% as.numeric(), rep(1000, length(1973:2019)))
  expect_equal(res_vpa_R2$input$dat$release.ratealive %>% as.numeric(), rep(0.1, length(1973:2019)))
  
  res_SRdat_R2 <- res_vpa_R2 %>% get.SRdata(weight.year = 1981:2018)
  expect_equal(names(res_SRdat_R2),
               c("year","SSB","R","weight","release_alive","release_all","release_ratealive"))  
  expect_equal(res_SRdat_R1$R,res_SRdat_R2$R)
  
  # normal vpa
  res_SRdat_R0 <- read.vpa(system.file("extdata","res_vpa_dummy.csv",package="frasyr")) %>%
      get.SRdata()
  expect_equal(names(res_SRdat_R0), c("year","SSB","R","weight"))

  # normal vpa with weight data
  res_SRdat_R0_8118 <- read.vpa(system.file("extdata","res_vpa_dummy.csv",package="frasyr")) %>%
      get.SRdata(weight.year = 1981:2018)  
  expect_equal(names(res_SRdat_R0_8118),
               c("year","SSB","R","weight"))
  
  # data.handlerによる放流データ利用
  load("res_vpa_files.rda")
  load("data_future_test.rda")  
  data_SR_release <- res_vpa_base0_nontune_release %>%
      get.SRdata(weight.year=1990:2100)
  res_SR_release <- fit.SR(SRdata=data_SR_release, AR=0)

  data_SR <- res_vpa_base0_nontune %>%
      get.SRdata(weight.year=1990:2100)
  res_SR <- fit.SR(SRdata=data_SR, AR=0)  

  expect_equal(names(data_SR_release),
               c("year","SSB","R","weight","release_alive"))
  expect_equal(res_SR_release$pars$a*res_SR_release$pars$b,2,tol=0.001)

  # 過去・将来に放流がある場合の将来予測
  res_future_release <- redo_future(data_future_test,
                                    list(res_SR=res_SR_release, recruit_intercept=2))
  res_future <- redo_future(data_future_test,list(res_SR=res_SR))
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
  
})
