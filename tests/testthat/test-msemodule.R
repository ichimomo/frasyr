test_that("VPA bootstrap procedure", {

  #  data("res_vpa_example2")
  res_vpa_example2 <-  load_data("~/Dropbox/git/frasyr_tool/data_SA2023/kataseto.rda")
  
  biopar <- derive_biopar(res_vpa_example2, derive_year=2017)
  res_ref <- ref.F(res_vpa_example2, Fcurrent=biopar$faa, waa=biopar$waa, maa=biopar$maa, M=biopar$M)
  F30per <- res_ref$summary$FpSPR.30.SPR[3] * biopar$faa
  
  vpa_boot1 <- boo.vpa(res_vpa_example2, B_ite=30, type="index", method="p", out="full") # index bootstrap
  vpa_boot2 <- boo.vpa(res_vpa_example2, B_ite=30, type="caa", B_cv=0, ess=500, out="full") # caa bootstrap (composition)
  vpa_boot3 <- boo.vpa(res_vpa_example2, B_ite=30, type="caa", B_cv=0.1, ess=10000, out="full") # caa bootstrap (catch number)
  # caaのブートストラップが一番影響が大きそう
  plot_boot(vpa_boot1)[[1]] + plot_boot(vpa_boot2)[[1]] + plot_boot(vpa_boot3)[[1]]

  # 残差プロット
  vpa_boot_resid1 <- purrr::map_dfr(vpa_boot1, function(x) plot_residual_vpa(x)$gg_data)
  # カタクチ瀬戸内の場合、最終年の親魚量はほぼ決定論的に決まっている
  vpa_boot_resid1 %>% ggplot() + geom_boxplot(aes(x=factor(year), y=resid))

  xx1 <-   unlist(res_vpa_example2$faa[3,])
  xx2 <- unlist(res_vpa_example2$input$dat$caa[3,]/res_vpa_example2$input$dat$caa[2,])
  plot(xx1~xx2)
  
  SR_boot1  <- do_bootSR(vpa_boot1, res_vpa_example2, type="fix", SR="HS",bio_par = biopar)
  SR_boot2  <- do_bootSR(vpa_boot2, res_vpa_example2, type="fix", SR="HS",bio_par = biopar)
  SR_boot3  <- do_bootSR(vpa_boot2, res_vpa_example2, type="auto",bio_par = biopar)
  #SR_boot4  <- do_bootSR(vpa_boot2, res_vpa_example2, type="regime", regime.year=2000)    
  
  #SR_boot1 %>% map_dfr(function(x) x$pars)
  #SR_boot2 %>% map_dfr(function(x) x$pars)
  #SR_boot3 %>% map_dfr(function(x) x$input$SR)

  tmpfunc <- function(vpa_boot, SR_boot, calc_det=FALSE){
    if(calc_det==FALSE){
      data_future <- make_future_data_boot(res_vpa=res_vpa_example2,nsim=1,
                                           vpa_boot=vpa_boot, SR_boot=SR_boot,
                                           M_year=2017, maa_year=2017, waa_year=2017, waa_catch_year=2017,
                                           plus_group=res_vpa_example2$input$plus.group, #faa_year=2017,
                                           futureF=F30per, currentF=F30per,
                                           start_ABC_year_name=2018)  
      res_future <- future_vpa(data_future$data, SPRtarget=30)
    }
    if(calc_det==TRUE){
      data_future <- make_future_data(res_vpa=res_vpa_example2,nsim=length(vpa_boot),
                                      res_SR=SR_boot[[1]],
                                      M_year=2017, maa_year=2017, waa_year=2017, waa_catch_year=2017,
                                      plus_group=res_vpa_example2$input$plus.group, #faa_year=2017,
                                      futureF=F30per, currentF=F30per,
                                      start_ABC_year_name=2018)
      res_future <- future_vpa(data_future$data, SPRtarget=30)      
    }
    return(res_future)
  }

  res_future <- list()
  res_future[[1]] <- tmpfunc(vpa_boot1, SR_boot1)
  res_future[[2]] <- tmpfunc(vpa_boot2, SR_boot2)
  res_future[[3]] <- tmpfunc(vpa_boot3, SR_boot3)
  res_future[[4]] <- tmpfunc(vpa_boot1, SR_boot1, calc_det=TRUE)
      
  plot_futures(res_vpa_example2, res_future, n_example=0)

  boxplot(t(res_future[[1]]$HCR_realized["2017",,"Fratio"]),)
  
#  hist(res_future1$HCR_realized["2017",,"Fratio"])
  
})
