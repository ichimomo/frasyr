context("test for graphic functions (level 0)")

test_that("plot_futures",{
  g1 <- plot_futures(vpares=res_vpa,future.list=list(res_future_0.8HCR,res_future_current))
  g2 <- plot_futures(vpares=NULL,future.list=list(res_future_0.8HCR,res_future_current))
  expect_equal(class(g1)[1],"gg")
  expect_equal(class(g2)[1],"gg")
})

test_that("plot.futures",{
  g1 <- plot.futures(fres.list=list(res_future_0.8HCR,res_future_current))
  expect_equal(class(g1)[1],"list")
})

test_that("plot.future",{
  g1 <- plot.future(res_future_0.8HCR)
  g2 <- plot.future(res_future_current)
  expect_equal(class(g1)[1],"list")
  expect_equal(class(g2)[1],"list")
})

test_that("plot_vpa",{
  g1 <- plot_vpa(res_vpa)
  g2 <- plot_vpa(list(res_vpa,res_vpa))
  expect_equal(class(g1)[1],"gg")
  expect_equal(class(g2)[1],"gg")
})

test_that("plot_Fref", {
  load(system.file("extdata","res_ref_f_pma.rda",package = "frasyr"))
  g1 <- plot_Fref(res_ref_f_pma)
  expect_equal(class(g1)[1],"data.frame")
})

test_that("plot_SRdata", {
  load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))
  g1 <- plot_SRdata(SRdata_pma, type="gg")
  g2 <- plot_SRdata(SRdata_pma)
  expect_equal(class(g1)[1],"gg")
  expect_equal(class(g2)[1],"NULL") 
})

test_that("SRplot_gg", {
  g1 <- SRplot_gg(res_sr_HSL1)
  g2 <- SRplot_gg(res_sr_HSL2)
  expect_equal(class(g1)[1],"gg")
  expect_equal(class(g2)[1],"gg") 
})

test_that("compare_SRfit",{
  data(res_sr_HSL1)
  data(res_sr_HSL2)
  (g1 <- compare_SRfit(list(HSL1=res_sr_HSL1, HSL2=res_sr_HSL2),
                      biomass.unit=1000, number.unit=1000))
  expect_equal(class(g1)[1],"gg")
})

test_that("SRregime_plot",{
  data(res_vpa)
  SRdata <- get.SRdata(res_vpa)
  resSRregime <- fit.SRregime(SRdata, SR="HS", method="L2",
                              regime.year=c(1994,2003), regime.key=c(0,1,0),
                              regime.par = c("a","b","sd")[2:3])
  g1 <- SRregime_plot(resSRregime, regime.name=c("Low","High"))
  # 本当の意味でのテストにはなっていない
  expect_equal(class(g1)[1],"gg")
  (g1 <- compare_SRfit(list(resSRregime, resSRregime),
                       biomass.unit=1000, number.unit=1000))
  expect_equal(class(g1)[1],"gg")  
})

test_that("SRregime_plot",{
    data(res_sr_HSL1)
    data(res_sr_HSL2)
    (g1 <- compare_SRfit(list(HSL1=res_sr_HSL1, HSL2=res_sr_HSL2),
                         biomass.unit=1000, number.unit=1000))
    expect_equal(class(g1)[1],"gg")      
    (g1 <- compare_SRfit(list(L1=list(res_sr_HSL1,res_sr_HSL1),
                              L2=list(res_sr_HSL1,res_sr_HSL2)),
                         biomass.unit=1000, number.unit=1000))    
    expect_equal(class(g1)[1],"gg")      

})

test_that("plot_waa", {
  g1 <- plot_waa(res_vpa)
  expect_equal(class(g1)[1],"list")
})

#test_that("plot_yield", {
#  refs.plot <- dplyr::filter(res_MSY$summary, RP.definition%in%c("Btarget0", "Blimit0", "Bban0"))
#  (graph_MSY$yield_curve_detail <- plot_yield(res_MSY$trace,
#                                              refs.plot,
#                                              refs.label=label_name_kobe,
#                                              future=list(res_future_0.8HCR),
#                                              past=res_vpa,label=FALSE,
#                                              refs.color=rep("black",3),
#                                              biomass.unit=1000,
#                                              AR_select=FALSE,
#                                              past_year_range=past_year_range_yieldcurve,
#                                              xlim.scale=0.7,ylim.scale=1.1
#                                              ) + theme_SH())
#  refs.plot <- dplyr::filter(res_MSY_HSL1$summary, RP.definition%in%c("Btarget0", "Blimit0", "Bban0"))
#  
#  load(system.file("extdata","refs_base_pma.rda",package = "frasyr"))
#  g1 <- plot_yield(MSY_obj=refs.plot, refs_base=refs_base_pma)
#  
#  g1 <- plot_yield(MSY_obj=res_MSY, refs_base=refs_base_pma)
#  g2 <- plot_yield(MSY_obj=res_MSY_pma, refs_base=refs_base_pma, biomass.unit=1, age.label.ratio = 0.9)
#  expect_equal(class(g1)[1],"gg")
#  expect_equal(class(g2)[1],"gg")
#})

test_that("plot_kobe_gg", {
  load(system.file("extdata","refs_base_pma.rda",package = "frasyr"))
  g1 <- plot_kobe_gg(vpares=res_vpa, refs_base=refs_base_pma)
  expect_equal(class(g1)[1],"gg")
})

test_that("plot_HCR", {
  g1 <- plot_HCR(SBtarget=1000000, SBlim=300000,SBban=100000, Ftarget=500000)
  expect_equal(class(g1)[1],"gg")
})

test_that("compare_eq_stat", {
  data(res_MSY_HSL1)
  data(res_MSY_HSL2)
  MSY_list <- tibble::lst(res_MSY_HSL1, res_MSY_HSL2)
  g1 <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="catch.mean")
  g2 <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="ssb.mean")
  g3 <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="rec.mean")
  gridExtra::grid.arrange(g1,g2,g3,ncol=1)
  g3.withCI <- compare_eq_stat(MSY_list,x_axis_name="fmulti",y_axis_name="rec.mean",plot_CI80=TRUE)
  expect_equal(class(g1)[1],"gg")
  expect_equal(class(g2)[1],"gg")
  expect_equal(class(g3)[1],"gg")
  expect_equal(class(g3.withCI)[1],"gg")
})

test_that("compare_MSY",{
  data(res_MSY)
  MSY_list <- tibble::lst(res_MSY_HSL1, res_MSY_HSL2)
  g1 <- compare_MSY(MSY_list)
  expect_equal(class(g1)[1],"gg")
})

#test_that("plot_bias_in_MSE",{
#  g1 <- plot_bias_in_MSE(fout=res_future_0.8HCR)
#  expect_equal(class(g1)[1],"gg")
#})

#test_that("plot_bias_in_MSE",{
#  fout <- format_to_old_future(res_future_0.8HCR)
#  g1 <- plot_bias_in_MSE(fout=fout)
#  expect_equal(class(g1)[1],"gg")
#})

test_that("compare_kobeII",{
  
})

test_that("plot_sprypr", {
  g1 <- plot_sprypr(result_vpa=res_vpa, type="perspr")
  expect_equal(class(g1)[1],"gg")
})
