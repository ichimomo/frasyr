context("test for graphic functions (level 0)")

test_that("plot_futures",{
  plot_futures(vpares=res_vpa_example,future.list=list(res_future_HSL2,res_future_HSL1),
               Btarget=100000, Blimit=50000, Bban=10000) %>% ggplot_build() %>%
    expect_error(NA)
  plot_futures(vpares=NULL,future.list=list(res_future_HSL2,res_future_HSL1)) %>%
    apply_minor_ticks() %>% ggplot_build() %>%
    expect_error(NA)
})

test_that("plot.futures",{
  g1 <- plot.futures(fres.list=list(res_future_HSL2,res_future_HSL1)) 
  expect_equal(class(g1)[1],"list")
})


test_that("plot.future",{
  g1 <- plot.future(res_future_HSL2)
  g2 <- plot.future(res_future_HSL1)
  expect_equal(class(g1)[1],"list")
  expect_equal(class(g2)[1],"list")
})

test_that("plot_vpa",{
  g1 <- plot_vpa(res_vpa_example) 
  g2 <- plot_vpa(list(res_vpa_example,res_vpa_example), is_minor_ticks=FALSE)
  expect_error(ggplot_build(g1), NA)
  expect_error(ggplot_build(g2), NA)
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
  expect_error(ggplot_build(g1), NA)
  expect_equal(class(g2)[1],"NULL")
})

test_that("SRplot_gg", {
    
  g1 <- plot_SR(res_sr_HSL1)
  g2 <- plot_SR(res_sr_HSL2, box.padding=1)
  g3 <- plot_SR(res_sr_HSL1, yscale=100000000, ylabel="億尾")

  expect_error(ggplot_build(g1), NA)
  expect_error(ggplot_build(g2), NA)
  expect_error(ggplot_build(g3), NA)

  # check the effect of tool_issue #529
  res_sr_tmp <- res_sr_HSL1
  res_sr_tmp$pars$rho <- 0.8
  g4 <- plot_SR(res_sr_tmp, plot_CI=TRUE)
  expect_error(ggplot_build(g4), NA)

})

test_that("compare_SRfit",{
  data(res_sr_HSL1)
  data(res_sr_HSL2)
  (g1 <- compare_SRfit(list(HSL1=res_sr_HSL1, HSL2=res_sr_HSL2),
                      biomass.unit=1000, number.unit=1000))
  expect_error(ggplot_build(g1), NA)
})

test_that("plot_SRregime",{
  data(res_vpa_example)
  SRdata <- get.SRdata(res_vpa_example)
  resSRregime <- fit.SRregime(SRdata, SR="HS", method="L2",
                              regime.year=c(1994,2003), regime.key=c(0,1,0),
                              regime.par = c("a","b","sd")[2:3])
  g1 <- plot_SRregime(resSRregime, regime.name=c("Low","High"))
  expect_error(ggplot_build(g1), NA)
  
  (g1 <- compare_SRfit(list(resSRregime, resSRregime),
                       biomass.unit=1000, number.unit=1000,newplot = T))
  expect_error(ggplot_build(g1), NA)
})

test_that("plot_waa", {
  g1 <- plot_waa(res_vpa_example)
  expect_equal(class(g1)[1],"list")
})

test_that("plot_yield", {
  refs.plot <- dplyr::filter(res_MSY_HSL2$summary, RP.definition%in%c("Btarget0", "Blimit0", "Bban0"))

  g1 <- purrr::map(c(TRUE,FALSE),
                   function(x)
                     plot_yield(res_MSY_HSL2$trace,
                                refs.plot,
                                refs.label=NULL,
                                future=list(res_future_HSL2),
                                past=res_vpa_example,label=FALSE,
                                refs.color=rep("black",3),
                                biomass.unit=1000,
                                AR_select=FALSE,
                                past_year_range=NULL,
                                plus_group=x,
                                xlim.scale=0.7,ylim.scale=1.1)
                   )

  expect_error(ggplot_build(g1[[1]]), NA)
  expect_error(ggplot_build(g1[[2]]), NA)
})

test_that("plot_kobe_gg", {
#  load(system.file("extdata","refs_base_pma.rda",package = "frasyr"))
  g1 <- plot_kobe_gg(vpares=res_vpa_example, refs_base=res_MSY_HSL1$summary, ylab_name="Uratio")
  expect_error(ggplot_build(g1), NA)

  FBdata <- tibble(year=1:10, Fratio=exp(rnorm(10,sd=0.1)), Uratio=exp(rnorm(10,sd=0.1)),
                   Bratio=exp(rnorm(10,sd=0.1)), DBratio=exp(rnorm(10,sd=0.1)))
  g1 <- plot_kobe_gg(FBdata=FBdata, refs_base=res_MSY_HSL1$summary)
  expect_error(ggplot_build(g1), NA)
  g1 <- plot_kobe_gg(FBdata=FBdata, refs_base=res_MSY_HSL1$summary, xlab_name="DBratio", ylab_name="Uratio")
  expect_error(ggplot_build(g1), NA)    
  
})

test_that("plot_HCR", {
  g1 <- plot_HCR(SBtarget=1000000, SBlim=300000,SBban=100000, Ftarget=500000)
  expect_error(ggplot_build(g1), NA)
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
  expect_error(ggplot_build(g1), NA)
  expect_error(ggplot_build(g2), NA)
  expect_error(ggplot_build(g3), NA)
  expect_error(ggplot_build(g3.withCI), NA)
})

test_that("compare_MSY",{
  data(res_MSY_HSL1); data(res_MSY_HSL2)
  MSY_list <- tibble::lst(res_MSY_HSL1, res_MSY_HSL2)
  g1 <- compare_MSY(MSY_list)
  expect_error(ggplot_build(g1), NA)
})

#test_that("plot_bias_in_MSE",{
#  g1 <- plot_bias_in_MSE(fout=res_future_HSL2)
#  expect_equal(class(g1)[1],"gg")
#})

#test_that("plot_bias_in_MSE",{
#  fout <- format_to_old_future(res_future_HSL2)
#  g1 <- plot_bias_in_MSE(fout=fout)
#  expect_equal(class(g1)[1],"gg")
#})

test_that("plot_sprypr", {
  g1 <- plot_sprypr(result_vpa=res_vpa_example, type="perspr")
  expect_error(ggplot_build(g1), NA)
})

test_that("plot_SR_AReffect", {
  SRdata <- get.SRdata(res_vpa_example)
  res_SR <- fit.SR(SRdata,SR="BH",AR=1,out.AR=FALSE)
  expect_equal(res_SR$sd.marginal, res_SR$pars$sd/sqrt(1-res_SR$pars$rho^2))  
  (g <- plot_SR_AReffect(res_SR))
  expect_error(ggplot_build(g), NA)
})

test_that("plot_Fcurrent", {

  g <- plot_Fcurrent(res_vpa_example)
  expect_error(ggplot_build(g), NA)

  res_vpa_tmp <- res_vpa_example
  res_vpa_tmp$faa <- rbind(res_vpa_example$faa,res_vpa_example$faa,res_vpa_example$faa,res_vpa_example$faa,res_vpa_example$faa)
  rownames(res_vpa_tmp$faa) <- 0:(nrow(res_vpa_tmp$faa)-1)
  g <- plot_Fcurrent(res_vpa_tmp, Fcurrent=seq(from=0,to=2,length=nrow(res_vpa_tmp$faa)))
  expect_error(ggplot_build(g), NA)
  
})
