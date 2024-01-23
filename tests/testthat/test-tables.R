context("Making tables")

#result_vpa    <- load_data("../../inst/extdata/res_vpa_pma.rda")
#result_msy    <- load_data("../../inst/extdata/res_MSY_pma_pre.rda")
#result_future <- load_data(".../../inst/extdata/res_future_Fcurrent_pma.rda")
yrs_bio       <- 2009:2011
#assertthat::assert_that(
#  head(colnames(result_vpa$ssb), 1) == "1982",
#  tail(colnames(result_vpa$ssb), 1) == "2011"
#)
yrs_pre_abc <- 2008:2011

#! make_stock_table have been moved to make_stock_table2
test_that("make_stock_table() works", {
##   tbl <- make_stock_table2(result_vpa        = res_vpa_example,
##                           result_msy        = res_MSY_HSL2,
##                           result_future     = generate_dummy_future_new_object(),
##                           yr_future_start   = 2012,
##                           yr_abc            = 2013,
##                           faa_pre_abc       = c(0.123, 0.234, 0.345, 0.345))

##   expect_df(tbl)
##   expect_equal(colnames(tbl),
##                c("Year", "Biomass", "SSB", "Catch", "F.Fmsy", "HarvestRate"))
##   expect_equal(tbl$Year,
##                2008:2013) # Six years until (yr_future_start + 1)
    ##   expect_true(all(!is.na(tbl$Biomass)))

    aa <- make_stock_table2(res_future_HSL1,
                      yr_future_ABC=2020,
                      0.1, 2015:2017)
    expect_equal(nrow(aa), 7)
    
})


