context("Making tables")

result_vpa <- load_data("../../inst/extdata/res_vpa_pma.rda")
result_msy <- load_data("../inst/extdata/res_MSY_pma_pre.rda")

assertthat::assert_that(
  head(colnames(result_vpa$ssb), 1) == "1982",
  tail(colnames(result_vpa$ssb), 1) == "2011"
)


test_that("make_stock_table() works", {
  recent_years <- 2008:2011
  tbl <- make_stock_table(result_vpa, result_msy,
                          yr_pre_abc = recent_years)
  expect_is(tbl, "data.frame")
  expect_equal(colnames(tbl),
               c("Year", "Biomass", "SSB", "Catch", "F.Fmsy", "HarvestRate"))
  expect_equal(tbl$Year,
               c(recent_years, 2012))
})
