context("Utils for selecting SR models")

data(res_vpa)

test_that("get.SRdata returns data frame", {
  expect_equal(class(get.SRdata(res_vpa, return.df = TRUE)),
               "data.frame")
  expect_equal(dim(get.SRdata(res_vpa, return.df = TRUE)),
               c(30, 3))
})

test_that("get.SRdata returns list", {
  expect_false(class(get.SRdata(res_vpa)) == "data.frame")
  expect_equal(length(get.SRdata(res_vpa)), 3)
  expect_null(dim(get.SRdata(res_vpa)))
})

test_that("get data frame from list object", {
  expect_equal(class(pull_df_from_list(list = res_vpa, dfname = "baa")),
               "data.frame")
  expect_equal(dim(pull_df_from_list(list = res_vpa, dfname = "baa")),
               c(4, 30))
})

test_that("Obtain SR data for plot quickly", {
  output <- make_SR_dframe("HS", "L1", res_vpa)
  expect_is(output, "data.frame")
  expect_equal(nrow(output), 100)
  expect_equal(ncol(output), 3)

  inputs_are <- function(SR, method, return_is) {
    expect_equal(make_SR_dframe(SR = SR, method = method, res_vpa) %>%
                   dplyr::pull(name) %>%
                   unique(),
                 return_is)
  }
  inputs_are("HS", "L1", return_is = "HS_L1")
  inputs_are("HS", "L2", return_is = "HS_L2")
  inputs_are("BH", "L1", return_is = "BH_L1")
  inputs_are("BH", "L2", return_is = "BH_L2")
  inputs_are("RI", "L1", return_is = "RI_L1")
  inputs_are("RI", "L2", return_is = "RI_L2")

})
