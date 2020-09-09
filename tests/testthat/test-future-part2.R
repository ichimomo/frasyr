library(frasyr)
library(testthat)

context("HCR_default") 

test_that("HCR_default",{
  HCR <- HCR_default(ssb=100000,Blimit=10000,Bban=1000,beta=0.8)
  expect_equal(HCR,0.8)
})

context("get_wcatch") 

test_that("get_wcatch",{
  expect_equal(get_wcatch(res_future_0.8HCR), apply(res_future_0.8HCR$wcaa,c(2,3),sum))
})
