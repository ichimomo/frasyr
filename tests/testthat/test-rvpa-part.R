library(frasyr)

context("data input")

test_that("caa, waa and maa data input error check",{
  ex1_caa <- read.csv(system.file("data-raw","ex1_caa.csv",package="frasyr"),row.names=1)
  ex1_waa <- read.csv(system.file("data-raw","ex1_waa.csv",package="frasyr"),row.names=1)
  ex1_maa <- read.csv(system.file("data-raw","ex1_maa.csv",package="frasyr"),row.names=1)
  expect_equal(nrow(ex1_caa),nrow(ex1_waa))
  expect_equal(nrow(ex1_caa),nrow(ex1_maa))
  ex2_caa <- read.csv(system.file("data-raw","ex2_caa.csv",package="frasyr"),row.names=1)
  ex2_waa <- read.csv(system.file("data-raw","ex2_waa.csv",package="frasyr"),row.names=1)
  ex2_maa <- read.csv(system.file("data-raw","ex2_maa.csv",package="frasyr"),row.names=1)
  expect_equal(nrow(ex2_caa),nrow(ex2_waa))
  expect_equal(nrow(ex2_caa),nrow(ex2_maa))
})

context("vpa")

test_that("output value check",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)

  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res_pma_check <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  #読み込んだ結果と照合
  expect_equal(res_pma_check,res_pma)

})


