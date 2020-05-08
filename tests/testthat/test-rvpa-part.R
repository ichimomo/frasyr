library(frasyr)

context("data input")

test_that("caa, waa and maa data input error check",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
  expect_equal(nrow(caa),nrow(waa))
  expect_equal(nrow(caa),nrow(maa))
})

context("vpa")

test_that("output value check",{
  caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
  waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
  maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)

  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
  res_vpa_pma_check <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                       term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  #照合内容

  testcontents <-c("term.f","np","minimum","gradient","hessian","code","Fc.at.age","ssb.coef","wcaa","naa","faa","baa","saa","ssb")

  #読み込んだ結果と照合

  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("res_vpa_pma$",testcontents[i]))),eval(parse(text=paste("res_vpa_pma_check$",testcontents[i]))))
  }

})


