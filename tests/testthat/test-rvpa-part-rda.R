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

  #読み込んだ結果と照合
  #expect_equal(res_pma_check,res_pma)
  expect_equal(res_vpa_pma$term.f, res_vpa_pma_check$term.f)
  expect_equal(res_vpa_pma$np, res_vpa_pma_check$np)
  expect_equal(res_vpa_pma$minimum, res_vpa_pma_check$minimum)
  expect_equal(res_vpa_pma$gradient, res_vpa_pma_check$gradient)
  expect_equal(res_vpa_pma$hessian, res_vpa_pma_check$hessian)
  expect_equal(res_vpa_pma$code, res_vpa_pma_check$code)
  for(i in 1:length(res_vpa_pma_check$Fc.at.age)){
    expect_equal(res_vpa_pma$Fc.at.age[i], res_vpa_pma_check$Fc.at.age[i])
  }
  expect_equal(res_vpa_pma$ssb.coef, res_vpa_pma_check$ssb.coef)
  for(i in 1:nrow(res_vpa_pma_check$wcaa)){
    for(j in 1:ncol(res_vpa_pma_check$wcaa)){
      expect_equal(res_vpa_pma$wcaa[i,j], res_vpa_pma_check$wcaa[i,j])
    }
  }
  for(i in 1:nrow(res_vpa_pma_check$naa)){
    for(j in 1:ncol(res_vpa_pma_check$naa)){
      expect_equal(res_vpa_pma$naa[i,j], res_vpa_pma_check$naa[i,j])
    }
  }
  for(i in 1:nrow(res_vpa_pma_check$faa)){
    for(j in 1:ncol(res_vpa_pma_check$faa)){
      expect_equal(res_vpa_pma$faa[i,j], res_vpa_pma_check$faa[i,j])
    }
  }
  for(i in 1:nrow(res_vpa_pma_check$baa)){
    for(j in 1:ncol(res_vpa_pma_check$baa)){
      expect_equal(res_vpa_pma$baa[i,j], res_vpa_pma_check$baa[i,j])
    }
  }
  for(i in 1:nrow(res_vpa_pma_check$saa)){
    for(j in 1:ncol(res_vpa_pma_check$saa)){
      expect_equal(res_vpa_pma$saa[i,j], res_vpa_pma_check$saa[i,j])
    }
  }
  for(i in 1:nrow(res_vpa_pma_check$ssb)){
    for(j in 1:ncol(res_vpa_pma_check$ssb)){
      expect_equal(res_vpa_pma$ssb[i,j], res_vpa_pma_check$ssb[i,j])
    }
  }

})


