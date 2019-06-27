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
  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)
  #上記引数での計算結果を読み込み
  termf_pma_check <- read.csv(system.file("extdata","termf_pma_check.csv",package="frasyr"),row.names=1)
  np_pma_check <- read.csv(system.file("extdata","np_pma_check.csv",package="frasyr"),row.names=1)
  minimum_pma_check <- read.csv(system.file("extdata","minimum_pma_check.csv",package="frasyr"),row.names=1)
  gradient_pma_check <- read.csv(system.file("extdata","gradient_pma_check.csv",package="frasyr"),row.names=1)
  hessian_pma_check <- read.csv(system.file("extdata","hessian_pma_check.csv",package="frasyr"),row.names=1)
  code_pma_check <- read.csv(system.file("extdata","code_pma_check.csv",package="frasyr"),row.names=1)
  Fc_at_age_pma_check <- read.csv(system.file("extdata","Fc_at_age_pma_check.csv",package="frasyr"),row.names=1)
  ssbcoef_pma_check <- read.csv(system.file("extdata","ssbcoef_pma_check.csv",package="frasyr"),row.names=1)
  wcaa_pma_check <- read.csv(system.file("extdata","wcaa_pma_check.csv",package="frasyr"),row.names=1)
  naa_pma_check <- read.csv(system.file("extdata","naa_pma_check.csv",package="frasyr"),row.names=1)
  faa_pma_check <- read.csv(system.file("extdata","faa_pma_check.csv",package="frasyr"),row.names=1)
  baa_pma_check <- read.csv(system.file("extdata","baa_pma_check.csv",package="frasyr"),row.names=1)
  saa_pma_check <- read.csv(system.file("extdata","saa_pma_check.csv",package="frasyr"),row.names=1)
  ssb_pma_check <- read.csv(system.file("extdata","ssb_pma_check.csv",package="frasyr"),row.names=1)

  #読み込んだ結果と照合（一つ一つの数値はしつこい？）

  expect_equal(res.pma$term.f, as.numeric(termf_pma_check[1]))
  expect_equal(res.pma$np, as.numeric(np_pma_check[1]))
  expect_equal(res.pma$minimum, as.numeric(minimum_pma_check[1]))
  expect_equal(res.pma$gradient, as.numeric(gradient_pma_check[1]))
  expect_equal(as.numeric(res.pma$hessian), as.numeric(hessian_pma_check[1]))
  expect_equal(res.pma$code, as.numeric(code_pma_check[1]))

  for(i in 1:nrow(Fc_at_age_pma_check)){
      expect_equal(as.numeric(res.pma$Fc.at.age[i]), Fc_at_age_pma_check[i,])
  }

  expect_equal(res.pma$ssb.coef, as.numeric(ssbcoef_pma_check[1]))

  for(i in 1:nrow(wcaa_pma_check)){
    for(j in 1:ncol(wcaa_pma_check)){
      expect_equal(res.pma$wcaa[i,j], wcaa_pma_check[i,j])
    }
  }
  for(i in 1:nrow(naa_pma_check)){
    for(j in 1:ncol(naa_pma_check)){
      expect_equal(res.pma$naa[i,j], naa_pma_check[i,j])
    }
  }
  for(i in 1:nrow(faa_pma_check)){
    for(j in 1:ncol(faa_pma_check)){
      expect_equal(res.pma$faa[i,j], faa_pma_check[i,j])
    }
  }
  for(i in 1:nrow(baa_pma_check)){
    for(j in 1:ncol(baa_pma_check)){
      expect_equal(res.pma$baa[i,j], baa_pma_check[i,j])
    }
  }
  for(i in 1:nrow(saa_pma_check)){
    for(j in 1:ncol(saa_pma_check)){
      expect_equal(res.pma$saa[i,j], saa_pma_check[i,j])
    }
  }
  for(i in 1:nrow(ssb_pma_check)){
    for(j in 1:ncol(ssb_pma_check)){
      expect_equal(res.pma$ssb[i,j], ssb_pma_check[i,j])
    }
  }
})
