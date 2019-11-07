library(devtools)
devtools::load_all() #library(frasyr)
# データの読み込みと資源量推定
caa   <- read.csv("data-raw/ex2_caa.csv",  row.names=1)
waa   <- read.csv("data-raw/ex2_waa.csv",  row.names=1)
maa   <- read.csv("data-raw/ex2_maa.csv",  row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)

# VPAによる資源量推定
res_vpa <- vpa(dat,fc.year=2015:2017,tf.year = 2015:2016,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=0.5)
devtools::use_data(res_vpa)

# SR関数のフィット
SRdata <- get.SRdata(res_vpa)
res_sr_HSL2 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
use_data(res_sr_HSL2)
