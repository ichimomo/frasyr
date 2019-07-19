library(frasyr)
# データの読み込みと資源量推定
caa   <- read.csv("data-raw/ex2_caa.csv",  row.names=1)
waa   <- read.csv("data-raw/ex2_waa.csv",  row.names=1)
maa   <- read.csv("data-raw/ex2_maa.csv",  row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)

# VPAによる資源量推定
res_vpa <- vpa(dat,fc.year=2015:2017,tf.year = 2015:2016,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=0.5)
devtools::use_data(res_vpa)

