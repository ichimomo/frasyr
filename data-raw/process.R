library(frasyr)
# データの読み込みと資源量推定
caa <- read.csv("../data-raw/caa_pma.csv",row.names=1)
waa <- read.csv("../data-raw/waa_pma.csv",row.names=1)
maa <- read.csv("../data-raw/maa_pma.csv",row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)

# VPAによる資源量推定
res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

devtools::use_data(res.pma)

