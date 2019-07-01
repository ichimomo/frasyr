# https://grasshoppermouse.github.io/2017/10/18/put-your-data-in-an-r-package/
devtools::github_install("ichimomo/frasyr")
library(frasyr)

## 例データの読み込み
caa   <- read.csv("http://cse.fra.affrc.go.jp/ichimomo/fish/caa.csv",  row.names=1)
waa   <- read.csv("http://cse.fra.affrc.go.jp/ichimomo/fish/caa.csv",  row.names=1)
maa   <- read.csv("http://cse.fra.affrc.go.jp/ichimomo/fish/caa.csv",  row.names=1)
M     <- read.csv("http://cse.fra.affrc.go.jp/ichimomo/fish/M.csv",    row.names=1)
index <- read.csv("http://cse.fra.affrc.go.jp/ichimomo/fish/index.csv",row.names=1)

## データの整形
dat <- data.handler(caa=caa, waa=waa, maa=maa, index=index, M=0.4)
## vout3; チューニング，選択率update
  # tf.yearは選択率の初期値として用いられる。
vout3 <- vpa(dat,tune=TRUE,sel.update=TRUE,Pope=FALSE,
             tf.year=1997:1999,sel.f=NULL, 
             abund=c("N"),min.age=c(0),max.age=c(7), # 資源量指数の設定
             alpha=1,p.init=0.5,max.dd = 0.00001,fc.year=1998:2000)
devtools::use_data(vout3)

