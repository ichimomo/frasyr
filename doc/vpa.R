## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.show='hold'----------------------------------------------------
library(frasyr)
caa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex2_caa.csv",  row.names=1)
waa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex2_waa.csv",  row.names=1)
maa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex2_maa.csv",  row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)

# VPAによる資源量推定
res_vpa <- vpa(dat,fc.year=2015:2017,tf.year = 2015:2016,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=0.5)
res_vpa$Fc.at.age # 将来予測やMSY計算で使うcurrent F (fc.yearのオプションでいつのFの平均かが指定される)
plot(res_vpa$Fc.at.age,type="b",xlab="Age",ylab="F",ylim=c(0,max(res_vpa$Fc.at.age)))

## ------------------------------------------------------------------------
#out.vpa(res_vpa) # vpa.csvというファイルが作成されます。VPAの結果のグラフ出力となるvpa.pdfも出力されます。
res_vpa2 <- read.vpa("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex2_vpa.csv") # vpa.csvを編集後、read.vpa関数で読み込みます

## ------------------------------------------------------------------------
# 読み込み
caa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_caa.csv",  row.names=1)
waa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_waa.csv",  row.names=1)
maa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_maa.csv",  row.names=1)
M     <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_M.csv",  row.names=1)
index     <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_index.csv",  row.names=1)

# 整形
dat <- data.handler(caa=caa, waa=waa, maa=maa, index=index, M=0.4)

## ------------------------------------------------------------------------
  # p.initは初期値。収束しない場合（Fがゼロになる場合など）は、初期値を変えてみる
  # fc.yearは、Fcurrentを計算する範囲。管理基準値と将来予測で使われる。
vout1 <- vpa(dat,tf.year=1997:1999,Pope=TRUE,fc.year=1998:2000,alpha=1,p.init=0.5) 
vout1a <- vpa(dat,tf.year=1997:1999,Pope=FALSE,fc.year=1998:2000,alpha=1,p.init=0.5)
vout1b <- vpa(dat,tf.year=1997:1999,Pope=TRUE,fc.year=1998:2000,alpha=0.5,p.init=0.5)
vout1c <- vpa(dat,tf.year=1999:1999,Pope=TRUE,fc.year=1998:2000,alpha=1,p.init=0.5) 

## ------------------------------------------------------------------------
## vout2; チューニング，選択率updateなし
vout2 <- vpa(dat,tune=TRUE,sel.update=FALSE,Pope=FALSE,
             tf.year=NULL,sel.f=vout1$saa$"2000", # 選択率の仮定
             abund=c("N"),min.age=c(0),max.age=c(6), # 資源量指数の設定
             alpha=1,p.init=0.5,max.dd = 0.00001,fc.year=1998:2000)

## vout3; チューニング，選択率update
  # tf.yearは選択率の初期値として用いられる。
vout3 <- vpa(dat,tune=TRUE,sel.update=TRUE,Pope=FALSE,
             tf.year=1997:1999,sel.f=NULL, 
             abund=c("N"),min.age=c(0),max.age=c(7), # 資源量指数の設定
             alpha=1,p.init=0.5,max.dd = 0.00001,fc.year=1998:2000)

## チューニング，選択率全推定
  # tf.yearも sel.fも必要ない
vout4 <- vpa(dat,tune=TRUE,sel.update=FALSE,term.F="all",
             tf.year=NULL,sel.f=NULL,
             abund=c("N"),min.age=c(0),max.age=c(6), # 資源量指数の設定
             alpha=1,p.init=0.5,max.dd = 0.00001,fc.year=1998:2000)

## ------------------------------------------------------------------------
ci0 <- profile_likelihood.vpa(vout3, method="ci",Alpha=0.80)$ci

## ------------------------------------------------------------------------
## ノンパラメトリックブートストラップ(method="n")
set.seed(1)
boot.sim1 <- boo.vpa(vout3,B=10,method="n")
  # boo.vpaはブートストラップ回数分のvpa関数の返り値のリストを返す
  # 値の取り出しは、リストの操作関数sapplyまたはlapplyを用いる
  # Bで繰り返し回数を指定します．ここでは10回ですが実際には1000回(B=1000)以上やってください
tf.dist1 <- sapply(boot.sim1,function(x) x$faa["2000"][7,])
ci1 <- quantile(tf.dist1,probs=c(0.1,0.9))

## パラメトリックブートストラップ(method="p")
set.seed(1)
boot.sim2 <- boo.vpa(vout3,B=10,method="p") # 実際には1000回以上(B=1000)やってください
tf.dist2 <- sapply(boot.sim2,function(x) x$faa["2000"][6,])
ci2 <- quantile(tf.dist2,probs=c(0.1,0.9))

## 平滑化ブートストラップ(method="r")
set.seed(1) 
boot.sim3 <- boo.vpa(vout3,B=10,method="r") # 実際には1000回以上(B=1000)やってください
tf.dist3 <- sapply(boot.sim3,function(x) x$faa["2000"][6,])
ci3 <- quantile(tf.dist3,probs=c(0.1,0.9))

## 4つの信頼区間の比較
rbind(ci0,ci1,ci2,ci3)

## ノンパラメトリックブートストラップ for vout4
set.seed(1)
boot.sim4 <- boo.vpa(vout4,B=10,method="n") # 実際には1000回以上(B=1000)やってください

tf.dist4 <- sapply(boot.sim4[boot.sim4!="try-error"],function(x) x$faa["2000"][7,])
ci4 <- quantile(tf.dist4,probs=c(0.1,0.9))

## 親魚量の信頼区間のプロット
Years <- colnames(dat$caa)
ssb.boot <- sapply(boot.sim1,function(x) colSums(x$ssb))
x <- t(apply(ssb.boot,1,quantile,probs=c(0.1,0.5,0.9)))
matplot(Years,x,ylim=c(0,max(x)),col=1,type=c("l","b","l"),
        pch=1,lty=c(2,1,2),ylab="Spawning biomass")

## ------------------------------------------------------------------------
resid <- as.numeric(log(vout3$pred.index) - log(dat$index))
plot(resid,type="b")
acf(resid) # 自己相関は特にない
# 正規性の検定
ks.test(resid,"pnorm",mean=mean(resid),sd=sd(resid))
# 分布が有意に正規分布から外れているわけではない
ks.test(c(resid,10),"pnorm",mean=mean(c(resid,10)),sd=sd(c(resid,10)))
# 大きな外れ値があると、p値が小さくなって、正規分布でない、となる。

