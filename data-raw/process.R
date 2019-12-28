library(devtools)
devtools::load_all() #library(frasyr)
# データの読み込みと資源量推定
caa   <- read.csv("ex2_caa.csv",  row.names=1)
waa   <- read.csv("ex2_waa.csv",  row.names=1)
maa   <- read.csv("ex2_maa.csv",  row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)

# VPAによる資源量推定
res_vpa <- vpa(dat,fc.year=2015:2017,tf.year = 2015:2016,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=0.5)
devtools::use_data(res_vpa)

# SR関数のフィット
SRdata <- get.SRdata(res_vpa)
res_sr_HSL2 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
use_data(res_sr_HSL2)

res_sr_HSL1 <- fit.SR(SRdata,SR="HS",method="L1",AR=0,hessian=FALSE)
use_data(res_sr_HSL1)

res_future_HSL2 <- future.vpa(res_vpa,
                      multi=1,
                      nyear=50, # 将来予測の年数
                      start.year=2011, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                      ABC.year=2019, # ABCを計算する年
                      waa.year=2015:2017, # 生物パラメータの参照年
                      maa.year=2015:2017,
                      M.year=2015:2017,
                      is.plot=TRUE, # 結果をプロットするかどうか
                      seed=1,
                      silent=FALSE,
                      recfunc=HS.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=res_sr_HSL2$pars$a,b=res_sr_HSL2$pars$b,
                                   rho=res_sr_HSL2$pars$rho, # ここではrho=0なので指定しなくてもOK
                                   sd=res_sr_HSL2$pars$sd,resid=res_sr_HSL2$resid))

res_MSY_HSL2 <- est.MSY(res_vpa, # VPAの計算結果
                 res_future_HSL2$input, # 将来予測で使用した引数
                 resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                 N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                 calc.yieldcurve=TRUE,
                 PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                 onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                 B0percent=c(0.2,0.3,0.4),
                 Bempirical=c(round(tail(colSums(res_vpa$ssb),n=1)),
                              round(max(colSums(res_vpa$ssb))),
                              24000, # 現行Blimit
                              res_sr_HSL2$pars$b) # HSの折れ点
                 ) # 計算したいB0%レベル

res_future_HSL1 <- future.vpa(res_vpa,
                      multi=1,
                      nyear=50, # 将来予測の年数
                      start.year=2011, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                      ABC.year=2019, # ABCを計算する年
                      waa.year=2015:2017, # 生物パラメータの参照年
                      maa.year=2015:2017,
                      M.year=2015:2017,
                      is.plot=TRUE, # 結果をプロットするかどうか
                      seed=1,
                      silent=FALSE,
                      recfunc=HS.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=res_sr_HSL1$pars$a,b=res_sr_HSL1$pars$b,
                                   rho=res_sr_HSL1$pars$rho, # ここではrho=0なので指定しなくてもOK
                                   sd=res_sr_HSL1$pars$sd,resid=res_sr_HSL1$resid))

res_MSY_HSL1 <- est.MSY(res_vpa, # VPAの計算結果
                 res_future_HSL1$input, # 将来予測で使用した引数
                 resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                 N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                 calc.yieldcurve=TRUE,
                 PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                 onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                 B0percent=c(0.2,0.3,0.4),
                 Bempirical=c(round(tail(colSums(res_vpa$ssb),n=1)),
                              round(max(colSums(res_vpa$ssb))),
                              24000, # 現行Blimit
                              res_sr_HSL1$pars$b) # HSの折れ点
                 ) # 計算したいB0%レベル

use_data(res_MSY_HSL1)
use_data(res_MSY_HSL2)
