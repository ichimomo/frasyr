## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5,
  fig.height=5
)

## ----SRdata--------------------------------------------------------------
# ライブラリとデータの読み出し
library(frasyr)
data(res_vpa)
# VPA結果を使って再生産データを作る
SRdata <- get.SRdata(res_vpa)
head(SRdata)

## ------------------------------------------------------------------------
# SSBとRのデータだけを持っている場合
SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
# 特定の期間のデータだけを使う場合
SRdata0 <- get.SRdata(res_vpa,years=1990:2000) 

## ----SRfit---------------------------------------------------------------
HS.par0 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
HS.par1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,hessian=FALSE)
BH.par0 <- fit.SR(SRdata,SR="BH",method="L2",AR=0,hessian=FALSE)
BH.par1 <- fit.SR(SRdata,SR="BH",method="L2",AR=1,hessian=FALSE)
RI.par0 <- fit.SR(SRdata,SR="RI",method="L2",AR=0,hessian=FALSE)
RI.par1 <- fit.SR(SRdata,SR="RI",method="L2",AR=1,hessian=FALSE)
c(HS.par0$AICc,HS.par1$AICc,BH.par0$AICc,BH.par1$AICc,RI.par0$AICc,RI.par1$AICc)

## ---- fig.cap="図：**観測値（○）に対する再生産関係式．plot=赤がHS，緑と青がBH, RIだが両者はほとんど重なっていて見えない**"----
plot_SRdata(SRdata)
lines(HS.par0$pred,col=2,type="l",lwd=3)
lines(BH.par0$pred,col=3,type="l",lwd=3)    
lines(RI.par0$pred,col=4,type="l",lwd=3)

## ---- eval=FALSE---------------------------------------------------------
#  # install.packages("TMB")　#TMBがインストールされてなければ
#  #library(TMB)
#  #compile("autoregressiveSR2.cpp")
#  #dyn.load(dynlib("autoregressiveSR2"))
#  #HS.par11 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,TMB=TRUE) #marginal likelihood

## ----future.vpa, fig.cap="**図：is.plot=TRUEで表示される図．資源量(Biomass)，親魚資源量(SSB), 漁獲量(Catch)の時系列．決定論的将来予測（Deterministic），平均値（Mean），中央値(Median)，80％信頼区間を表示**"----
fres.HS <- future.vpa(res_vpa,
                      multi=1, # res_vpa$Fc.at.ageに掛ける乗数
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2015:2017, # 生物パラメータの参照年
                      maa.year=2015:2017,
                      M.year=2015:2017,
                      is.plot=FALSE, # 結果をプロットするかどうか
                      seed=1,
                      silent=TRUE,
                      recfunc=HS.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                   rho=HS.par0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                   sd=HS.par0$pars$sd,resid=HS.par0$resid))

## ----future.vpaAR--------------------------------------------------------
fres.HS.AR <- future.vpa(res_vpa,
                      multi=1,
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2015:2017, # 生物パラメータの参照年
                      maa.year=2015:2017,
                      M.year=2015:2017,is.plot=FALSE, 
                      seed=1, silent=TRUE,recfunc=HS.recAR, 
                      # recfuncに対する引数 => 自己相関ありのオプションで計算した結果を入れる
                      rec.arg=list(a=HS.par1$pars$a,b=HS.par1$pars$b,
                                   rho=HS.par1$pars$rho, # 自己相関が高い場合、この値が>0となる
                                   sd=HS.par1$pars$sd,
                                   resid=HS.par1$resid, # 再生産関係における残差の時系列
                                   resid.year=NULL # 近年の残差を何年分平均して将来予測に使うか？NULLの場合は、最後の年の残差を使う
                                   ) 
                      )

## ----future.vpa2, fig.cap="**図：is.plot=TRUEで表示される図．資源量(Biomass)，親魚資源量(SSB), 漁獲量(Catch)の時系列．決定論的将来予測（Deterministic），平均値（Mean），中央値(Median)，80％信頼区間を表示**"----
fres.BH <- future.vpa(res_vpa,
                      multi=1,
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2015:2017, # 生物パラメータの参照年
                      maa.year=2015:2017,
                      M.year=2015:2017,
                      is.plot=FALSE, # 結果をプロットするかどうか
                      seed=1,
                      silent=TRUE,
                      recfunc=BH.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=BH.par0$pars$a,b=BH.par0$pars$b,rho=BH.par0$rho,
                                   sd=BH.par0$pars$sd,resid=BH.par0$resid))

## ------------------------------------------------------------------------
fres.HS2 <- do.call(future.vpa,fres.HS$input)

## ------------------------------------------------------------------------
# 引数をinput.tmpに代入．
input.tmp <- fres.HS2$input
# 引数の一部を変える
input.tmp$multi <- 0.5 # current Fの1/2で漁獲
fres.HS3 <- do.call(future.vpa,input.tmp)

## ---- fig.cap="図：plot.futures関数の結果"-------------------------------
par(mfrow=c(2,2))
#plot.futures(list(fres.HS,fres.HS3),legend.text=c("F=Fcurrent","F=0.5Fcurrent"),target="SSB")
#plot.futures(list(fres.HS,fres.HS3),legend.text=c("F=Fcurrent","F=0.5Fcurrent"),target="Catch")
#plot.futures(list(fres.HS,fres.HS3),legend.text=c("F=Fcurrent","F=0.5Fcurrent"),target="Biomass") 

## ---- fig.cap="図：plot.futures関数の結果"-------------------------------
byear <- 2015:2017 # 生物パラメータを平均する期間を2009年から2011年とする
rps.year <- 2001:2011
fout.rps <- future.vpa(res_vpa,currentF=NULL, multi=1, 
                    nyear=15,start.year=2012,N=10000,ABC.year=2013, 
                    waa.year=byear,maa.year=byear,M.year=byear,
                    rec.new=NULL,is.plot=FALSE,
                    recfunc=RPS.simple.rec,
                    rec.arg=list(rps.year=rps.year,
                      upper.ssb=Inf,bias.corrected=TRUE,rpsmean=FALSE,
                      upper.recruit=Inf))
                    

## ------------------------------------------------------------------------
# 残差リサンプリングによる将来予測
fres.HS4 <- future.vpa(res_vpa,
                          multi=1,
                          nyear=50, # 将来予測の年数
                          start.year=2012, # 将来予測の開始年
                          N=100, # 確率的計算の繰り返し回数
                          ABC.year=2013, # ABCを計算する年
                          waa.year=2015:2017, # 生物パラメータの参照年
                          maa.year=2015:2017,
                          M.year=2015:2017,
                          is.plot=FALSE, # 結果をプロットするかどうか
                          seed=1,
                          recfunc=HS.recAR, # 再生産関係の関数（HS.rec=Hockey-stick)                                
                          rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                       rho=HS.par0$pars$rho,
                                       sd=HS.par0$pars$sd,bias.correction=TRUE,
                                       resample=TRUE,resid=HS.par0$resid))

## ----eval=FALSE----------------------------------------------------------
#  par(mfrow=c(2,2))
#  plot(fres.HS$vssb[,-1],fres.HS$naa[1,,-1],xlab="SSB",ylab="Recruits")
#  plot(fres.HS4$vssb[,-1],fres.HS4$naa[1,,-1],xlab="SSB",ylab="Recruits")
#  #plot.futures(list(fres.HS,fres.HS4)) # 両者の比較

## ---- fig.cap="Frecオプションを使った場合は、結果の図に目的とする年・資源量のところに赤線が入ります。これが将来予測の結果と一致しているか確かめてください。もし一致していない場合、multi（初期値）かFrecのオプションのFrangeを指定してやり直してください"----
# たとえば現状の資源量に維持するシナリオ
fres.currentSSB <- future.vpa(res_vpa,
                      multi=0.8,
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2015:2017, # 生物パラメータの参照年
                      maa.year=2015:2017,
                      M.year=2015:2017,seed=1,
                      is.plot=FALSE, # 結果をプロットするかどうか
                      Frec=list(stochastic=TRUE,future.year=2023,Blimit=rev(colSums(res_vpa$ssb))[1],scenario="blimit",target.probs=50),
                      recfunc=HS.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                   rho=HS.par0$pars$rho,                                    
                                   sd=HS.par0$pars$sd,bias.corrected=TRUE))

## ------------------------------------------------------------------------
lm.res <- plot_waa(res_vpa) # weight at ageが資源尾数の関数になっているかどうか，確認してみる．この例の場合は特に有意な関係はない
# lm.resの中に回帰した結果が年齢分だけ入っています
fres.HS6 <- fres.HS
fres.HS6$input$waa.fun <- TRUE
fres.HS6$input$N <- 1000
fres.HS6 <- do.call(future.vpa, fres.HS6$input)

## ----options-------------------------------------------------------------
fres.HS5 <- future.vpa(res_vpa,
                       multi=1,
                       nyear=50, # 将来予測の年数
                       start.year=2012, # 将来予測の開始年
                       N=100, # 確率的計算の繰り返し回数
                       ABC.year=2013, # ABCを計算する年
                       waa.year=2015:2017, # 生物パラメータの参照年
                       maa.year=2015:2017,
                       M.year=2015:2017,is.plot=FALSE,
                       recfunc=HS.recAR, 
                       rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                       rho=HS.par0$pars$rho,
                                       sd=HS.par0$pars$sd,bias.correction=TRUE,
                                    resample=TRUE,resid=HS.par0$resid),
                       rec.new=list(year=2012,rec=100),
                       pre.catch=list(year=2013,wcatch=100)
                       )

## ----ref.F, fig.cap="図：YPR, SPR曲線"-----------------------------------
byear <- 2015:2017 # 生物パラメータを平均する期間を2009年から2011年とする
res_bref <- ref.F(res_vpa, # VPAの計算結果
                  waa.year=byear, maa.year=byear, M.year=byear, 
                  rps.year=2000:2017, # Fmedを計算するときに用いるRPSの範囲
                  max.age=Inf, # SPR計算で仮定する年齢の最大値 
                  pSPR=c(10,20,30,35,40)) # F_%SPRを計算するときに，何パーセントのSPRを計算するか

## ----ref.F2--------------------------------------------------------------
res_bref$summary

## ----ref.F3, fig.cap="図：YPR, SPR曲線 (x軸などを変更した場合)"----------
# 横軸や縦線で示す管理基準値を調整する場合、plot_Fref関数を使う
# x.labelは res_bref$summaryの行名、vline.textは res_bref$summaryの列名前に対応させて指定する
plot_Fref(res_bref,xlabel="Fref/Fcur", vline.text=c("FpSPR.20.SPR","FpSPR.30.SPR","FpSPR.40.SPR"))

