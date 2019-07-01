## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.widths=100,
  fig.heights=100
)

## ----data-read-----------------------------------------------------------

library(frasyr)
# データの読み込み
caa <- read.csv("../data-raw/caa_pma.csv",row.names=1)
waa <- read.csv("../data-raw/waa_pma.csv",row.names=1)
maa <- read.csv("../data-raw/maa_pma.csv",row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
names(dat)


## ----vpa-----------------------------------------------------------------
# VPAによる資源量推定
res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

## ------------------------------------------------------------------------
res.pma$Fc.at.age # 将来予測やMSY計算で使うcurrent F (fc.yearのオプションでいつのFの平均かが指定される)
plot(res.pma$Fc.at.age,type="b",xlab="Age",ylab="F",ylim=c(0,max(res.pma$Fc.at.age)))

## ------------------------------------------------------------------------
#out.vpa(res.pma) # vpa.csvというファイルが作成されます。VPAの結果のグラフ出力となるvpa.pdfも出力されます。
res.pma2 <- read.vpa("../data-raw/vpa.csv") # vpa.csvを編集後、read.vpa関数で読み込みます

## ----ref.F, fig.cap="図：YPR, SPR曲線"-----------------------------------
byear <- 2009:2011 # 生物パラメータを平均する期間を2009年から2011年とする
rres.pma <- ref.F(res.pma, # VPAの計算結果
                  waa.year=byear, maa.year=byear, M.year=byear, # weight at age, maturity at age, Mは2009から2011年までの平均とする
                  rps.year=2000:2011, # Fmedを計算するときに用いるRPSの範囲
                  max.age=Inf, # SPR計算で仮定する年齢の最大値 
                  pSPR=c(10,20,30,35,40), # F_%SPRを計算するときに，何パーセントのSPRを計算するか
                  Fspr.init=1)

## ----ref.F2--------------------------------------------------------------
rres.pma$summary

## ----ref.F3, fig.cap="図：YPR, SPR曲線 (x軸などを変更した場合)"----------
# 横軸や縦線で示す管理基準値を調整する場合、plot_Fref関数を使う
# x.labelは res.pma$summaryの行名、vline.textは res.pma$summaryの列名前に対応させて指定する
plot_Fref(rres.pma,xlabel="Fref/Fcur", vline.text=c("FpSPR.20.SPR","FpSPR.30.SPR","FpSPR.40.SPR"))

## ----SRdata--------------------------------------------------------------
# VPA結果を使って再生産データを作る
SRdata <- get.SRdata(res.pma)
head(SRdata)

## ------------------------------------------------------------------------
# SSBとRのデータだけを持っている場合
SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
# 特定の期間のデータだけを使う場合
SRdata0 <- get.SRdata(res.pma,years=1990:2000) 

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
fres.HS <- future.vpa(res.pma,
                      multi=1, # res.pma$Fc.at.ageに掛ける乗数
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2009:2011, # 生物パラメータの参照年
                      maa.year=2009:2011,
                      M.year=2009:2011,
                      is.plot=FALSE, # 結果をプロットするかどうか
                      seed=1,
                      silent=TRUE,
                      recfunc=HS.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                   rho=HS.par0$pars$rho, # ここではrho=0なので指定しなくてもOK
                                   sd=HS.par0$pars$sd,resid=HS.par0$resid))

## ----future.vpaAR--------------------------------------------------------
fres.HS.AR <- future.vpa(res.pma,
                      multi=1,
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2009:2011, # 生物パラメータの参照年
                      maa.year=2009:2011,
                      M.year=2009:2011,is.plot=FALSE, 
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
fres.BH <- future.vpa(res.pma,
                      multi=1,
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2009:2011, # 生物パラメータの参照年
                      maa.year=2009:2011,
                      M.year=2009:2011,
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
rps.year <- 2001:2011
fout.rps <- future.vpa(res.pma,currentF=NULL, multi=1, 
                    nyear=15,start.year=2012,N=10000,ABC.year=2013, 
                    waa.year=byear,maa.year=byear,M.year=byear,
                    rec.new=NULL,is.plot=FALSE,
                    recfunc=RPS.simple.rec,
                    rec.arg=list(rps.year=rps.year,
                      upper.ssb=Inf,bias.corrected=TRUE,rpsmean=FALSE,
                      upper.recruit=Inf))
                    

## ---- fig.cap="Frecオプションを使った場合は、結果の図に目的とする年・資源量のところに赤線が入ります。これが将来予測の結果と一致しているか確かめてください。もし一致していない場合、multi（初期値）かFrecのオプションのFrangeを指定してやり直してください"----
# たとえば現状の資源量に維持するシナリオ
fres.currentSSB <- future.vpa(res.pma,
                      multi=0.8,
                      nyear=50, # 将来予測の年数
                      start.year=2012, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数
                      ABC.year=2013, # ABCを計算する年
                      waa.year=2009:2011, # 生物パラメータの参照年
                      maa.year=2009:2011,
                      M.year=2009:2011,seed=1,
                      is.plot=FALSE, # 結果をプロットするかどうか
                      Frec=list(stochastic=TRUE,future.year=2023,Blimit=rev(colSums(res.pma$ssb))[1],scenario="blimit",target.probs=50),
                      recfunc=HS.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                   rho=HS.par0$pars$rho,                                    
                                   sd=HS.par0$pars$sd,bias.corrected=TRUE))

## ------------------------------------------------------------------------
# 残差リサンプリングによる将来予測
fres.HS4 <- future.vpa(res.pma,
                          multi=1,
                          nyear=50, # 将来予測の年数
                          start.year=2012, # 将来予測の開始年
                          N=100, # 確率的計算の繰り返し回数
                          ABC.year=2013, # ABCを計算する年
                          waa.year=2009:2011, # 生物パラメータの参照年
                          maa.year=2009:2011,
                          M.year=2009:2011,
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

## ------------------------------------------------------------------------
lm.res <- plot_waa(res.pma) # weight at ageが資源尾数の関数になっているかどうか，確認してみる．この例の場合は特に有意な関係はない
# lm.resの中に回帰した結果が年齢分だけ入っています
fres.HS6 <- fres.HS
fres.HS6$input$waa.fun <- TRUE
fres.HS6$input$N <- 1000
fres.HS6 <- do.call(future.vpa, fres.HS6$input)

## ----options-------------------------------------------------------------
# 残差リサンプリングによる将来予測
fres.HS5 <- future.vpa(res.pma,
                       multi=1,
                       nyear=50, # 将来予測の年数
                       start.year=2012, # 将来予測の開始年
                       N=100, # 確率的計算の繰り返し回数
                       ABC.year=2013, # ABCを計算する年
                       waa.year=2009:2011, # 生物パラメータの参照年
                       maa.year=2009:2011,
                       M.year=2009:2011,is.plot=FALSE,
                       recfunc=HS.recAR, 
                       rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
                                       rho=HS.par0$pars$rho,
                                       sd=HS.par0$pars$sd,bias.correction=TRUE,
                                    resample=TRUE,resid=HS.par0$resid),
                       rec.new=list(year=2012,rec=100),
                       pre.catch=list(year=2013,wcatch=100)
                       )

## ----msy, fig.cap="**図：est.MSYのis.plot=TRUEで計算完了時に表示される図．Fの強さに対する平衡状態の親魚資源量（左）と漁獲量（右）．推定された管理基準値も表示．**", fig.height=5----

# MSY管理基準値の計算
MSY.HS <- est.MSY(res.pma, # VPAの計算結果
                 fres.HS$input, # 将来予測で使用した引数
#                 nyear=NULL, # 何年計算するかは、指定しなければ関数内部で世代時間の20倍の年数を計算し、それを平衡状態とする
                 resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                 N=100, # 将来予測の年数，繰り返し回数
                 PGY=c(0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                 onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                 B0percent=c(0.3,0.4)) # 計算したいB0%レベル

## ----summary-------------------------------------------------------------
# 結果の表示(tibbleという形式で表示され、最初の10行以外は省略されます)
(refs.all <- MSY.HS$summary_tb)

# 全データをじっくり見たい場合
View(refs.all)

# のちの使用のために、Bmsy, Blimit, Bban, Fmsyを定義しておく
refs <- list(BmsyAR=as.numeric(MSY.HS$summaryAR$SSB[1]),
             BlimAR=as.numeric(MSY.HS$summaryAR$SSB[6]),
             BbanAR=as.numeric(MSY.HS$summaryAR$SSB[8]),
             Bmsy=as.numeric(MSY.HS$summary$SSB[1]),
             Blim=as.numeric(MSY.HS$summary$SSB[6]),
             Bban=as.numeric(MSY.HS$summary$SSB[8]),
             Fmsy=as.numeric(MSY.HS$summary$"Fref/Fcur"[1]),
             MSY=as.numeric(MSY.HS$summary$Catch[1]),
             Umsy=as.numeric(MSY.HS$summary$Catch[1])/as.numeric(MSY.HS$all.stat$biom.mean[1]))

# また、各管理基準値に対して以下のようなラベルをつける
# どの管理基準値をどのように定義するか、ここで指定します
refs.all$RP.definition <- NA 
refs.all$RP.definition[refs.all$RP_name=="MSY" & refs.all$AR==FALSE] <- "Btarget0"  # RP_nameがMSYでARがなしのものをBtargetとする
refs.all$RP.definition[refs.all$RP_name=="B0-20%" & refs.all$AR==FALSE] <- "Btarget1"  # たとえばBtargetの代替値をいちおう示す場合
refs.all$RP.definition[refs.all$RP_name=="PGY_0.95_lower" & refs.all$AR==FALSE] <- "Btarget2" 
refs.all$RP.definition[refs.all$RP_name=="PGY_0.9_lower" & refs.all$AR==FALSE] <- "Blow0"
refs.all$RP.definition[refs.all$RP_name=="PGY_0.6_lower" & refs.all$AR==FALSE] <- "Blimit0"
refs.all$RP.definition[refs.all$RP_name=="PGY_0.1_lower" & refs.all$AR==FALSE] <- "Bban0"
refs.all$RP.definition[refs.all$RP_name=="Ben-19431" & refs.all$AR==FALSE] <- "Bcurrent"
refs.all$RP.definition[refs.all$RP_name=="Ben-63967" & refs.all$AR==FALSE] <- "Bmax"
refs.all$RP.definition[refs.all$RP_name=="Ben-24000" & refs.all$AR==FALSE] <- "Blimit1"
refs.all$RP.definition[refs.all$RP_name=="Ben-51882" & refs.all$AR==FALSE] <- "B_HS"

# 定義した結果を見る
refs.all %>% select(RP_name,RP.definition)

# refs.allの中からRP.definitionで指定された行だけを抜き出す
(refs.base <- refs.all %>%
     dplyr::filter(!is.na(RP.definition)) %>% # RP.definitionがNAでないものを抽出
     arrange(desc(SSB)) %>% # SSBを大きい順に並び替え
     select(RP.definition,RP_name,SSB,Catch,U,Fref2Fcurrent)) #　列を並び替え


## ----abc-----------------------------------------------------------------
input.abc <- MSY.HS$input$msy # MSY計算で使った引数を使う
input.abc$N <- 1000 # 実際に計算するときは10000以上を使ってください
input.abc$HCR <- list(Blim=refs$Blim,
                      Bban=refs$Bban,
                      beta=0.8) # 算定ルールのデフォルトは0.8
input.abc$nyear <- 30 # ABC計算時には長期間計算する必要はない
input.abc$ABC.year <- 2013 # ここでABC.yearを設定しなおしてください
input.abc$is.plot <- FALSE
fres.abc1 <- do.call(future.vpa,input.abc)

par(mfrow=c(1,1))
hist(fres.abc1$ABC,main="distribution of ABC") # ABCの分布
ABC <- mean(fres.abc1$ABC) # 平均値をABCとする

## SSBの将来予測結果
par(mfrow=c(1,1))
#plot.future(fres.abc1,what=c(FALSE,TRUE,FALSE),is.legend=TRUE,lwd=2,
#            col="darkblue",N=5,label=rep(NA,3))
#draw.refline(cbind(unlist(refs[c(1,1,2,3)+3]),unlist(refs[c(1,1,2,3)])),horiz=TRUE,lwd=1,scale=1)

## 漁獲量の将来予測結果
par(mfrow=c(1,1))
#plot.future(fres.abc1,what=c(FALSE,FALSE,TRUE),is.legend=TRUE,lwd=2,
#            col="darkblue",N=5,label=rep(NA,3))
#points(fres.abc1$input$ABC.year,ABC,pch=20,col=2,cex=3)
#text(fres.abc1$input$ABC.year+1,ABC,"ABC",col=2)

## 実際に、どんなFが将来予測で使われているか
boxplot(t(fres.abc1$faa[1,,]/fres.abc1$faa[1,1,]),ylab="multiplier to current F")

## ----HCR-----------------------------------------------------------------
# どんなHCRなのか書いてみる
ssb.abc <- mean(fres.abc1$vssb[rownames(fres.abc1$vssb)%in%fres.abc1$input$ABC.year,]) # ABC計算年のssbをとる
#plot.HCR(beta=beta,bban=refs$Bban,blimit=refs$Blim,btarget=refs$Bmsy,lwd=2,
#         xlim=c(0,refs$Bmsy*2),ssb.cur=ssb.abc,Fmsy=refs$Fmsy,yscale=0.7,scale=1000)

## ----probability---------------------------------------------------------
# 将来の親魚資源量がBMSYやBlimitを上回る確率の表示
plot(apply(fres.abc1$vssb>refs$Bmsy,1,mean)*100,type="b",ylab="Probability",ylim=c(0,100))
points(apply(fres.abc1$vssb>refs$BmsyAR,1,mean)*100,pch=2,type="b")
points(apply(fres.abc1$vssb>refs$Blim,1,mean)*100,pch=1,col=2,type="b")
points(apply(fres.abc1$vssb>refs$BlimAR,1,mean)*100,pch=2,col=2,type="b")
abline(h=c(50,90),col=c(1,2))
legend("bottomright",col=c(1,1,2,2),title="Probs",pch=c(1,2,1,2),legend=c(">Btarget_Eq",">Btarget_AR",">Blimit_Eq",">Blimit_AR"))

## ----ref.label='data-read', eval=FALSE-----------------------------------
#  
#  library(frasyr)
#  # データの読み込み
#  caa <- read.csv("../data-raw/caa_pma.csv",row.names=1)
#  waa <- read.csv("../data-raw/waa_pma.csv",row.names=1)
#  maa <- read.csv("../data-raw/maa_pma.csv",row.names=1)
#  dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
#  names(dat)
#  

## ----ref.label='vpa',  eval=FALSE----------------------------------------
#  # VPAによる資源量推定
#  res.pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
#                 term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

## ----ref.label='SRdata', eval=FALSE--------------------------------------
#  # VPA結果を使って再生産データを作る
#  SRdata <- get.SRdata(res.pma)
#  head(SRdata)

## ----ref.label='SRfit', eval=FALSE---------------------------------------
#  HS.par0 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
#  HS.par1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1,hessian=FALSE)
#  BH.par0 <- fit.SR(SRdata,SR="BH",method="L2",AR=0,hessian=FALSE)
#  BH.par1 <- fit.SR(SRdata,SR="BH",method="L2",AR=1,hessian=FALSE)
#  RI.par0 <- fit.SR(SRdata,SR="RI",method="L2",AR=0,hessian=FALSE)
#  RI.par1 <- fit.SR(SRdata,SR="RI",method="L2",AR=1,hessian=FALSE)
#  c(HS.par0$AICc,HS.par1$AICc,BH.par0$AICc,BH.par1$AICc,RI.par0$AICc,RI.par1$AICc)

## ----ref.label='future.vpa', eval=FALSE----------------------------------
#  fres.HS <- future.vpa(res.pma,
#                        multi=1, # res.pma$Fc.at.ageに掛ける乗数
#                        nyear=50, # 将来予測の年数
#                        start.year=2012, # 将来予測の開始年
#                        N=100, # 確率的計算の繰り返し回数
#                        ABC.year=2013, # ABCを計算する年
#                        waa.year=2009:2011, # 生物パラメータの参照年
#                        maa.year=2009:2011,
#                        M.year=2009:2011,
#                        is.plot=FALSE, # 結果をプロットするかどうか
#                        seed=1,
#                        silent=TRUE,
#                        recfunc=HS.recAR, # 再生産関係の関数
#                        # recfuncに対する引数
#                        rec.arg=list(a=HS.par0$pars$a,b=HS.par0$pars$b,
#                                     rho=HS.par0$pars$rho, # ここではrho=0なので指定しなくてもOK
#                                     sd=HS.par0$pars$sd,resid=HS.par0$resid))

## ----ref.label='msy', eval=FALSE-----------------------------------------
#  
#  # MSY管理基準値の計算
#  MSY.HS <- est.MSY(res.pma, # VPAの計算結果
#                   fres.HS$input, # 将来予測で使用した引数
#  #                 nyear=NULL, # 何年計算するかは、指定しなければ関数内部で世代時間の20倍の年数を計算し、それを平衡状態とする
#                   resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
#                   N=100, # 将来予測の年数，繰り返し回数
#                   PGY=c(0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
#                   onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
#                   B0percent=c(0.3,0.4)) # 計算したいB0%レベル

## ----ref.label='summary', eval=FALSE-------------------------------------
#  # 結果の表示(tibbleという形式で表示され、最初の10行以外は省略されます)
#  (refs.all <- MSY.HS$summary_tb)
#  
#  # 全データをじっくり見たい場合
#  View(refs.all)
#  
#  # のちの使用のために、Bmsy, Blimit, Bban, Fmsyを定義しておく
#  refs <- list(BmsyAR=as.numeric(MSY.HS$summaryAR$SSB[1]),
#               BlimAR=as.numeric(MSY.HS$summaryAR$SSB[6]),
#               BbanAR=as.numeric(MSY.HS$summaryAR$SSB[8]),
#               Bmsy=as.numeric(MSY.HS$summary$SSB[1]),
#               Blim=as.numeric(MSY.HS$summary$SSB[6]),
#               Bban=as.numeric(MSY.HS$summary$SSB[8]),
#               Fmsy=as.numeric(MSY.HS$summary$"Fref/Fcur"[1]),
#               MSY=as.numeric(MSY.HS$summary$Catch[1]),
#               Umsy=as.numeric(MSY.HS$summary$Catch[1])/as.numeric(MSY.HS$all.stat$biom.mean[1]))
#  
#  # また、各管理基準値に対して以下のようなラベルをつける
#  # どの管理基準値をどのように定義するか、ここで指定します
#  refs.all$RP.definition <- NA
#  refs.all$RP.definition[refs.all$RP_name=="MSY" & refs.all$AR==FALSE] <- "Btarget0"  # RP_nameがMSYでARがなしのものをBtargetとする
#  refs.all$RP.definition[refs.all$RP_name=="B0-20%" & refs.all$AR==FALSE] <- "Btarget1"  # たとえばBtargetの代替値をいちおう示す場合
#  refs.all$RP.definition[refs.all$RP_name=="PGY_0.95_lower" & refs.all$AR==FALSE] <- "Btarget2"
#  refs.all$RP.definition[refs.all$RP_name=="PGY_0.9_lower" & refs.all$AR==FALSE] <- "Blow0"
#  refs.all$RP.definition[refs.all$RP_name=="PGY_0.6_lower" & refs.all$AR==FALSE] <- "Blimit0"
#  refs.all$RP.definition[refs.all$RP_name=="PGY_0.1_lower" & refs.all$AR==FALSE] <- "Bban0"
#  refs.all$RP.definition[refs.all$RP_name=="Ben-19431" & refs.all$AR==FALSE] <- "Bcurrent"
#  refs.all$RP.definition[refs.all$RP_name=="Ben-63967" & refs.all$AR==FALSE] <- "Bmax"
#  refs.all$RP.definition[refs.all$RP_name=="Ben-24000" & refs.all$AR==FALSE] <- "Blimit1"
#  refs.all$RP.definition[refs.all$RP_name=="Ben-51882" & refs.all$AR==FALSE] <- "B_HS"
#  
#  # 定義した結果を見る
#  refs.all %>% select(RP_name,RP.definition)
#  
#  # refs.allの中からRP.definitionで指定された行だけを抜き出す
#  (refs.base <- refs.all %>%
#       dplyr::filter(!is.na(RP.definition)) %>% # RP.definitionがNAでないものを抽出
#       arrange(desc(SSB)) %>% # SSBを大きい順に並び替え
#       select(RP.definition,RP_name,SSB,Catch,U,Fref2Fcurrent)) #　列を並び替え
#  

## ----ref.label='abc', eval=FALSE-----------------------------------------
#  input.abc <- MSY.HS$input$msy # MSY計算で使った引数を使う
#  input.abc$N <- 1000 # 実際に計算するときは10000以上を使ってください
#  input.abc$HCR <- list(Blim=refs$Blim,
#                        Bban=refs$Bban,
#                        beta=0.8) # 算定ルールのデフォルトは0.8
#  input.abc$nyear <- 30 # ABC計算時には長期間計算する必要はない
#  input.abc$ABC.year <- 2013 # ここでABC.yearを設定しなおしてください
#  input.abc$is.plot <- FALSE
#  fres.abc1 <- do.call(future.vpa,input.abc)
#  
#  par(mfrow=c(1,1))
#  hist(fres.abc1$ABC,main="distribution of ABC") # ABCの分布
#  ABC <- mean(fres.abc1$ABC) # 平均値をABCとする
#  
#  ## SSBの将来予測結果
#  par(mfrow=c(1,1))
#  #plot.future(fres.abc1,what=c(FALSE,TRUE,FALSE),is.legend=TRUE,lwd=2,
#  #            col="darkblue",N=5,label=rep(NA,3))
#  #draw.refline(cbind(unlist(refs[c(1,1,2,3)+3]),unlist(refs[c(1,1,2,3)])),horiz=TRUE,lwd=1,scale=1)
#  
#  ## 漁獲量の将来予測結果
#  par(mfrow=c(1,1))
#  #plot.future(fres.abc1,what=c(FALSE,FALSE,TRUE),is.legend=TRUE,lwd=2,
#  #            col="darkblue",N=5,label=rep(NA,3))
#  #points(fres.abc1$input$ABC.year,ABC,pch=20,col=2,cex=3)
#  #text(fres.abc1$input$ABC.year+1,ABC,"ABC",col=2)
#  
#  ## 実際に、どんなFが将来予測で使われているか
#  boxplot(t(fres.abc1$faa[1,,]/fres.abc1$faa[1,1,]),ylab="multiplier to current F")

## ----ref.label='HCR', eval=FALSE-----------------------------------------
#  # どんなHCRなのか書いてみる
#  ssb.abc <- mean(fres.abc1$vssb[rownames(fres.abc1$vssb)%in%fres.abc1$input$ABC.year,]) # ABC計算年のssbをとる
#  #plot.HCR(beta=beta,bban=refs$Bban,blimit=refs$Blim,btarget=refs$Bmsy,lwd=2,
#  #         xlim=c(0,refs$Bmsy*2),ssb.cur=ssb.abc,Fmsy=refs$Fmsy,yscale=0.7,scale=1000)

## ----ref.label='probability', eval=FALSE---------------------------------
#  # 将来の親魚資源量がBMSYやBlimitを上回る確率の表示
#  plot(apply(fres.abc1$vssb>refs$Bmsy,1,mean)*100,type="b",ylab="Probability",ylim=c(0,100))
#  points(apply(fres.abc1$vssb>refs$BmsyAR,1,mean)*100,pch=2,type="b")
#  points(apply(fres.abc1$vssb>refs$Blim,1,mean)*100,pch=1,col=2,type="b")
#  points(apply(fres.abc1$vssb>refs$BlimAR,1,mean)*100,pch=2,col=2,type="b")
#  abline(h=c(50,90),col=c(1,2))
#  legend("bottomright",col=c(1,1,2,2),title="Probs",pch=c(1,2,1,2),legend=c(">Btarget_Eq",">Btarget_AR",">Blimit_Eq",">Blimit_AR"))

