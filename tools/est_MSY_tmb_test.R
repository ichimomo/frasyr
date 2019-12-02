library(frasyr)
data(res_vpa)
SRdata <- get.SRdata(res_vpa)

## モデルのフィット(網羅的に試しています)
# 網羅的なパラメータ設定
SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), L.type = c("L1", "L2"))
SR.list <- list()
for (i in 1:nrow(SRmodel.list)) {
    SR.list[[i]] <- fit.SR(SRdata, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i], 
        AR = SRmodel.list$AR.type[i], hessian = FALSE)
}

SRmodel.list$AICc <- sapply(SR.list, function(x) x$AICc)
SRmodel.list$delta.AIC <- SRmodel.list$AICc - min(SRmodel.list$AICc)
SR.list <- SR.list[order(SRmodel.list$AICc)]  # AICの小さい順に並べたもの
(SRmodel.list <- SRmodel.list[order(SRmodel.list$AICc), ]) # 結果   

plot_SRdata(SRdata)
for(i in 1:nrow(SRmodel.list)){
  lines(SR.list[[i]]$pred,col=i)
}

## 将来予測の実施
SRmodel.base <- SR.list[[1]] # AIC最小モデルを今後使っていく
res_future_Fcurrent <- future.vpa(res_vpa,
                      multi=1,
                      nyear=100, # 将来予測の年数
                      start.year=2017, # 将来予測の開始年
                      N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                      ABC.year=2016, # ABCを計算する年
                      waa.year=2015:2017, # 生物パラメータの参照年
                      maa.year=2015:2017,
                      M.year=2015:2017,
                      is.plot=TRUE, # 結果をプロットするかどうか
                      seed=1,
                      silent=FALSE,
                      recfunc=HS.recAR, # 再生産関係の関数
                      # recfuncに対する引数
                      rec.arg=list(a=SRmodel.base$pars$a,b=SRmodel.base$pars$b,
                                   rho=SRmodel.base$pars$rho, # ここではrho=0なので指定しなくてもOK
                                   sd=SRmodel.base$pars$sd,
                                   resid=SRmodel.base$resid))

## MSY管理基準値の計算
R_time <- system.time(res_MSY <- est.MSY(res_vpa, # VPAの計算結果
                 res_future_Fcurrent$input, # 将来予測で使用した引数
                 resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                 N=1000, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                 calc.yieldcurve=TRUE,
                 nyear=100,
                 PGY=0.6,#NULL,#c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                 onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                 B0percent=NULL,#c(0.2,0.3,0.4),
                 Bempirical=NULL,#c(round(tail(colSums(res_vpa$ssb),n=1)),
#                              round(max(colSums(res_vpa$ssb))),
#                              24000, # 現行Blimit
#                              SRmodel.base$pars$b) # HSの折れ点
                 )) # 計算したいB0%レベル
            
#------------ TMB

devtools::load_all()

res1 <- tmb_future(res_vpa, future_initial_year_name=2017,nsim=1000,nyear=30)

res2 <- tmb_future(res_vpa, future_initial_year_name=2000,nsim=1000,nyear=60,
                   start_F_year_name=2001, start_random_rec_year_name=2017)

res3 <- tmb_future(res_vpa, future_initial_year_name=1995,
                   start_F_year_name=2017,nsim=1000,nyear=60,compile=FALSE,
                   start_random_rec_year_name=2017)

plot(rownames(res1$ssb), rowMeans(res1$ssb),type="l")
points(rowMeans(res2$ssb),type="l",col=2)
points(rownames(res3$ssb), rowMeans(res3$ssb),type="l",col=3)
