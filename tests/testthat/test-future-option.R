library(frasyr)

context("future MSE option")

## test_that("MSE option check",{
##     data(res_vpa)

##     SRdata <- get.SRdata(res_vpa, years=1988:2016) 
##     head(SRdata)

##     if(0){
##     SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), L.type = c("L1", "L2"))
##     SR.list <- list()
##     for (i in 1:nrow(SRmodel.list)) {
##         SR.list[[i]] <- fit.SR(SRdata, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i], 
##                                AR = SRmodel.list$AR.type[i], hessian = FALSE)
##     }
##     SRmodel.list$AICc <- sapply(SR.list, function(x) x$AICc)
##     SRmodel.list$delta.AIC <- SRmodel.list$AICc - min(SRmodel.list$AICc)
##     SR.list <- SR.list[order(SRmodel.list$AICc)]  # AICの小さい順に並べたもの
##     (SRmodel.list <- SRmodel.list[order(SRmodel.list$AICc), ]) # 結果
##     SRmodel.base <- SR.list[[1]] # AIC最小モデルを今後使っていく
##     SRmodel.R1 <- SR.list[[12]] # 候補となる別の加入シナリオ
##     future.Fcurrent <- future.vpa(res_vpa,
##                                   multi=1,
##                                   nyear=50, # 将来予測の年数
##                                   start.year=2018, # 将来予測の開始年
##                                   N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
##                       ABC.year=2019, # ABCを計算する年
##                       waa.year=2015:2017, # 生物パラメータの参照年
##                       maa.year=2015:2017,
##                       M.year=2015:2017,
##                       is.plot=TRUE, # 結果をプロットするかどうか
##                       seed=1,
##                       silent=FALSE,
##                       recfunc=HS.recAR, # 再生産関係の関数
##                       # recfuncに対する引数
##                       rec.arg=list(a=SRmodel.base$pars$a,b=SRmodel.base$pars$b,
##                                    rho=SRmodel.base$pars$rho, # ここではrho=0なので指定しなくてもOK
##                                    sd=SRmodel.base$pars$sd,resid=SRmodel.base$resid))
## # MSY管理基準値の計算; base caseのシナリオをもとにした管理基準値
## MSY.base <- est.MSY(res_vpa, # VPAの計算結果
##                     future.Fcurrent$input, # 将来予測で使用した引数
##                  resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
##                  N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
##                  calc.yieldcurve=TRUE,
##                  PGY=c(0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
##                  onlylower.pgy=TRUE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
##                  B0percent=NULL,
##                  Bempirical=NULL
##                  ) # 計算したいB0%レベル
## (refs.all <- MSY.base$summary_tb)

## # refs.allの中からRP.definitionで指定された行だけを抜き出す
## (refs.base <- refs.all %>%
##     dplyr::filter(!is.na(RP.definition)) %>% # RP.definitionがNAでないものを抽出
##     arrange(desc(SSB)) %>% # SSBを大きい順に並び替え
##     select(RP.definition,RP_name,SSB,SSB2SSB0,Catch,Catch.CV,U,Fref2Fcurrent)) #　列を並び替え

## # HCRによる将来予測
## input.abc <- future.Fcurrent$input # Fcurrentにおける将来予測の引数をベースに将来予測します
## input.abc$multi <- derive_RP_value(refs.base,"Btarget0")$Fref2Fcurrent # currentFへの乗数を"Btarget0"で指定した値に
## input.abc$silent <- TRUE
## input.abc$HCR <- list(Blim=derive_RP_value(refs.base,"Blimit0")$SSB,
##                       Bban=derive_RP_value(refs.base,"Bban0")$SSB,
##                       beta=0.8,year.lag=0) # BlimitはBlimit0, BbanはBban0の値
## input.abc$N <- 1000
## future.default <- do.call(future.vpa,input.abc) # デフォルトルールの結果→図示などに使う

## input.mse <- input.abc
## input.mse$N <- 100
## input.mse$use.MSE <- TRUE # MSE仕様での将来予測
## input.mse$is.plot <- FALSE
## future.mse <- do.call(future.vpa,input.mse)

## # use.MSEオプションで、ARありバージョンには十分対応していない→今後の課題
## input.mse_R1 <- input.mse
## input.mse_R1$N <- 100
## input.mse_R1$MSE.options$recfunc <- future.Fcurrent_R1$recfunc
## future.Fcurrent_R1$rec.arg$sd2 <- future.Fcurrent_R1$rec.arg$sd2[1:input.mse_R1$N]
## input.mse_R1$MSE.options$rec.arg <- future.Fcurrent_R1$rec.arg
## future.mse_R1 <- do.call(future.vpa,input.mse_R1)  

##     all.table <- purrr::map_dfr(list(future.mse,future.default,future.mse_R1),convert_future_table,.id="scenario")
## all.table %>% dplyr::filter(stat=="catch",year<2025,year>2018) %>%
##     ggplot() +
##     geom_boxplot(aes(x=factor(year),y=value,fill=scenario))
   
##     }
## }

