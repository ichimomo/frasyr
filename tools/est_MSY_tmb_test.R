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
            

###res_vpaの構造を保ったまま将来予測する

#### --- age-structured
nsim  <- 10
nage  <- nrow(res_vpa$naa)
nyear <- 30 # number of future projection
vpa_nyear <- ncol(res_vpa$naa)
future_initial_year <- vpa_nyear ## ここはvpa_nyear以下、任意
total_nyear <- future_initial_year + nyear
allyear_name <- min(as.numeric(colnames(res_vpa$naa)))+c(0:(total_nyear-1))
allyear_label <- c(rep("VPA",future_initial_year),rep("future",nyear))
print(tibble(allyear_name, allyear_label))

recruit_age <- min(as.numeric(rownames(res_vpa$naa)))

waa_mat <- M_mat <- maa_mat <- naa_mat <- faa_mat <- array(0, dim=c(nage, total_nyear, nsim)) 
resid_mat <- array(0, dim=c(total_nyear, nsim))

# input VPA part parameter (VPA設定を使うか将来予測設定を使うか、２つのオプションあり。ここではデータがある年のものはそれを使うやりかた)
waa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa)
maa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$maa)
M_mat  [,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$M)
naa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$naa)
faa_mat  [,1:vpa_nyear,] <- as.matrix(res_vpa$faa)
resid_mat[1:future_initial_year,] <- SRmodel.base$resid[1:future_initial_year]

# assume future parameter（いろいろオプションがあるが、最終年のものを使う）
waa_mat[,(vpa_nyear+1):total_nyear,] <- waa_mat[,vpa_nyear,]
maa_mat[,(vpa_nyear+1):total_nyear,] <- maa_mat[,vpa_nyear,]
M_mat  [,(vpa_nyear+1):total_nyear,] <- M_mat  [,vpa_nyear,]
faa_mat[,(vpa_nyear+1):total_nyear,] <- faa_mat[,vpa_nyear,]
resid_mat[(future_initial_year+1):total_nyear,] <-
    rnorm(nsim*nyear, mean=-0.5 * (SRmodel.base$pars$sd)^2, sd=SRmodel.base$pars$sd)

tmb_data_msy <- list(naa_mat=naa_mat,
#                     deviance_init=0,
                     SR = 0,
                     rec_par_a = SRmodel.base$pars$a,
                     rec_par_b = SRmodel.base$pars$b,
                     rec_par_rho = 0,
                     bias_corrected_mean = -0.5 * (SRmodel.base$pars$sd)^2, # bias correction factorを入れる
                     rec_resid_mat = resid_mat,
                     waa_mat = waa_mat,
                     maa_mat = maa_mat,                     
                     M_mat = M_mat,
                     faa_mat = faa_mat,                     
                     Pope = ifelse(isTRUE(res_vpa$input$Pope),1,0),
                     total_nyear = total_nyear,                     
                     future_initial_year = future_initial_year,                     
                     nsim = nsim,
                     nage = nage,
                     recruit_age = recruit_age,                                          
                     obj_catch = 0, # 0: mean, 1:geomean
                     objective = 0, # 0: MSY, 1: PGY, 2: percentB0 or Bempirical
                     objective_value = 12000,
                     num_to_mass_scale = 1
                     )
x_init <- 0; x_lower <- -3; x_upper <- 4

#system("rm est_MSY_tmb_AS2.so")
#system("rm est_MSY_tmb_AS2.o")



devtools::load_all()
# comple & load cpp file
use_rvpa_tmb(TmbFile = "est_MSY_tmb",
             CppDir = system.file("executable",package="frasyr"),
             RunDir = getwd(), overwrite=TRUE) 

objAD <- TMB::MakeADFun(tmb_data_msy, list(x=x_init), DLL="est_MSY_tmb")

msy_optim <- nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                    lower=list(x=x_lower), upper=list(x=x_upper))#,contol=nlminb_control)

(multi_msy <- as.numeric(exp(msy_optim$par)))

ssb_msy <- objAD$report()$spawner_mat
c(rev(apply(ssb_msy,1,mean))[1],res_MSY$summary$SSB[1])
matplot(ssb_msy, type = "b")

#objAD$report()

c(msy <- exp(-msy_optim$objective),res_MSY$summary$Catch[1])

objAD$report()$N_mat[,,1]

matplot(t(objAD$report()$N_mat[,,1]))

round(objAD$report()$N_mat[,1:8,1],2)
round(res_future_Fcurrent$naa[,1:8,1],2)
tmb_data_msy$number_init

# F=0, sd=0の場合の平衡状態のnumbers at age= 1486.19  901.42  546.74  842.79 一致
# F=0, sd=0の場合の平衡状態のssb= 491083.1 一致
#objAD$report()$spawner_mat[100,1] 
#objAD$report()$spawner_mat[100,1]

# sd=0のときのFmultiplier　frasyr; F multiplier= 0.5268391 , tmb= 0.5268342 (6桁目でずれるがまあまあOK)
