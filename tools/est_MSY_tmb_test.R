#library(frasyr)
devtools::load_all()
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
# 5万回、100年で43秒
a3 <- system.time(res_future_Fcurrent <- future.vpa(res_vpa,
                      multi=1,
                      nyear=100, # 将来予測の年数
                      start.year=2017, # 将来予測の開始年
                      N=50000, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
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
                                   resid=SRmodel.base$resid)))

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

# 2017年をinitial valueにして、2018年からF, 加入について将来予測する場合
devtools::load_all()
future_data1 <- make_future_data(res_vpa,
                      nsim = 1000, # number of simulation
                      nyear = 50, # number of future year
                      future_initial_year_name = 2017,
                              start_F_year_name = 2018,
                              start_biopar_year_name=2018,
                              start_random_rec_year_name = 2018,                                # biopar setting
                          waa_year=2016:2018, waa=NULL,
                          waa.catch_year=2016:2018, waa.catch=NULL,
                          maa_year=2016:2018, maa=NULL,
                          M_year=2016:2018, M=NULL,
                          # faa setting
                          faa_year=2016:2018, faa=NULL,
                          # HCR setting
                          HCR_beta=1,
                          HCR_Blimit=-1,
                          HCR_Bban=-1,
                          HCR_year_lag=0,
                          # SR setting
                          res_SR=SRmodel.base,                       
                          seed_number=1
                      )

res_future_tmb <- future_vpa(future_data1$data, optim_method="tmb",
                             x_init = 0,
                             x_lower = -3,
                             x_upper = 3,
                             trace.multi=c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100),
                             compile=TRUE)

future_data1 <- make_future_data(res_vpa,
                      nsim = 1000, # number of simulation
                      nyear = 50, # number of future year
                      future_initial_year_name = 2017,
                              start_F_year_name = 2018,
                              start_biopar_year_name=2018,
                              start_random_rec_year_name = 2018,                                # biopar setting
                          waa_year=2016:2018, waa=NULL,
                          waa.catch_year=2016:2018, waa.catch=NULL,
                          maa_year=2016:2018, maa=NULL,
                          M_year=2016:2018, M=NULL,
                          # faa setting
                          faa_year=2016:2018, faa=NULL,
                          # HCR setting
                          HCR_beta=1,
                          HCR_Blimit=20000,
                          HCR_Bban=-1,
                          HCR_year_lag=0,
                          # SR setting
                          res_SR=SRmodel.base,                       
                          seed_number=1
                      )

res_future_tmb <- future_vpa(future_data1$data,
                             optim_method="none",
                             x_init = 0,
                             x_lower = -3,
                             x_upper = 3,
                             trace.multi=c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100),
                             compile=TRUE)

#> res_future_tmb$multi
#[1] 0.5402367

## 以下、もう動かない
# 同じシミュレーションをもう一回できるかどうか=>完全に再現できる
res1_replicate <- tmb_future(res_vpa,skip_setting=TRUE,tmb_data=res1$tmb_data)

# 一部の設定を変えてできるか？
tmb_data_dummy <- res1$tmb_data
tmb_data_dummy$naa_mat[] <- res1$naa[,,1]
tmb_data_dummy$rec_resid_mat[] <- res1$tmb_data$rec_resid_mat[,1]
tmb_data_dummy$rec_resid_mat[35:nrow(tmb_data_dummy$rec_resid_mat),] <- rnorm(100)
tmb_data_dummy$future_initial_year <- 35
tmb_data_dummy$start_F_year <- 35
res1_replicate2 <- tmb_future(res_vpa,skip_setting=TRUE,tmb_data=tmb_data_dummy)

# R関数と同じ結果が出るか？
a1 <- system.time(res1.n10 <- tmb_future(res_vpa,nsim=50000,nyear=100,
                   SRmodel=SRmodel.base,                   
                   future_initial_year_name=2017,
                   start_F_year_name=2018,
                   start_random_rec_year_name=2018,
                   x_init=0,x_upper=0,x_lower=0)) # 5万回30年で73秒

tmb_data_dummy <- res1.n10$tmb_data # 5万回30年で61秒 => future.vpaのほうが全然速い。ベクトル化か
tmb_data_dummy$x <- 0
tmb_data_dummy$what_return <- "stat"
a2 <- system.time(res1_replicate3 <- do.call(est_MSY_R,tmb_data_dummy))

# オプションuse_tmbがうまく動くか
res1_R <- tmb_future(res_vpa,nsim=1000,nyear=30,
                     SRmodel=SRmodel.base,
                     optim_method="both",
                     future_initial_year_name=2017,
                     start_F_year_name=2018,
                     start_random_rec_year_name=2018,
                     x_init=0,x_upper=3,x_lower=-5) 
round(res1_R$tmb$multi_msy,4)==round(res1_R$R$multi_msy,4)

# n=1でうまく動くかどうか
res1_n1 <- tmb_future(res_vpa,nsim=1,nyear=30,
                      SRmodel=SRmodel.base,
                      future_initial_year_name=2017,
                      start_F_year_name=2018,
                      start_random_rec_year_name=2018)

# 2000年をinitial valueにして過去のFに従って漁獲, 加入も将来予測=>同じになった
res2 <- tmb_future(res_vpa,nsim=1000,nyear=47,
                   SRmodel=SRmodel.base,                                      
                   future_initial_year_name=2000,
                   start_F_year_name=2018,
                   start_random_rec_year_name=2018)

# 2000年をinitial valueにして過去のFは変えるが、残差は固定=>うまくいっていない。
res3 <- tmb_future(res_vpa,nsim=1000,nyear=47,
                   SRmodel=SRmodel.base,                                      
                   future_initial_year_name=1989,
                   start_F_year_name=1989,
                   start_random_rec_year_name=2018)

# 2000年をinitial valueにして、過去のFを変えずに、残差を変動。
res4 <- tmb_future(res_vpa,nsim=1000,nyear=47,
                   SRmodel=SRmodel.base,                                      
                   future_initial_year_name=2000,
                   start_F_year_name=2018,
                   start_random_rec_year_name=2001)

# res1と同じ設定だが、パラメータの範囲を固定する(current Fでの予測)
res1_fix <- tmb_future(res_vpa,nsim=1000,nyear=30,
                   SRmodel=SRmodel.base,                                          
                   future_initial_year_name=2017,
                   start_F_year_name=2018,
                   start_random_rec_year_name=2018,
                   x_lower=0, x_upper=0, x_init=0)


plot(rownames(res1$ssb), rowMeans(res1$ssb),type="l")
points(rowMeans(res2$ssb),rowMeans(res1$ssb), type="l",col=2)
points(rownames(res3$ssb), rowMeans(res3$ssb),type="l",col=3)
points(rownames(res4$ssb), rowMeans(res4$ssb),type="l",col=4)
points(rownames(res1_fix$ssb), rowMeans(res1_fix$ssb),type="l",col=5)

all_ssb <- bind_rows(convert_2d_future(res1$ssb, "SSB", label="res1"),
                     convert_2d_future(res2$ssb, "SSB", label="res2"),
                     convert_2d_future(res3$ssb, "SSB", label="res3"),
                     convert_2d_future(res4$ssb, "SSB", label="res4"))


all_ssb %>% ggplot() +
    geom_boxplot(aes(x=factor(year),y=value,fill=label)) +
    facet_wrap(.~label)
          

#--------------

