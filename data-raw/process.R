library(devtools)
devtools::load_all() #library(frasyr)
# データの読み込みと資源量推定
caa   <- read.csv("ex2_caa.csv",  row.names=1)
waa   <- read.csv("ex2_waa.csv",  row.names=1)
maa   <- read.csv("ex2_maa.csv",  row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)

# VPAによる資源量推定
res_vpa_example <- vpa(dat,fc.year=2015:2017,tf.year = 2015:2016,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=0.5)
use_data(res_vpa_example, overwrite=TRUE)

# SR関数のフィット
SRdata <- get.SRdata(res_vpa)
res_sr_HSL2 <- fit.SR(SRdata,SR="HS",method="L2",AR=0,hessian=FALSE)
use_data(res_sr_HSL2, overwrite=TRUE)

res_sr_HSL1 <- fit.SR(SRdata,SR="HS",method="L1",AR=0,hessian=FALSE)
use_data(res_sr_HSL1, overwrite=TRUE)

# normal lognormal ----
data_future_HSL2 <- make_future_data(res_vpa_example, # VPAの結果
                                     nsim = 100, # シミュレーション回数
                                     nyear = 20, # 将来予測の年数
                                     future_initial_year_name = 2017, 
                                     start_F_year_name = 2018, 
                                     start_biopar_year_name=2018, 
                                     start_random_rec_year_name = 2018,
                                     # biopar setting
                                     waa_year=2015:2017, waa=NULL, 
                                     waa_catch_year=2015:2017, waa_catch=NULL,
                                     maa_year=2015:2017, maa=NULL,
                                     M_year=2015:2017, M=NULL,
                                     # faa setting
                                     faa_year=2015:2017, 
                                     currentF=NULL,futureF=NULL, 
                                     # HCR setting (not work when using TMB)
                                     start_ABC_year_name=2019, # HCRを適用する最初の年
                                     HCR_beta=1, # HCRのbeta
                                     HCR_Blimit=-1, # HCRのBlimit
                                     HCR_Bban=-1, # HCRのBban
                                     HCR_year_lag=0, # HCRで何年遅れにするか
                                     HCR_function_name = "HCR_default",
                                     # SR setting
                                     res_SR=res_sr_HSL2, 
                                     seed_number=1, 
                                     resid_type="lognormal", 
                                     resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                                     bias_correction=TRUE, # バイアス補正をするかどうか
                                     recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                                     # Other
                                     Pope=res_vpa$input$Pope,
                                     fix_recruit=list(year=c(2020,2021),rec=c(1000,2000)),
                                     fix_wcatch=list(year=c(2020,2021),wcatch=c(1000,2000))
                                     )

res_future_HSL2 <- future_vpa(tmb_data=data_future_HSL2$data,
                              optim_method="none")

res_MSY_HSL2 <- est_MSYRP(data_future=data_future_HSL2, candidate_PGY=c(0.1,0.6),
                          candidate_B0=c(0.2), candidate_Babs=20000)

data_future_HSL1 <- redo_future(data_future_HSL2, list(res_SR=res_sr_HSL1), only_data=TRUE)
res_future_HSL1  <- future_vpa(tmb_data=data_future_HSL1$data,
                              optim_method="none")
res_MSY_HSL1 <- est_MSYRP(data_future=data_future_HSL1, candidate_PGY=c(0.1,0.6),
                          candidate_B0=c(0.2), candidate_Babs=20000)

use_data(res_MSY_HSL1, overwrite=TRUE)
use_data(res_MSY_HSL2, overwrite=TRUE)
use_data(res_future_HSL1, overwrite=TRUE)
use_data(res_future_HSL2, overwrite=TRUE)
use_data(data_future_HSL1, overwrite=TRUE)
use_data(data_future_HSL2, overwrite=TRUE)

# VPA dummy data
data_base <- readr::read_csv("all_dummy_data_base.csv") 
vpadat_base0 <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
                             waa=to_vpa_data(data_base, label_name="waa"),
                             maa=to_vpa_data(data_base, label_name="maa"),
                             M  = 0.2,
                             index = to_vpa_data(data_base, label_name="abund"),
                             maa.tune = NULL,
                             waa.catch = NULL,
                             catch.prop = NULL)
res_vpa_base0_tune <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, abund = c("B", "B"),
                            Pope = TRUE, p.init = 0.5, tune=TRUE, sel.update=TRUE)
#res_vpa_example <- res_vpa_base0_tune
save(vpadat_base0,   file="../data/vpadat_base0.rda")
#save(res_vpa_example,file="../data/res_vpa_example.rda")

data_estb <-readr::read_csv("../inst/extdata/all_dummy_data_estb.csv")
vpadat_estb <- data.handler(caa=to_vpa_data(data_estb, label_name="caa"),
                            waa=to_vpa_data(data_estb, label_name="waa"),
                            maa=to_vpa_data(data_estb, label_name="maa"),
                            M  = 0.4,
                            index = to_vpa_data(data_estb, label_name="abund"),
                            maa.tune = NULL,
                            waa.catch = NULL,
                            catch.prop = NULL)

res_vpa_estb_tune4l_b <- vpa(vpadat_estb,last.catch.zero = FALSE, min.age=c(0,0,0,0,0,0),max.age=c(3,3,0,0,3,3),
                             Pope = TRUE,  tune=TRUE, term.F="all", est.method="ml" ,b.est=TRUE,p.init=c(0.2,0.3,0.6,0.6),abund=c("N","N","N","N","N","N"),fc.year=1998:2000,link=c("id","id","id","id","id","id"))

res_vpa_estb <- res_vpa_estb_tune4l_b
use_data(vpadat_estb, overwrite=TRUE)
use_data(res_vpa_estb, overwrite=TRUE)

# 2022/02/14 res_vpaはバグを多く誘発するためdataから削除。そのかわりに昔のres_vpaをres_vpa_orgとして、dataフォルダにおいたができるだけ使わないこと
