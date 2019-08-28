#- HCRに従った将来予測の実施
devtools::load_all()
#library(frasyr)
library(tidyverse)

#-- 1) 読み込みファイルの設定
#--- VPA結果が保存されているファイルの名前
vpa_file_path <- "data/res_vpa.rda"
#--- 将来予測で仮定する再生産関係の推定結果が保存されているファイルの名前
SR_file_path <- "data/res_sr_HSL2.rda"
#--- MSY推定結果が保存されているファイルの名前
MSY_file_path <- "data/res_MSY.rda"
#--- 将来予測結果を保存するファイルの名前
future_file_path <- "data/res_futures.rda"

#-- 2) 将来予測の基本設定
#--- 漁獲量の計算方法（1:VPAと同じ設定を使う, 2:Popeの近似式を使う, 3:Bavanovの式を使う）
is_pope <- 1
#--- 乱数のシード
future_seed <- 1
#--- MSY計算時のシミュレーション回数(1000回以上推奨)
future_nsim <- 1000
#--- 計算した結果の簡単な図を示す（1:示す,1以外:しない）
future_est_plot <- 1

#-- 生物パラメータ
#--- 将来予測で仮定する年齢別体重(資源量計算用)の設定(1:年で指定, 2:直接指定, 3: MSY計算と同じ)
set_waa_in_future <- 3
if(set_waa_in_future==1){ # 1の場合にはこちらを設定
    waa_year_in_future <- 2015:2017
}
if(set_waa_in_future==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    waa_in_future <- c(100,200,300,400)
}
#--- 将来予測で仮定する年齢別体重(漁獲量計算用)の設定(0:資源量計算用と同じ, 1:年で指定, 2:直接指定, 3: MSY計算と同じ)
set_waa.catch_in_future <- 3
if(set_waa.catch_in_future==1){ # 1の場合にはこちらを設定
    waa.catch_year_in_future <- 2015:2017
}
if(set_waa.catch_in_future==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    waa.catch_in_future <- c(100,200,300,400)
}
#---- 年齢別体重を資源尾数の関数として計算する (0:しない, 1:する, 2: MSY計算と同じ)
waa_fun_future <- 2

#--- 将来予測で仮定する年齢別成熟率の設定 (1:年で指定, 2:直接指定, 3: MSY計算と同じ)
set_maa_in_future <- 3
if(set_maa_in_future==1){ # 1の場合にはこちらを設定
    maa_year_in_future <- 2015:2017
}
if(set_maa_in_future==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    maa_in_future <- c(0,0,0.5,1)
}

#--- 将来予測で仮定する自然死亡係数の設定 (1: 年で指定, 2: 直接指定, 3:MSY計算と同じ))
set_M_in_future <- 3
if(set_M_in_future==1){ # 1の場合にはこちらを設定
    M_year_in_future <- 2015:2017
}
if(set_M_in_future==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    M_in_future <- c(0,0,0.5,1)
}

#-- 再生産関係
#--- MSY計算とすべて同じ仮定を使うか (1: 使う, 0: 使わずすべて手動で設定する)
use_MSY_SR <- 1
if(use_MSY_SR==0){ # すべて手動で計算する場合、以下のオプションを設定
    #---  バイアス補正（シミュレーションの平均と決定論的予測を一致させる）(1:する, 1以外: しない)
    bias_correction_future <- 1
    #--- 加入変動の誤差分布 (1: 対数正規誤差, 2: 残差リサンプリング)
    SR_error_future <- 1
    #---- SR_error_future=1（対数正規分布の誤差）の場合の設定
    if(SR_error_future==1){ # 対数正規分布の場合は自己相関のオプションも選ぶ
        #--- 自己相関の仮定 (-1: 推定結果どおりに設定する, 0: 推定結果にかかわらず自己相関を「なし」にする,
        #---               0以上の数字：推定結果にかかわらず自己相関係数をここで設定した値にする)
        set_AR_future <- -1
    }
    #---- SR_error_future=2（リサンプリング誤差）の場合の設定
    if(SR_error_future==2){ 
        #--- リサンプリングの年の範囲(0: 全年, それ以外の数字: 指定した加入年の残差のみリサンプリング）
        set_resample_year <- 0 # or 1990:2000
    }
}

#-- 年数や回数などの設定
#--- 将来予測開始年
future_start_year <- 2018
#--- ABC計算年（この年からHCRに沿った漁獲を開始する）
future_ABC_year   <- 2020
#--- 将来予測の実施年数
future_nyear   <- 30
#--- シミュレーションの一回目は決定論的予測の結果とする (1: する, 0: しない)
det_run <- 1
#-- 直近の加入や漁獲の仮定
#--- 特定の年の加入を外部から与える (0: 設定しない, 1: 設定する)
set_specific_recruit <- 0
#---- 加入を外部から与える場合の設定
if(set_specific_recruit==1){
    recruit_year <- c(2018,2019)
    recruit_number <- c(10000,20000)
}
#--- 特定の年の漁獲量を外部から与える (0: 設定しない, 1: 設定する)
set_specific_catch <- 0
#---- 漁獲量を外部から与える場合の設定
if(set_specific_catch==1){
    catch_year <- c(2018,2019)
    catch_weight <- c(10000,10000)
}

#-- 漁獲シナリオの設定
#--- HCRを実施するときのSSB参照年 (0: ABC算定年とSSBの参照年が同じ。1以上の整数：時間遅れの年数。デフォルトはゼロ)
HCR_year_lag <- 0
#--- HCRが有効な場合の管理基準値 (0: MSY_resで指定されたBtarget0, Blimit0, Bban0をそのまま使う。
#---                         1: MSY_resの設定を上書きする)
overwrite_RP <- 0
if(overwrite_RP==1){ # MSY_resの設定を上書きする場合
    #--- 目標管理基準値の選択 (0: MSY,
    #---                   1以上の数字: MSY_res$summaryの行数,
    #---                   負の数字: インタラクティブに決定)
    set_Btarget <- 0
    #--- 限界管理基準値の選択 (0: 60%MSY,
    #---                   1以上の数字: MSY_res$summaryの行数,
    #---                   負の数字: インタラクティブに決定)
    set_Blimit  <- 0
    #--- 禁漁水準の選択      (0: 10%MSY,
    #---                   1以上の数字: MSY_res$summaryの行数,
    #---                   負の数字: インタラクティブに決定)
    set_Bban  <- 0
}
#---- HCRの将来予測におけるデフォルトのベータ
beta_default <- 0.8
#---- ベータをいろいろ変える将来予測におけるベータの範囲
beta_table <- c(0.5,0.6,0.7,0.8,0.9,1)

#--- ABC.year以前のF at age (1: vpaの結果の中のFc.at.ageをそのまま使う, 2:手動でF at ageを設定する)
set_Fc_at_age_preABC <- 1
#---- 上で2を選んだ場合はこちらも設定する
if(set_Fc_at_age_preABC==2){ 
    Fc_at_age_preABC <- c(0.2,0.3,0.3,0.44) 
}
#--- ABC.year以降のF at age (1: ABC.year以前のF at ageと共通, 
#---                        2: 手動でF at ageを設定する,
#---                        3: Fmsy at age (=MSY_resで"Btarget0"の管理基準値に対応するF at age),
#---                        4: ABC.year以前のF at ageと同じ選択率で%SPRがFmsy相当になるようなF(未実装())
set_Fc_at_age_afterABC <- 3
#---- 上で2を選んだ場合はこちらも設定する
if(set_Fc_at_age_afterABC==2){ 
    Fc_at_age_afterABC <- c(0.2,0.3,0.3,0.44) # ベクトルとして直接指定する場合
}
#-- 出力の調整：下記の各項目について表として出力したい年数を入れるか、表が必要ない場合はマイナス値を入れる
#--- 将来の平均漁獲量
year_catch_average <- c(2020:2030,2040,2050)
#--- 目標管理基準値を上回る確率
year_catch_average <- c(2020:2030,2040,2050)
#--- 目標管理基準値を上回る確率
year_ssbtarget_prob <- c(2020:2030,2040,2050)
#--- 限界管理基準値を上回る確率
year_ssblimit_prob  <- c(2020:2030,2040,2050)
#--- 禁漁水準を上回る確率
year_ssbban_prob    <- -1
#--- 過去の最低親魚量を上回る確率
year_ssbmin_prob    <- -1
#--- 過去の最高親魚量を上回る確率
year_ssbmax_prob    <- -1
#--- 漁獲量のAAV
year_catch_aav      <- -1
#--- Fの削減率の平均
year_Fsakugen_mean  <- -1


####################################################
### 以下は基本的には編集しないこと
####################################################
res_vpa_name <- load(vpa_file_path)
res_vpa <- get(res_vpa_name)

res_SR_name <- load(SR_file_path)
res_SR <- get(res_SR_name)

res_MSY_name <- load(MSY_file_path)
res_MSY <- get(res_MSY_name)

input_MSY <- res_MSY$input.list[[1]]

{if(use_MSY_SR==0){
    if(SR_error_future==1){# lognormal
        if(res_SR_future$input$SR=="HS") SRfun_future <- HS.recAR
        if(res_SR_future$input$SR=="BH") SRfun_future <- BH.recAR
        if(res_SR_future$input$SR=="RI") SRfun_future <- RI.recAR

        rho_future <- switch(as.character(set_AR_future),
                          "-1"=res_SR_future$pars$rho,
                          "0"=0,
                          set_AR_future)
        opt_SR_future <- list(a=res_SR_future$pars$a,
                           b=res_SR_future$pars$b,
                           rho=rho_future,
                           sd=res_SR_future$pars$sd,
                           bias.correction=ifelse(bias_correction_future==1,TRUE,FALSE),
                           resample=FALSE,resid=res_SR_future$pars$resid,
                           resid.year=NULL)
        if(rho_future>0){
            opt_SR_future$resid.year <- AR_average_year_future 
        }
    }

    if(SR_error_future==2){# lognormal
        SRfun_future <- resample.rec
        set_resample_year <- ifelse(set_resample_year==0,
                                    TRUE,
                                    set_resample_year)
        resid_selected <- res_SR_future$resid[res_SR_future$input$SRdata$year%in%
                                              set_resample_year]
        opt_SR_future <- list(resample=TRUE,
                              SR=res_SR_future$input$SR,
                              resid=resid_selected)    
    }
    opt_SR_future$bias.correction <- ifelse(bias_correction_future==1,TRUE,FALSE)
 }
 else{
     opt_SR_future <- input_MSY$rec.arg
     SRfun_future <- input_MSY$recfunc
}}

if(set_Fc_at_age_afterABC==3){
    Fc_at_age_preABC <- res_MSY$Fvector %>%
        slice(which(res_MSY$summary$RP.definition=="Btarget0")) %>%
        as.numeric()
    cat("Set future F vector as: ")
    Fc_at_age_preABC  %>% round(2) %>% print()
}
if(set_Fc_at_age_afterABC==4){
    stop("error : not implemented yet, sorry")
}


if(overwrite_RP==1){
    options(tibble.width=100)
    # define RP.definition for Btarget
    if(set_Btarget!=0){
        if(set_Btarget<0){
            print(select(res_MSY$summary,-AR,-Catch.CV))
            set_Btarget <- readline("Enter row number to be Btarget: ")
            set_Btarget <- as.integer(set_Btarget)
        }
        res_MSY$summary$RP.definition[1] <- NA    
        res_MSY$summary$RP.definition[set_Btarget] <- "Btarget0"
    }
    # define RP.definition for Blimit
    if(set_Blimit!=0){
        if(set_Blimit<0){
            print(select(res_MSY$summary,-AR,-Catch.CV))
            set_Blimit <- readline("Enter row number to be Blimit: ")
            set_Blimit <- as.integer(set_Blimit)
        }
        res_MSY$summary$RP.definition[which(res_MSY$summary$RP.definition=="Blimit0")] <- NA
        res_MSY$summary$RP.definition[set_Blimit] <- "Blimit0"
    }
    # define RP.definition for Bban
    if(set_Bban!=0){
        if(set_Bban<0){
            print(select(res_MSY$summary,-AR,-Catch.CV,-RP.definition))
            set_Bban <- readline("Enter row number to be Bban: ")
            set_Bban <- as.integer(set_Bban)
        }
        res_MSY$summary$RP.definition[which(res_MSY$summary$RP.definition=="Bban0")] <- NA    
        res_MSY$summary$RP.definition[set_Bban] <- "Bban0"
    }
}

HCR.future <- list(Blim    = derive_RP_value(res_MSY$summary,"Blimit0")$SSB,
                   Bban    = derive_RP_value(res_MSY$summary,"Bban0")$SSB,
                   beta    = beta_default,
                   year.lag= HCR_year_lag)

input_future_0.8HCR <- list(
    res0     =res_vpa,
    currentF =switch(as.character(set_Fc_at_age_preABC),"1"=NULL,Fc_at_age_preABC),
    multi    =1,
    N        =future_nsim,    
    futureF  =switch(set_Fc_at_age_afterABC,
                     NULL,
                     Fc_at_age_preABC,
                     Fc_at_age_preABC,
                     stop("Set appropriate number (1-3) in set_Fc_at_age_afterABC")),
    nyear    =future_nyear,
    Pope     =switch(is_pope,
                     res_vpa$input$Pope,
                     TRUE,
                     FALSE,
                     stop("Set appropriate number (1-3) in is_pope")),
    outtype  ="FULL",
    multi.year=1,
    start.year=future_start_year,
    ABC.year  =future_ABC_year,
    waa.year      =switch(set_waa_in_future,
                          waa_year_in_future,
                          input_MSY$waa.year,
                          NULL,
                          stop("Set appropriate number (1-3) in set_waa_in_future")),
    waa.catch.year=switch(set_waa.catch_in_future,
                          waa.catch_year_in_future,
                          input_MSY$waa.catch.year,
                          NULL,
                          stop("Set appropriate number (1-3) in set_waa.catch_in_future")),
    maa.year      =switch(set_maa_in_future,
                          maa_year_in_future,
                          input_MSY$maa.year,                          
                          NULL,
                          stop("Set appropriate number (1-3) in set_maa_in_future")),
    M.year        =switch(set_M_in_future,
                          M_year_in_future,
                          input_MSY$M.year,
                          NULL,
                          stop("Set appropriate number (1-3) in set_M_in_future")),            
    waa           =switch(set_waa_in_future,
                          NULL,
                          waa_in_future,
                          input_MSY$waa,                          
                          stop("Set appropriate number (1-3) in set_waa_in_future")),
    waa.catch     =switch(set_waa.catch_in_future,
                          NULL,                          
                          waa.catch_in_future,
                          input_MSY$waa.catch,                          
                          stop("Set appropriate number (1-3) in set_waa.catch_in_future")),            
    maa           =switch(set_maa_in_future,
                          NULL,
                          maa_in_future,
                          input_MSY$maa,
                          stop("Set appropriate number (1-3) in set_maa_in_future")),    
    M             =switch(set_M_in_future,
                          NULL,
                          M_in_future,
                          input_MSY$M,                          
                          stop("Set appropriate number (1-3) in set_M_in_future")),
    seed       = future_seed,
    strategy   = "F", 
    HCR        = HCR.future,
    use.MSE    = FALSE,
    MSE.options= NULL,
    beta       = NULL,
    delta      = NULL,
    Blim       = 0,
    Bban       = 0,
    plus.group = res_vpa$input$plus.group,
    silent     = TRUE,
    is.plot    = future_est_plot, 
    random.select=NULL, 
    recfunc    = SRfun_future, 
    rec.arg    = opt_SR_future,
    rec.new    = (if(set_specific_recruit==0) NULL else list(year=recruit_year,rec=recruit_number)),
    pre.catch  = (if(set_specific_catch==0) NULL else list(year=catch_year,wcatch=catch_weight)),
    waa.fun    = ifelse(waa_fun_future==0,
                        FALSE,
                 ifelse(waa_fun_MSY==1,TRUE,input_MSY$waa.fun)),
    det.run    = as.logical(det_run),    
    naa0=NULL,eaa0=NULL,ssb0=NULL,faa0=NULL,
    add.year=0)

# Default HCR Run
res_future_0.8HCR <- do.call(future.vpa,input_future_0.8HCR)

# current F run
input_future_current <- input_future_0.8HCR
input_future_current$futureF <- NULL
res_future_current <- do.call(future.vpa,input_future_current)

plot_futures(res_vpa,list(res_future_0.8HCR,res_future_current))

# kobe II table
kobeII.data <- calc_kobeII_matrix(res_future_0.8HCR,
                                  res_MSY$summary,
                                  Btarget=c("Btarget0"), 
                                  Blimit=c("Blimit0"),
                                  beta=beta_table)
kobeII.table <- make_kobeII_table(kobeII.data,
                                  res_vpa        = res_vpa,
                                  year.catch     = year_catch_average,
                                  year.Fsakugen  = year_Fsakugen_mean,
                                  year.ssbtarget = year_ssbtarget_prob,
                                  year.ssblimit  = year_ssblimit_prob,
                                  year.ssbban    = year_ssbban_prob,
                                  year.ssbmin    = year_ssbmin_prob,
                                  year.ssbmax    = year_ssbmax_prob,                                                       year.aav       = year_catch_aav) 


# save results
cat("\n***** Summary results *****\n")
save(res_future_0.8HCR,res_future_current,file=future_file_path)
cat(paste("create output file of",future_file_path,": Future simualtion results\n"))


