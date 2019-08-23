#- MSY管理基準値を計算する
devtools::load_all()

#-- 1) 読み込みファイルの設定
#--- VPA結果が保存されているファイルの名前
vpa_file_path <- "data/res_vpa.rda"
#--- MSY推定で仮定する再生産関係の推定結果が保存されているファイルの名前
SR_file_path_MSY <- "data/res_sr_HSL2.rda"
#--- MSY推定結果を保存するファイル/オブジェクトの名前
MSY_file_path <- "data/res_MSY.rda"

#-- 2) MSY推定の設定（F一定の条件下での将来予測をもとにする）
#--- Fcurrent（1: vpaの結果の中のFc.at.ageをそのまま使う, 2:手動でFcurrentを設定する）
use_vpa_Fc_at_age <- 1
#---- 上で2を選んんだ場合はこちらも設定
if(use_vpa_Fc_at_age==2){ 
    Fc_at_age <- c(0.2,0.3,0.3,0.44) # ベクトルとして与えるFcurrentの指定
}

#-- 各種計算方法
#--- 漁獲量の計算方法（1:VPAと同じ設定を使う, 2:Popeの近似式を使う, 3:Bavanovの式を使う）
is_pope <- 1
#--- 乱数のシード
MSY_seed <- 1
#--- MSY計算時のシミュレーション回数(1000回以上推奨)
MSY_nsim <- 10
#--- 計算した結果の簡単な図を示す（1:示す,1以外:しない）
MSY_est_plot <- 1

#-- 生物パラメータ
#--- MSY計算時の年齢別体重(資源量計算用)の設定(1:年で指定する、2:直接指定する)
set_waa_in_MSY <- 1
if(set_waa_in_MSY==1){ # 1の場合にはこちらを設定
    waa_year_in_MSY <- 2015:2017
}
if(set_waa_in_MSY==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    waa_in_MSY <- c(100,200,300,400)
}
#--- MSY計算時の年齢別体重(漁獲量計算用)の設定(0:資源量計算用と同じ、1:年数で指定、2:直接指定)
set_waa.catch_in_MSY <- 0
if(set_waa.catch_in_MSY==1){ # 1の場合にはこちらを設定
    waa.catch_year_in_MSY <- 2015:2017
}
if(set_waa.catch_in_MSY==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    waa.catch_in_MSY <- c(100,200,300,400)
}
#---- 年齢別体重を資源尾数の関数として計算する(0:しない, 1:する)
waa_fun_MSY <- 0

#--- MSY計算時の年齢別成熟率の設定(1:年数で指定、2:直接指定)
set_maa_in_MSY <- 1
if(set_maa_in_MSY==1){ # 1の場合にはこちらを設定
    maa_year_in_MSY <- 2015:2017
}
if(set_maa_in_MSY==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    maa_in_MSY <- c(0,0,0.5,1)
}

#--- MSY計算時の自然死亡係数の設定(1:年数で指定、2:直接指定)
set_M_in_MSY <- 1
if(set_M_in_MSY==1){ # 1の場合にはこちらを設定
    M_year_in_MSY <- 2015:2017
}
if(set_M_in_MSY==2){ # 2の場合にはこちらを設定。年毎に異る場合は年齢×年の行列を入れる
    M_in_MSY <- c(0,0,0.5,1)
}

#-- 再生産関係
#---  バイアス補正（シミュレーションの平均と決定論的予測が一致するようにする）(1:する、1以外: しない)
bias_correction_MSY <- 1
#--- 加入変動の誤差分布(1: 対数正規誤差, 2: 残差リサンプリング)
SR_error_MSY <- 1
#---- SR_error_MSY=1（対数正規分布の誤差）の場合の設定
if(SR_error_MSY==1){ # 対数正規分布の場合は自己相関のオプションも選ぶ
    #--- 自己相関の仮定 (-1: 推定結果どおりに設定する, 0: 推定結果にかかわらず自己相関を「なし」にする,
    #---               0以上の数字：推定結果にかかわらず自己相関係数をここで設定した値にする)
    set_AR_MSY <- -1
}
#---- SR_error_MSY=2（リサンプリング誤差）の場合の設定
if(SR_error_MSY==2){ 
    #--- リサンプリングの年の範囲(0: 全年, それ以外の数字: 指定した加入年の残差のみリサンプリング）
    set_resample_year <- 0 # or 1990:2000
}

#-- MSY推定時の設定
#--- 漁獲量曲線を細かく書くか（0: 書かない（計算時間短縮), 1: 書く）
calc_yieldcurve <- 1
#--- 管理基準値を計算するPGYレベル（-1: 計算しない, 0から1までの数字のベクトル: MSYx割合に対応する親魚レベル）
set_PGY <- c(0.6,0.1)
#--- PGYの下側のみの管理基準値を計算する（1: 下側のみ（計算時間短縮）, 2:両側）
only_lowerPGY <- 1
#--- 管理基準値を計算するB0レベル（-1: 計算しない, 0から1までの数字のベクトル: B0x割合に対応する親魚レベル）
set_B0 <- c(0.2,0.4)
#--- 特定の親魚量レベルを直接指定する（-1: 計算しない, 親魚量のベクトル: その親魚量に維持するときの管理基準値）
set_Babs <- -1
#--- 平衡状態にあると仮定する年の決め方（1: 世代時間から計算する, 2: 具体的な年数を与える）
set_nyear_MSY <- 1
#---- 世代時間から計算する場合
if(set_nyear_MSY==1){
    #--- 世代時間の推定方法（0: 自動的に計算, 1以上の数：ここで指定された値を世代時間（年）とする）
    set_GT <- 0
    #--- 世代時間の何倍を平衡状態に置くか（デフォルトは20)
    GT_mutipiler <- 20
}
#---- 具体的な年数を与える場合
if(set_nyear_MSY==2){
    nyear_MSY <- 100
}
#--- 複数シミュレーションの結果をあわせて最大化するときの統計量 (1: 平均（デフォルト）, 2: 中央値)
stat_maximize <- 1
#--- 自己相関を考慮した管理基準値を計算するか (0:しない, 1:する)
calc_RP_with_AR <- 0
#---- 自己相関を考慮した管理基準値を計算する場合(ARありの場合だけ有効)
if(calc_RP_with_AR==1){
    #--- 何年分の残差の平均から平衡状態からスタートさせるか
    AR_average_year_MSY <- 5    
    #--- 平衡状態から何年進めたときの値を管理基準値とするか（デフォルトは5年）
    forward_year_ARRP <- 5
    #--- 残差を手動で設定する場合
    current.resid <- 0
}
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


####################################################
### 以下は基本的には編集しないこと
####################################################
res_vpa_name <- load(vpa_file_path)
res_vpa <- get(res_vpa_name)

res_SR_name_MSY <- load(SR_file_path_MSY)
res_SR_MSY <- get(res_SR_name_MSY)

future_MSY_year <- as.numeric(rev(colnames(res_vpa$naa))[1])

if(SR_error_MSY==1){# lognormal
    if(res_SR_MSY$input$SR=="HS") SRfun_MSY <- HS.recAR
    if(res_SR_MSY$input$SR=="BH") SRfun_MSY <- BH.recAR
    if(res_SR_MSY$input$SR=="RI") SRfun_MSY <- RI.recAR

    rho_MSY <- switch(as.character(set_AR_MSY),
                      "-1"=res_SR_MSY$pars$rho,
                      "0"=0,
                      set_AR_MSY)
    opt_SR_MSY <- list(a=res_SR_MSY$pars$a,
                       b=res_SR_MSY$pars$b,
                       rho=rho_MSY,
                       sd=res_SR_MSY$pars$sd,
                       bias.correction=ifelse(bias_correction_MSY==1,TRUE,FALSE),
                       resample=FALSE,resid=res_SR_MSY$pars$resid,
                       resid.year=NULL)
    if(rho_MSY>0){
        opt_SR_MSY$resid.year <- AR_average_year_MSY 
    }
}

if(SR_error_MSY==2){# lognormal
    SRfun_MSY <- resample.rec
    set_resample_year <- ifelse(set_resample_year==0,
                                    TRUE,
                                    set_resample_year)
    resid_selected <- res_SR_MSY$resid[res_SR_MSY$input$SRdata$year%in%
                                       set_resample_year]
    opt_SR_MSY <- list(resample=TRUE,
                       SR=res_SR_MSY$input$SR,
                       resid=resid_selected)    
}

opt_SR_MSY$bias.correction <- ifelse(bias_correction_MSY==1,TRUE,FALSE)

input_future_MSY <- list(
    res0     =res_vpa,
    currentF =switch(as.character(use_vpa_Fc_at_age),"1"=NULL,Fc_at_age),
    multi    =1,
    N        =10,    
    futureF  =NULL,
    nyear    =10,
    Pope     =switch(is_pope,
                     res_vpa$input$Pope,
                     TRUE,
                     FALSE,
                     stop("Set appropriate number (1-3) in is_pope")),
    outtype  ="FULL",
    multi.year=1,
    start.year=future_MSY_year+1,
    ABC.year  =future_MSY_year+1,
    waa.year      =switch(set_waa_in_MSY,
                          waa_year_in_MSY,
                          NULL,
                          stop("Set appropriate number (1-2) in set_waa_in_MSY")),
    waa.catch.year=switch(set_waa.catch_in_MSY,
                          waa.catch_year_in_MSY,
                          NULL,
                          stop("Set appropriate number (1-2) in set_waa.catch_in_MSY")),
    maa.year      =switch(set_maa_in_MSY,
                          maa_year_in_MSY,
                          NULL,
                          stop("Set appropriate number (1-2) in set_maa_in_MSY")),
    M.year        =switch(set_M_in_MSY,
                          M_year_in_MSY,
                          NULL,
                          stop("Set appropriate number (1-2) in set_M_in_MSY")),            
    waa           =switch(set_waa_in_MSY,
                          NULL,
                          waa_in_MSY,
                          stop("Set appropriate number (1-2) in set_waa_in_MSY")),
    waa.catch     =switch(set_waa.catch_in_MSY,
                          NULL,                          
                          "2"=waa.catch_in_MSY,
                          waa.catch_in_MSY,
                          stop("Set appropriate number (1-2) in set_waa.catch_in_MSY")),            
    maa           =switch(set_maa_in_MSY,
                          NULL,
                          maa_in_MSY,
                          stop("Set appropriate number (1-2) in set_maa_in_MSY")),    
    M             =switch(set_M_in_MSY,
                          NULL,
                          M_in_MSY,
                          stop("Set appropriate number (1-2) in set_M_in_MSY")),
    seed       =MSY_seed,
    strategy   ="F", 
    HCR        =NULL,
    use.MSE    =FALSE,
    MSE.options=NULL,
    beta       =NULL,
    delta      =NULL,
    Blim       =0,
    Bban       =0,
    plus.group =res_vpa$input$plus.group,
    silent     =TRUE,
    is.plot    =FALSE, 
    random.select=NULL, 
    recfunc    =SRfun_MSY, 
    rec.arg    =opt_SR_MSY,
    rec.new    =NULL,
    pre.catch  =NULL,
    waa.fun    =ifelse(waa_fun_MSY==0,FALSE,TRUE), 
    naa0=NULL,eaa0=NULL,ssb0=NULL,faa0=NULL,
    add.year=0, det.run=TRUE)

# test run
res_future_preMSY <- do.call(future.vpa,input_future_MSY)

is.estAR.RP <- as.logical(calc_RP_with_AR) && as.logical(rho_MSY>0)

input_est_MSY <-
    list(vpares         =res_vpa,
         farg           =input_future_MSY,
         N              =MSY_nsim,
         is.plot        =ifelse(MSY_est_plot==1,TRUE,FALSE),
         calc.yieldcurve=ifelse(calc_yieldcurve==1,TRUE,FALSE),
         onlylower.pgy  =ifelse(only_lowerPGY==1,TRUE,FALSE),         
         PGY            =(if(set_PGY[1] <0) NULL else set_PGY),
         B0percent      =(if(set_B0[1]  <0) NULL else set_B0),
         Bempirical     =(if(set_Babs[1]<0) NULL else set_Babs),
         seed           =MSY_seed,
         eyear          =0, # 将来予測の最後のeyear+1年分を平衡状態とする
         nyear          =switch(set_nyear_MSY,
                                NULL,
                                nyear_MSY,
                                stop("Set appropriate number (1-2) in set_nyear_MSY")),
         long.term      =GT_mutipiler,
         GT             =(if(set_GT==0) NULL else set_GT),
         FUN            =switch(stat_maximize,
                                mean,
                                median,
                                stop("Set appropriate number (1-2) in stat_maximize")),
         optim.method   ="optimize",
         max.target     ="catch.mean", # この設定は意味ない
         trace.multi    =c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100),
         estAR.RP       =is.estAR.RP,
         resid.year     =ifelse(isTRUE(is.estAR.RP),AR_average_year_MSY,0),
         mY             =ifelse(isTRUE(is.estAR.RP),forward_year_ARRP,5),
         current.resid  =NULL         
         ) 

# do estimation
res_MSY <- do.call(est.MSY,input_est_MSY)

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

# save results
cat("\n***** Summary results *****\n")
print(select(res_MSY$summary,-AR,-Catch.CV))
save(res_MSY,file=MSY_file_path)
options(tibble.width=NULL)
cat(paste("create output file of",MSY_file_path,": MSY estimation results\n"))

if(0){
# summarize MSY reference point
(refs.all <- input_est_MSY$summary_tb)

# refs.allの中からRP.definitionで指定された行だけを抜き出す
(refs.base <- refs.all %>%
    dplyr::filter(!is.na(RP.definition)) %>% # RP.definitionがNAでないものを抽出
    arrange(desc(SSB)) %>% # SSBを大きい順に並び替え
    select(RP.definition,RP_name,SSB,SSB2SSB0,Catch,Catch.CV,U,Fref2Fcurrent)) #　列を並び替え

# HCRによる将来予測
input.abc <- future.Fcurrent$input # Fcurrentにおける将来予測の引数をベースに将来予測します
input.abc$multi <- derive_RP_value(refs.base,"Btarget0")$Fref2Fcurrent # currentFへの乗数を"Btarget0"で指定した値に
input.abc$silent <- TRUE
input.abc$HCR <- list(Blim=derive_RP_value(refs.base,"Blimit0")$SSB,
                      Bban=derive_RP_value(refs.base,"Bban0")$SSB,
                      beta=0.8,year.lag=0) # BlimitはBlimit0, BbanはBban0の値
input.abc$N <- 1000
future.default <- do.call(future.vpa,input.abc) # デフォルトルールの結果→図示などに使う
}
