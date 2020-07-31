# basic dynamics ----

#'
#' future_vpaにインプットとして入れる将来予測の空のarrayを生成する関数
#'
#' @param res_vpa vpaの結果 (vpa関数の返り値)
#' @param nsim シミュレーションの繰り返し回数
#' @param nyear 将来予測の実施年数
#' @param plus_age プラスグループとして計算する行（年齢ではないことに注意）。デフォルト値（NULL）ならfuture_initial_year_name年にNA以外の数値が入っている一番大きい年齢をプラスグループの年齢とする
#' @param future_initial_year_name 将来予測の「初期値となる」年齢別資源尾数を参照する年。この年の年齢別資源尾数を使って翌年の個体群動態が将来予測で決定される
#' @param start_F_year_name 将来予測でF全体にmultiplierを乗じる場合、multiplierを乗じる最初の年
#' @param start_biopar_year_name 生物パラメータを将来の生物パラメータとして設定された値に置き換える年の最初の年
#' @param start_random_rec_year_name 将来の加入を再生産関係から予測する最初の年
#' @param waa_year 将来の年齢別体重を過去の平均値とする場合、過去のパラメータを平均する期間, maa_year, M_yearも同様
#' @param waa 将来の年齢別体重を直接与える場合, maa, M_yearも同様
#' @param waa_fun log(weight)~log(number)の回帰式から将来のweightを予測する
#' @param start_waafun_year_name 上記の設定がスタートする最初の年。それ以外の年は上で設定されたパラメータが使われる
#' @param faa_year 将来のFを過去の平均値とする場合、平均をとる年を指定する。下のcurrentF, futureFが指定されている場合にはこの設定は無視される
#' @param currentF start_ABC_yar_name以前に使うFのベクトル（いわゆるcurrent F）
#' @param futureF  start_ABC_yar_name以降に使うFのベクトル（いわゆるFmsy）
#' @param start_ABC_year_name HCRを有効にする年
#' @param HCR_beta HCRのbeta
#' @param HCR_Blimi HCRのBlimi
#' @param HCR_Bban HCRのBban
#' @param HCR_year_lag HCRするときにいつのタイミングのssbを参照するか.0の場合、ABC計算年のSSBを参照する。正の値1を入れると1年前のssbを参照する
#' @param HCR_beta_year betaを年によって変える場合。tibble(year=2020:2024, beta=c(1.3,1.2,1.1,1,0.9))　のようにtibble形式で与える
#' @param Pope 漁獲方程式にPopeの近似式を使うかどうか。与えない場合には、VPAのオプションが引き継がれる
#' @param fix_recruit 将来予測において再生産関係を無視して加入量を一定値で与える場合、その加入の値。list(year=2020, rec=1000)のように与える。
#' @param fix_wcatch 将来予測において漁獲量をあらかじめ決める場合
#' @param res_SR 再生産関係の推定関数 (fit.SR　of fit.SRregime) の返り値
#' @param seed_number 乱数のシードの数
#' @param resid_type 加入量の残差の発生方法。"lognormal":対数正規分布, "resample": リサンプリング, "backward": backward resampling
#' @param bias_correction 将来予測でバイアス補正をするかどうか
#' @param resample_year_range "resampling", "backward"で有効。0の場合、推定に使ったデータから計算される残差を用いる。年の範囲を入れると、対象とした年の範囲で計算される残差を用いる。
#' @param backward_duration "backward"の場合、何年で１ブロックとするか。"backward"で有効。デフォルトは5 。
#' @param recruit_intercept 将来の加入の切片。将来の加入は R=f(ssb) + intercept となる。
#' @param model_average_option model averagingをする場合のオプション. SR_matのlistとweightをlist形式で入れる(list(SR_list=list(res_SR1,res_SR2),weight=c(0.5,0.5)))
#' @param regime_shift_option res_SRにfit.SRregimeの返り値を入れた場合に指定する。将来予測で再生産関係のどのフェーズがおこるかを指定する。list(future_regime=将来のregimeの仮定。keyで指定された番号を入れる)
#' @param special_setting list形式で与えるmake_future_dataの返り値のdataと同じ名前の要素について、最後にデータをここで示されたarrayのシミュレーション1回めの値で上書きする。arrayのデータに対してのみ有効。
#' 
#' @export
#' @encoding UTF-8

make_future_data <- function(res_vpa,
                             nsim = 1000, # number of simulation
                             nyear = 50, # number of future year
                             plus_age  = NULL, # if null, equal to row number as plus group
                             future_initial_year_name = 2017,
                             start_F_year_name = 2018,
                             start_biopar_year_name=2018,
                             start_random_rec_year_name = 2018,                          
                             # biopar setting
                             waa_year, waa=NULL,
                             waa_catch_year, waa_catch=NULL,
                             waa_fun = FALSE,
                             start_waafun_year_name = start_biopar_year_name, 
                             maa_year, maa=NULL,
                             M_year, M=NULL,
                             # faa setting
                             faa_year=NULL,
                             currentF=NULL,futureF=NULL,
                             # HCR setting (not work when using TMB)
                             start_ABC_year_name=2019,
                             HCR_beta=1,
                             HCR_Blimit=-1,
                             HCR_Bban=-1,
                             HCR_year_lag=0,
                             HCR_beta_year=NULL, # tibble(year=2020:2024, beta=c(1.3,1.2,1.1,1,0.9))
                             # Other
                             Pope=res_vpa$input$Pope,
                             fix_recruit=NULL, # list(year=2020, rec=1000)
                             fix_wcatch=NULL, # list(year=2020, wcatch=2000)                          
                             # SR setting
                             res_SR=NULL,                       
                             seed_number=1,
                             resid_type="lognormal", # or resample or backward
                             bias_correction=TRUE,                          
                             resample_year_range=0, # only when "resample" or backward
                             backward_duration=5, # only when backward
                             recruit_intercept=0, # number of additional recruitment (immigration or enhancement)
                             model_average_option=NULL,
                             regime_shift_option =NULL,
                             silent=FALSE,
                             # special
                             special_setting=NULL
) 
{
  
  argname <- ls()
  input <- lapply(argname,function(x) eval(parse(text=x)))
  names(input) <- argname
  
  # define age and year
  nage <- nrow(res_vpa$naa)
  age_name    <- as.numeric(rownames(res_vpa$naa))
  recruit_age <- min(as.numeric(rownames(res_vpa$naa)))
  
  vpa_nyear           <- ncol(res_vpa$naa)
  future_initial_year <- which(colnames(res_vpa$naa)==future_initial_year_name)
  if(length(future_initial_year)  ==0) stop("future_initial_year_name is invalid.")
  total_nyear         <- future_initial_year + nyear
  allyear_name        <- min(as.numeric(colnames(res_vpa$naa)))+c(0:(total_nyear-1))
  allyear_label       <- c(rep("VPA",future_initial_year),rep("future",nyear))
  start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
  if(length(start_random_rec_year)==0) stop("start_random_rec_year_name is invalid.")
  
  tmpdata <- tibble(allyear_name, allyear_label) %>%
    group_by(allyear_label) %>%
    summarize(start=min(allyear_name),end=max(allyear_name))
  if(silent==FALSE) print(tmpdata)
  
  if(is.null(plus_age)) plus_age <- max(which(!is.na(res_vpa$naa[,future_initial_year])))
  
  # define empty array
  waa_mat <- waa_catch_mat <- M_mat <- maa_mat <- naa_mat <- faa_mat <- caa_mat <- waa_catch_mat <-
    array(0, dim=c(nage, total_nyear, nsim),
          dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim))
  class(waa_mat) <- class(M_mat) <- class(maa_mat) <- class(naa_mat) <- class(faa_mat) <- class(caa_mat) <- class(waa_catch_mat) <- "myarray"                                                                                  
  SR_mat <- array(0, dim=c(total_nyear, nsim, 15),
                  dimnames=list(year=allyear_name, nsim=1:nsim,
                                par=c("a","b","rho", #1-3
                                      "SR_type", # 4
                                      "rand_resid", # 5
                                      "deviance", #6
                                      "recruit","ssb",
                                      "intercept","sd",#9-10
                                      "bias_factor", #11
                                      "blank2","blank3","blank4","blank5")))  
  HCR_mat <- array(0, dim=c(total_nyear, nsim, 8),
                   dimnames=list(year=allyear_name, nsim=1:nsim,
                                 par=c("beta","Blimit","Bban","gamma","year_lag", #1-5
                                       "beta_gamma","wcatch","Fratio")))  # 6-8
  class(SR_mat)  <- "myarray"
  class(HCR_mat) <- "myarray"
  
  HCR_mat[,,"Blimit"] <- HCR_mat[,,"Bban"] <- -1
  HCR_mat[,,"beta"] <- HCR_mat[,,"beta_gamma"] <- 1
  
  
  # fill vpa data 
  waa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa)
  maa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$maa)
  M_mat  [,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$M)
  faa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$faa)
  naa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$naa)
  caa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$caa)
  
  waa_mat <- make_array(waa_mat, waa, waa_year, start_biopar_year_name)
  maa_mat <- make_array(maa_mat, maa, maa_year, start_biopar_year_name)
  M_mat   <- make_array(M_mat  , M  , M_year  , start_biopar_year_name)
  
  if(is.null(res_vpa$input$dat$waa.catch)){
    waa_catch_mat <- waa_mat
  }
  else{
    waa_catch_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa.catch)
    waa_catch_mat <- make_array(waa_catch_mat, waa_catch, waa_catch_year, start_biopar_year_name)
  }    
  
  # set SR parameter
  SR_mat <- set_SR_mat(res_vpa=res_vpa,
                       res_SR=res_SR,
                       SR_mat=SR_mat, seed_number=seed_number,
                       start_random_rec_year_name=start_random_rec_year_name, resid_type=resid_type,
                       resample_year_range=resample_year_range,
                       bias_correction=bias_correction,
                       recruit_intercept=recruit_intercept,
                       recruit_age=recruit_age,
                       backward_duration=backward_duration,
                       model_average_option=model_average_option,
                       regime_shift_option=regime_shift_option)
  
  # when fix recruitment
  if(!is.null(fix_recruit)) naa_mat[1,as.character(fix_recruit$year),] <- fix_recruit$rec           
  
  # set F & HCR parameter
  start_F_year <- which(allyear_name==start_F_year_name)
  start_ABC_year <- which(allyear_name==start_ABC_year_name)    
  faa_mat <- make_array(faa_mat, currentF, faa_year, start_F_year_name)
  faa_mat <- make_array(faa_mat, futureF,  faa_year, start_ABC_year_name)
  
  HCR_mat[start_ABC_year:total_nyear,,"beta"    ] <- HCR_beta
  HCR_mat[start_ABC_year:total_nyear,,"Blimit"  ] <- HCR_Blimit
  HCR_mat[start_ABC_year:total_nyear,,"Bban"    ] <- HCR_Bban    
  HCR_mat[start_ABC_year:total_nyear,,"year_lag"] <- HCR_year_lag
  
  if(!is.null(HCR_beta_year)){
    HCR_mat[as.character(HCR_beta_year$year),,"beta"] <- HCR_beta_year$beta
  }
  
  
  # when fix wcatch
  # VPA期間中のwcatchをfixするかどうか？？
  #    if(isTRUE(Pope)){
  #        tmp <- naa_mat*(1-exp(-faa_mat))*exp(-M_mat/2) * waa_catch_mat
  #    }
  #    else{
  #        tmp <- naa_mat*(1-exp(-faa_mat-M))*faa_mat/(faa_mat+M) * waa_catch_mat
  #    }
  #    HCR_mat[,,"wcatch"] <- apply(tmp,c(2,3),sum)
  HCR_mat[as.character(fix_wcatch$year), ,"wcatch"] <- fix_wcatch$wcatch
  
  waa_mat[is.na(waa_mat)] <- 0
  waa_catch_mat[is.na(waa_catch_mat)] <- 0            
  maa_mat[is.na(maa_mat)] <- 0    
  M_mat  [is.na(M_mat)  ] <- 0
  faa_mat[is.na(faa_mat)] <- 0
  naa_mat[is.na(naa_mat)] <- 0        
  caa_mat[is.na(caa_mat)] <- 0
  
  
  # set data and parameter for TMB
  tmb_data <- list(naa_mat=naa_mat,
                   caa_mat=caa_mat,
                   SR_mat=SR_mat,
                   waa_mat = waa_mat,
                   waa_catch_mat = waa_catch_mat,
                   maa_mat = maa_mat,                     
                   M_mat = M_mat,
                   faa_mat = faa_mat,
                   Pope = as.numeric(Pope),
                   total_nyear = total_nyear,                     
                   future_initial_year = future_initial_year,
                   start_ABC_year=start_ABC_year,
                   start_random_rec_year=start_random_rec_year,
                   nsim = nsim,
                   nage = nage,
                   plus_age = plus_age,
                   recruit_age = recruit_age,
                   HCR_mat = HCR_mat,
                   obj_stat = 0, # 0: mean, 1:geomean
                   objective = 0, # 0: MSY, 1: PGY, 2: percentB0 or Bempirical
                   obj_value = -1
  )
  
  if(isTRUE(waa_fun)){
    waa_rand_mat <- array(0,dim=c(nage,total_nyear,nsim),
                          dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim))
    waa_par_mat <- array(0,dim=c(nage,nsim,3),
                         dimnames=list(age=age_name, nsim=1:nsim, pars=c("sd", "b0", "b1")))
    class(waa_rand_mat) <- class(waa_par_mat) <- "myarray"
    waa_fun_year <- which(allyear_name %in% start_waafun_year_name:max(allyear_name))
    
    for(a in 1:nage){        
      for(i in 1:nsim){
        log.w <- as.numeric(log(waa_mat[a,,i]))
        log.n <- as.numeric(log(naa_mat[a,,i]))
        observed <- log.n>-Inf
        log.w <- log.w[observed]
        log.n <- log.n[observed]
        tmp <- lm(log.w~log.n)
        waa_par_mat[a,i,c("b0","b1")] <- as.numeric(tmp$coef[1:2])
        waa_par_mat[a,i,c("sd")] <- sqrt(mean(tmp$residual^2))
        waa_rand_mat[a,observed,i] <- tmp$residual
        waa_rand_mat[a,waa_fun_year,i] <- rnorm(length(waa_fun_year),-0.5*waa_par_mat[a,i,c("sd")]^2,waa_par_mat[a,i,c("sd")])
      }}
    tmb_data$waa_rand_mat <- waa_rand_mat
    tmb_data$waa_par_mat <- waa_par_mat
    tmb_data$waa_mat[,waa_fun_year,] <- 0
  }
  
  if(!is.null(special_setting)){
    set_name <- names(special_setting)
    for(i in 1:length(set_name)){
      tmb_data[[which(set_name[[i]]==tmb_data)[[1]]]][] <- special_setting[[i]][]
    }}
  
  return(tibble::lst(data=tmb_data,input=input))
}

#' 将来予測の実施関数
#'
#' @param tmb_data make_future_dataの返り値
#' @param SPR_target 目標とする\%SPR。NULL以外の値の場合、過去〜将来のそれぞれの年・シミュレーションが、目標とするF\%SPRに対して何倍にあたるか(F/Ftarget)を計算して、HCR_matの"Fratio"に入れる。HCRが生きている年については"beta_gamma"と一致するはず。
#'
#' @export
#' @encoding UTF-8

future_vpa <- function(tmb_data,
                       optim_method="none", # or "R" or "none"
                       multi_init = 1,
                       multi_lower = 0.0001,
                       multi_upper = 10,
                       objective ="MSY", # or PGY, percentB0, Bempirical
                       obj_value = 0,                         
                       obj_stat  ="mean",
                       do_MSE=NULL,
                       MSE_input_data=NULL,
                       MSE_nsim=NULL,
                       compile=FALSE,
                       output_format="new",
                       attach_input=TRUE,
                       SPRtarget=NULL,
                       calc_SPR_year_name=NULL){
  
  argname <- ls()
  input <- lapply(argname,function(x) eval(parse(text=x)))
  names(input) <- argname
  
  # conversion of estimation option
  tmb_data$objective <-
    dplyr::case_when(objective=="MSY"   ~ 0,
                     objective=="PGY"   ~ 1,
                     objective=="SSB"   ~ 2)        
  tmb_data$obj_stat <-
    dplyr::case_when(obj_stat=="mean"   ~ 0,
                     obj_stat=="median" ~ 1)
  tmb_data$obj_value <- obj_value
  
  if(optim_method=="tmb" && !is.null(tmb_data$waa_par_mat)){
    cat("Warning: waa_fun option cannot be used in TMB. The optimization is conducted by R..\n")
    optim_method <- "R"
  }
  
  if(optim_method=="tmb"){
    
    #        if(!is.null(rec_new) | !is.null(wcatch_fix)){
    #            stop("rec_new or wcatch_fix option cannot be used in \"tmb\" option\n")
    #        }
    
    # comple & load cpp file
    use_rvpa_tmb(TmbFile = "est_MSY_tmb",
                 CppDir = system.file("executable",package="frasyr"),
                 RunDir = getwd(), overwrite=compile) 
    
    objAD <- TMB::MakeADFun(tmb_data, list(x=log(multi_init)), DLL="est_MSY_tmb")
    msy_optim <- nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                        lower=list(x=log(multi_lower)),
                        upper=list(x=log(multi_upper)))#,contol=nlminb_control)
    
    multi <- as.numeric(exp(msy_optim$par))
    msy <- exp(-msy_optim$objective)
    ad_report <- objAD$report()
    
    res_future <- naming_adreport(tmb_data, ad_report)
    res_future$multi <- multi
  }
  
  #--- R
  if(optim_method=="R" | optim_method=="none"){
    
    tmb_data$do_MSE <- do_MSE
    tmb_data$MSE_input_data <- MSE_input_data
    tmb_data$MSE_nsim <- MSE_nsim
    
    R_obj_fun <- function(x, tmb_data, what_return="obj"){
      tmb_data$x <- x
      tmb_data$what_return <- what_return
      obj <- safe_call(future_vpa_R, tmb_data)
      return(obj)
    }
    
    if(optim_method=="R"){
      # nlminbを使うとrangeによって壁に当たることがある。optimizeを使ったほうがよさそうなので、optimizeに変更
      # 一方、tmbはnlminbでないとうまくいかないらしい
      #            msy_optim <- nlminb(start=0, objective=R_obj_fun, tmb_data=tmb_data,
      #                                lower=list(x=log(multi_lower)), upper=list(x=log(multi_upper)))
      #            tmb_data$x <- msy_optim$par            
      msy_optim <- optimize(R_obj_fun,lower=log(multi_lower),
                            upper=log(multi_upper), tmb_data=tmb_data)            
      tmb_data$x <- msy_optim$minimum
      
    }
    else{
      tmb_data$x <- log(multi_init)            
    }
    tmb_data$what_return <- "stat"
    res_future <- safe_call(future_vpa_R, tmb_data)
  }
  
  res_future$input <- input
  
  if(output_format!="new"){
    res_future <- format_to_old_future(res_future)
  }
  else{
    class(res_future) <- "future_new"    
  }
  if(!isTRUE(attach_input)) res_future$input <- NULL
  
  # modify results
  # remove HCR control parameter before start_ABC_year
  res_future$HCR_mat[1:(tmb_data$start_ABC_year-1),,] <- NA
  # add Fratio if needed
  if(!is.null(SPRtarget)){
    if(is.null(calc_SPR_year_name)){
       calc_SPR_year <- 1:dim(res_future$faa)[[2]]
    }
    else{
       calc_SPR_year <- which(dimnames(res_future$faa)[[2]] %in% calc_SPR_year_name)
    }
    
    for(i in calc_SPR_year){
      for(j in 1:dim(res_future$faa)[[3]]){
        if(j>2 &&
           all(res_future$faa[,i,j]==res_future$faa[,i,j-1]) &&
           all(res_future$waa[,i,j]==res_future$waa[,i,j-1]) &&
           all(res_future$waa_catch_mat[,i,j:(j-1)]==res_future$waa_catch_mat[,i,j:(j)]) &&
           all(res_future$maa[,i,j]==res_future$maa[,i,j-1]) &&
           all(res_future$M  [,i,j]==res_future$M  [,i,j-1])){
          res_future$HCR_mat[i,j,"Fratio"] <- res_future$HCR_mat[i,j-1,"Fratio"]
        }
        else{
          tmp <- res_future$naa[,i,j]>0
          res_future$HCR_mat[i,j,"Fratio"] <-
            calc_Fratio(faa=res_future$faa[tmp,i,j],
                        waa=res_future$waa[tmp,i,j],
                        maa=res_future$input$tmb_data$maa_mat[tmp,i,j],
                        M  =res_future$input$tmb_data$M_mat[tmp,i,j],
                        waa.catch=res_future$waa_catch_mat[tmp,i,j],
                        SPRtarget=SPRtarget)
        }
      }}}
  
  return(res_future)
  
  # 足りないもの
  # 推定結果の簡易的グラフ(MSY_est_plot)
}

#' @export

future_vpa_R <- function(naa_mat,
                         caa_mat,
                         SR_mat,
                         waa_mat,
                         waa_catch_mat,
                         M_mat,
                         maa_mat,                      
                         faa_mat,
                         Pope,
                         total_nyear,
                         future_initial_year,
                         start_ABC_year, # the year for estimating multiplier F & conduct HCR
                         start_random_rec_year,
                         nsim,
                         nage,
                         plus_age,
                         recruit_age,
                         obj_stat,
                         objective,
                         obj_value,
                         x,
                         what_return="obj",
                         HCR_mat,
                         do_MSE=NULL,
                         MSE_input_data=NULL,
                         MSE_nsim = NULL,
                         waa_par_mat  = NULL, # option for waa_fun
                         waa_rand_mat = NULL
){
  
  options(deparse.max.lines=10)
  
  argname <- ls()
  tmb_data <- lapply(argname,function(x) eval(parse(text=x)))
  names(tmb_data) <- argname
  
  if(isTRUE(do_MSE)){
    MSE_seed <- MSE_input_data$input$seed_number + 1        
    if(!is.null(MSE_nsim)) MSE_input_data$input$nsim <- MSE_nsim
    if( is.null(MSE_nsim)) MSE_nsim <- MSE_input_data$input$nsim
    SR_MSE <- SR_mat
    SR_MSE[,,"recruit"] <- SR_MSE[,,"ssb"] <- 0
    dimnames(SR_MSE)$par[12] <- "real_true_catch"
    dimnames(SR_MSE)$par[13] <- "pseudo_true_catch"        
  }                
  
  F_mat <- N_mat <-  naa_mat
  F_mat[] <- N_mat[] <-  0
  class(F_mat) <- class(N_mat) <- "myarray"
  
  N_mat <- naa_mat
  spawner_mat <- apply(N_mat * waa_mat * maa_mat, c(2,3) , sum)
  
  F_mat[,1:(start_ABC_year-1),] <- faa_mat[,1:(start_ABC_year-1),]
  F_mat[,start_ABC_year:total_nyear,] <- faa_mat[,start_ABC_year:total_nyear,] * exp(x)
  
  ## この時点でもwaa_funを入れる必要がある
  
  for(t in future_initial_year:total_nyear){
    
    if(!is.null(waa_par_mat)) waa_mat[,t,] <- waa_catch_mat[,t,] <- update_waa_mat(waa=waa_mat[,t,],rand=waa_rand_mat[,t,],naa=N_mat[,t,],pars_b0=waa_par_mat[,,"b0"],pars_b1=waa_par_mat[,,"b1"])
    spawner_mat[t,] <- colSums(N_mat[,t,] * waa_mat[,t,] * maa_mat[,t,])
    
    if(t>=start_random_rec_year & all(N_mat[1,t,]==0)){
      spawn_t <- t-recruit_age
      N_mat[1,t,] <- purrr::pmap_dbl(tibble(x=SR_mat[t,,"SR_type"],
                                            ssb=spawner_mat[spawn_t,],
                                            a=SR_mat[t,,"a"],b=SR_mat[t,,"b"]),
                                     function(x,ssb,a,b){
                                       fun <- list(SRF_HS,SRF_BH,SRF_RI)[[x]];
                                       fun(ssb,a,b)
                                     })
      N_mat[1,t,] <- N_mat[1,t,]*exp(SR_mat[t,,"deviance"]) + 
        SR_mat[t,,"intercept"]
      if(is.na(N_mat[1,t,1])) stop("Error: Recruitment cannot be estimated correctly...")
      if(!is.null(waa_par_mat)) waa_mat[1,t,] <- waa_catch_mat[1,t,] <- update_waa_mat(waa=waa_mat[1,t,],rand=waa_rand_mat[1,t,],naa=N_mat[1,t,],pars_b0=waa_par_mat[1,,"b0"],pars_b1=waa_par_mat[1,,"b1"])
    }
    
    if(t>=start_ABC_year){
      # harvest control rule
      ssb_tmp <- spawner_mat[cbind(t-HCR_mat[t,,"year_lag"],1:nsim)]
      HCR_mat[t,,"beta_gamma"] <- HCR_default(ssb_tmp, HCR_mat[t,,"Blimit"],
                                              HCR_mat[t,,"Bban"], HCR_mat[t,,"beta"])
      F_mat[,t,] <- sweep(F_mat[,t,],2,HCR_mat[t,,"beta_gamma"],FUN="*")
    }
    
    if(isTRUE(do_MSE) && t>=start_ABC_year){
      MSE_dummy_data <- safe_call(make_future_data,MSE_input_data$input)$data
      MSE_dummy_data <- MSE_dummy_data %>%
        purrr::list_modify(future_initial_year   = t-2,
                           start_random_rec_year = t-1,
                           start_ABC_year        = t,
                           total_nyear           = t,
                           x                     = 0, # = 1 in normal scale
                           what_return           = "stat",
                           do_MSE                = FALSE)
      for(i in 1:nsim){
        #                cat(t,"-",i,"\n")                
        MSE_dummy_data$naa_mat[] <-  N_mat[,,i] # true dynamics
        MSE_dummy_data$naa_mat[,(t-1):t,] <-  0 # estiamted as future
        MSE_dummy_data$faa_mat[,1:(t-1),] <-  F_mat[,1:(t-1),i] # we know true F even in future
        MSE_dummy_data$faa_mat[,t,] <-  MSE_input_data$data$faa[,t,i] # alpha in ABC year is depends on future SSB
        MSE_dummy_data$waa_mat[] <-  waa_mat[,,i] # in case
        MSE_dummy_data$waa_catch_mat[] <-  waa_catch_mat[,,i] # in case
        MSE_dummy_data$maa_mat[] <-  maa_mat[,,i] # in case
        MSE_dummy_data$M_mat[]   <-  M_mat[,,i] # in case
        for(k in 1:MSE_nsim){
          MSE_dummy_data$SR_mat[,k,]  <- SR_mat[,i,]
          MSE_dummy_data$SR_mat[,k,"ssb"]  <- spawner_mat[,i] # true ssb
          MSE_dummy_data$SR_mat[,k,"recruit"]  <- N_mat[1,,i] # true recruit         
        }
        # re-calculate past deviance and produce random residual in future
        MSE_dummy_data$SR_mat <- 
          set_SR_mat(res_vpa   = NULL, # past deviande is calculated by true ssb
                     res_SR    = MSE_input_data$input$res_SR,
                     SR_mat    = MSE_dummy_data$SR_mat,
                     seed_number=MSE_seed,
                     start_random_rec_year_name = dimnames(naa_mat)[[2]][t-1],
                     recruit_age = recruit_age,
                     resid_type                 = MSE_input_data$input$resid_type,
                     resample_year_range        = dimnames(naa_mat)[[2]][1]:dimnames(naa_mat)[[2]][t-2],
                     backward_duration          = MSE_input_data$input$backward_duration,
                     bias_correction            = MSE_input_data$input$bias_correction,
                     recruit_intercept          = MSE_input_data$input$recruit_intercept,
                     model_average_option       = MSE_input_data$input$model_average_option,
                     regime_shift_option        = MSE_input_data$input$regime_shift_option)
        
        res_tmp <- safe_call(future_vpa_R,MSE_dummy_data) # do future projection
        #                if(t>55) browser()
        HCR_mat[t,i,"wcatch"] <- mean(apply(res_tmp$wcaa[,t,],2,sum)) # determine ABC in year t
        SR_MSE[t,i,"recruit"] <- mean(res_tmp$naa[1,t,])
        SR_MSE[t,i,"ssb"]     <- mean(res_tmp$SR_mat[t,,"ssb"])
        if(Pope==1){
          SR_MSE[t,i,"real_true_catch"] <- sum(N_mat[,t,i]*(1-exp(-F_mat[,t,i]))*exp(-M_mat[,t,i]/2) * waa_catch_mat[,t,i])
        }
        else{
          SR_MSE[t,i,"real_true_catch"] <- sum(N_mat[,t,i]*(1-exp(-F_mat[,t,i]-M_mat[,t,i]))*F_mat[,t,i]/(F_mat[,t,i]+M_mat[,t,i]) * waa_catch_mat[,t,i])
        }
        
        MSE_seed <- MSE_seed+1
      }
    }        
    
    if(sum(HCR_mat[t,,"wcatch"])>0){
      F_max_tmp <- apply(F_mat[,t,],2,max)
      #            saa.tmp <- sweep(F_mat[,t,],2,F_max_tmp,FUN="/")
      fix_catch_multiplier <- purrr::map_dbl(which(F_max_tmp>0),
                                             function(x) caa.est.mat(N_mat[,t,x],F_mat[,t,x],#saa.tmp[,x],
                                                                     waa_catch_mat[,t,x],M_mat[,t,x],
                                                                     HCR_mat[t,x,"wcatch"],
                                                                     set_max1=FALSE,
                                                                     Pope=as.logical(Pope))$x)
      F_mat[,t,which(F_max_tmp>0)] <- sweep(F_mat[,t,which(F_max_tmp>0)],#saa.tmp[,which(F_max_tmp>0)],
                                            2, fix_catch_multiplier, FUN="*")
      HCR_mat[t,which(F_max_tmp>0),"beta_gamma"] <- HCR_mat[t,which(F_max_tmp>0),"beta_gamma"]*fix_catch_multiplier
    }
    
    if(t<total_nyear){
      # forward calculation                 
      for(iage in 1:(plus_age-1)) {
        N_mat[iage+1,t+1,] <- N_mat[iage,t,]*exp(-M_mat[iage,t,]-F_mat[iage,t,])
      }
      N_mat[plus_age,t+1,] <- N_mat[plus_age,t+1,] + N_mat[plus_age,t,]*exp(-M_mat[plus_age,t,]-F_mat[plus_age,t,])
      if(!is.null(waa_par_mat)) waa_mat[,t,] <- waa_catch_mat[,t,] <- update_waa_mat(waa=waa_mat[,t,],rand=waa_rand_mat[,t,],naa=N_mat[,t,],pars_b0=waa_par_mat[,,"b0"],pars_b1=waa_par_mat[,,"b1"])                        
    }
  }
  
  if(Pope==1){
    wcaa_mat <- N_mat*(1-exp(-F_mat))*exp(-M_mat/2) * waa_catch_mat
  }
  else{
    wcaa_mat <- N_mat*(1-exp(-F_mat-M_mat))*F_mat/(F_mat+M_mat) * waa_catch_mat
  }
  
  if(isTRUE(do_MSE)){
    F_pseudo_mat <- MSE_input_data$data$faa
    beta_gamma <- HCR_default(spawner_mat,
                              MSE_input_data$data$HCR_mat[,,"Blimit"],
                              MSE_input_data$data$HCR_mat[,,"Bban"],
                              MSE_input_data$data$HCR_mat[,,"beta"])
    F_pseudo_mat[] <- sweep(F_pseudo_mat,c(2,3),beta_gamma,FUN="*")
    
    if(Pope==1){
      wcaa_tmp <- N_mat*(1-exp(-F_pseudo_mat))*exp(-M_mat/2) * waa_catch_mat
    }
    else{
      wcaa_tmp <- N_mat*(1-exp(-F_pseudo_mat-M_mat))*F_pseudo_mat/
        (F_pseudo_mat+M_mat) * waa_catch_mat
    }
    SR_MSE[,,"pseudo_true_catch"] <- apply(wcaa_tmp, c(2,3), sum)
  }
  
  if(objective<2){
    last_catch <- colSums(wcaa_mat[,total_nyear,])
    if(obj_stat==0) obj <- mean(last_catch)
    if(obj_stat==1) obj <- geomean(last_catch)
  }
  else{
    if(obj_stat==0) obj <- mean(spawner_mat[total_nyear,])
    if(obj_stat==1) obj <- geomean(spawner_mat[total_nyear,])
  }
  
  
  {if(objective==0) obj <- -log(obj)
    else{
      obj <- (log(obj/obj_value))^2
    }
  }
  
  if(what_return=="obj")  return(obj)
  if(what_return=="stat"){
    tmb_data$SR_mat[,,"ssb"]  <- spawner_mat
    tmb_data$SR_mat[,,"recruit"]  <- N_mat[1,,]
    res <- list(naa=N_mat, wcaa=wcaa_mat, faa=F_mat, SR_mat=tmb_data$SR_mat,
                HCR_mat=HCR_mat,multi=exp(x),waa=waa_mat, waa_catch_mat=waa_catch_mat)
    if(isTRUE(do_MSE)) res$SR_MSE <- SR_MSE
    return(res)
  }
  
}

# utility functions for stock recruitment relationship ----

#'
#' 将来予測用の再生産関係の設定を行う関数
#'
#' 再生産関係をres_SRで与えると、res_vpaを見ながら残差を再計算したのち、start_random_rec_year_name以降の加入のdeviationを計算しSR_mat[,,"deviance"]に入れる。
#'
#' @param res_vpa VPAの推定結果
#' @param res_SR 再生産関係の推定結果
#' @param SR_mat 将来予測用の再生産関係パラメータが格納する３次元行列
#' @param seed_number シード番号
#' @param start_random_rec_year_name ランダム加入を仮定する最初の年
#' @param resid_type 残差の発生パターン；対数正規分布は"lognormal"、単純リサンプリングは"resampling"、backward-resamplingは"backward"
#' @param resample_year_range "resampling", "backward"で有効。0の場合、推定に使ったデータから計算される残差を用いる。年の範囲を入れると、対象とした年の範囲で計算される残差を用いる。
#' @param backward_duration "backward"の場合、何年で１ブロックとするか。"backward"で有効。デフォルトは5 。
#' @param model_average_option model averagingをする場合のオプション. SR_matのlistとweightをlist形式で入れる(list(SR_list=list(res_SR1,res_SR2),weight=c(0.5,0.5))). 上で設定されたres_SRは使われない.
#' @param regime_shift_option レジームシフトを仮定する場合のオプション. この場合, res_SRにはfit.SRregimeの結果オブジェクトを入れる. オプションの設定は list(future_regime=将来のregimeの仮定。keyで指定された番号を入れる)
#' 
#' @export
#' @encoding UTF-8

set_SR_mat <- function(res_vpa=NULL,
                       start_random_rec_year_name,                       
                       SR_mat,                       
                       res_SR,
                       seed_number,
                       resid_type="lognormal",
                       bias_correction=TRUE,                       
                       resample_year_range=0,
                       backward_duration=5,
                       recruit_intercept=0,
                       recruit_age=0,
                       model_average_option=NULL,
                       regime_shift_option=NULL
){
  
  allyear_name <- dimnames(SR_mat)[[1]]
  start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
  random_rec_year_period <- (start_random_rec_year):length(allyear_name)
  
  if(!resid_type%in%c("lognormal","resample","backward")) stop("resid_type is invalid.")
  
  # define SR function
  if(res_SR$input$SR=="HS"){
    SR_mat[,,"SR_type"] <- 1
    SRF <- SRF_HS
  }
  if(res_SR$input$SR=="BH"){
    SR_mat[,,"SR_type"] <- 2
    SRF <- SRF_BH        
  }
  if(res_SR$input$SR=="RI"){
    SR_mat[,,"SR_type"] <- 3
    SRF <- SRF_RI                
  }
  
  # define SR parameter
  if(is.null(regime_shift_option)){
    SR_mat[,,"a"] <- res_SR$pars$a
    SR_mat[,,"b"] <- res_SR$pars$b
    SR_mat[,,"sd"] <- res_SR$pars$sd        
  }
  else{
    regime_data <- res_SR$regime_resid %>%
      left_join(res_SR$regime_pars, by="regime") %>%
      bind_cols(res_SR$input$SRdata)
    SR_mat[as.character(regime_data$year),,"a"] <- regime_data$a
    SR_mat[as.character(regime_data$year),,"b"] <- regime_data$b
    SR_mat[as.character(regime_data$year),,"sd"] <- regime_data$sd
    future_regime_par <- res_SR$regime_pars %>% dplyr::filter(regime==regime_shift_option$future_regime)
    SR_mat[random_rec_year_period,,"a"] <- future_regime_par$a 
    SR_mat[random_rec_year_period,,"b"] <- future_regime_par$b 
    SR_mat[random_rec_year_period,,"sd"] <- future_regime_par$sd

    # どのレジームに属するか明示的に指定されないケースが出てくる。
    # そういう場合は将来予測のレジームのパラメータを使うようにする→要改善
    missing_year <- which(!(1:dim(SR_mat)[[1]] %in% c(which(allyear_name %in% as.character(regime_data$year)),random_rec_year_period)))

    if(length(missing_year)>0){
        SR_mat[missing_year,,"a"] <- future_regime_par$a 
        SR_mat[missing_year,,"b"] <- future_regime_par$b 
        SR_mat[missing_year,,"sd"] <- future_regime_par$sd        
    }

    res_SR$pars$rho <- 0
  }
  SR_mat[,,"rho"] <- res_SR$pars$rho            
  SR_mat[,,"intercept"] <- recruit_intercept
  
  if(!is.null(res_vpa)){
    SR_mat[1:(start_random_rec_year-1),,"ssb"] <- as.numeric(colSums(res_vpa$ssb,na.rm=T))[1:(start_random_rec_year-1)]
    SR_mat[1:(start_random_rec_year-1),,"recruit"] <- as.numeric(res_vpa$naa[1,1:(start_random_rec_year-1)])
  }
  
  recruit_range <- (recruit_age+1):(start_random_rec_year-1)
  ssb_range     <- 1:(start_random_rec_year-1-recruit_age)    
  
  # re-culcurate recruitment deviation    
  SR_mat[recruit_range,,"deviance"] <- SR_mat[recruit_range,,"rand_resid"] <- 
    log(SR_mat[recruit_range,,"recruit"]) -
    log(SRF(SR_mat[ssb_range,,"ssb"],SR_mat[recruit_range,,"a"],SR_mat[recruit_range,,"b"]))
  
  # define future recruitment deviation
  set.seed(seed_number)
  nsim <- dim(SR_mat)[[2]]
  
  if(resid_type=="lognormal"){
    if(isTRUE(bias_correction)){
      #            sd_with_AR <- sqrt(res_SR$pars$sd^2/(1-res_SR$pars$rho^2))
      #            bias_factor <- 0.5* sd_with_AR^2
      sd_with_AR <- sqrt(SR_mat[,,"sd"]^2/(1-SR_mat[,,"rho"]^2))
      SR_mat[,,"bias_factor"] <- 0.5 * sd_with_AR^2
      SR_mat[-random_rec_year_period,,"bias_factor"] <- 0
    }
    else{
      #            bias_factor <- 0
      SR_mat[,,"bias_factor"] <- 0
    }
    tmp_SR <- t(SR_mat[random_rec_year_period,,"rand_resid"])
    tmp_SR[] <- rnorm(nsim*length(random_rec_year_period), mean=0, sd=res_SR$pars$sd)
    SR_mat[random_rec_year_period,,"rand_resid"] <- t(tmp_SR)
    
    for(t in random_rec_year_period){
      SR_mat[t, ,"deviance"] <- SR_mat[t-1, ,"deviance"]*SR_mat[t,,"rho"] + SR_mat[t, ,"rand_resid"] 
    }
    SR_mat[random_rec_year_period,,"deviance"] <- SR_mat[random_rec_year_period,,"deviance"] - SR_mat[random_rec_year_period,,"bias_factor"]
  }
  
  if(resid_type=="resample" | resid_type=="backward"){
    # 推定された残差をそのまま使う
    if(resample_year_range==0){
      #            sampled_residual <- res_SR$resid[res_SR$input$w==1]
      #            if(isTRUE(bias_correction)) bias_factor <- log(mean(exp(sampled_residual))) else bias_factor <- 0
      #            SR_mat[random_rec_year_period,,"rand_resid"] <- sample(sampled_residual, nsim*length(random_rec_year_period), replace=TRUE)
      #            SR_mat[random_rec_year_period,,"deviance"] <- SR_mat[random_rec_year_period,,"rand_resid"]-bias_factor
      resample_year_range <- sort(res_SR$input$SRdata$year[res_SR$input$w==1])
    }
    
    sampled_residual <- SR_mat[as.character(resample_year_range),,"rand_resid"]
    if(isTRUE(bias_correction)){
      #            bias_factor <- log(colMeans(exp(sampled_residual)))
      SR_mat[random_rec_year_period,,"bias_factor"] <- rep(log(colMeans(exp(sampled_residual))),
                                                           each=length(random_rec_year_period))
    }
    else{
      #            bias_factor <- rep(0,ncol(sampled_residual))
      SR_mat[random_rec_year_period,,"bias_factor"] <- 0            
    }
    for(i in 1:ncol(sampled_residual)){
      if(resid_type=="resample"){
        SR_mat[random_rec_year_period,i,"rand_resid"] <- sample(sampled_residual[,i], length(random_rec_year_period), replace=TRUE)
      }
      if(resid_type=="backward"){
        SR_mat[random_rec_year_period,i,"rand_resid"] <- sample_backward(sampled_residual[,i], length(random_rec_year_period), backward_duration)
      }            
      SR_mat[random_rec_year_period,i,"deviance"] <- SR_mat[random_rec_year_period,i,"rand_resid"]-SR_mat[random_rec_year_period,i,"bias_factor"]
    }
  }
  
  if(!is.null(model_average_option)){
    weight <- arrange_weight(model_average_option$weight,nsim)
    SR_mat <- average_SR_mat(res_vpa     = res_vpa,
                             res_SR_list = model_average_option$SR_list,
                             range_list  = weight,
                             SR_mat      = SR_mat,
                             seed_number = seed_number+1,
                             start_random_rec_year_name=start_random_rec_year_name,
                             resid_type  = resid_type,
                             recruit_intercept=recruit_intercept,
                             resample_year_range=resample_year_range,                                 
                             recruit_age = recruit_age,
                             regime_shift_option = regime_shift_option,
                             bias_correction = bias_correction
    )
  }
  return(SR_mat)
}

#' @export
SRF_HS <- function(x,a,b) ifelse(x>b,b*a,x*a)

#' @export
SRF_BH <- function(x,a,b) a*x/(1+b*x)

#' @export
SRF_RI <- function(x,a,b) a*x*exp(-b*x)

#'
#' 将来予測用の三次元行列（年齢×年×シミュレーション）を与えられたら, pars.yearで指定された期間のパラメータを平均するか、parで指定されたパラメータを、year_replace_future以降の年で置き換える
#'
#' @param d3_mat 将来予測用の３次元行列
#' @param pars 置き換えるべき生物パラメータ
#' @param pars.year この期間の生物パラメータを平均して、将来のパラメータとする
#' @param year_replace_future 生物パラメータを置き換える最初の年
#' @encoding UTF-8
#' @export
#' 

make_array <- function(d3_mat, pars, pars.year, year_replace_future){
  if(length(dim(pars))==3){
    return(pars)
  }
  else{
    years <- dimnames(d3_mat)[[2]]
    if(is.null(pars)){
      pars.future <- rowMeans(d3_mat[,years%in%pars.year,1])
    }
    else{
      if(length(pars)==dim(d3_mat)[[1]]) pars.future <- pars
      else stop("length of parameter is different from what is expected.")
    }
    d3_mat[,which(year_replace_future==years):length(years),] <- pars.future
    if(sum(year_replace_future==years)==0) stop("year_replace_future is invalid.")
    
    return(d3_mat)
  }
}


#' @param weight
#' @param nsim
#' @export
arrange_weight <- function(weight, nsim){
  weight <- weight / sum(weight)
  weight <- round(cumsum(weight) * nsim)
  weight2 <- c(1,weight[-length(weight)]+1)
  purrr::map(1:length(weight),function(x) weight2[x]:weight[x])
}

#'
#' モデル平均的な再生産関係を与える
#'
#' @param res_vpa VPAの推定結果
#' @param res_SR_list 再生産関係の推定結果のリスト
#' @param range_list 
#' @param SR_mat 将来予測用の再生産関係パラメータを格納する３次元行列
#' @param seed_number シード番号
#' @param start_random_rec_year_name ランダム加入を仮定する最初の年
#' @param resid_type 残差の発生パターン；対数正規分布は"lognormal"、単純リサンプリングは"resampling"
#' @param resample_year_range 0の場合、推定に使ったデータから計算される残差を用いる。年の範囲を入れると、対象とした年の範囲で計算される残差を用いる
#' 
#' @export
#' @encoding UTF-8

average_SR_mat <- function(res_vpa,
                           res_SR_list,
                           range_list=list(1:500,501:1000),
                           SR_mat,
                           seed_number,
                           start_random_rec_year_name,
                           recruit_age,
                           resid_type="lognormal",
                           resample_year_range=0,
                           regime_shift_option=NULL,
                           recruit_intercept=0,
                           bias_correction=TRUE){
  
  allyear_name <- dimnames(SR_mat)[[1]]
  start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
  random_rec_year_period <- (start_random_rec_year):length(allyear_name)
  
  for(i in 1:length(res_SR_list)){
    SR_mat_tmp <- set_SR_mat(res_vpa=res_vpa,
                             start_random_rec_year_name, 
                             SR_mat=SR_mat,
                             res_SR=res_SR_list[[i]],
                             seed_number=seed_number+i,
                             resid_type=resid_type,        
                             resample_year_range=resample_year_range,
                             recruit_age=recruit_age,
                             recruit_intercept=recruit_intercept,
                             bias_correction=bias_correction,
                             regime_shift_option=regime_shift_option)
    SR_mat[,as.character(range_list[[i]]),] <-
      SR_mat_tmp[,range_list[[i]],]
  }
  
  return(SR_mat)
}

#' 過去にさかのぼってブロックサンプリングをおこなう
#' 
#' @param residual リサンプリングする残差
#' @param n 将来にわたって何年分のリサンプリング残差を作るか
#' @param duration 1ブロックの年の長さ
#' 
#' @examples
#'
#' set.seed(1)
#' res <- sample_backward(rep(1:5,each=5), 30, 5)
#' apply(matrix(res,5,6),2,min)
#' 
#' @export
#' @encoding UTF-8
#' 

sample_backward <- function(residual, n, duration){
  residual_rev <- rev(residual)
  nblock <- floor(length(residual_rev)/duration)
  block <- (1:(nblock))*(duration)
  block2 <- block-(duration)+1
  block[length(block)] <- length(residual_rev)
  block.list <- purrr::map(1:length(block),function(x) residual_rev[block2[x]:block[x]])
  # calculate sampling probability in the case of tail block (different length)
  block.probability <- sapply(block.list,length)
  block.probability <- block.probability/sum(block.probability)
  resid_future <- ceiling(1:n/duration)
  resid_future <- ifelse(resid_future>nblock,nblock,resid_future)
  
  for(i in 1:n){
    # block is resampled according to the data length of each block
    if(i%%duration==1){
      block_select <- sample(x=1:resid_future[i],size=1,prob=block.probability[1:resid_future[i]])
    }
    # iごとに指定されたブロックをリサンプリング＝重複あり
    resid_future[i] <- sample(x=block.list[[block_select]],size=1)        
  }
  return(resid_future)
}


#' @export
print.myarray <- function(x) cat("array :", dim(x),"\n")


# utility functions for other parts -----

#'
#' @export
#' @encoding UTF-8

naming_adreport <- function(tmb_data, ad_report){
  #    ssb <- ad_report$spawner_mat
  wcaa <- ad_report$catch_mat
  naa <- ad_report$N_mat
  faa <- ad_report$F_mat
  dimnames(wcaa) <- dimnames(naa) <- dimnames(faa) <-
    dimnames(tmb_data$naa)
  class(wcaa) <- class(naa) <- class(faa) <- "myarray"
  
  tmb_data$SR_mat[,,"ssb"] <- ad_report$spawner_mat
  tmb_data$SR_mat[,,"recruit"] <- ad_report$N_mat[1,,]
  
  return(list(wcaa=wcaa, naa=naa, faa=faa,
              SR_mat        = tmb_data$SR_mat,
              HCR_mat       = tmb_data$HCR_mat,
              waa           = tmb_data$waa_mat,
              waa_catch_mat = tmb_data$waa_catch_mat))
}


#'
#' @export
#' @encoding UTF-8

trace_future <- function(tmb_data,
                         trace.multi=c(seq(from=0,to=0.9,by=0.1),1,
                                       seq(from=1.1,to=2,by=0.1),3:5,7,20,100),
                         ncore=0){
  #    R_obj_fun <- function(x, tmb_data, what_return="obj"){
  #        tmb_data$x <- x
  #        tmb_data$what_return <- what_return
  #        obj <- do.call(future_vpa_R, tmb_data)
  #        return(obj)
  #    }    
  
  if(ncore==0){
    trace <- purrr::map_dfr(trace.multi, function(x)
      #       R_obj_fun(x=x, tmb_data = tmb_data, what_return="stat") %>%
      future_vpa(tmb_data,
                 optim_method="none",
                 multi_init=x,
                 multi_lower=x,
                 multi_upper=x) %>%
        get.stat(use_new_output=TRUE) %>% as_tibble())
  }
  else{
    library(foreach)
    #        library(doParallel)
    cl <- parallel::makeCluster(ncore, type="FORK")
    doParallel::registerDoParallel(cl)
    trace <- foreach::foreach(x=trace.multi,.combine="rbind")%dopar%{
      future_vpa(tmb_data,
                 optim_method="none",
                 multi_init=x,
                 multi_lower=x,
                 multi_upper=x) %>%
        get.stat(use_new_output=TRUE) %>% as_tibble()
    }
    parallel::stopCluster(cl)
  }
  
  return(trace)
}

#'
#' @export
#' @encoding UTF-8

get_summary_stat <- function(all.stat){
  
  sumvalue <- all.stat %>% as_tibble %>%
    mutate(SSB2SSB0=all.stat$ssb.mean/all.stat$ssb.mean[2]) %>%
    select(RP_name,ssb.mean,SSB2SSB0,biom.mean,U.mean,catch.mean,catch.CV,Fref2Fcurrent)
  colnames(sumvalue) <- c("RP_name","SSB","SSB2SSB0","B","U","Catch","Catch.CV","Fref/Fcur")
  
  sumvalue <- bind_cols(sumvalue,all.stat[,substr(colnames(all.stat),1,1)=="F"])
  
  sumvalue$RP.definition <- NA
  sumvalue$RP.definition[sumvalue$RP_name=="MSY"]           <- "Btarget0"
  sumvalue$RP.definition[sumvalue$RP_name=="PGY_0.6_lower"] <- "Blimit0"    
  sumvalue$RP.definition[sumvalue$RP_name=="PGY_0.1_lower"] <- "Bban0"
  sumvalue <- sumvalue %>% select(1,ncol(sumvalue),2:(ncol(sumvalue)-1))    
  
  Fvector <- select(sumvalue,num_range("F",0:40))
  
  #    sumvalue$perSPR <- NA
  #    for(i in 1:nrow(Fvector)){
  #        sumvalue$perSPR[i] <- calc_perspr(input.list[[1]],Fvector[i,])
  #    }
  
  tibble::lst(sumvalue,Fvector)
  
}

#'
#' @export
#' @encoding UTF-8

format_to_old_future <- function(fout){
  fout_old <- fout[c("naa","faa","multi","input","waa")]
  #    fout_old$waa       <- fout$input$tmb_data$waa_mat
  fout_old$waa.catch <- fout$waa_catch_mat        
  fout_old$maa       <- fout$input$tmb_data$maa_mat
  fout_old$M         <- fout$input$tmb_data$M_mat        
  fout_old$vssb      <- apply(fout$naa * fout_old$waa * fout_old$maa, c(2,3), sum, na.rm=T)
  fout_old$vbiom     <- apply(fout$naa * fout_old$waa, c(2,3),sum, na.rm=T)
  fout_old$vwcaa     <- apply(fout$wcaa,c(2,3),sum, na.rm=T)
  fout_old$currentF  <- fout$faa[,fout$input$tmb_data$start_ABC_year-1,1]
  fout_old$futureF   <- fout$faa[,fout$input$tmb_data$start_ABC_year,1]
  fout_old$finalmeanF<- fout$faa[,dim(fout$faa)[[2]],] %>% apply(1,mean) # newly define
  fout_old$caa       <- fout$wcaa/fout_old$waa
  fout_old$multi     <- fout$multi
  fout_old$recruit   <- fout$SR_mat[,,"recruit"]
  fout_old$beta_gamma     <- fout$HCR_mat[,,"beta_gamma"]
  fout_old$alpha     <- fout$HCR_mat[,,"beta_gamma"]
  fout_old$Fratio     <- fout$HCR_mat[,,"Fratio"]        
  return(fout_old)
}


#' 
#' do.callのsafe版
#'
#' do.callで与えたリストの中にfuncで定義されていないものが混じっていた場合に、実際にdo.callを呼び出す前にerorrを吐いて関数をストップさせる。非常に大きいオブジェクトを与えていながらdo.callで上記の場面でエラーが出ると、長時間Rがフリーズするのを避けるため。force=TRUEにすると、func内で定義されていない引数はリストから除外してdo.callを実行する.
#'
#' @export
#' @encoding UTF-8

safe_call <- function(func,args,force=FALSE,...){
  argname <- names(formals(func))
  check_argument <- names(args) %in% argname
  if(sum(check_argument==FALSE)>0){
    if(force==FALSE){
      stop(paste(names(args)[check_argument==FALSE]), " is not used in func\n")
    }
    else{
      args <- args[check_argument==TRUE]
    }
  }
  return(do.call(func,args,...))
}

#'
#' @export
#' @encoding UTF-8

HCR_default <- function(ssb, Blimit, Bban, beta){
  beta_gamma <- beta
  tmp <- ssb < Blimit
  beta_gamma[tmp] <- beta[tmp]*(ssb[tmp]-Bban[tmp])/(Blimit[tmp]-Bban[tmp])
  beta_gamma[beta_gamma < 0] <- 0
  return(beta_gamma)
}


#'
#' 
#'  

if(0){
  update_waa_mat(waa=matrix(c(rep(0,5),10),2,3),
                 naa=matrix(1,2,3),
                 rand=matrix(0,2,3),
                 pars=array(c(0,0,log(1),log(1),2,2),dim=c(2,3),                                                            dimnames=list(age=1:2, pars=c("sd", "b0", "b1"))))
  #   pars
  #age sd b0 b1
  #  1  0  0  2
  #  2  0  0  2
}

update_waa_mat <- function(waa,rand,naa,pars_b0,pars_b1){
  waa_tmp <- exp(pars_b0+pars_b1*log(naa)+rand)
  waa[waa==0 & naa>0] <- waa_tmp[waa==0 & naa>0]
  waa
}



