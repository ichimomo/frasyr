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
#'
#' @export
#' 

make_array <- function(d3_mat, pars, pars.year, year_replace_future){
    years <- dimnames(d3_mat)[[2]]
    if(is.null(pars)){
        pars.future <- rowMeans(d3_mat[,years%in%pars.year,1])
    }
    else{
        if(length(pars)==dim(d3_mat)[[1]]) pars.future <- pars
        else stop("length of parameter is bad..")
    }
    d3_mat[,which(year_replace_future==years):length(years),] <- pars.future
    d3_mat
}

#'
#' 対数正規分布の残差分布を作る関数
#'
#' 再生産関係をres_SRで与えると、res_vpaを見ながら残差を再計算したのち、start_random_rec_year_name以降の加入のdeviationを計算しSR_mat[,,"deviance"]に入れる。
#'
#' @param res_vpa VPAの推定結果
#' @param res_SR 再生産関係の推定結果
#' @param SR_mat 将来予測用の再生産関係パラメータが格納する３次元行列
#' @param seed_number シード番号
#' @param start_random_rec_year_name ランダム加入を仮定する最初の年
#' @param resid_type 残差の発生パターン；対数正規分布は"lognormal"、単純リサンプリングは"resampling"
#' @param resample_year_range 0の場合、推定に使ったデータから計算される残差を用いる。年の範囲を入れると、対象とした年の範囲で計算される残差を用いる。
#' 
#' @export
#' 

set_SR_mat <- function(res_vpa,
                       res_SR,
                       SR_mat,
                       seed_number,
                       start_random_rec_year_name,
                       resid_type="lognormal",
                       resample_year_range=0,
                       intercept=0,
                       bias_correction=TRUE){
    
    allyear_name <- dimnames(SR_mat)[[1]]
    start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
    random_rec_year_period <- (start_random_rec_year):length(allyear_name)
    
  
    # define SR function
    if(res_SR$input$SR=="HS"){
        SR_mat[random_rec_year_period,,"SR_type"] <- 1
        SRF <- SRF_HS
    }
    if(res_SR$input$SR=="BH"){
        SR_mat[random_rec_year_period,,"SR_type"] <- 2
        SRF <- SRF_BH        
    }
    if(res_SR$input$SR=="RI"){
        SR_mat[random_rec_year_period,,"SR_type"] <- 3
        SRF <- SRF_RI                
    }

    # define SR parameter
    SR_mat[,,"a"] <- res_SR$pars$a
    SR_mat[,,"b"] <- res_SR$pars$b
    SR_mat[,,"rho"] <- res_SR$pars$rho
    SR_mat[,,"intercept"] <- intercept

    # re-culcurate recruitment deviation
    SR_mat[1:(start_random_rec_year-1),,"ssb"] <- as.numeric(colSums(res_vpa$ssb))[1:(start_random_rec_year-1)]
    SR_mat[1:(start_random_rec_year-1),,"recruit"] <- as.numeric(res_vpa$naa[1,1:(start_random_rec_year-1)])
    SR_mat[1:(start_random_rec_year-1),,"deviance"] <- SR_mat[1:(start_random_rec_year-1),,"rand_resid"] <- 
        log(SR_mat[1:(start_random_rec_year-1),,"recruit"]) -
        log(SRF(SR_mat[1:(start_random_rec_year-1),,"ssb"],
                SR_mat[1:(start_random_rec_year-1),,"a"],
                SR_mat[1:(start_random_rec_year-1),,"b"]))

    # define future recruitment deviation
    set.seed(seed_number)
    if(resid_type=="lognormal"){
        if(isTRUE(bias_correction)){
            sd_with_AR <- sqrt(res_SR$pars$sd^2/(1-res_SR$pars$rho^2))
            bias_factor <- 0.5* sd_with_AR^2
        }
        else{
            bias_factor <- 0
        }
        tmp_SR <- t(SR_mat[random_rec_year_period,,"rand_resid"])
        tmp_SR[] <- rnorm(dim(SR_mat)[[2]]*length(random_rec_year_period), mean=0, sd=res_SR$pars$sd)
        SR_mat[random_rec_year_period,,"rand_resid"] <- t(tmp_SR)
   
        for(t in random_rec_year_period){
            SR_mat[t, ,"deviance"] <- SR_mat[t-1, ,"deviance"]*SR_mat[t,,"rho"] + SR_mat[t, ,"rand_resid"] 
        }
        SR_mat[random_rec_year_period,,"deviance"] <- SR_mat[random_rec_year_period,,"deviance"] - bias_factor
    }
    if(resid_type=="resample"){
        # 推定された残差をそのまま使う
        if(resample_year_range==0){
            sampled_residual <- res_SR$resid[res_SR$input$w==1]
            if(isTRUE(bias_correction)) bias_factor <- log(mean(exp(sampled_residual))) else bias_factor <- 0
            SR_mat[random_rec_year_period,,"rand_resid"] <- sample(sampled_residual, dim(SR_mat)[[2]]*length(random_rec_year_period), replace=TRUE)
            SR_mat[random_rec_year_period,,"deviance"] <- SR_mat[random_rec_year_period,,"rand_resid"]-bias_factor
        }
        else{ # VPA結果から計算しなおした残差を使う
            sampled_residual <- SR_mat[as.character(resample_year_range),,"rand_resid"]
            if(isTRUE(bias_correction)){
                bias_factor <- log(colMeans(exp(sampled_residual)))
            }
            else{
                bias_factor <- rep(0,ncol(sampled_residual))
                }
            for(i in 1:ncol(sampled_residual)){
                SR_mat[random_rec_year_period,i,"rand_resid"] <- sample(sampled_residual[,i], length(random_rec_year_period), replace=TRUE)
                SR_mat[random_rec_year_period,i,"deviance"] <- SR_mat[random_rec_year_period,,"rand_resid"]-bias_factor[i]
            }
        }
    }
    if(resid_type=="backward"){
        ## coming soon
    }
    
    return(SR_mat)
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
#' 

average_SR_mat <- function(res_vpa,
                           res_SR_list,
                           range_list=list(1:500,501:1000),
                           SR_mat,
                           seed_number,
                           start_random_rec_year_name,
                           resid_type="lognormal",
                           resample_year_range=0,
                           bias_correction=TRUE){
    
    allyear_name <- dimnames(SR_mat)[[1]]
    start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
    random_rec_year_period <- (start_random_rec_year):length(allyear_name)

    for(i in 1:length(res_SR_list)){
        SR_mat_tmp <- set_SR_mat(res_vpa, res_SR_list[[i]], SR_mat, seed_number+i,
                                 start_random_rec_year_name, resid_type=resid_type,
                                 resample_year_range=resample_year_range,
                                 bias_correction=bias_correction)
        SR_mat[,as.character(range_list[[i]]),] <-
            SR_mat_tmp[,range_list[[i]],]
    }
    
    return(SR_mat)
}


#' @export
print.myarray <- function(x) cat("array :", dim(x),"\n")


#'
#' future_vpaにインプットとして入れる将来予測の空のarrayを生成する関数
#' 
#' @param nsim シミュレーションの繰り返し回数
#' @param nyear 将来予測の実施年数
#' @param future_initial_year_name 将来予測の「初期値となる」年齢別資源尾数を参照する年。この年のFを使って翌年の個体群動態が将来予測で決定される
#' @param start_biopar_year_name 将来の生物パラメータとして生物パラメータを置き換える年の最初の年
#' @param start_random_rec_year_name 将来の加入を再生産関係から予測する最初の年
#' @param waa_year 将来の年齢別体重を過去の平均値とする場合、過去のパラメータを平均する期間
#' @param waa 将来の年齢別体重を直接与える場合
#' 
#' @export

make_future_data <- function(res_vpa,
                          nsim = 1000, # number of simulation
                          nyear = 50, # number of future year
                          future_initial_year_name = 2017,
                          start_F_year_name = 2018,
                          start_biopar_year_name=2018,
                          start_random_rec_year_name = 2018,                          
                          # biopar setting
                          waa_year, waa=NULL,
                          waa_catch_year, waa_catch=NULL,
                          maa_year, maa=NULL,
                          M_year, M=NULL,
                          # faa setting
                          faa_year,
                          currentF=NULL,futureF=NULL,
                          # HCR setting (not work when using TMB)
                          start_ABC_year_name=2019,
                          HCR_beta=1,
                          HCR_Blimit=-1,
                          HCR_Bban=-1,
                          HCR_year_lag=0,
                          # SR setting
                          res_SR=NULL,                       
                          seed_number=1,
                          resid_type="lognormal", # or resample
                          resample_year_range=0, # only when "resample"
                          bias_correction=TRUE,
                          recruit_intercept=0, # number of additional recruitment (immigration or enhancement)
                          # Other
                          Pope=res_vpa$input$Pope
                          ) 
{

    argname <- ls()
    input <- lapply(argname,function(x) eval(parse(text=x)))
    names(input) <- argname
    
    # define age and year
    nage        <- nrow(res_vpa$naa)
    age_name    <- as.numeric(rownames(res_vpa$naa))
    recruit_age <- min(as.numeric(rownames(res_vpa$naa)))
   
    vpa_nyear           <- ncol(res_vpa$naa)
    future_initial_year <- which(colnames(res_vpa$naa)==future_initial_year_name)
    total_nyear         <- future_initial_year + nyear
    allyear_name        <- min(as.numeric(colnames(res_vpa$naa)))+c(0:(total_nyear-1))
    allyear_label       <- c(rep("VPA",future_initial_year),rep("future",nyear))
    print(tibble(allyear_name, allyear_label))

    # define empty array
    waa_mat <- waa_catch_mat <- M_mat <- maa_mat <- naa_mat <- faa_mat <- caa_mat <- 
        array(0, dim=c(nage, total_nyear, nsim),
              dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim))
    class(waa_mat) <- class(M_mat) <- class(maa_mat) <- class(naa_mat) <- class(faa_mat) <- class(caa_mat) <- "myarray"                                                                                  
    SR_mat <- array(0, dim=c(total_nyear, nsim, 9),
                    dimnames=list(year=allyear_name, nsim=1:nsim,
                                  par=c("a","b","rho", #1-3
                                        "SR_type", # 4
                                        "rand_resid", # 5
                                        "deviance", #6
                                        "recruit","ssb",
                                        "intercept")))  #9
    HCR_mat <- array(0, dim=c(total_nyear, nsim, 8),
                    dimnames=list(year=allyear_name, nsim=1:nsim,
                                  par=c("beta","Blimit","Bban","gamma","year_lag", #1-5
                                        "alpha","par2","par3")))  # 6-8
    class(SR_mat)  <- "myarray"
    class(HCR_mat) <- "myarray"

    HCR_mat[,,"Blimit"] <- HCR_mat[,,"Bban"] <- -1
    HCR_mat[,,"beta"] <- HCR_mat[,,"alpha"] <- 1
    
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

    if(is.null(res_vpa$input$dat$waa_catch)){
        waa_catch_mat <- waa_mat
    }
    else{
        waa_catch_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa_catch)
        waa_catch_mat <- make_array(waa_catch_mat, waa_catch, waa_catch_year)
    }    

    # set SR parameter
    SR_mat <- set_SR_mat(res_vpa, res_SR, SR_mat, seed_number,
                         start_random_rec_year_name, resid_type=resid_type,
                         resample_year_range=resample_year_range,
                         bias_correction=bias_correction)

    SR_mat[,,"intercept"] <- recruit_intercept

    # set F & HCR parameter
    start_F_year <- which(allyear_name==start_F_year_name)
    start_ABC_year <- which(allyear_name==start_ABC_year_name)    
    faa_mat <- make_array(faa_mat, currentF, faa_year, start_F_year_name)
    faa_mat <- make_array(faa_mat, futureF,  faa_year, start_ABC_year_name)
    
    HCR_mat[start_ABC_year:total_nyear,,"beta"    ] <- HCR_beta
    HCR_mat[start_ABC_year:total_nyear,,"Blimit"  ] <- HCR_Blimit
    HCR_mat[start_ABC_year:total_nyear,,"Bban"    ] <- HCR_Bban    
    HCR_mat[start_ABC_year:total_nyear,,"year_lag"] <- HCR_year_lag
    
    
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
                     nsim = nsim,
                     nage = nage,
                     recruit_age = recruit_age,
                     HCR_mat = HCR_mat,
                     obj_stat = 0, # 0: mean, 1:geomean
                     objective = 0, # 0: MSY, 1: PGY, 2: percentB0 or Bempirical
                     obj_value = -1
                     )

    return(tibble::lst(data=tmb_data,input=input))
}

#'
#' @export
#' 


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
    
    return(list(wcaa=wcaa, naa=naa, faa=faa, SR_mat=tmb_data$SR_mat))
}

#' 将来予測の実施関数
#'
#' @param tmb_data make_future_dataの返り値
#'
#' @export
#'

future_vpa <- function(tmb_data,
                       optim_method="tmb", # or "R" or "none"
                       multi_init = 1,
                       multi_lower = 0,
                       multi_upper = 10,
                       objective ="MSY", # or PGY, percentB0, Bempirical
                       obj_value = 0,                         
                       obj_stat  ="mean", 
                       compile=FALSE){

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
    
    if(optim_method=="tmb" | optim_method=="both"){
        
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

    #--- R ----
    if(optim_method=="R" | optim_method=="none"){
        R_obj_fun <- function(x, tmb_data, what_return="obj"){
            tmb_data$x <- x
            tmb_data$what_return <- what_return
            obj <- do.call(future_vpa_R, tmb_data)
            return(obj)
        }
        
        if(optim_method=="R"){
            msy_optim <- nlminb(start=0, objective=R_obj_fun, tmb_data=tmb_data,
                                lower=list(x=log(multi_lower)), upper=list(x=log(multi_upper)))
            msy <- exp(-msy_optim$objective)
            multi <- exp(msy_optim$par)
            tmb_data$x <- log(multi)
        }
        else{
            tmb_data$x <- log(multi_init)            
        }
        tmb_data$what_return <- "stat"
        res_future <- do.call(future_vpa_R, tmb_data)
    }

    res_future$input <- input

    return(res_future)

    # 足りないもの
    # 推定結果の簡易的グラフ(MSY_est_plot)
    # waa.fun
    # 生産関係の細かい設定
    # FperSPRの計算結果を出力する
    # yield curveが変
    # do_MSE option
    # 2do_futureで、SR関数を置き換えるオプション
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
                      nsim,
                      nage,
                      recruit_age,
                      obj_stat,
                      objective,
                      obj_value,
                      x,
                      what_return="obj",
                      HCR_mat                      
                      ){

    argname <- ls()
    tmb_data <- lapply(argname,function(x) eval(parse(text=x)))
    names(tmb_data) <- argname    

    F_mat <- N_mat <-  naa_mat
    F_mat[] <- N_mat[] <-  0
    class(F_mat) <- class(N_mat) <- "myarray"
    
    N_mat[,1:future_initial_year,] <- naa_mat[,1:future_initial_year,]
    spawner_mat <- apply(N_mat * waa_mat * maa_mat, c(2,3) , sum)

    F_mat[,1:(start_ABC_year-1),] <- faa_mat[,1:(start_ABC_year-1),]
    F_mat[,start_ABC_year:total_nyear,] <- faa_mat[,start_ABC_year:total_nyear,] * exp(x)

    for(t in future_initial_year:total_nyear){
        spawner_mat[t,] <- colSums(N_mat[,t,] * waa_mat[,t,] * maa_mat[,t,])
        if(t>future_initial_year){
            spawn_t <- t-recruit_age
            N_mat[1,t,] <- purrr::pmap_dbl(tibble(x=SR_mat[t,,"SR_type"],
                             ssb=spawner_mat[spawn_t,],
                             a=SR_mat[t,,"a"],b=SR_mat[t,,"b"]),
                      function(x,ssb,a,b){
                          fun <- list(SRF_HS,SRF_BH,SRF_RI)[[x]];
                          fun(ssb,a,b)
                      })
            N_mat[1,t,] <- N_mat[1,t,]*exp(SR_mat[t,,"deviance"]) + 
                SR_mat[t,,"intercept"];            
        }

        # harvest control rule
        ssb_tmp <- spawner_mat[cbind(t-HCR_mat[t,,"year_lag"],1:nsim)]
        Blimit <- HCR_mat[t,,"Blimit"]
        HCR_mat[t, ,"alpha"] <- HCR_mat[t, ,"beta"]
        tmp <- ssb_tmp < Blimit
        HCR_mat[t, tmp ,"alpha"] <- HCR_mat[t,tmp,"beta"]*(ssb_tmp[tmp]-HCR_mat[t,tmp,"Bban"])/
            (HCR_mat[t,tmp,"Blimit"]-HCR_mat[t,tmp,"Bban"])
        HCR_mat[t, HCR_mat[t,,"alpha"]<0 ,"alpha"] <- 0

        F_mat[,t,] <- F_mat[,t,]*HCR_mat[t,,"alpha"]

        if(t<total_nyear){
            # forward calculation                 
            for(iage in 1:(nage-1)) {
                N_mat[iage+1,t+1,] <- N_mat[iage,t,]*exp(-M_mat[iage,t,]-F_mat[iage,t,])
            }
            N_mat[nage,t+1,] <- N_mat[nage,t+1,] + N_mat[nage,t,]*exp(-M_mat[nage,t,]-F_mat[nage,t,])
        }
    }

    if(Pope==1){
        wcaa_mat <- N_mat*(1-exp(-F_mat))*exp(-M_mat/2) * waa_catch_mat
    }
    else{
        wcaa_mat <- N_mat*(1-exp(-F_mat-M))*F_mat/(F_mat+M) * waa_catch_mat
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
        list(naa=N_mat, wcaa=wcaa_mat, faa=F_mat, SR_mat=tmb_data$SR_mat,
             alpha=HCR_mat[,,"alpha"],multi=exp(x))
    }

}

#'
#' @export
#' 


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
#' 

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
#' 


format_to_old_future <- function(fout){
    fout_old <- fout[c("naa","faa","multi","alpha")]
    fout_old$waa       <- fout$input$tmb_data$waa_mat
    fout_old$waa.catch <- fout$input$tmb_data$waa_catch_mat        
    fout_old$maa       <- fout$input$tmb_data$maa_mat
    fout_old$M         <- fout$input$tmb_data$M_mat        
    fout_old$vssb      <- apply(fout$naa * fout_old$waa * fout_old$maa, c(2,3), sum)
    fout_old$vbiom     <- apply(fout$naa * fout_old$waa, c(2,3),sum)
    fout_old$vwcaa     <- apply(fout$wcaa,c(2,3),sum)
    fout_old$currentF  <- fout$faa[,fout$input$tmb_data$start_ABC_year-1,1]
    fout_old$futureF   <- fout$faa[,fout$input$tmb_data$start_ABC_year,1]    
    fout_old$caa       <- fout$wcaa/fout_old$waa
    fout_old$multi     <- fout$multi
    fout_old$recruit   <- fout$SR_mat[,,"recruit"]
    return(fout_old)
}
