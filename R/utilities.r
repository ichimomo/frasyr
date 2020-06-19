#'
#' @import ggplot2
#' @import magrittr          
#' @import dplyr             
#' @import tidyr             
#' @import tibble 
#' @import readr
#' @import stringr           
#' @importFrom magrittr %>%
#' @importFrom magrittr %T>%
#' @importFrom dplyr filter
#' 
NULL

# 将来予測や管理基準値計算に使われる関数群 ----

# 個体群動態に関する基本的関数 ----

#'
#' 年齢別パラメータを与えて、１年分前進計算する
#'
#' @export
#' @encoding UTF-8

calc.rel.abund <- function(sel,Fr,na,M,waa,waa.catch=NULL,maa,min.age=0,max.age=Inf,Pope=TRUE,ssb.coef=0){
  if(is.null(waa.catch)) waa.catch <- waa
  rel.abund <- rep(NA, na)
  rel.abund[1] <- 1
  for (i in 2:(na-1)) {
    rel.abund[i] <- rel.abund[i-1]*exp(-M[i-1]-sel[i-1]*Fr)
  }
  rel.abund[na] <- rel.abund[na-1]*exp(-M[na-1]-sel[na-1]*Fr)*(1-exp(-((max.age-min.age)-(na-2))*(M[na]+sel[na]*Fr)))/(1-exp(-M[na]-sel[na]*Fr))
  
  if(isTRUE(Pope)){
    ypr1 <- rel.abund*waa.catch[1:na]*(1-exp(-sel[1:na]*Fr))*exp(-M[1:na]/2)
  }
  else{
    # use Baranov catch equation
    ypr1 <- rel.abund*(1-exp(-sel[1:na]*Fr-M[1:na]))*sel[1:na]*Fr/
      (sel[1:na]*Fr+M[1:na])*waa.catch[1:na]
  }
  spr <- rel.abund*waa[1:na]*maa[1:na]*exp(-ssb.coef*(sel[1:na]*Fr+M[1:na])) 
  return(list(rel.abund=rel.abund,ypr=ypr1,spr=spr))
}

#'
#' 年齢別生物パラメータとFと漁獲量を与えると与えた漁獲量と一致するFへの乗数を返す
#'  
#' @export
#' @encoding UTF-8

caa.est.mat <- function(naa,saa,waa,M,catch.obs,Pope,set_max1=TRUE){
  if(set_max1==TRUE) saa <- saa/max(saa)
  tmpfunc <- function(logx,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,out=FALSE,Pope=Pope){
    x <- exp(logx)
    if(isTRUE(Pope)){
      caa <- naa*(1-exp(-saa*x))*exp(-M/2)
    }
    else{
      caa <- naa*(1-exp(-saa*x-M))*saa*x/(saa*x+M)
    }
    wcaa <- caa*waa
    if(out==FALSE){
      return(log((sum(wcaa,na.rm=T)-catch.obs)^2))
    }
    else{
      return(caa)
    }
  }
  tmp <- optimize(tmpfunc,log(c(0.000001,10)),catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=FALSE)#,tol=.Machine$double.eps)
  tmp2 <- tmpfunc(logx=tmp$minimum,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=TRUE)
  return(list(x=exp(tmp$minimum),caa=tmp2))
}

# #　上の関数とどっちが使われているか不明,,,多分使われていないのでコメントアウトする
# caa.est <- function(naa,saa,waa,M,catch.obs,Pope){
#   saa <- saa/max(saa)
#   tmpfunc <- function(x,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,out=FALSE,Pope=Pope){
#     if(isTRUE(Pope)){
#       caa <- naa*(1-exp(-saa*x))*exp(-M/2)
#     }
#     else{
#       caa <- naa*(1-exp(-saa*x-M))*saa*x/(saa*x+M)
#     }
#     wcaa <- caa*waa
#     if(out==FALSE){
#       return((sum(wcaa,na.rm=T)-catch.obs)^2)
#     }
#     else{
#       return(caa)
#     }
#   }
#   tmp <- optimize(tmpfunc,c(0,5),catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=FALSE)
#   tmp2 <- tmpfunc(x=tmp$minimum,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=TRUE)
#   return(list(x=tmp$minimum,caa=tmp2))
# }

#' @export
#' @encoding UTF-8

solv.Feq <- function(cvec,nvec,mvec){
  Fres <- rep(0,length(cvec))
  # cat(nvec," ")
  for(i in 1:length(cvec)){
    F0 <- cvec[i]/nvec[i]
    F1 <- cvec[i]*(F0+mvec[i])/nvec[i]/(1-exp(-F0-mvec[i]))
    if(round(cvec[i],6)<round(nvec[i],6)){
      while(abs(F0-F1)>0.0001 ){
        F0 <- F1
        F1 <- cvec[i]*(F0+mvec[i])/nvec[i]/(1-exp(-F0-mvec[i]))
        if(F0-F1==-Inf) cat("\n",cvec[i]," ",nvec[i]," \n")
      }
      Fres[i] <- F1
    }
    else{
      Fres[i] <- 10
      cat("Warning: catch exceeded tot_num at: ",i," ",
          round(cvec[i],6)," ",round(nvec[i],6),"\n")
    }
  }
  Fres
}

#' VPAの結果に格納されている生物パラメータから世代時間を計算
#'
#' @export
#' @encoding UTF-8

Generation.Time <- function(vpares,
                            maa.year=2014:2015,
                            maa=NULL,
                            M.year=2014:2015,
                            M=NULL,
                            Plus = 19
){
  
  if(is.null(maa)){
    maa <- vpares$input$dat$maa
    maa <- rowMeans(maa[,colnames(maa) %in% maa.year,drop=F],na.rm=T)
    maa <- maa[!is.na(maa)]
  }
  if(is.null(M)){
    M <- vpares$input$dat$M
    M <- rowMeans(M[,colnames(M) %in% M.year,drop=F],na.rm=T)
    M <- M[!is.na(M)]
  }
  
  age <- as.numeric(names(maa))
  maa <- c(maa, rep(1,Plus))
  M <- c(M, rep(M[length(M)],Plus))
  age <- c(age, max(age)+1:Plus)
  A <- length(M)
  L <- c(1,exp(-cumsum(M[-A])))
  G <- sum(age*L*maa)/sum(L*maa)
  
  return(G)
}

### dynamics MSYを計算してみる                                                                                  
dyn.msy <- function(naa.past,naa.init=NULL,fmsy,a,b,resid,resid.year,waa,maa,M,SR=TRUE){
  nyear <- length(resid)
  if(is.null(naa.init)) nage <- nrow(naa.past) else nage <- length(naa.init)
  naa <- matrix(0,nage,nyear)
  ssb <- numeric()
  if(is.null(naa.init)) naa[,1] <- naa.past[,colnames(naa.past)==min(resid.year)]
  else naa[,1] <- naa.init
  colnames(naa) <- resid.year
  if(is.null(naa.init)){
    waa <- waa[,colnames(naa.past)%in%resid.year]
    maa <- maa[,colnames(naa.past)%in%resid.year]
    M <- M[,colnames(naa.past)%in%resid.year]
  }
  for(i in 2:nyear){
    ssb[i-1] <- sum(naa[,i-1]*waa[,i-1]*maa[,i-1],na.rm=T)
    if(SR==TRUE){
      naa[1,i] <- HS(ssb[i-1],a,b)*exp(resid[i])
    }
    else{
      naa[1,i] <- naa.past[1,i]
    }
    for(j in 2:(nage-1)) naa[j,i] <- naa[j-1,i-1] * exp(-fmsy[j-1]-M[j-1,i-1])
    naa[nage,i] <- naa[nage-1,i-1] * exp(-fmsy[j-1]-M[j-1,i-1]) + naa[nage,i-1] * exp(-fmsy[nage]-M[nage,i-1])
  }
  i <- nyear ; ssb[i] <- sum(naa[,i]*waa[,i]*maa[,i])
  list(naa=naa,ssb=ssb)
}



# 再生産関係を仮定しない管理基準値計算 ----

#'
#' 再生産関係を仮定しない管理基準値計算(SPR,YPR,F0.1,Fmax)のための関数
#'
#' @param res VPAの出力結果
#' @param sel 仮定する選択率．NULLの場合，res$Fc.at.ageが使われる
#' @param waa 仮定する年齢別体重。直接の値を入れるか，waa.yearで年を指定するやり方のどちらでも動く。直接指定するほうが優先。
#' @param maa 仮定する年齢別成熟率。直接の値を入れるか，waa.yearで年を指定するやり方のどちらでも動く。直接指定するほうが優先。
#' @param M 仮定する年齢別死亡率。直接の値を入れるか，waa.yearで年を指定するやり方のどちらでも動く。直接指定するほうが優先。
#' @param waa.catch　仮定する年齢別体重（漁獲量計算用）。直接の値を入れるか，waa.yearで年を指定するやり方のどちらでも動く。直接指定するほうが優先。
#' @param M.year 年を指定して生物パラメータを仮定する場合．年の範囲の平均値が用いられる．NULLの場合，VPA最終年の値が使われる
#' @param waa.year 年を指定して生物パラメータを仮定する場合．年の範囲の平均値が用いられる．NULLの場合，VPA最終年の値が使われる
#' @param maa.year 年を指定して生物パラメータを仮定する場合．年の範囲の平均値が用いられる．NULLの場合，VPA最終年の値が使われる
#' @param rps.year Fmedの計算に使うRPSの年の範囲．NULLの場合，全範囲が用いられる
#' @param max.age 加入年齢を０歳としたときに、SPR計算で考慮される最大の年齢（年齢の数ではないことに注意, デフォルトはInf）。加入年齢が１歳以上のときは、SPR計算で考慮したい年齢-加入年齢を入力する、またはmin.ageの引数に加入年齢を設定する。
#' @param min.age  加入年齢が0歳でないときに指定できる(デフォルトは0)
#' @param  pSPR = seq(10,90,by=10), # F\%SPRを計算するときの％SPR
#' @param d 0.001
#' @param  Fem.init 経験的管理基準値(Fmed, Fmean, Fhigh, Flow)の初期値 (default=0.5)
#' @param  Fmax.init Fmaxの初期値 (default=1.5)
#' @param  F0.1.init F0.1の初期値 (default=0.7)
#' @param  iterlim 
#' @param  plot 結果のプロットを表示するかどうか
#' @param  Pope Popeの式を使うか
#' @param  F.range YPR, SPR曲線を書くときのFの範囲（Fの最大値のスケール）、かつ、F\%SPRを計算するときの初期値を決めるために利用される。F\%SPRの推定がうまくいかない場合はこの範囲を調整してください。
#' @encoding UTF-8
#'
#' @note F_SPRのF管理基準値の初期値は　与えられたFのもとでのSPR/目的のSPR　を初期値とするように調整されるので不要。
#'
#' @export
#' @import tibble 
#' @encoding UTF-8
#' 

# ref.F
ref.F <- function(
  res, # VPAの結果のオブジェクト
  #  sel=NULL, # 仮定する選択率．NULLの場合，res$Fc.at.ageが使われる
  Fcurrent=NULL, # Fcurrentの仮定．NULLの場合，res$Fc.at.ageが使われる  
  waa=NULL, # 仮定する生物パラメータ．直接の値を入れるか，年を指定するやり方のどちらでも動く。直接指定するほうが優先。
  maa=NULL,
  M=NULL,
  waa.catch=NULL,
  M.year=NULL, 
  waa.year=NULL, # 年を指定して生物パラメータを仮定する場合．年の範囲の平均値が用いられる．NULLの場合，VPA最終年の値が使われる
  maa.year=NULL,
  rps.year = NULL, # Fmedの計算に使うRPSの年の範囲．NULLの場合，全範囲が用いられる
  max.age = Inf, # 加入年齢を０歳としたときに、SPR計算で考慮される最大の年齢（年齢の数ではないことに注意）。加入年齢が１歳以上のときは、SPR計算で考慮したい年齢-加入年齢を入力する、またはmin.ageの引数に加入年齢を設定する。
  min.age = 0, # 加入年齢が0歳でないときに指定できる
  d = 0.001,
  Fem.init = 0.5, 
  Fmax.init = 1.5, # Fmaxの初期値
  F0.1.init = 0.7, # F0.1の初期値
  pSPR = seq(10,90,by=10), # F%SPRを計算するときの％SPR
  iterlim=1000,
  plot=TRUE,
  Pope=NULL, # 2014.7.4追加
  F.range = seq(from=0,to=2,length=101)  # YPR, SPR曲線を書くときのFの範囲
){
  
  argname <- ls()
  arglist <- lapply(argname,function(x) eval(parse(text=x)))
  names(arglist) <- argname
  
  if(is.null(Pope)) Pope <- res$input$Pope
  
  naa <- res$naa
  ssb <- res$ssb
  ny <- ncol(naa)
  years <- dimnames(naa)[[2]]
  ages <- dimnames(naa)[[1]]
  
  if(is.null(Fcurrent)){
    Fcurrent <- res$Fc.at.age
  }
  sel <- Fcurrent/max(Fcurrent,na.rm=TRUE)
  na <- sum(!is.na(Fcurrent))
  
  if(is.null(waa.year)) waa.year <- rev(years)[1]
  if(is.null(maa.year)) maa.year <- rev(years)[1]
  if(is.null(M.year)) M.year <- rev(years)[1]
  if(is.null(rps.year)) rps.year <- as.numeric(colnames(res$naa))
  
  if(is.null(waa))  waa <- apply_year_colum(res$input$dat$waa,waa.year)
  if(is.null(M))    M   <- apply_year_colum(res$input$dat$M,M.year)
  if(is.null(maa))  maa <- apply_year_colum(res$input$dat$maa,maa.year)
  
  if(is.null(waa.catch)){
    if(is.null(res$input$dat$waa.catch)){
      waa.catch <- waa
    }
    else{
      waa.catch <- apply_year_colum(res$input$dat$waa.catch,waa.year)
    }
  }
  
  ssb.coef <- ifelse(is.null(res$ssb.coef),0,res$ssb.coef)
  
  min.age <- min(as.numeric(rownames(res$naa)))
  if(min.age==0) slide.tmp <- TRUE else slide.tmp <- -1:-min.age
  rps.data <- data.frame(year=as.numeric(names(colSums(ssb,na.rm=T))),
                         ssb=as.numeric(colSums(ssb,na.rm=T)),
                         recruit=as.numeric(c(naa[1,slide.tmp],rep(NA,min.age))))
  if (sum(is.na(rps.data$year))>0) rps.data <- rps.data[-which(is.na(rps.data$year)),]
  rps.data$rps <- rps <- rps.data$recruit/rps.data$ssb
  #  rps <- as.numeric(naa[1,]/colSums(ssb,na.rm=TRUE))
  
  #  if (is.null(rps.year)) rps.year <- years
  
  tmp <- rps.data$year %in% rps.year
  rps.q <- quantile(rps[tmp], na.rm=TRUE, probs=c(0.1,0.5,0.9))
  rps.q <- c(rps.q,mean(rps[tmp], na.rm=TRUE))  
  names(rps.q)[4] <- "mean"
  spr.q <- 1/rps.q
  
  original.spr <- calc.rel.abund(Fcurrent,1,na,M,waa,waa.catch,maa,min.age=min.age,
                                 max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
  original.spr0 <- calc.rel.abund(Fcurrent,0,na,M,waa,waa.catch,maa,min.age=min.age,
                                  max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
  original.perspr <- sum(original.spr$spr,na.rm=T)/sum(original.spr0$spr,na.rm=T)
  
  
  # Fcurrent
  Fcurrent_max_mean <- c(max(Fcurrent,na.rm=T), mean(Fcurrent,na.rm=T))
  
  # grid search
  Fcurrent_max <- Fcurrent_max_mean[1]
  F.range <- sort(c(F.range,  Fcurrent_max))
  spr0 <- sum(calc.rel.abund(Fcurrent,0,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)$spr,na.rm=T)  
  tmp <- lapply(F.range, function(x) calc.rel.abund(sel,x,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef))
  ypr <- sapply(tmp,function(x) sum(x$ypr,na.rm=T))
  pspr <- sapply(tmp,function(x) sum(x$spr,na.rm=T))/spr0*100
  ypr.spr <- data.frame(F.range=F.range,ypr=ypr,pspr=pspr)
  ypr.spr$Frange2Fcurrent  <- ypr.spr$F.range/Fcurrent_max
  
  # F.spr
  
  spr.f.est <- function(log.p, out=FALSE, sub="med", spr0=NULL){
    Fr <- exp(log.p)
    
    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    spr <- sum(tmp$spr,na.rm=T)
    if (isTRUE(out)) obj <- spr
    else{
      if(sub=="mean") obj <- (spr-spr.q[4])^2       
      if(sub=="low") obj <- (spr-spr.q[3])^2 
      if(sub=="med") obj <- (spr-spr.q[2])^2
      if(sub=="high") obj <- (spr-spr.q[1])^2
      if(is.numeric(sub)) obj <- (spr/spr0-sub/100)^2
      
    }
    
    return(obj)
  }
  
  spr0 <- spr.f.est(-Inf, out=TRUE)
  
  Fmed.res <- nlm(spr.f.est, Fem.init, out=FALSE, sub="med", iterlim = iterlim)
  Fmean.res <- nlm(spr.f.est, Fem.init, out=FALSE, sub="mean", iterlim = iterlim)
  Flow.res <- nlm(spr.f.est, Fem.init, out=FALSE, sub="low", iterlim = iterlim)
  Fhigh.res <- nlm(spr.f.est, Fem.init, out=FALSE, sub="high", iterlim = iterlim)
  
  Fmean <- exp(Fmean.res$estimate)  
  Fmed <- exp(Fmed.res$estimate)
  Flow <- exp(Flow.res$estimate)
  Fhigh <- exp(Fhigh.res$estimate)
  
  if (!is.null(pSPR)){
    FpSPR <- NULL
    
    for (i in pSPR){
      tmp <- which.min(abs(ypr.spr$pspr-i))+c(-1,1)
      tmp <- ifelse(tmp<=0,1,tmp)
      tmp <- ifelse(tmp>length(ypr.spr$F.range),length(ypr.spr$F.range),tmp)        
      Fspr.init <- log(ypr.spr$F.range[tmp]) #original.perspr/i*100
      FpSPR.res <- optimize(spr.f.est, Fspr.init,out=FALSE,sub=i,spr0=spr0)
      #nlm(spr.f.est, Fspr.init, out=FALSE, sub=i, spr0=spr0, iterlim=iterlim)
      #      cat("Estimate F%spr: initial value=", Fspr.init," : estimated value=",exp(FpSPR.res$estimate),"\n")
      FpSPR <- c(FpSPR, exp(FpSPR.res$minimum))
    }
    names(FpSPR) <- paste(pSPR,"%SPR",sep="")
  }
  
  # Fmax
  
  ypr.f.est <- function(log.p, out=FALSE){
    Fr <- exp(log.p)
    
    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    ypr <- sum(tmp$ypr,na.rm=T)    
    
    if (isTRUE(out)) obj <- ypr else obj <- -ypr
    
    return(obj)
  }
  
  Fmax.res <- nlm(ypr.f.est, log(Fmax.init), out=FALSE)
  
  Fmax <- exp(Fmax.res$estimate)
  
  # F0.1
  
  Fp <- function(log.p, out=FALSE){
    Fr <- exp(log.p)
    
    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    ypr <- sum(tmp$ypr,na.rm=T)
    if (isTRUE(out)) obj <- ypr else obj <- -ypr
    
    return(obj)
  }
  
  F0.1.est <- function(log.p){
    p <- exp(log.p)
    ref.trend <- (Fp(log(d))-Fp(log(0)))/d
    trend <- (Fp(log(p+d)) - Fp(log(p)))/d
    
    obj <- (ref.trend/10 - trend)^2
    
    obj
  }
  
  F0.1.res <- nlm(F0.1.est,log(F0.1.init))
  
  F0.1 <- exp(F0.1.res$estimate)
  
  
  # output
  f.mean <- function(x) mean(x*sel, na.rm=T)
  
  Fmean <- c(Fmean, f.mean(Fmean))  
  Fmed <- c(Fmed, f.mean(Fmed))
  Flow <- c(Flow, f.mean(Flow))
  Fhigh <- c(Fhigh, f.mean(Fhigh))
  Fmax <- c(Fmax, f.mean(Fmax))
  F0.1 <- c(F0.1, f.mean(F0.1))
  
  names(Fcurrent_max_mean) <- names(Fmed) <- names(Fmean) <- names(Flow) <- names(Fhigh) <- names(Fmax) <- names(F0.1) <- c("max","mean")
  Res <- list(min.age=min.age, max.age=max.age, rps.q=rps.q, spr.q=spr.q, Fcurrent=Fcurrent_max_mean,
              sel=sel, Fmed=Fmed, Flow=Flow, Fhigh=Fhigh, Fmax=Fmax, F0.1=F0.1, Fmean=Fmean,rps.data=rps.data)
  
  if (!is.null(pSPR)){
    FpSPR <- rbind(FpSPR, sapply(FpSPR, f.mean))
    rownames(FpSPR) <- c("max","mean")
    Res$FpSPR <- FpSPR
  }
  
  #---- make summary
  Res$summary <- as.data.frame(Res[substr(names(Res),1,1)=="F"])
  Res$summary <- rbind(Res$summary,Res$summary[1,]/Res$summary[1,1])
  dimnames(Res$summary)[[1]][3] <- "Fref/Fcur"
  
  Res$currentSPR <- list(SPR=sum(original.spr$spr,na.rm=T),
                         perSPR=original.perspr,
                         YPR=sum(original.spr$spr,na.rm=T),
                         Fcurrent=Fcurrent)
  
  Res$ypr.spr  <- ypr.spr #data.frame(F.range=F.range,ypr=ypr,spr=spr)
  Res$waa <- waa
  Res$waa.catch <- waa.catch  
  Res$maa <- maa
  Res$arglist <- arglist
  Res$spr0 <- spr0
  
  class(Res) <- "ref"
  
  if(isTRUE(plot)){
    plot_Fref(Res)
  }    
  return(Res)
}

#' 毎年のFの\%SPRやターゲットした\%SPRに相当するFの大きさを計算する
#' 
#' VPA計算結果を使って毎年のF at ageがどのくらいのSPR, YPRに相当するかを計算する。また、各年のFが、目標としたSPR（target.SPR）を達成するためのF(Ftarget)の何倍(F/Ftarget)に相当するかも計算する。F/Ftargetは数値的に探索するが、そのときのF/Ftargetの上限をFmaxにて指定する。十分大きい値（デフォルトは１０）を与えておけば大丈夫だが、Ftargetが非常に小さい数字になりそうな場合にはこの値をもっと大きくとるなどする。また、SPRの計算は、デフォルトでは等比級数の和の公式を使って、無限大の年齢までSPRを足しているが、max.ageを指定することで、有限の年齢までの和も計算できる。
#'
#' @param dres vpa関数の返り値
#' @param target.SPR 目標とするSPR。この値を入れると、結果の$ysdata$"F/Ftarget"で、その年のFが目標としたSPR(％)を達成するためのF（Ftarget）の何倍になっているかを返す。デフォルトは30が入っている。このとき、SPRを計算するための生物パラメータ（年齢別体重・成熟率・死亡率）はそれぞれの年で仮定されているものを用いる。
#' @param Fmax F/Ftargetを推定するときに探索するFの乗数の最大値
#' @param max.age SPRやYPRの計算をするときに最大何歳まで考慮するか（デフォルトは無限大)。値の指定の仕方はhelp(ref.F)を参照のこと
#' @encoding UTF-8
#'
#' @examples
#' data(res_vpa)
#' Fratio <- get.SPR(res_vpa,target.SPR=12)$ysdata$"F/Ftarget"
#' plot(colnames(res_vpa$naa),Fratio,type="b")
#' 
#'
#' @export
#' @encoding UTF-8

get.SPR <- function(dres,target.SPR=30,Fmax=10,max.age=Inf){
  dres$ysdata <- matrix(0,ncol(dres$faa),5)
  dimnames(dres$ysdata) <- list(colnames(dres$faa),c("perSPR","YPR","SPR","SPR0","F/Ftarget"))
  for(i in 1:ncol(dres$faa)){
    dres$Fc.at.age <- dres$faa[,i] # Fc.at.ageに対象年のFAAを入れる
    if(all(dres$Fc.at.age>0, na.rm=T)){
      byear <- colnames(dres$faa)[i] # 何年の生物パラメータを使うか
      
      a <- ref.F(dres,waa.year=byear,maa.year=byear,M.year=byear,rps.year=2000:2011,
                 pSPR=round(target.SPR),
                 F.range=c(seq(from=0,to=ceiling(max(dres$Fc.at.age,na.rm=T)*Fmax),
                               length=301),max(dres$Fc.at.age,na.rm=T)),plot=FALSE,max.age=max.age)
      # YPRと%SPR
      dres$ysdata[i,1:2] <- (as.numeric(rev(a$ypr.spr[which(a$ypr.spr$Frange2Fcurrent==1)[1],2:3])))
      # SPR                                                                                               
      dres$ysdata[i,3] <- a$spr0*dres$ysdata[i,1]/100
      # SPR0                                                                                              
      dres$ysdata[i,4] <- a$spr0
      # relative F
      dres$ysdata[i,5] <- 1/a$summary[3,grep("SPR",colnames(a$summary))][1]
    }
    else{
      break;
    }
  }
  dres$ysdata <- as.data.frame(dres$ysdata)
  dres$target.SPR <- target.SPR
  return(dres)
}

#' vpaデータ, 選択率の参照年1, 漁獲圧の参照年2を与えて、参照年2の漁獲圧のもとで参照年1の選択率となるF at ageを作成する関数。漁獲圧の変換は％SPR換算
#'
#' sel_year 選択率の参照年
#' faa_year 漁獲圧の参照年
#'
#' @encoding UTF-8
#' @export

convert_faa_perSPR <- function(res_vpa, sel_year, faa_year, Fcurrent_MSY=NULL, Fsel_MSY=NULL){
  if(is.null(Fcurrent_MSY)){
    Fcurrent_MSY <- apply_year_colum(res_vpa$faa,target_year=faa_year)
  }
  Fcurrent.per.spr <- ref.F(res_vpa,Fcurrent=Fcurrent_MSY,waa.year=faa_year,M.year=faa_year,
                            maa.year=faa_year,plot=FALSE,pSPR=NULL)$currentSPR$perSPR
  cat("%SPR in Fcurrent ",Fcurrent.per.spr,"\n")
  if(is.null(Fsel_MSY)){
    Fsel_MSY <- apply_year_colum(res_vpa$faa,target_year=sel_year)
  }
  
  Fmultiplier <- ref.F(res_vpa,Fcurrent=Fsel_MSY,waa.year=faa_year,M.year=faa_year,
                       maa.year=faa_year,pSPR=Fcurrent.per.spr*100,plot=FALSE)
  cat("----------------------------- \n ")                
  cat("%SPR in Fsel_MSY",Fmultiplier$currentSPR$perSPR," = ")        
  Fmultiplier <- rev(Fmultiplier$summary)[3,1] %>% unlist()
  cat("(F multiplier=",Fmultiplier,")\n")
  Fcurrent_MSY_new <- Fsel_MSY * Fmultiplier
  Fref_tmp <- ref.F(res_vpa,Fcurrent=Fcurrent_MSY_new,waa.year=faa_year,
                    M.year=faa_year,maa.year=faa_year,
                    pSPR=Fcurrent.per.spr*100,plot=FALSE)
  cat("%SPR in new-Fcurrent",Fref_tmp$currentSPR$perSPR,"\n")
  cat("----------------------------- \n ")                        
  matplot(cbind(Fcurrent_MSY_new,Fcurrent_MSY,Fsel_MSY),type="b",pch=1:3,ylab="F",xlab="Age (recruit age=1)")
  legend("bottomright",pch=1:3,col=1:3,c("Fcurrent","Fsel_MSY","new-Fcurrent"))        
  return(Fcurrent_MSY_new)
}




#' @export
#' @encoding UTF-8
make_summary_table <- function(mat_data,side=1,probs=c(0.1,0.5,0.8)){
  res_mat <- cbind(apply(mat_data,side,mean),
                   t(apply(mat_data,side,quantile,probs=probs)))
  colnames(res_mat)[1] <- "average"
  res_mat %>% as_tibble()
}


# 結果の入出力&サマリーの作成 ----

#'
#' VPA結果をcsvファイルに出力する
#' 
#' @param res  VPAの結果
#' @param srres fit.SRの結果
#' @param msyres est.MSYの結果
#' @param fres_current future.vpaの結果(Fcurrent)
#' @param fres_HCR future.vpaの結果(F with HCR)
#' @param kobeII kobeII.matrixの結果
#' @param filename csvファイルとpdfファイルの両方のファイル名を指定する場合（拡張子なしで指定）
#' @param csvname csvファイルのファイル名
#' @param pdfname pdfファイルのファイル名
#' @encoding UTF-8
#' @export

out.vpa <- function(res=NULL,    # VPA result
                    srres=NULL,  # fit.SR result
                    msyres=NULL, # est.MSY result
                    fres_current=NULL,   # future projection result
                    fres_HCR=NULL,   # future projection result                    
                    kobeII=NULL, # kobeII result
                    filename="vpa", # filename without extension
                    csvname=NULL,
                    pdfname=NULL,
                    ci.future=c(0.1,0.5,0.9)
){
  old.par <- par()  
  exit.func <- function(){
    dev.off()
    options(warn=0)      
  }
  on.exit(exit.func())
  
  if(!is.null(filename)){
    csvname <- paste(filename,".csv",sep="")
    pdfname <- paste(filename,".pdf",sep="")
  }
  
  pdf(pdfname)
  par(mfrow=c(3,2),mar=c(3,3,2,1))  
  options(warn=-1)
  
  write.table2 <- function(x,title.tmp="",is.plot=TRUE,...){
    if(is.plot){
      if(!is.null(dim(x))){
        matplot(colnames(x),t(x),type="b",ylim=c(0,max(x,na.rm=T)),pch=substr(rownames(x),1,1))
      }
      else{
        barplot(x)
      }
      title(title.tmp)
    }
    if(!is.null(dim(x))){
      tmp <- matrix("",nrow(x)+1,ncol(x)+1)
      tmp[-1,-1] <- as.character(unlist(x))
      tmp[-1,1] <- rownames(x)
      tmp[1,-1] <- colnames(x)
    }
    else{
      tmp <- x
    }
    write.table(tmp,append=T,sep=",",quote=FALSE,file=csvname,col.names=F,row.names=F,...)
  }
  
  write(paste("# frasyr outputs at ",date()," & ",getwd()),file=csvname)  
  
  if(!is.null(res)){
    write("# VPA results",file=csvname, append=T)
    
    write("\n# catch at age",file=csvname,append=T)    
    write.table2(res$input$dat$caa,title.tmp="Catch at age")
    
    write("\n# maturity at age",file=csvname,append=T)    
    write.table2(res$input$dat$maa,title.tmp="Maturity at age")
    
    write("\n# weight at age for biomass calculation",file=csvname,append=T)    
    write.table2(res$input$dat$waa,title.tmp="Weight at age (for biomass)")
    
    if(!is.null(res$input$dat$waa.catch)){
      write("\n# weight at age for catch calculation",file=csvname,append=T)    
      write.table2(res$input$dat$waa.catch,title.tmp="Weight at age (for catch)")    
    }
    
    write("\n# M at age",file=csvname,append=T)    
    write.table2(res$input$dat$M,title.tmp="M at age")          
    
    write("\n# fishing mortality at age",file=csvname,append=T)    
    write.table2(res$faa,title.tmp="F at age")
    
    write("\n# Current F",file=csvname,append=T)    
    write.table2(res$Fc.at.age,title.tmp="Current F")
    
    write("\n# numbers at age",file=csvname,append=T)    
    write.table2(res$naa,title.tmp="Numbers at age")
    
    write("\n# total and spawning biomass ",file=csvname,append=T)
    x <- rbind(colSums(res$ssb, na.rm=T),colSums(res$baa, na.rm=T),colSums(res$wcaa, na.rm=T))
    rownames(x) <- c("Spawning biomass","Total biomass","Catch biomass")
    write.table2(x,title.tmp="Total and spawning biomass")
    
    write("\n# YPR & SPR history ",file=csvname,append=T)    
    get.SPR(res)$ysdata %>% as_tibble() %>% select(-"F/Ftarget") %>%
      write_csv(path=csvname,append=T, col_names=TRUE)                   
  }
  
  if(!is.null(srres)){
    sr_summary <- 
      as_tibble(srres$pars) %>% mutate(AICc   =srres$AICc,
                                       AIC    =srres$AIC,
                                       method=srres$input$method,
                                       type  =srres$input$SR)      
    write("\n# SR fit resutls",file=csvname,append=T)
    write_csv(sr_summary,path=csvname,append=T,
              col_names=TRUE)
  }  
  
  if(!is.null(msyres)){
    write("\n# MSY Reference points",file=csvname,append=T)
    write_csv(msyres$summary,path=csvname,append=T,
              col_names=TRUE)
  }
  
  tmpfunc <- function(fres){
    if(class(fres)=="future_new"){
      fres <- format_to_old_future(fres)
    }      
    
    write("\n# future F at age",file=csvname,append=T)
    write.table2(apply(fres$faa,c(1,2),mean),title.tmp="Average future F at age")
    
    write("\n# future numbers at age",file=csvname,append=T)
    write.table2(apply(fres$naa,c(1,2),mean),title.tmp="Average future numbers at age")
    
    write("\n# future maturity at age",file=csvname,append=T)
    write.table2(apply(fres$maa,c(1,2),mean),title.tmp="Average future numbers at age")
    
    write("\n# future weight at age",file=csvname,append=T)
    write.table2(apply(fres$waa,c(1,2),mean),title.tmp="Average future numbers at age")
    
    write("\n# future total spawning biomass",file=csvname,append=T)    
    make_summary_table(fres$vssb,1,probs=ci.future) %>%
      write_csv(path=csvname,append=TRUE, col_names = TRUE)
    
    write("\n# future total biomass",file=csvname,append=T)    
    make_summary_table(fres$vbiom,1,probs=ci.future) %>%
      write_csv(path=csvname,append=TRUE, col_names = TRUE)
    
    write("\n# future total catch",file=csvname,append=T)    
    make_summary_table(fres$vwcaa,1,probs=ci.future) %>%
      write_csv(path=csvname,append=TRUE, col_names = TRUE)
  }  
  
  if(!is.null(fres_current)){
    write("\n# future projection under F current",file=csvname,append=T)  
    tmpfunc(fres_current)
  }
  
  if(!is.null(fres_HCR)){
    write("\n# future projection under HCR",file=csvname,append=T)  
    tmpfunc(fres_HCR)
  }
  
  if(!is.null(kobeII)){
    write("\n# Kobe II table",file=csvname,append=T)  
    kobeII.table_name <- names(kobeII)
    for(i in 1:length(kobeII.table_name)){
      write(str_c("\n# ",kobeII.table_name[i]),file=csvname,append=T)        
      write_csv(kobeII[kobeII.table_name[i]][[1]],path=csvname,append=TRUE,
                col_names = TRUE)
    }
  }
  
  
}

#' csvファイルとしてまとめられた資源計算結果を読み込んでRのオブジェクトにする
#' @param tfile 資源計算結果がまとめられたcsvファイルの名前
#' @encoding UTF-8
#'
#' @export

read.vpa <- function(tfile,
                     caa.label="catch at age",
                     maa.label="maturity at age",
                     waa.label="weight at age",
                     waa.biomass.label="weight at age for biomass calculation",
                     waa.catch.label="weight at age for catch calculation",                     
                     M.label="M at age",
                     faa.label="fishing mortality at age",
                     Fc.label="Current F",
                     naa.label="numbers at age",
                     Blimit=NULL,
                     Pope=NULL, # VPA計算時にどっちを使っているか入れる（TRUE or FALSE）。デフォルトはNULLでcaa,faa,naaの関係から自動判別するが、自動判別の結果はcatで出力されるので、それをみて正しく判断されているか確認してください。
                     fc.year=NULL){
  
  tmpdata <- read.csv(tfile,header=F,as.is=F,colClasses="character")
  
  tmpfunc <- function(tmpdata,label,type=NULL){
    flags <- which(substr(tmpdata[,1],1,1)=="#")
    flag.name <- tmpdata[flags,1]
    flag.name <- gsub("#","",flag.name)
    flag.name <- gsub(" ","",flag.name)
    get.flag <- which(flag.name==gsub(" ","",label))
    {if(length(get.flag)>0){
      tdat <- tmpdata[(flags[get.flag]+1):(flags[min(which(flags[get.flag]<flags))]-1),]
      if(!is.null(type)){
        tdat <- tdat[,!apply(tdat=="",2,all)]
        tdat <- as.numeric(tdat)
      }
      else{
        tmp.names <- list(tdat[-1,1],tdat[1,-1])
        tdat <- tdat[,!apply(tdat=="",2,all)]
        tdat <- tdat[!apply(tdat=="",1,all),]
        tdat <- tdat[,!apply(is.na(tdat),2,all)]
        tdat <- tdat[!apply(is.na(tdat),1,all),]
        tdat <- sapply((tdat[-1,-1]),as.numeric)
        tmp.names <- lapply(tmp.names,function(x) x[x!=""])
        tmp.names <- lapply(tmp.names,function(x) x[!is.na(x)])                
        dimnames(tdat) <- tmp.names                
        tdat <- as.data.frame(tdat)
      }
    }
      else{
        tdat <- NULL
      }}
    tdat
  }
  
  dres <- list()
  dres$naa <- tmpfunc(tmpdata,naa.label)
  dres$faa <- tmpfunc(tmpdata,faa.label)
  dres$Fc.at.age <- tmpfunc(tmpdata,Fc.label,type="Fc")
  #  dres$Fc.at.age <- dres$Fc.at.age[!is.na(dres$Fc.at.age)]
  dres$Fc.at.age <- dres$Fc.at.age[1:nrow(dres$naa)]    
  #  if(length(dres$Fc.at.age)!=nrow(dres$naa)) stop("Dimension of Fc.at.age and numbers at age is differerent.")
  
  dres$input <- list()
  dres$input$dat <- list()
  dres$input$dat$maa <- tmpfunc(tmpdata,maa.label)
  dres$input$dat$caa <- tmpfunc(tmpdata,caa.label)    
  dres$input$dat$M <- tmpfunc(tmpdata,M.label)
  dres$input$dat$waa <- tmpfunc(tmpdata,waa.label)
  if(is.null(dres$input$dat$waa)) dres$input$dat$waa <- tmpfunc(tmpdata,waa.biomass.label)        
  dres$input$dat$waa.catch <- tmpfunc(tmpdata,waa.catch.label)      
  if(is.null(dres$input$dat$waa.catch)) dres$input$dat$waa.catch <- dres$input$dat$waa
  
  dres$ssb <- dres$input$dat$waa * dres$input$dat$maa * dres$naa
  dres$ssb <- as.data.frame(dres$ssb)
  
  dres$baa <- dres$input$dat$waa * dres$naa
  dres$baa <- as.data.frame(dres$baa)
  
  # setting total catch
  dres$wcaa <- dres$input$dat$waa.catch * dres$input$dat$caa
  dres$wcaa <- as.data.frame(dres$wcaa)
  
  dres$Blimit <- Blimit
  
  ## catch at ageの計算時にpopeの近似式を使っているかどうか、通常は外から情報として与えてほしいところだが、与えられない場合、入力されたcaa,faa,naaの関係を見て、Popeで計算されているのかそうでないのかを判断してdres$input$Popeに入れる
  if(is.null(Pope)){
    caa.pope  <- dres$naa*(1-exp(-dres$faa))*exp(-dres$input$dat$M/2)
    diff.pope <- mean(unlist(dres$input$dat$caa/caa.pope),na.rm=T)
    
    faa <- dres$faa
    M <- dres$input$dat$M
    caa.bara <- dres$naa*faa/(faa+M)*(1-exp(-faa-M))
    diff.bara <- mean(unlist(dres$input$dat$caa/caa.bara),na.rm=T)
    
    if(abs(1-mean(diff.bara))>abs(1-mean(diff.pope))){
      dres$input$Pope <- TRUE
      
      cat("Pope is TRUE... OK? (mean difference=", 1-mean(diff.pope),")\n")
    }
    else{
      dres$input$Pope <- FALSE
      cat("Pope is FALSE... OK? (mean difference=", 1-mean(diff.bara),")\n")            
    }
  }
  else{
    dres$input$Pope <- Pope
  }    
  if(is.null(dres$Fc.at.age) && !is.null(fc.year)) dres$Fc.at.age <- apply(dres$faa[,colnames(dres$faa)%in%fc.year],1,mean)
  
  # その他、他関数で必要になるVPAへのインプット
  dres$input$last.catch.zero <- FALSE
  
  return(dres)
}

#' tidy形式のVPAへのインプットデータをdata_dandlerへ渡す形式に変換する（暫定版）
#'
#' @export
#' @encoding UTF-8
#' 

to_vpa_data <- function(x, label_name){
  vpa_data <- x %>% dplyr::filter(label==label_name) %>% as.data.frame()
  if(label_name!="abund"){  
    rownames(vpa_data) <- str_c("X",vpa_data$year)
    vpa_data <- vpa_data %>% select(-year, -label) 
    vpa_data$value <- vpa_data$fishery <- NULL    
    vpa_data <- as.data.frame(t(vpa_data))
    rownames(vpa_data) <- str_sub(rownames(vpa_data),5,5)
  }
  else{
    vpa_data <- vpa_data %>% 
      select(-starts_with("age_")) %>%
      pivot_wider(names_from=year, values_from=value) 
    abund_name <- vpa_data$fishery
    vpa_data <- as.data.frame(vpa_data)
    vpa_data$label <- vpa_data$fishery <- NULL
    rownames(vpa_data) <- abund_name
    colnames(vpa_data) <- str_c("X",colnames(vpa_data))
  }
  vpa_data  
}


#' @export
#' @encoding UTF-8
#' 

get.stat <- function(fout,eyear=0,tmp.year=NULL, use_new_output=FALSE){
  
  if(isTRUE(use_new_output)){
    fout <- format_to_old_future(fout)
    col.target <- TRUE
  }
  else{
    col.target <- ifelse(fout$input$N==0,1,-1)
  }
  tc <- fout$caa * fout$waa.catch     
  tmp <- as.numeric(fout$vssb[(nrow(fout$vssb)-eyear):nrow(fout$vssb),col.target])
  if(is.null(tmp.year)) tmp.year <- (nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa)
  a <- data.frame("catch.mean"=mean(fout$vwcaa[tmp.year,col.target]),
                  "catch.sd"=sd(fout$vwcaa[tmp.year,col.target]),
                  "catch.geomean"=geomean(fout$vwcaa[tmp.year,col.target]),
                  "catch.median"=median(fout$vwcaa[tmp.year,col.target],na.rm=T),
                  "catch.L10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.1),
                  "catch.H10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.9),
                  "ssb.mean"=mean(fout$vssb[tmp.year,col.target]),
                  "ssb.sd"=sd(fout$vssb[tmp.year,col.target]),                        
                  "ssb.geomean"=geomean(fout$vssb[tmp.year,col.target]),
                  "ssb.median"=median(fout$vssb[tmp.year,col.target],na.rm=T),
                  "ssb.L10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.1),
                  "ssb.H10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.9),
                  "biom.mean"=mean(fout$vbiom[tmp.year,col.target]),
                  "biom.sd"=sd(fout$vbiom[tmp.year,col.target]),                        
                  "biom.geomean"=geomean(fout$vbiom[tmp.year,col.target]),
                  "biom.median"=median(fout$vbiom[tmp.year,col.target],na.rm=T),
                  "biom.L10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.1),
                  "biom.H10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.9),
                  "rec.mean"=mean(unlist(fout$naa[1,,])[tmp.year,col.target]),
                  "rec.sd"=sd(unlist(fout$naa[1,,])[tmp.year,col.target]),
                  "rec.geomean"=geomean(unlist(fout$naa[1,,])[tmp.year,col.target]),
                  "rec.median"=median(unlist(fout$naa[1,,])[tmp.year,col.target],na.rm=T),
                  "rec.L10"=quantile(unlist(fout$naa[1,,])[tmp.year,col.target],na.rm=T,probs=0.1),
                  "rec.H10"=quantile(unlist(fout$naa[1,,])[tmp.year,col.target],na.rm=T,probs=0.9),
                  #                    "lower.HSpoint"=lhs,
                  "Fref2Fcurrent"=fout$multi,
                  fmulti=fout$multi
  )
  a$U.mean <- a$catch.mean/a$biom.mean
  a$U.median <- a$catch.median/a$biom.median
  a$U.geomean <- a$catch.geomean/a$biom.geomean
  
  a$catch.CV <- a$catch.sd/a$catch.mean
  a$ssb.CV <- a$ssb.sd/a$ssb.mean
  a$biom.CV <- a$biom.sd/a$biom.mean
  a$rec.CV <- a$rec.sd/a$rec.mean
  
  #    Faa <- as.data.frame(t(fout$multi * fout$input$res0$Fc.at.age))
  # Faa <- as.data.frame(t(fout$multi * fout$currentF))
  if("finalmeanF" %in% names(fout)) Faa <- as.data.frame(t(fout$finalmeanF)) else Faa <- as.data.frame(t(fout$multi * fout$futureF))
  colnames(Faa) <- paste("F",dimnames(fout$naa)[[1]],sep="")
  res.stat1 <- cbind(a,Faa) # ここまで、get.stat
  
  agename <- dimnames(fout$naa)[[1]]
  nage <- dim(fout$naa)[[1]]    
  tb <- fout$naa * fout$waa 
  if(is.null(fout$waa.catch)) fout$waa.catch <- fout$waa
  ssb <- fout$naa * fout$waa *fout$maa  
  tb.mat <- tc.mat <- ssb.mat <- matrix(0,nage,6)
  for(i in 1:nage){
    tb.mat[i,1] <- mean(tb[i,tmp.year,col.target])
    tb.mat[i,2] <- median(tb[i,tmp.year,col.target])
    tb.mat[i,3] <- geomean(tb[i,tmp.year,col.target])
    tb.mat[i,4] <- mean(tb[i,tmp.year,1])
    tb.mat[i,5:6] <- quantile(tb[i,tmp.year,col.target],probs=c(0.1,0.9),na.rm=T)
    
    tc.mat[i,1] <- mean(tc[i,tmp.year,col.target])
    tc.mat[i,2] <- median(tc[i,tmp.year,col.target])
    tc.mat[i,3] <- geomean(tc[i,tmp.year,col.target])
    tc.mat[i,4] <- mean(tc[i,tmp.year,1])
    tc.mat[i,5:6] <- quantile(tc[i,tmp.year,col.target],probs=c(0.1,0.9),na.rm=T)            
    
    ssb.mat[i,1] <- mean(ssb[i,tmp.year,col.target])
    ssb.mat[i,2] <- median(ssb[i,tmp.year,col.target])
    ssb.mat[i,3] <- geomean(ssb[i,tmp.year,col.target])
    ssb.mat[i,4] <- mean(ssb[i,tmp.year,1])
    ssb.mat[i,5:6] <- quantile(ssb[i,tmp.year,col.target],probs=c(0.1,0.9),na.rm=T)                        
  }
  tc.mat <- as.numeric(tc.mat)
  tb.mat <- as.numeric(tb.mat)
  ssb.mat <- as.numeric(ssb.mat)        
  
  # mean, ME; median, geomean; geometric mean
  names(tc.mat) <- c(paste("TC-mean-A",agename,sep=""),paste("TC-median-A",agename,sep=""),
                     paste("TC-geomean-A",agename,sep=""),paste("TC-det-A",agename,sep=""),
                     paste("TC-L10-A",agename,sep=""),paste("TC-H10-A",agename,sep=""))
  names(tb.mat) <- c(paste("TB-mean-A",agename,sep=""),paste("TB-median-A",agename,sep=""),
                     paste("TB-geomean-A",agename,sep=""),paste("TB-det-A",agename,sep=""),
                     paste("TB-L10-A",agename,sep=""),paste("TB-H10-A",agename,sep=""))
  names(ssb.mat) <- c(paste("SSB-GA-A",agename,sep=""),paste("SSB-median-A",agename,sep=""),
                      paste("SSB-geomean-A",agename,sep=""),paste("SSB-det-A",agename,sep=""),
                      paste("SSB-L10-A",agename,sep=""),paste("SSB-H10-A",agename,sep=""))
  res.stat2 <- as.data.frame(t(c(tb.mat,tc.mat,ssb.mat)))
  res.stat <- cbind(res.stat1,res.stat2) %>% as_tibble()
  return(res.stat)    
}

get.stat3 <- get.stat


get.stat4 <- function(fout,Brefs,
                      refyear=c(2019:2023,2028,2038)){
  col.target <- ifelse(fout$input$N==0,1,-1)
  years <- as.numeric(rownames(fout$vwcaa))
  
  if(is.null(refyear)){
    refyear <- c(seq(from=min(years),to=min(years)+5),
                 c(min(years)+seq(from=10,to=20,by=5)))
  }
  
  catch.mean <- rowMeans(fout$vwcaa[years%in%refyear,col.target])
  names(catch.mean) <- str_c("Catch",names(catch.mean))
  catch.mean <- as_tibble(t(catch.mean))
  
  Btarget.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]>Brefs$Btarget) %>%
    t() %>% as_tibble()
  names(Btarget.prob) <- str_c("Btarget_prob",names(Btarget.prob))
  
  #    Blow.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]>Brefs$Blow) %>%
  #        t() %>% as_tibble()
  #    names(Blow.prob) <- str_c("Blow_prob",names(Blow.prob))
  
  Blimit.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]<Brefs$Blimit) %>%
    t() %>% as_tibble()
  names(Blimit.prob) <- str_c("Blimit_prob",names(Blimit.prob))
  
  Bban.prob <- rowMeans(fout$vssb[years%in%refyear,col.target]<Brefs$Bban) %>%
    t() %>% as_tibble()
  names(Bban.prob) <- str_c("Bban_prob",names(Bban.prob))
  
  return(bind_cols(catch.mean,Btarget.prob,Blimit.prob,Bban.prob))
}

geomean <- function(x)
{
  ifelse(all(x > 0), exp(mean(log(x))), NA)
}

get.trace <- function(trace){
  trace <- trace  %>% as_tibble() %>%
    select(starts_with("TC-mean"),ssb.mean,fmulti,catch.CV) %>%
    mutate(label=as.character(1:nrow(.)))
  
  trace <- trace %>% gather(value=value,key=age,-label,-fmulti,-ssb.mean,-catch.CV) %>%
    mutate(age=str_extract(age, "(\\d)+")) %>%
    mutate(age=factor(age)) %>%
    mutate(age=fct_reorder(age,length(age):1))
  return(trace)
}


#' 列が年である行列に対して、年を指定するとその年のあいだの平均値（等）を返す関数
#'
#' @encoding UTF-8
#' @export
#'
#' 

apply_year_colum <- function(mat,target_year,stat="mean"){
  if(target_year[1]<0){
    target_year <- rev(colnames(mat))[-target_year] %>%
      as.numeric() %>% sort()
  }
  mat <- as.data.frame(mat)
  apply(mat[as.character(target_year)],1,get(stat))
}


# VPAや将来予測の結果をtidyに ----
convert_df <- function(df,name){
  df %>%
    as_tibble %>%
    mutate(age = as.numeric(rownames(df))) %>%
    gather(key=year, value=value, -age, convert=TRUE) %>%
    group_by(year) %>%
    #        summarise(value=sum(value)) %>%
    mutate(type="VPA",sim="s0",stat=name)
}

#'
#' @export
#' @encoding UTF-8
convert_2d_future <- function(df, name, label="tmp"){
  df %>%
    as_tibble %>%
    mutate(year=rownames(df)) %>%
    gather(key=sim, value=value, -year, convert=TRUE) %>%
    mutate(year=as.numeric(year), stat=name, label=label)
}

#' future.vpaの結果オブジェクトをtibble形式に変換する関数
#'
#' @param fout future.vpaの結果のオブジェクト
#'
#' @encoding UTF-8
#' @export
#'

convert_future_table <- function(fout,label="tmp"){
  
  U_table <- fout$vwcaa/fout$vbiom
  if(is.null(fout$Fsakugen)) fout$Fsakugen <- -(1-fout$faa[1,,]/fout$currentF[1])
  if(is.null(fout$recruit))  fout$recruit <- fout$naa[1,,]
  
  ssb      <- convert_2d_future(df=fout$vssb,   name="SSB",     label=label)
  catch    <- convert_2d_future(df=fout$vwcaa,  name="catch",   label=label)
  biomass  <- convert_2d_future(df=fout$vbiom,  name="biomass", label=label)
  U_table  <- convert_2d_future(df=U_table,     name="U",       label=label)
  beta_gamma    <- convert_2d_future(df=fout$alpha,  name="beta_gamma",   label=label)
  Fsakugen <- convert_2d_future(df=fout$Fsakugen, name="Fsakugen",   label=label)
  recruit  <- convert_2d_future(df=fout$recruit, name="Recruitment",   label=label)
  if(!is.null(fout$Fratio)){
    Fratio <- convert_2d_future(df=fout$Fratio, name="Fratio",   label=label)
  }
  else{
    Fratio <- NULL
  }
  
  Fsakugen_ratio <- Fsakugen %>%
    mutate(value=value+1)
  Fsakugen_ratio$stat <- "Fsakugen_ratio"
  
  bind_rows(ssb,catch,biomass,beta_gamma,Fsakugen,Fsakugen_ratio,recruit, U_table, Fratio)
}


convert_vector <- function(vector,name){
  vector %>%
    as_tibble %>%
    mutate(year = as.integer(names(vector))) %>%
    mutate(type="VPA",sim="s0",stat=name,age=NA)
}

#' VPAの結果オブジェクトをtibble形式に変換する関数
#'
#' @param vpares vpaの結果のオブジェクト
#' @encoding UTF-8
#'
#'
#' @export

convert_vpa_tibble <- function(vpares,SPRtarget=NULL){
  
  if (is.null(vpares$input$dat$waa.catch)) {
    total.catch <- colSums(vpares$input$dat$caa*vpares$input$dat$waa,na.rm=T)
  } else {
    total.catch <- colSums(vpares$input$dat$caa*vpares$input$dat$waa.catch,na.rm=T)
  }
  U <- total.catch/colSums(vpares$baa, na.rm=T)
  
  SSB <- convert_vector(colSums(vpares$ssb,na.rm=T),"SSB") %>%
    dplyr::filter(value>0&!is.na(value))
  Biomass <- convert_vector(colSums(vpares$baa,na.rm=T),"biomass") %>%
    dplyr::filter(value>0&!is.na(value))
  FAA <- convert_df(vpares$faa,"fishing_mortality") %>%
    dplyr::filter(value>0&!is.na(value))
  Recruitment <- convert_vector(colSums(vpares$naa[1,,drop=F]),"Recruitment") %>%
    dplyr::filter(value>0&!is.na(value))
  
  if(!is.null(SPRtarget)){
    if(is.null(vpares$input$dat$waa.catch)) waa.catch <- vpares$input$dat$waa
    else waa.catch <- vpares$input$dat$waa.catch
    Fratio <- purrr::map_dfc(1:ncol(vpares$naa),
                             function(i){
                               tmp <- !is.na(vpares$faa[,i])
                               calc_Fratio(faa=vpares$faa[tmp,i],
                                           maa=vpares$input$dat$maa[tmp,i],
                                           waa=vpares$input$dat$waa[tmp,i],
                                           M  =vpares$input$dat$M[tmp,i],
                                           waa.catch=waa.catch[tmp],
                                           SPRtarget=SPRtarget)
                             })
    colnames(Fratio) <- colnames(vpares$naa)
    Fratio <- convert_df(Fratio,"Fratio")
  }
  else{
    Fratio <- NULL
  }
  
  all_table <- bind_rows(SSB,
                         Biomass,
                         convert_vector(U[U>0],"U"),
                         convert_vector(total.catch[total.catch>0],"catch"),
                         convert_df(vpares$naa,"fish_number"),
                         FAA,
                         convert_df(vpares$input$dat$waa,"weight"),
                         convert_df(vpares$input$dat$maa,"maturity"),
                         convert_df(vpares$input$dat$caa,"catch_number"),
                         convert_df(vpares$input$dat$M,  "natural_mortality"),                         
                         Recruitment,
                         Fratio)
}

#' fit.SRの結果をtibble形式に治す
#'
#' @param SR_result fit.SRの結果のオブジェクト
#' @encoding UTF-8
#'
#' @export
#'

convert_SR_tibble <- function(res_SR){
  bind_rows(tibble(value=as.numeric(res_SR$pars),type="parameter",name=names(res_SR$pars)),
            res_SR$pred %>% mutate(type="prediction",name="prediction"),
            res_SR$input$SRdata %>% as_tibble() %>%
              mutate(type="observed",name="observed",residual=res_SR$resid))
}

#' 管理基準値の表を作成する
#'
#' @param refs_base est.MSYから得られる管理基準値の表
#' @encoding UTF-8
#'
#' @export
#'

make_RP_table <- function(refs_base){
  #    require(formattable)
  #    require(tidyverse,quietly=TRUE)
  table_output <- refs_base %>%
    select(-RP_name) %>% # どの列を表示させるか選択する
    # 各列の有効数字を指定
    mutate(SSB=round(SSB,-floor(log10(min(SSB)))),
           SSB2SSB0=round(SSB2SSB0,2),
           Catch=round(Catch,-floor(log10(min(Catch)))),
           Catch.CV=round(Catch.CV,2),
           U=round(U,2),
           Fref2Fcurrent=round(Fref2Fcurrent,2)) %>%
    rename("管理基準値"=RP.definition,"親魚資源量"=SSB,"B0に対する比"=SSB2SSB0,
           "漁獲量"=Catch,"漁獲量の変動係数"=Catch.CV,"漁獲率"=U,"努力量の乗数"=Fref2Fcurrent)
  
  table_output  %>%
    # 表をhtmlで出力
    formattable::formattable(list(親魚資源量=color_bar("olivedrab"),
                                       漁獲量=color_bar("steelblue"),
                                       漁獲率=color_bar("orange"),
                                       努力量の乗数=color_bar("tomato")))
  
  #    return(table_output)
  
}

#' 管理基準値表から目的の管理基準値を取り出す関数
#'
#' @param refs_base est.MSYから得られる管理基準値の表
#' @param RP_name 取り出したい管理基準値の名前
#' @encoding UTF-8
#'
#' @export
#'

derive_RP_value <- function(refs_base,RP_name){
  refs_base[refs_base$RP.definition%in%RP_name,]
}


# kobe II matrix など、パフォーマンスを計算する関数 ----

#'
#' @export
#'
#' @encoding UTF-8

make_kobeII_table <- function(kobeII_data,
                              res_vpa,
                              year.catch=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.ssb=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.Fsakugen=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.ssbtarget=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.ssblimit=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.ssbban=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.ssbmin=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.ssbmax=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.aav=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              year.catchdiff=(max(as.numeric(colnames(res_vpa$naa)))+1:10),
                              Btarget=0,
                              Blimit=0,
                              Bban=0){
  # 平均漁獲量
  (catch.mean <- kobeII_data %>%
     dplyr::filter(year%in%year.catch,stat=="catch") %>% # 取り出す年とラベル("catch")を選ぶ
     group_by(HCR_name,beta,year) %>%
     summarise(catch.mean=mean(value)) %>%  # 値の計算方法を指定（漁獲量の平均ならmean(value)）
     # "-3"とかの値で桁数を指定
     spread(key=year,value=catch.mean) %>% ungroup() %>%
     arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
     mutate(stat_name="catch.mean"))
  
  # 平均親魚
  (ssb.mean <- kobeII_data %>%
      dplyr::filter(year%in%year.ssb,stat=="SSB") %>%
      group_by(HCR_name,beta,year) %>%
      summarise(ssb.mean=mean(value)) %>%
      spread(key=year,value=ssb.mean) %>% ungroup() %>%
      arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
      mutate(stat_name="ssb.mean"))
  
  # 親魚, 下10%
  (ssb.ci10 <- kobeII_data %>%
      dplyr::filter(year%in%year.ssb,stat=="SSB") %>%
      group_by(HCR_name,beta,year) %>%
      summarise(ssb.ci10=quantile(value,probs=0.1)) %>%
      spread(key=year,value=ssb.ci10) %>% ungroup() %>%
      arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
      mutate(stat_name="ssb.ci10"))
  
  # 親魚, 上10%
  (ssb.ci90 <- kobeII_data %>%
      dplyr::filter(year%in%year.ssb,stat=="SSB") %>%
      group_by(HCR_name,beta,year) %>%
      summarise(ssb.ci90=quantile(value,probs=0.9)) %>%
      spread(key=year,value=ssb.ci90) %>% ungroup() %>%
      arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
      mutate(stat_name="ssb.ci90"))
  
  # 1-currentFに乗じる値=currentFからの努力量の削減率の平均値（実際には確率分布になっている）
  (Fsakugen.table <- kobeII_data %>%
      dplyr::filter(year%in%year.Fsakugen,stat=="Fsakugen") %>% # 取り出す年とラベル("catch")を選ぶ
      group_by(HCR_name,beta,year) %>%
      summarise(Fsakugen=round(mean(value),2)) %>%
      spread(key=year,value=Fsakugen) %>% ungroup() %>%
      arrange(HCR_name,desc(beta)) %>% # HCR_nameとbetaの順に並び替え
      mutate(stat_name="Fsakugen"))
  
  # SSB>SSBtargetとなる確率
  ssbtarget.table <- kobeII_data %>%
    dplyr::filter(year%in%year.ssbtarget,stat=="SSB") %>%
    group_by(HCR_name,beta,year) %>%
    summarise(ssb.over=round(100*mean(value>Btarget))) %>%
    spread(key=year,value=ssb.over) %>%
    ungroup() %>%
    arrange(HCR_name,desc(beta))%>%
    mutate(stat_name="Pr(SSB>SSBtarget)")
  
  # SSB>SSBlimとなる確率
  ssblimit.table <- kobeII_data %>%
    dplyr::filter(year%in%year.ssblimit,stat=="SSB") %>%
    group_by(HCR_name,beta,year) %>%
    summarise(ssb.over=round(100*mean(value>Blimit))) %>%
    spread(key=year,value=ssb.over)%>%
    ungroup() %>%
    arrange(HCR_name,desc(beta))%>%
    mutate(stat_name="Pr(SSB>SSBlim)")
  
  # SSB>SSBbanとなる確率
  ssbban.table <- kobeII_data %>%
    dplyr::filter(year%in%year.ssbban,stat=="SSB") %>%
    group_by(HCR_name,beta,year) %>%
    summarise(ssb.over=round(100*mean(value>Bban))) %>%
    spread(key=year,value=ssb.over)%>%
    ungroup() %>%
    arrange(HCR_name,desc(beta))%>%
    mutate(stat_name="Pr(SSB>SSBban)")
  
  # SSB>SSBmin(過去最低親魚量を上回る確率)
  ssb.min <- min(unlist(colSums(res_vpa$ssb, na.rm=T)))
  ssbmin.table <- kobeII_data %>%
    dplyr::filter(year%in%year.ssbmin,stat=="SSB") %>%
    group_by(HCR_name,beta,year) %>%
    summarise(ssb.over=round(100*mean(value>ssb.min))) %>%
    spread(key=year,value=ssb.over)%>%
    ungroup() %>%
    arrange(HCR_name,desc(beta))%>%
    mutate(stat_name="Pr(SSB>SSBmin)")
  
  # SSB>SSBmax(過去最低親魚量を上回る確率)
  ssb.max <- max(unlist(colSums(res_vpa$ssb, na.rm=T)))
  ssbmax.table <- kobeII_data %>%
    dplyr::filter(year%in%year.ssbmax,stat=="SSB") %>%
    group_by(HCR_name,beta,year) %>%
    summarise(ssb.over=round(100*mean(value>ssb.max))) %>%
    spread(key=year,value=ssb.over)%>%
    ungroup() %>%
    arrange(HCR_name,desc(beta))%>%
    mutate(stat_name="Pr(SSB>SSBmax)")
  
  # オプション: Catch AAV mean
  calc.aav <- function(x)sum(abs(diff(x)))/sum(x[-1])
  catch.aav.table <- kobeII_data %>%
    dplyr::filter(year%in%year.aav,stat=="catch") %>%
    group_by(HCR_name,beta,sim) %>%
    dplyr::summarise(catch.aav=(calc.aav(value))) %>%
    group_by(HCR_name,beta) %>%
    summarise(catch.aav.mean=mean(catch.aav)) %>%
    arrange(HCR_name,desc(beta))%>%
    mutate(stat_name="catch.aav")
  
  res_list <- list(catch.mean   = catch.mean,
                   ssb.mean         = ssb.mean,
                   ssb.lower10percent            = ssb.ci10,
                   ssb.upper90percent            = ssb.ci90,
                   prob.over.ssbtarget  = ssbtarget.table,
                   prob.over.ssblimit   = ssblimit.table,
                   prob.over.ssbban     = ssbban.table,
                   prob.over.ssbmin     = ssbmin.table,
                   prob.over.ssbmax     = ssbmax.table,
                   catch.aav       = catch.aav.table)
  return(res_list)
  
}


#' kobeII matrixの簡易版（Btarget, Blimitは決め打ちでβのみ変える)
#'
#' @encoding UTF-8
#' @export
#'
#'

beta.simulation <- function(finput,beta_vector,year.lag=0,type="old"){
  
  tb <- NULL
  
  for(i in 1:length(beta_vector)){
    if(type=="old"){
      finput$HCR$beta <- beta_vector[i]
      finput$is.plot <- FALSE
      finput$silent <- TRUE
      fres_base <- do.call(future.vpa,finput) # デフォルトルールの結果→図示などに使う
    }
    else{
      finput$tmb_data$HCR_mat[,,"beta"] <- beta_vector[i]
      if(!is.null(finput$MSE_input_data)) finput$MSE_input_data$input$HCR_beta <- beta_vector[i]
      fres_base <- do.call(future_vpa,finput) # デフォルトルールの結果→図示などに使う
      fres_base <- format_to_old_future(fres_base)
    }
    tmp <- convert_future_table(fres_base,label=beta_vector[i]) %>%
      rename(HCR_name=label)  %>% mutate(beta=beta_vector[i])
    tb <- bind_rows(tb,tmp)
  }
  return(tb)
}


#' MSYを達成するときの\%SPRを計算する
#'
#' @param finput 将来予測インプット
#' @param fout 将来予測のアウトプット（finputがない場合)
#' @param Fvector Fのベクトル
#' @encoding UTF-8
#' @export
calc_perspr <- function(finput=NULL,
                        fout=NULL,
                        res_vpa=NULL,
                        Fvector,
                        Fmax=10,
                        max.age=Inf,
                        target.col=NULL
){
  if(!is.null(finput)){
    # MSYにおける将来予測計算をやりなおし
    finput$outtype <- "FULL"
    fout.tmp <- do.call(future.vpa,finput)
    res_vpa <- finput$res0
  }
  else{
    fout.tmp <- fout
  }
  # 生物パラメータはその将来予測で使われているものを使う
  if(is.null(target.col)){
    waa.tmp           <- fout.tmp$waa[,dim(fout.tmp$waa)[[2]],1]
    waa.catch.tmp <- fout.tmp$waa.catch[,dim(fout.tmp$waa.catch)[[2]],1]
    maa.tmp           <- fout.tmp$maa[,dim(fout.tmp$maa)[[2]],1]
    M.tmp                <- fout.tmp$M[,dim(fout.tmp$M)[[2]],1]
  }
  else{
    waa.tmp           <- fout.tmp$waa[,target.col,1]
    waa.catch.tmp <- fout.tmp$waa.catch[,target.col,1]
    maa.tmp           <- fout.tmp$maa[,target.col,1]
    M.tmp               <- fout.tmp$M[,target.col,1]
  }
  
  # 緊急措置。本来ならどこをプラスグループとして与えるかを引数として与えないといけない
  allsumpars <- waa.tmp+waa.catch.tmp+maa.tmp+M.tmp
  waa.tmp <- waa.tmp[allsumpars!=0]
  waa.catch.tmp <- waa.catch.tmp[allsumpars!=0]
  maa.tmp <- maa.tmp[allsumpars!=0]
  M.tmp <- M.tmp[ allsumpars!=0]     
  Fvector <- Fvector %>%  as.numeric()
  Fvector <- Fvector[allsumpars!=0]
  ## ここまで緊急措置
  
  # SPRを計算
  spr.current <- ref.F(res_vpa,Fcurrent=Fvector,
                       waa=waa.tmp,
                       waa.catch=waa.catch.tmp,pSPR=NULL,
                       maa=maa.tmp,M=M.tmp,rps.year=as.numeric(colnames(res_vpa$naa)),
                       F.range=c(seq(from=0,to=ceiling(max(res_vpa$Fc.at.age,na.rm=T)*Fmax),
                                     length=101),max(res_vpa$Fc.at.age,na.rm=T)),
                       plot=FALSE,max.age=max.age)$currentSPR$perSPR
  spr.current
}

#' kobeIItable から任意の表を指名して取り出す
#'
#' @param kobeII_table \code{make_kobeII_table}の出力
#' @param name \code{kobeII_table}の要素名
#'
#' @encoding UTF-8
pull_var_from_kobeII_table <- function(kobeII_table, name) {
  table <- kobeII.table[[name]]
  table %>%
    dplyr::arrange(desc(beta)) %>%
    dplyr::select(-HCR_name, -stat_name)
}

#' kobeIItableから取り出した表を整形
#'
#' - 報告書に不要な列を除去する
#' - 単位を千トンに変換
#' @param beta_table \code{pull_var_from_kobeII_table}で取得した表
#' @param divide_by 表の値をこの値で除する．トンを千トンにする場合には1000
#' @param round TRUEなら値を丸める．漁獲量は現状整数表示なのでデフォルトはTRUE
format_beta_table <- function(beta_table, divide_by = 1, round = TRUE) {
  beta   <- beta_table %>%
    dplyr::select(beta) %>%
    magrittr::set_colnames("\u03B2") # greek beta in unicode
  values <- beta_table %>%
    dplyr::select(-beta) / divide_by
  if (round == TRUE) return(cbind(beta, round(values)))
  cbind(beta, values)
}

#' 値の大小に応じて表の背景にグラデーションをつける
#' @param beta_table \code{format_beta_table}で整形したβの表
#' @param color 表の背景となる任意の色
colorize_table <- function(beta_table, color) {
  beta_table %>%
    formattable::formattable(list(formattable::area(col = -1) ~
                                    formattable::color_tile("white", color)))
}

#' 表を画像として保存
#'
#' # @inheritParams \code{\link{formattable::as.htmlwidget}}
#' # @inheritParams \code{\link{htmltools::html_print}}
#' # @inheritParams \code{\link{webshot::webshot}}
#' @param table ファイルとして保存したい表
#' @examples
#' \dontrun{
#' your_table %>%
#'  export_formattable(file = "foo.png")
#' }
#' @export
export_formattable <- function(table, file, width = "100%", height = NULL,
                               background = "white", delay = 0.1) {
  widget <- formattable::as.htmlwidget(table, width = width, height = height)
  path   <- htmltools::html_print(widget, background = background, viewer = NULL)
  url    <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot::webshot(url,
                   file = file,
                   selector = ".formattable_widget",
                   delay = delay)
}

#' kobeIItableから任意の表を取得し画像として保存
#'
#' @inheritParams \code{\link{pull_var_from_kobeII_table}}
#' @inheritParams \code{\link{format_beta_table}}
#' @inheritParams \code{\link{colorize_table}}
#' @inheritParams \code{\link{export_formattable}}
export_kobeII_table <- function(name, divide_by, color, fname, kobeII_table) {
  kobeII_table %>%
    pull_var_from_kobeII_table(name) %>%
    format_beta_table(divide_by = divide_by) %>%
    colorize_table(color) %>%
    export_formattable(fname)
}

#' β調整による管理効果を比較する表を画像として一括保存
#'
#' @inheritParams \code{\link{pull_var_from_kobeII_table}}
#' @param fname_ssb 「平均親魚量」の保存先ファイル名
#' @param fname_catch 「平均漁獲量」の保存先ファイル名
#' @param fname_ssb_above_target 「親魚量が目標管理基準値を上回る確率」の保存先ファイル名
#' @param fname_ssb_above_limit 「親魚量が限界管理基準値を上回る確率」の保存先ファイル名
#' @examples
#' \dontrun{
#' export_kobeII_tables(kobeII.table)
#' }
#' @export
export_kobeII_tables <- function(kobeII_table,
                                 fname_ssb = "tbl_ssb.png",
                                 fname_catch = "tbl_catch.png",
                                 fname_ssb_above_target = "tbl_ssb>target.png",
                                 fname_ssb_above_limit = "tbl_ssb>limit.png") {
  blue   <- "#96A9D8"
  green  <- "#B3CE94"
  yellow <- "#F1C040"
  
  purrr::pmap(list(name = c("ssb.mean", "catch.mean",
                            "prob.over.ssbtarget", "prob.over.ssblimit"),
                   divide_by = c(1000, 1000, 1, 1),
                   color = c(blue, green, yellow, yellow),
                   fname = c(fname_ssb, fname_catch,
                             fname_ssb_above_target, fname_ssb_above_limit)),
              .f = export_kobeII_table,
              kobeII_table = kobeII_table)
}




#'
#' 将来予測の結果のリストを入れると、代表的なパフォーマンス指標をピックアップする
#'
#' @param future_list future_vpaまたはfuture.vpaの返り値のリスト
#' @param res_vpa vpaの返り値
#' @param ABC_year 特にABCに注目したい年
#' @param is_MSE MSEの結果を使うかどうか。MSEの結果の場合には、ABCや加入尾数の真の値との誤差も出力する
#' @param indicator_year ABC_yearから相対的に何年後を指標として取り出すか
#' @param Btarget 目標管理基準値の値
#' @param Blimit 限界管理基準値の値
#' @param Bban 禁漁水準の値
#' @param type 出力の形式。"long"は縦長（ggplotに渡すときに便利）、"wide"は横長（数字を直接比較するときに便利）
#' @param biomass.unit 資源量の単位
#'
#' @export
#' @encoding UTF-8 
#'

get_performance <- function(future_list,res_vpa,ABC_year=2021,
                            is_MSE=FALSE,
                            indicator_year=c(0,5,10),Btarget=0, Blimit=0, Bban=0,
                            type=c("long","wide")[1],biomass.unit=10000,...){
  
  future_list_original <- future_list
  future_list <- purrr::map(future_list,
                            function(x) if(class(x)=="future_new")
                              format_to_old_future(x) else x)
  
  if(is.null(names(future_list))) names(future_list) <- 1:length(future_list)
  
  future_tibble <- purrr::map_dfr(1:length(future_list),
                                  function(i) convert_future_table(future_list[[i]],
                                                                   label=names(future_list)[i]) %>%
                                    rename(HCR_name=label) %>% mutate(beta=NA))
  
  kobe_res <- make_kobeII_table(future_tibble,res_vpa,
                                year.ssb   = ABC_year+indicator_year,
                                year.catch = ABC_year+indicator_year,
                                year.ssbtarget = ABC_year+indicator_year,
                                year.ssblimit  = ABC_year+indicator_year,
                                year.ssbban=NULL, year.ssbmin=NULL, year.ssbmax=NULL,
                                year.aav = c(ABC_year,ABC_year-1),
                                Btarget= Btarget,
                                Blimit = Blimit,
                                Bban   = Bban)
  
  if(isTRUE(is_MSE)){
    error_table <- purrr::map_dfr(1:length(future_list_original), function(i){
      if("SR_MSE" %in% names(future_list_original[[i]])){
        plot_bias_in_MSE(future_list_original[[i]], out="stat") %>%
          dplyr::filter(year %in% (ABC_year + indicator_year)) %>%
          group_by(year, stat) %>%
          summarise(mean_error=mean(Relative_error_normal)) %>%
          mutate(HCR_name=names(future_list)[i])
      }
      else{
        NULL
      }
    })
    error_table <- error_table %>%
      gather(key=stat_name,value=value,-HCR_name,-year,-stat) %>%
      mutate(unit="",stat_category="推定バイアス") %>%
      mutate(stat_year_name=str_c(stat_category,year)) %>%
      ungroup(year) %>%
      mutate(year=as.character(year)) %>%
      mutate(stat_name=stat) %>%
      select(-stat)
  }
  else{
    error_table <- NULL
  }
  
  junit <- c("","十","百","千","万")[log10(biomass.unit)+1]
  
  stat_data <- tibble(stat_name=c("ssb.mean","catch.mean","Pr(SSB>SSBtarget)","Pr(SSB>SSBlim)",
                                  "catch.aav"),
                      stat_category=c("平均親魚量 ", "平均漁獲量 ", "目標上回る確率 ", "限界上回る確率 ",
                                      "漁獲量変動"))
  
  kobe_res <- purrr::map_dfr(kobe_res[c("ssb.mean", "catch.mean", "prob.over.ssbtarget",
                                        "prob.over.ssblimit", "catch.aav")],
                             function(x) x %>% select(-beta) %>%
                               gather(key=year,value=value,-HCR_name,-stat_name)) %>%
    mutate(value=ifelse(stat_name %in% c("ssb.mean", "catch.mean"), value/biomass.unit, value)) %>%
    mutate(unit =ifelse(stat_name %in% c("ssb.mean", "catch.mean"), str_c(junit, "トン"), "%")) %>%
    mutate(unit =ifelse(stat_name %in% c("catch.aav"), "", unit)) %>%
    left_join(stat_data) %>%
    mutate(stat_year_name=str_c(stat_category,year))
  
  kobe_res <- bind_rows(kobe_res,error_table)
  
  if(type=="wide"){
    kobe_res <- kobe_res  %>%
      select(-year, -stat_name, -stat_category) %>%
      spread(key=HCR_name,value=value) %>%
      select(2:ncol(.),1)
  }
  
  return(tibble::lst(kobe_res,error_table))
}

#'
#' 短期的将来予測における複数の管理方策のパフォーマンスを比較する表を出力する
#'
#' @param future_list 将来予測の結果のリスト
#' @param res_vpa VPAの結果
#' @param ... get_performanceで必要な引数
#'
#' @encoding UTF-8
#' @export
#'

compare_future_performance <- function(future_list,res_vpa,res_MSY,
                                       biomass.unit=1000,is_MSE=FALSE,...){
  perform_res <- get_performance(future_list=future_list, res_vpa=res_vpa,
                                 Btarget=derive_RP_value(res_MSY$summary,"Btarget0")$SSB,
                                 Blimit =derive_RP_value(res_MSY$summary,"Blimit0")$SSB,
                                 Bban   =derive_RP_value(res_MSY$summary,"Bban0")$SSB,
                                 type="long",is_MSE=is_MSE,biomass.unit=biomass.unit,...)
  
  g1_ssb0 <- perform_res$kobe_res %>% dplyr::filter(stat_name=="ssb.mean") %>%
    ggplot() +
    geom_bar(aes(x=HCR_name,y=value,fill=stat_category),stat="identity") +
    facet_wrap(stat_category~year,ncol=1)+coord_flip()+
    geom_label(aes(x=HCR_name,y=max(value)/2,label=str_c(round(value),unit)),
               alpha=0.5)+
    theme_SH() + theme(legend.position="top")+xlab("Senario") +
    guides(fill=guide_legend(title=""))+
    scale_fill_manual(values=c("lightblue"))
  
  g1_ssb <- g1_ssb0 +
    geom_hline(yintercept=derive_RP_value(res_MSY$summary,"Btarget0")$SSB/biomass.unit,
               col="#00533E")+
    geom_hline(yintercept=derive_RP_value(res_MSY$summary,"Blimit0")$SSB/biomass.unit,
               col="#edb918")
  
  g1_catch0 <- g1_ssb0 %+% dplyr::filter(perform_res$kobe_res, stat_name=="catch.mean")+
    scale_x_discrete(labels=rep("",4))+
    scale_fill_manual(values=c("lightgreen"))
  
  g1_catch <- g1_catch0 +
    geom_hline(yintercept=derive_RP_value(res_MSY$summary,"Btarget0")$Catch/10000,
               col="#00533E")+
    geom_hline(yintercept=rev(colSums(res_vpa$wcaa,na.rm=T))[1]/biomass.unit,
               col="gray",lty=2)
  
  g1_probtar <- g1_catch0 %+% dplyr::filter(perform_res$kobe_res, stat_name=="Pr(SSB>SSBtarget)")+
    scale_fill_manual(values=c("lightblue"))+
    geom_hline(yintercept=c(0,50,100),col="gray",lty=2)
  
  g1_problim <- g1_catch0 %+% dplyr::filter(perform_res$kobe_res, stat_name=="Pr(SSB>SSBlim)")+
    scale_fill_manual(values=c("gray"))+
    geom_hline(yintercept=c(0,50,100),col="gray",lty=2)
  
  g1_error_table <- perform_res$error_table  %>%
    ggplot() +
    geom_bar(aes(x=HCR_name,y=value,fill=stat_category),stat="identity") +
    facet_grid(year~stat_name)+coord_flip()+
    geom_label(aes(x=HCR_name,y=max(value)/2,label=str_c(round(value,2),unit)),
               alpha=0.5)+
    theme_SH() + theme(legend.position="top")+xlab("Senario") +
    guides(fill=guide_legend(title=""))+
    scale_fill_manual(values=c("lightpink"))
  
  g1_performance <- gridExtra::marrangeGrob(list(g1_ssb,g1_probtar,g1_catch,g1_problim),
                                            widths=c(1.3,1,1,1),nrow=1,ncol=4)
  
  list(g1_performance, g1_error_table, perform_res)
  #ggsave(g1_performance, filename="g1_performance.png" ,path=output_folder,
  #           width=20, height= 10)
  
  #    ggsave(g1_error_table, filename="g1_error_table.png" ,path=output_folder,
  #       width=15, height= 8)
  
}


#'
#' calculate F/Ftarget based on F_\%SPR multiplier
#'
#' @param faa F at age
#' @param waa weight at age
#' @param maa maturity at age
#' @param M natural morality at age
#' @param SPRtarget target SPR
#'
#' @export
#' @encoding UTF-8 
#'


calc_Fratio <- function(faa, waa, maa, M, SPRtarget=30, waa.catch=NULL,Pope=TRUE){
  tmpfunc <- function(x,SPR0=0,...){
    SPR_tmp <- calc.rel.abund(sel=faa,Fr=exp(x),na=length(faa),M=M, waa=waa, waa.catch=waa.catch,
                              min.age=0,max.age=Inf,Pope=Pope,ssb.coef=0,maa=maa)$spr %>% sum()
    sum(((SPR_tmp/SPR0*100)-SPRtarget)^2)
  }
  if(sum(faa)==0){ return(NA) }
  else{
    tmp <- !is.na(faa)
    SPR0 <- calc.rel.abund(sel=faa,Fr=0,na=length(faa),M=M, waa=waa, waa.catch=waa.catch,maa=maa,
                           min.age=0,max.age=Inf,Pope=Pope,ssb.coef=0)$spr %>% sum()
    opt_res <- optimize(tmpfunc,interval=c(-10,10),SPR0=SPR0)
    SPR_est <- calc.rel.abund(sel=faa,Fr=exp(opt_res$minimum),na=length(faa),
                              M=M, waa=waa, waa.catch=waa.catch,maa=maa,
                              min.age=0,max.age=Inf,Pope=Pope,ssb.coef=0)$spr %>% sum()
    SPR_est <- SPR_est/SPR0 * 100
    if(abs(SPR_est-SPRtarget)>0.01) {browser(); return(NA)}
    else return(1/exp(opt_res$minimum))
  }
}

#'
#' @export
#'
#'

calc_akaike_weight <- function(AIC) exp(-AIC/2)/sum(exp(-AIC/2))

#'
#' 使うフォルダ名を与えると一連の結果の関数を読み込む関数
#' 
#' @export
#' @encoding UTF-8 
#'

load_folder <- function(folder_name){

  tmpfunc <- function(folder_name, file_name){
    if(isTRUE(file.exists(str_c(folder_name,"/",file_name)))){
      if(isTRUE(str_detect(file_name, pattern=".rda"))){
        a <- load(str_c(folder_name,"/",file_name))
        a <- get(a)
      }
      if(isTRUE(str_detect(file_name, pattern=".csv"))){
        a <- read_csv(str_c(folder_name,"/",file_name))
      }      
    }
    else{
      a <- NA
    }
    return(a)
  }
    
  file_name <- c("res_MSY.rda","res_SR.rda","res_future_0.8HCR.rda","kobeII.table.rda","model_selection.csv")

  res_all <- list()
  for(i in 1:length(file_name)){
    res_all[[i]] <- purrr::map(folder_name, function(x) tmpfunc(x,file_name[i]))
    res_all[[i]] <- res_all[[i]][!is.na(res_all[[i]])]
  }
  names(res_all) <- file_name

  res_all$res_vpa <- purrr::map(res_all$res_MSY.rda, function(x) if(!is.na(x)) x$res_vpa else NA)
  res_all$res_vpa <- res_all$res_vpa[!is.na(res_all$res_vpa)]
  invisible(res_all)
}


#' Make Kobe ratio (?) for frasyr_tool
#'
#' @param result_vpa object createb dy vpa()
#' @param result_msy object created by script 1do_MSYest.R of SC meeting
#' @return A tibble object
#' @export
make_kobe_ratio <- function(result_vpa, result_msy) {

  assertthat::assert_that(
    assertthat::has_name(result_vpa, c("ssb")),
    assertthat::has_name(result_msy, c("summary")),
    assertthat::has_name(result_msy$summary, c("perSPR"))
  )

  return_kobe_ratio <- function() {
    tibble::tibble(year   = get_year(),
                   Fratio = get_f_ratio(),
                   Bratio = calc_b_ratio()) %>%
      ad_hoc_filter()
  }

  get_year <- function() {
    colnames(result_vpa$ssb)
  }

  get_f_ratio <- function() {
    target_spr  <- derive_RP_value(result_msy$summary,"Btarget0")$perSPR * 100
    spr_history <- get.SPR(result_vpa,
                           target_spr,
                           max.age = Inf, Fmax = 1)

    assertthat::assert_that(
      assertthat::validate_that(is.list(spr_history)),
      assertthat::has_name(spr_history, "ysdata"),
      assertthat::has_name(spr_history$ysdata, "F/Ftarget"),
      assertthat::are_equal(rownames(spr_history$ysdata),
                            colnames(result_vpa$ssb))
    )

    force(spr_history$ysdata$`F/Ftarget`)
  }


  calc_b_ratio <- function() {
    ssb        <- colSums(result_vpa$ssb, na.rm=T)
    target_ssb <- derive_RP_value(result_msy$summary,"Btarget0")$SSB

    force(ssb / target_ssb)
  }

  ad_hoc_filter <- function(koberatio) {
    dplyr::filter(koberatio, !is.na(Bratio)) # When does Bratio become 'NA'?
  }

  return_kobe_ratio()
}
