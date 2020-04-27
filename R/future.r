#'
#' @import ggplot2
#' @import magrittr          
#' @import dplyr             
#' @import tidyr             
#' @import tibble            
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' 
NULL

# 個体群動態の記述に使われる基本的関数 ----

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

# 結果のプロットのための関数 ----

#' ref.Fの出力をプロットするための関数
#'
#' @param rres ref.Fの出力結果
#' @encoding UTF-8
#'
#' @export

plot_Fref <- function(rres,xlabel="max", # or, "mean","Fref/Fcur"
                      vline.text=c("F0.1","Fmax","Fcurrent","Fmed") # and "FpSPR.20.SPR" etc..
){
  old.par <- par()
  par(mar=c(4,4,1,4))
  F.range <- rres$ypr.spr$F.range
  if(xlabel=="Fref/Fcur") F.range <- F.range/rres$summary$Fcurrent[1]*rres$summary$Fcurrent[3]
  if(xlabel=="mean") F.range <- F.range/rres$summary$Fcurrent[1]*rres$summary$Fcurrent[2]    
  spr <- rres$ypr.spr$pspr
  ypr <- rres$ypr.spr$ypr
  plot(F.range,spr,xlab=xlabel,ylab="%SPR",type="l",ylim=c(0,max(spr)))
  par(new=T)
  plot(F.range,ypr,axes=F,xlab="",ylab="",lty=2,type="l",ylim=c(0,max(ypr)))
  axis(side=4)
  mtext("YPR",side=4,line=2)
  n.line <- which(rownames(rres$summary) %in% xlabel)
  abline(v=xx <- c(rres$summary[vline.text][n.line,]))
  text(xx,max(ypr)*seq(from=0.5,to=0.3,length=length(vline.text)),vline.text)
  legend("topright",lty=1:2,legend=c("SPR","YPR"))
  old.par[c("cin","cxy","csi","cra","csy","din","page")] <- NULL
  par(old.par)
}




#' 複数の将来予測の結果をプロットする（ggplotは使わず）
#'
#' @param fres.list future.vpaからの出力結果をリストで並べたもの
#' @encoding UTF-8
#' 
#' 
#'
#' @export

plot.futures <- function(fres.list,conf=c(0.1,0.5,0.9),target="SSB",legend.text="",xlim.tmp=NULL,y.scale=1,det.run=TRUE){
  
  if(legend.text=="") legend.text <- names(fres.list)
  if(is.null(legend.text)) legend.text <- 1:length(fres.list)
  
  for(i in 1:length(fres.list)){
    if(class(fres.list[[i]])=="future_new")
      fres.list[[i]] <- format_to_old_future(fres.list[[i]])
    det.run <- FALSE
  }
  
  if(isTRUE(det.run)) select_col <- -1  else select_col <- TRUE    
  
  if(target=="SSB")  aa <- lapply(fres.list,function(x) apply(x$vssb[,select_col],1,quantile,probs=conf))
  if(target=="Biomass") aa <- lapply(fres.list,function(x) apply(x$vbiom[,select_col],1,quantile,probs=conf))
  if(target=="Catch") aa <- lapply(fres.list,function(x) apply(x$vwcaa[,select_col],1,quantile,probs=conf))
  if(target=="Recruit"){
    if(is.null(x$recruit)) x$recruit <- x$naa
    aa <- lapply(fres.list,function(x) apply(x$recruit[,select_col],1,quantile,probs=conf))
  }
  
  if(is.null(xlim.tmp)) xlim.tmp <- as.numeric(range(unlist(sapply(aa,function(x) colnames(x)))))
  plot(0,max(unlist(aa)),type="n",xlim=xlim.tmp,
       ylim=y.scale*c(0,max(unlist(aa))),xlab="Year",ylab=target)
  lapply(1:length(aa),function(i) matpoints(colnames(aa[[i]]),t(aa[[i]]),col=i,type="l",lty=c(2,1,2)))
  legend("bottomright",col=1:length(fres.list),legend=legend.text,lty=1)
  invisible(aa)
}

#' 一つの将来予測の結果をプロットする（ggplotは使わず）
#'
#' @param fres0 future.vpaからの出力結果
#' @encoding UTF-8
#' 
#' 
#'
#' @export

plot.future <- function(fres0,ylim.tmp=NULL,xlim.tmp=NULL,vpares=NULL,what=c(TRUE,TRUE,TRUE),conf=0.1,N.line=0,det.run=TRUE,
                        label=c("Biomass","SSB","Catch"),is.legend=TRUE,add=FALSE,col=NULL,...){
  ## 暗黙に、vssbなどのmatrixの1列目は決定論的なランの結果と仮定している
  
  if(class(fres0)=="future_new"){
    fres0 <- format_to_old_future(fres0)
    det.run <- FALSE
  }
  
  if(is.null(col)) col <- 1                        
  matplot2 <- function(x,add=FALSE,...){
    if(add==FALSE) matplot(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)
    if(add==TRUE) matpoints(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)
  }
  
  if(is.null(xlim.tmp)) xlim.tmp <- range(as.numeric(rownames(fres0$vssb)))
  
  if(isTRUE(det.run)) select_col <- -1  else select_col <- TRUE
  
  if(what[1]){
    matplot2(x <- t(apply(fres0$vbiom[,select_col],1,
                          quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[1])),
             xlim=xlim.tmp,
             ylab=label[1],main=label[1],add=add,...)
    points(rownames(fres0$vbiom),apply(fres0$vbiom[,select_col],1,mean),type="b",pch=1)
    if(isTRUE(det.run)) points(rownames(fres0$vbiom),as.numeric(fres0$vbiom[,1]),type="b",pch=3)
    if(!is.null(vpares)){
      points(colnames(vpares$baa),colSums(vpares$baa),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vbiom),fres0$vbiom[,2:(N.line+1)],col="gray",type="l",lty=1)
  }
  
  if(what[2]){
    matplot2(x <- t(apply(fres0$vssb[,select_col],1,quantile,
                          probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[2])),
             xlim=xlim.tmp,           
             ylab=label[2],main=label[2],add=add,...)
    points(rownames(fres0$vssb),apply(fres0$vssb[,select_col],1,mean),type="b",pch=1)    
    if(isTRUE(det.run)) points(rownames(fres0$vssb),as.numeric(fres0$vssb[,1]),type="b",pch=3)
    if(!is.null(fres0$input$Frec))
      if(!is.null(fres0$input$Frec$scenario))
        if(fres0$input$Frec$scenario!="catch.mean"){
          abline(h=fres0$input$Frec$Blimit,col=2)
          abline(v=fres0$input$Frec$future.year,col=2)            
        }
    if(!is.null(vpares)){
      points(colnames(vpares$ssb),colSums(vpares$ssb),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vssb),fres0$vssb[,2:(N.line+1)],col="gray",type="l",lty=1)
  }
  
  if(what[3]){
    matplot2(x <- t(apply(fres0$vwcaa[,select_col],1,
                          quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[3])),
             xlim=xlim.tmp,           
             ylab=label[3],main=label[3],add=add,...)
    points(rownames(fres0$vwcaa),apply(fres0$vwcaa[,select_col],1,mean),type="b",pch=1)        
    if(isTRUE(det.run)) points(rownames(fres0$vwcaa),as.numeric(fres0$vwcaa[,1]),type="b",pch=3)
    if(!is.null(fres0$input$Frec))
      if(fres0$input$Frec$scenario=="catch.mean"){
        abline(h=fres0$input$Frec$Blimit,col=2)
        abline(v=fres0$input$Frec$future.year,col=2)                    
      }    
    if(!is.null(vpares)){
      points(colnames(vpares$baa),colSums(vpares$input$dat$caa*vpares$input$dat$waa),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vwcaa),fres0$vwcaa[,2:(N.line+1)],col="gray",type="l",lty=1)    
  }
  if(is.legend){
    if(sum(what)>1) plot(1:10,type = "n",ylab = "", xlab = "", axes = F)
    legend("topleft",lty=c(NA,NA,1,2),legend=c("Deterministic","Mean","Median",paste(100-(conf*2)*100,"%conf")),pch=c(3,1,NA,NA))
  }
  
}

#' @export
make_summary_table <- function(mat_data,side=1,probs=c(0.1,0.5,0.8)){
  res_mat <- cbind(apply(mat_data,side,mean),
                   t(apply(mat_data,side,quantile,probs=probs)))
  colnames(res_mat)[1] <- "average"
  res_mat %>% as_tibble()
}


# 結果の入出力 ----

#'
#' VPA結果をcsvファイルに出力する
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

#' @export
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


#' SRdataをプロットする
#'
#' @param vpares VPAの結果のオブジェクト
#' @param frag   MSY計算時の将来予測で用いる引数のリスト
#' @param type "classic"(通常プロット) or "gg"(ggplot)
#' 
#' @encoding UTF-8
#' 
#'
#' @export
#'
#' 

plot_SRdata <- function(SRdata, type=c("classic","gg")[1]){
  if(type=="classic") plot(SRdata$SSB,SRdata$R,xlab="SSB",ylab="R",xlim=c(0,max(SRdata$SSB)),ylim=c(0,max(SRdata$R)))
  if(type=="gg") ggplot2::qplot(y=R,x=SSB,data=as_tibble(SRdata),
                                xlab="SSB",ylab="R",xlim=c(0,max(SRdata$SSB)),ylim=c(0,max(SRdata$R))) + theme_SH()
}



plotRadial <- function(index,base=1,col.tmp=NULL,lwd=2,...){
  old.par <- par()
  layout(matrix(c(1,2),2,1),heights=c(2,1))
  
  index2 <- sweep(matrix(unlist(index),nrow(index),ncol(index)),2,as.numeric(unlist(index[base,])),FUN="/")
  
  if(is.null(col.tmp)) col.tmp <- brewer.pal(nrow(index2-1),"Set1")
  
  radial.plot(index2,rp.type="p",lwd=lwd,show.grid.labels=FALSE,
              labels=colnames(index),
              radial.lim=c(0,1.5),clockwise=TRUE,start=1,
              line.col=c(NA,col.tmp),
              poly.col=c(rgb(40/255,96/255,163/255,0.2),rep(NA,nrow(index2)-1)), # MSYだけ色で塗る
              ...
  )
  refname <- rownames(index)
  par(mar=c(1,0,1,0))
  plot(0,10,type="n",axes=FALSE,ylab="")
  legend("topleft",legend=refname,
         col=c(rgb(40/255,96/255,163/255,0.2),col.tmp),
         ncol=2,lwd=c(10,rep(lwd,length(refname)-1)))
  layout(matrix(c(1),1,1),heights=c(1))
  par(old.par)
  invisible(index2)
}


## 管理基準値を取り出す関数
get.Bref <- function(res,SRfunc="hs",B0=c(0.3),SPR0=c(0.3),HS=c(1,1.3),PGY=c("PGY_0.9_upper_hs","PGY_0.9_lower_hs")){
  sumref <- res$summary[rownames(res$summary)==SRfunc,]
  refpoints <- list()
  ## MSY管理基準値をピックアップ
  refpoints$BMSY <- sumref$"SSB_MSY"
  
  ## B0基準の管理基準値はB0×％
  ## B0の値はmout$summary$"B0(SSB)"にある。１番目がHSの結果
  refpoints$B0per <- sumref$"B0(SSB)"[1] * B0 # B0_10,20,30,35,40%の値
  names(refpoints$B0per) <- paste(B0*100,"%",sep="")
  
  ## B_HS関連の管理基準値
  refpoints$BHS <- sumref$b[1] *  HS
  names(refpoints$BHS) <- paste("B_HSx",HS,sep="")
  
  ## B_PGY関連の管理基準値(HSをもとにしたものはPGY.biom.hsにあります)
  x <- res$PGY.biom.hs["ssb.mean"]
  refpoints$BPGY <- x[match(PGY,rownames(x)),1]
  names(refpoints$BPGY) <- PGY
  
  ## SSB_current
  refpoints$SSBcur <- rev(as.numeric(res$vpares$ssb))[1]
  
  ## SSB_max
  refpoints$SSBmax <- max(as.numeric(res$vpares$ssb))
  return(unlist(refpoints))
}

plot.RP <- function(rdata,RP=NULL,biomass.scale=1,ymax=1,is.text=TRUE){
  n <- length(rdata)
  rdata <- sort(rdata)
  if(is.null(RP)) RP <- names(rdata)
  ymax <- ymax * seq(from=0.5,to=1,length=n)
  for(j in 1:n){
    abline(v=rdata[j]/biomass.scale,lty=1,lwd=2,col=rgb(40/255,96/255,40/255,0.5))
    if(isTRUE(is.text)){
      text(rdata[j]/biomass.scale,ymax[j],
           paste(RP[j],"=\n",format(round(rdata[j]/biomass.scale),big.mark=","),"",sep=""),adj=0)
    }
  }
}

#'
#' @export
#' 

plot_waa <- function(vres){
  lm.list <- list()
  nage <- nrow(vres$naa)
  col.tmp <- rainbow(nage)    
  logx <- log(unlist(vres$naa))
  logy <- log(unlist(vres$input$dat$waa))
  ages <- as.numeric(rep(rownames(vres$naa),ncol(vres$naa)))
  u.age <- unique(ages)
  plot(logx,logy,col=col.tmp[1+ages],xlab="log(N)",ylab="log(weight)")
  for(i in 1:length(u.age)){
    tmp <- ages==u.age[i] & logy>-Inf & logx>-Inf
    if(sum(tmp,na.rm=TRUE)>0){
      lm.list[[i]] <- lm(logy[tmp]~logx[tmp])
      l.type <- ifelse(summary(lm.list[[i]])$coeff[2,4]<0.05,1,2)
      if(!is.na(l.type)) abline(lm.list[[i]],col=col.tmp[1+ages[i]],lty=l.type)
    }
  }
  title(vres$stockid,line=0.2)
  legend("bottomleft",lty=c(1:2,rep(1,nage)),
         col=c(1,1,col.tmp),
         legend=c("p<0.05","p>0.05",paste("Age",u.age)))    
  return(lm.list)
}


#'
#' @export
#' 

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


draw.refline <- function(reftable,horiz=TRUE,scale=1000,lwd=3){
  if(horiz==FALSE){
    abline(v=reftable[1:4,2]/scale,
           col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty="22")
    abline(v=reftable[1:4,1]/scale,
           col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty=1)
  }
  else{
    abline(h=reftable[1:4,2]/scale,
           col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty="22")
    abline(h=reftable[1:4,1]/scale,
           col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty=1)
  }
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
