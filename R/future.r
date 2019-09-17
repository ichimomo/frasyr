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
#' @param  pSPR = seq(10,90,by=10), # F%SPRを計算するときの％SPR
#' @param d 0.001
#' @param  Fem.init 経験的管理基準値(Fmed, Fmean, Fhigh, Flow)の初期値 (default=0.5)
#' @param  Fmax.init Fmaxの初期値 (default=1.5)
#' @param  F0.1.init F0.1の初期値 (default=0.7)
#' @param  iterlim 
#' @param  plot 結果のプロットを表示するかどうか
#' @param  Pope Popeの式を使うか
#' @param  F.range YPR, SPR曲線を書くときのFの範囲（Fの最大値のスケール）、かつ、F%SPRを計算するときの初期値を決めるために利用される。F%SPRの推定がうまくいかない場合はこの範囲を調整してください。
#'
#' @note F_SPRのF管理基準値の初期値は　与えられたFのもとでのSPR/目的のSPR　を初期値とするように調整されるので不要。
#'
#' @export
#' @import tibble 

# ref.F
ref.F <- function(
  res, # VPAの結果のオブジェクト
  sel=NULL, # 仮定する選択率．NULLの場合，res$Fc.at.ageが使われる
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

  if(is.null(sel)){
    Fc.at.age <- res$Fc.at.age
    sel <- Fc.at.age/max(Fc.at.age,na.rm=TRUE)
  }
  else{
    Fc.at.age <- sel
  }
  sel <- sel/max(sel,na.rm=T)
    
  na <- sum(!is.na(sel))

  if(is.null(waa.year)) waa.year <- rev(years)[1]
  if(is.null(maa.year)) maa.year <- rev(years)[1]
  if(is.null(M.year)) M.year <- rev(years)[1]
  if(is.null(rps.year)) rps.year <- as.numeric(colnames(res$naa))
  
  if(is.null(waa))  waa <- apply(as.matrix(as.data.frame(res$input$dat$waa)[as.character(waa.year)]),1,mean)
  if(is.null(M))  M <- apply(as.matrix(as.data.frame(res$input$dat$M)[as.character(M.year)]),1,mean)
  if(is.null(maa))  maa <- apply(as.matrix(as.data.frame(res$input$dat$maa)[as.character(maa.year)]),1,mean)

  if(is.null(waa.catch)){
      if(is.null(res$input$dat$waa.catch)){
          waa.catch <- waa
      }
      else{
          waa.catch <- apply(as.matrix(as.data.frame(res$input$dat$waa.catch)[as.character(waa.year)]),1,mean)
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

  original.spr <- calc.rel.abund(sel,1,na,M,waa,waa.catch,maa,min.age=min.age,
                                 max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
  original.spr0 <- calc.rel.abund(sel,0,na,M,waa,waa.catch,maa,min.age=min.age,
                                  max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
  original.perspr <- sum(original.spr$spr)/sum(original.spr0$spr)
    

    # Fcurrent
    Fcurrent <- c(max(Fc.at.age,na.rm=T), mean(Fc.at.age,na.rm=T))
    
    # grid search
    F_current <- Fcurrent[1]
    F.range <- sort(c(F.range,  F_current))
    spr0 <- sum(calc.rel.abund(sel,0,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)$spr)  
    tmp <- lapply(F.range, function(x) calc.rel.abund(sel,x,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef))
    ypr <- sapply(tmp,function(x) sum(x$ypr))
    pspr <- sapply(tmp,function(x) sum(x$spr))/spr0*100
    ypr.spr <- data.frame(F.range=F.range,ypr=ypr,pspr=pspr)
    ypr.spr$Frange2Fcurrent  <- ypr.spr$F.range/F_current    
    
  # F.spr

  spr.f.est <- function(log.p, out=FALSE, sub="med", spr0=NULL){
    Fr <- exp(log.p)

    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    spr <- sum(tmp$spr)
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
      Fspr.init <- ypr.spr$F.range[which.min(abs(ypr.spr$pspr-i))] #original.perspr/i*100
      FpSPR.res <- nlm(spr.f.est, Fspr.init, out=FALSE, sub=i, spr0=spr0, iterlim=iterlim)
      cat("Estimate F%spr: initial value=", Fspr.init," : estimated value=",exp(FpSPR.res$estimate),"\n")
      FpSPR <- c(FpSPR, exp(FpSPR.res$estimate))
    }
    names(FpSPR) <- paste(pSPR,"%SPR",sep="")
  }

  # Fmax

  ypr.f.est <- function(log.p, out=FALSE){
    Fr <- exp(log.p)
  
    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    ypr <- sum(tmp$ypr)    

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
    ypr <- sum(tmp$ypr)
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

  names(Fcurrent) <- names(Fmed) <- names(Fmean) <- names(Flow) <- names(Fhigh) <- names(Fmax) <- names(F0.1) <- c("max","mean")

  Res <- list(sel=sel, min.age=min.age, max.age=max.age, rps.q=rps.q, spr.q=spr.q, Fcurrent=Fcurrent, Fmed=Fmed, Flow=Flow, Fhigh=Fhigh, Fmax=Fmax, F0.1=F0.1, Fmean=Fmean,rps.data=rps.data)
  
  if (!is.null(pSPR)){
    FpSPR <- rbind(FpSPR, sapply(FpSPR, f.mean))
    rownames(FpSPR) <- c("max","mean")
    Res$FpSPR <- FpSPR
  }

  #---- make summary
  Res$summary <- as.data.frame(Res[substr(names(Res),1,1)=="F"])
  Res$summary <- rbind(Res$summary,Res$summary[1,]/Res$summary[1,1])
  dimnames(Res$summary)[[1]][3] <- "Fref/Fcur"

  Res$currentSPR <- list(SPR=sum(original.spr$spr),
                         perSPR=original.perspr,
                       YPR=sum(original.spr$spr))
 
  Res$ypr.spr  <- ypr.spr #data.frame(F.range=F.range,ypr=ypr,spr=spr)
  Res$waa <- waa
  Res$waa.catch <- waa.catch  
  Res$maa <- maa
  #------------------------------

  Res$arglist <- arglist
  Res$spr0 <- spr0

  class(Res) <- "ref"
    
  if(isTRUE(plot)){
      plot_Fref(Res)
  }    
  return(Res)
}

#' ref.Fの出力をプロットするための関数
#'
#' @param rres ref.Fの出力結果
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


#' 将来予測のための関数
#'
#' @param res0 VPAの出力結果
#' @param currentF ABC算定年-1年目までに用いるF at age。NULLの場合、VPA結果のres0$Fc.at.ageを用いる
#' @param multi ABC算定年以降にcurrentFに乗じる係数
#' @param futureF ABC算定年以降に用いられるF at age。NULLの場合にはcurrentFをそのまま用いる。currentF, multi, futureFをあわせて、将来のF at ageはABC算定年-1年まではcurrentF, ABC算定年以降はmulti x futureF (NULLの場合、currentF)が用いられる
#' @param use.MSE 簡易MSEを実施するかどうか
#' @param MSE.option 簡易MSEのoption
#' 
#'
#' @export

##----------------------- 将来予測関数 ----------------------------
## multiのオプションは管理後のFのmultiplier（管理前後でselectivityが同じ）
future.vpa <-
    function(res0,
             currentF=NULL, # 管理前のF
             multi=1, # 管理後（ABC.yearから）のF (current F x multi)
             futureF=NULL,              
             nyear=10,Pope=res0$input$Pope,
             outtype="FULL",
             multi.year=1,#ある特定の年だけFを変えたい場合。デフォルトは1。変える場合は、指定した年またはタイムステップの要素数のベクトルで指定。
             # 年数の指定
             start.year=NULL, # 将来予測の開始年，NULLの場合はVPA計算の最終年の次の年
             ABC.year=NULL, # ABC yearを計算する年。NULLの場合はVPA計算の最終年の次の次の年
             waa.year=NULL, # VPA結果から生物パラメータをもってきて平均する期間
             waa.catch.year=NULL, # VPA結果から生物パラメータをもってきて平均する期間             
             # NULLの場合，VPAの最終年のパラメータを持ってくる
             maa.year=NULL, # VPA結果から生物パラメータをもってきて平均する期間
             M.year=NULL, # VPA結果から生物パラメータをもってきて平均する期間
             seed=NULL,
             strategy="F", # F: 漁獲係数一定, E: 漁獲割合一定、C: 漁獲量一定（pre.catchで漁獲量を指定）
             HCR=NULL,# HCRを使う場合、list(Blim=154500, Bban=49400,beta=1,year.lag=0)のように指定するか、以下の引数をセットする,year.lag=0で将来予測年の予測SSBを使う。-2の場合は２年遅れのSSBを使う
             use.MSE=FALSE,MSE.options=NULL,
             beta=NULL,delta=NULL,Blim=0,Bban=0,
             plus.group=res0$input$plus.group,
             N=1000,# 確率的なシミュレーションをする場合の繰り返し回数。
             # N+1の結果が返され、1列目に決定論的な結果が                       
             # 0を与えると決定論的な結果のみを出力
             silent=FALSE, is.plot=TRUE, # 計算条件を出力、プロットするか
             random.select=NULL, # 選択率をランダムリサンプリングする場合、ランダムリサンプリングする年を入れる
             # strategy="C"または"E"のときのみ有効
             pre.catch=NULL, # list(year=2012,wcatch=13000), 漁獲重量をgivenで与える場合
             # list(year=2012:2017,E=rep(0.5,6)), 漁獲割合をgivenで与える場合                       
             ##-------- 加入に関する設定 -----------------
             rec.new=NULL, # 指定した年の加入量
             # 年を指定しないで与える場合は、自動的にスタート年の加入になる。
             # list(year=, rec=)で与える場合は、対応する年の加入を置き換える。
             ##--- 加入関数
             recfunc=HS.recAR, # 再生産関係の関数
             rec.arg=list(a=1,b=1,rho=0,sd=0,c=1,bias.correction=TRUE,
                          resample=FALSE,resid=0,resid.year=NULL), # 加入の各種設定
             ##--- Frecオプション；Frec計算のための設定リストを与えると、指定された設定でのFrecに対応するFで将来予測を行う
             Frec=NULL,
             # list(stochastic=TRUE, # TRUEの場合、stochastic simulationで50%の確率でBlimitを越す(PMS, TMI)
             # FALSEの場合、RPS固定のprojectionがBilmitと一致する(NSK)
             #      future.year=2018, # 何年の資源量を見るか？
             #      Blimit=450*1000,  # Blimit (xトン)
             #      scenario="catch.mean" or "blimit" (デフォルトはblimit; "catch.mean"とするとstochastic simulationにおける平均漁獲量がBlimitで指定した値と一致するようになる)
             #      Frange=c(0.01,2*mult)) # Fの探索範囲
             waa=NULL,waa.catch=NULL,maa=NULL,M=NULL, # 季節毎の生物パラメータ、または、生物パラメータを外から与える場合
             waa.fun=FALSE, #waaをnaaのfunctionとするか(waa.catchが別に与えられているときには機能しない)
             naa0=NULL,eaa0=NULL,ssb0=NULL,faa0=NULL,
             add.year=0, # 岡村オプションに対応。=1で1年分余計に計算する
             det.run=TRUE # 1回めのランは決定論的将来予測をする（完璧には対応していない）
             ){

        
        argname <- ls()
        arglist <- lapply(argname,function(x) eval(parse(text=x)))
        names(arglist) <- argname
        
        if(is.null(res0$input$unit.waa)) res0$input$unit.waa <- 1
        if(is.null(res0$input$unit.caa)) res0$input$unit.caa <- 1
        if(is.null(res0$input$unit.biom)) res0$input$unit.biom <- 1  
        if(is.null(plus.group)) plus.group <- TRUE
        if(is.null(Pope)) Pope <- FALSE
        
        ##--------------------------------------------------
        if(isTRUE(det.run)) N <- N + 1
        years <- as.numeric(dimnames(res0$naa)[[2]])
        
        ##------------- set default options
        if(is.null(currentF)) currentF <- res0$Fc.at.age
        if(is.null(waa.year)) waa.year <- rev(years)[1]
        if(is.null(waa.catch.year)) waa.catch.year <- rev(years)[1]        
        if(is.null(maa.year)) maa.year <- rev(years)[1]
        if(is.null(M.year)) M.year <- rev(years)[1]
        if(is.null(start.year)) start.year <- rev(years)[1]+1
        if(is.null(ABC.year)) ABC.year <- rev(years)[1]+1
        arglist$ABC.year <- ABC.year

        ##------------- set SR options
        rec.arg.org <- rec.arg
        rec.arg <- set_SR_options(rec.arg,N=N,silent=silent,eaa0=eaa0,det.run=det.run)

        ##------------- set HCR options
        
        if(!is.null(HCR) && is.null(HCR$year.lag)) HCR$year.lag <- 0
        if(!is.null(beta)){
            HCR$beta <- beta
            HCR$Blim <- Blim
            HCR$Bban <- Bban
        }

        ##------------- set options for MSE
        if(isTRUE(use.MSE)){
            if(is.null(MSE.options$N)) MSE.options$N <- N            
            if(is.null(MSE.options$recfunc)) MSE.options$recfunc <- recfunc
            if(is.null(MSE.options$rec.arg)){
                MSE.options$rec.arg <- set_SR_options(rec.arg.org,
                                                      N=MSE.options$N,silent=silent,eaa0=eaa0)
            }
            else{
                MSE.options$rec.arg <- set_SR_options(MSE.options$rec.arg,
                                                      N=MSE.options$N,silent=silent,eaa0=eaa0)
            }
            if(is.null(MSE.options$max.ER)) MSE.options$max.ER <- 0.7            
        }
        ##-------------        
        
        #  fyears <- seq(from=start.year,to=start.year+nyear-1,by=1/ts)
        fyears <- seq(from=start.year,to=start.year+nyear+add.year,by=1)
        
        fyear.year <- floor(fyears)
        ntime <- length(fyears)
        ages <- as.numeric(dimnames(res0$naa)[[1]]) # ages:VPAで考慮される最大年齢数
        min.age <- min(as.numeric(ages))

        year.overlap <- years %in% start.year   
                 {if(sum(year.overlap)==0){
                      nage <- sum(!is.na(res0$naa[,ncol(res0$naa)])) # nage:将来予測で考慮すべき年の数
                  }
                  else{
                      nage <- sum(!is.na(res0$naa[,year.overlap])) 
                  }}
        
        if(!silent){
            arglist.tmp <-  arglist
            arglist.tmp$res0 <- NULL
            arglist.tmp$Bban <- arglist.tmp$Bblim <- arglist.tmp$beta <- arglist.tmp$ssb0 <- arglist.tmp$strategy <- NULL
            print(arglist.tmp)
        }
        
        # シードの設定
        if(is.null(seed)) arglist$seed <- as.numeric(Sys.time())
        
        #------------Frecオプションの場合 -------------
        if(!is.null(Frec)){
            multi.org <- multi
            if(is.null(Frec$stochastic)) Frec$stochastice <- TRUE
            if(is.null(Frec$target.probs)) Frec$target.probs <- 50
            if(is.null(Frec$scenario)) Frec$scenario <- "blimit" # 2017/12/25追記 
            if(is.null(Frec$Frange)) Frec$Frange <- c(0.01,multi.org*2)   # 2017/12/25追記(探索するFの範囲の指定)
            if(is.null(Frec$future.year)) Frec$future.year <- fyears[length(fyears)]-1
            #      arglist$Frec <- Frec
            
            getFrec <- function(x,arglist){
                set.seed(arglist$seed)
                arglist.tmp <- arglist
                arglist.tmp$multi <- x
                arglist.tmp$silent <- TRUE      
                arglist.tmp$Frec <- NULL
                arglist.tmp$is.plot <- FALSE
                if(Frec$stochastic==FALSE){
                    arglist.tmp$N <- 0
                }      
                fres.tmp <- do.call(future.vpa,arglist.tmp)
                tmp <- rownames(fres.tmp$vssb)==Frec$future.year
                if(all(tmp==FALSE)) stop("nyear should be longer than Frec$future.year.")
                if(Frec$stochastic==TRUE){
                    if(Frec$scenario=="blimit"){          
                        is.lower.ssb <- fres.tmp$vssb<Frec$Blimit
                        probs <- (sum(is.lower.ssb[tmp,-1],na.rm=T)-1)/
                            (length(is.lower.ssb[tmp,-1])-1)*100
                        return.obj <- probs-Frec$target.probs
                    }
                    # stochastic projectionにおける平均漁獲量を目的の値に一致させる 
                    if(Frec$scenario=="catch.mean"){
                        return.obj <- (log(Frec$Blimit)-log(mean(fres.tmp$vwcaa[tmp,-1])))^2
                    }
                    # stochastic projectionにおける平均親魚資源量を目的の値に一致させる 
                    if(Frec$scenario=="ssb.mean"){
                        return.obj <- (log(Frec$Blimit)-log(mean(fres.tmp$vssb[tmp,-1])))^2
                    }                
                }
                else{
                    return.obj <- Frec$Blimit-fres.tmp$vssb[tmp,1]
                }
                #        return(ifelse(Frec$method=="nibun",return.obj,return.obj^2))
                return(return.obj^2)                
            }
            
            res <- optimize(getFrec,interval=Frec$Frange,arglist=arglist)        
            multi <- res$minimum
            cat("F multiplier=",multi,"\n")
        }
        
        #-------------- main function ---------------------
        waa.org <- waa
        waa.catch.org <- waa.catch
        maa.org <- maa
        M.org <- M
        
        if(strategy=="C"|strategy=="E") multi.catch <- multi else multi.catch <- 1
        
        faa <- naa <- waa <- waa.catch <- maa <- M <- caa <- 
            array(NA,dim=c(length(ages),ntime,N),dimnames=list(age=ages,year=fyears,nsim=1:N))
        
        allyears <- sort(unique(c(fyears,years)))

        # 全部のデータを記録したフレーム  
        naa_all <- waa_all <- waa_catch_all <- maa_all <- faa_all <- 
            array(NA,dim=c(length(ages),length(allyears),N),dimnames=list(age=ages,year=allyears,nsim=1:N))
        naa_all[,1:length(years),] <- unlist(res0$naa)
        faa_all[,1:length(years),] <- unlist(res0$faa)        
        waa_all[,1:length(years),] <- unlist(res0$input$dat$waa)
        if(is.null(res0$input$dat$waa.catch)){
            waa_catch_all[,1:length(years),] <- unlist(res0$input$dat$waa)
        }else{
            waa_catch_all[,1:length(years),] <- unlist(res0$input$dat$waa.catch)
        }
        maa_all[,1:length(years),] <- unlist(res0$input$dat$maa)      
        i_all <- which(allyears%in%start.year)
        
        alpha <- thisyear.ssb <- array(1,dim=c(ntime,N),dimnames=list(year=fyears,nsim=1:N))
        
        # future biological patameter
        if(is.null(  M.org))   M.org <- apply(as.matrix(res0$input$dat$M[,years %in% M.year]),1,mean)
        if(is.null(waa.org)) waa.org <- apply(as.matrix(res0$input$dat$waa[,years %in% waa.year]),1,mean)
        if(is.null(maa.org)) maa.org <- apply(as.matrix(res0$input$dat$maa[,years %in% maa.year]),1,mean)
        if(is.null(waa.catch.org)){
            if(!is.null(res0$input$dat$waa.catch)) waa.catch.org <- apply(as.matrix(res0$input$dat$waa.catch[,years %in% waa.catch.year]),1,mean)
            else waa.catch.org <- waa.org
        }
        
        M[] <- M.org
        waa[] <- waa.org
        waa_all[,(length(years)+1):dim(waa_all)[[2]],] <- waa.org
        maa[] <- maa.org
        maa_all[,(length(years)+1):dim(maa_all)[[2]],] <- maa.org
        waa.catch[] <- waa.catch.org
        waa_catch_all[,(length(years)+1):dim(maa_all)[[2]],] <- waa.catch.org        
        
        if(is.null(futureF)) futureF <- currentF         
        # future F matrix
        faa[] <- futureF*multi # *exp(rnorm(length(faa),0,F.sigma))
        faa_all[is.na(faa_all)] <- futureF*multi
        # ABCyear以前はcurrent Fを使う。
        faa[,fyears<min(ABC.year),] <- currentF
        faa_all[,allyears%in%fyears[fyears<min(ABC.year)],] <- currentF
        
        ## VPA期間と将来予測期間が被っている場合、VPA期間のFはVPAの結果を使う
        overlapped.years <- list(future=which(fyear.year %in% years),vpa=which(years %in% fyear.year))
        if(length(overlapped.years$future)>0){  
            #          for(jj in 1:length(vpayears.overlapped)){
            for(j in 1:length(overlapped.years$future)){
                if(any(res0$faa[,overlapped.years$vpa[j]]>0) && !is.null(res0$input$dat$waa[,overlapped.years$vpa[j]])){ # もしfaaがゼロでないなら（PMIの場合、2012までデータが入っているが、faaはゼロになっているので
                    faa[,overlapped.years$future[j],] <- res0$faa[,overlapped.years$vpa[j]]
                    waa[,overlapped.years$future[j],] <- res0$input$dat$waa[,overlapped.years$vpa[j]]
                    if(!is.null(res0$input$dat$waa.catch)){
                        waa.catch[,overlapped.years$future[j],] <- res0$input$dat$waa.catch[,overlapped.years$vpa[j]]
                    }
                    else{
                        waa.catch[,overlapped.years$future[j],] <- res0$input$dat$waa[,overlapped.years$vpa[j]]
                    }
                }
            }}
        #}
        
        tmp <- aperm(faa,c(2,1,3))
        tmp <- tmp*multi.year
        faa <- aperm(tmp,c(2,1,3))
        
        #  vpa.multi <- ifelse(is.null(vpa.mode),1,vpa.mode$multi)
        # rps assumption
        rps.mat <- array(NA,dim=c(ntime,N),dimnames=list(fyears,1:N))
        eaa <- ABC.mat <- matrix(0,ntime,N)
        rec.tmp <- list(rec.resample=NULL,tmparg=NULL)
        
        if (waa.fun){ #年齢別体重の予測関数
            WAA <- res0$input$dat$waa
            NAA <- res0$naa
            #      nage <- nrow(WAA)
            WAA.res <- lapply(1:nage, function(i) {
                log.w <- as.numeric(log(WAA[i,]))
                log.n <- as.numeric(log(NAA[i,]))
                lm(log.w~log.n)
            })
            WAA.cv <- sapply(1:nage, function(i) sqrt(mean(WAA.res[[i]]$residuals^2)))
            WAA.b0 <- sapply(1:nage, function(i) as.numeric(WAA.res[[i]]$coef[1]))
            WAA.b1 <- sapply(1:nage, function(i) as.numeric(WAA.res[[i]]$coef[2]))
            ##      waa.rand <- array(0,dim=c(al,nyear+1-min.age,N))
            set.seed(0)      
            cv.vec <- rep(WAA.cv,N*ntime)
            waa.rand <- array(rnorm(length(cv.vec),-0.5*cv.vec^2,cv.vec),dim=c(nage,ntime,N))
            waa.rand[,,1] <- 0
        }
        
        set.seed(arglist$seed)        

        # 将来予測の最初の年の設定；バリエーションがありややこしいのでここで設定される
        if(!start.year%in%years){
            # VPA結果が2011年までで、将来予測の開始年が2012年の場合      
            if(start.year==(max(years)+1)){
            {if(is.null(res0$input$dat$M)){
                 M.lastyear <- M.org
             }
             else{
                 M.lastyear <- res0$input$dat$M[,length(years)]
             }}
            # 1年分forwardさせた年齢構成を初期値とする
            tmp <- forward.calc.simple(res0$faa[1:nage,length(years)],
                                       res0$naa[1:nage,length(years)],
                                       M.lastyear[1:nage],
                                       plus.group=plus.group)
            naa[1:nage,1,] <- naa_all[1:nage,i_all,] <- tmp

            
            if(fyears[1]-min.age < start.year){
                thisyear.ssb[1,] <- sum(res0$ssb[,as.character(fyears[1]-min.age)],na.rm=T)
                #                thisyear.ssb <- rep(thisyear.ssb,N)
            }
            else{
                if(waa.fun){
                    waa[2:nage,1,] <- waa_all[2:nage,i_all,] <-
                        t(sapply(2:nage, function(ii) as.numeric(exp(WAA.b0[ii]+WAA.b1[ii]*log(naa[ii,1,])+waa.rand[ii,1,]))))
                }
                thisyear.ssb[1,] <- colSums(naa[,1,]*waa[,1,]*maa[,1,],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom                           }
            
            thisyear.ssb[1,] <- thisyear.ssb[1,]+(1e-10)
            
            if(!is.null(ssb0)) thisyear.ssb[1,] <- colSums(ssb0)
            
            rec.tmp <- recfunc(thisyear.ssb[1,],res0,
                               rec.resample=rec.tmp$rec.resample,
                               rec.arg=rec.arg)
            eaa[1,] <- rec.tmp$rec.resample[1:N]
            rec.arg$resid <- rec.tmp$rec.resample # ARオプションに対応
            
            if(!is.null(rec.tmp$rec.arg)) rec.arg <- rec.tmp$rec.arg
            naa[1,1,] <- naa_all[1,i_all,] <- rec.tmp$rec
            if (waa.fun) {
                waa[1,1,] <- waa_all[1,i_all,] <-
                    as.numeric(exp(WAA.b0[1]+WAA.b1[1]*log(naa[1,1,])+waa.rand[1,1,])) 
            }
            rps.mat[1,] <- naa[1,1,]/thisyear.ssb[1,]          
            }
            else{
                stop("ERROR Set appropriate year to start projection\n")
            }
        }
        else{
            # VPA期間と将来予測期間が被っている場合にはVPAの結果を初期値として入れる
            naa[,1,] <- naa_all[,i_all,] <- res0$naa[,start.year==years]
        }

        # もし引数naa0が与えられている場合にはそれを用いる
        if(!is.null(naa0)){
            naa[,1,] <- naa_all[,i_all,] <- naa0
            if(is.null(faa0)) faa0 <- res0$Fc.at.age
            faa[] <- faa0*multi
        }      
        
        if(!is.null(rec.new)){
            if(!is.list(rec.new)){
                naa[1,1,] <- naa_all[1,i_all,] <- rec.new
            }
            else{ # rec.newがlistの場合
                naa[1,fyears%in%rec.new$year,] <- naa_all[,allyears%in%rec.new$year,] <- rec.new$rec
            }}

        # 2年目以降の将来予測
        for(i in 1:(ntime-1)){
            
            #漁獲量がgivenの場合
            if(!is.null(pre.catch) && fyears[i]%in%pre.catch$year){
                if(!is.null(pre.catch$wcatch)){
                    if(fyears[i]<ABC.year){
                        tmpcatch <- as.numeric(pre.catch$wcatch[pre.catch$year==fyears[i]]) 
                    }
                    else{
                        tmpcatch <- as.numeric(pre.catch$wcatch[pre.catch$year==fyears[i]]) * multi.catch                  
                    }
                }
                if(!is.null(pre.catch$E)){
                    biom <- sum(naa[,i,]*waa[,i,]*res0$input$unit.waa/res0$input$unit.biom)
                    if(fyears[i]<ABC.year){
                        tmpcatch <- as.numeric(pre.catch$E[pre.catch$year==fyears[i]])  * biom
                    }
                    else{
                        tmpcatch <- as.numeric(pre.catch$E[pre.catch$year==fyears[i]]) * biom * multi.catch                  
                    }
                }
                
                saa.tmp <- sweep(faa[,i,],2,apply(faa[,i,],2,max),FUN="/")
                tmp <- lapply(1:dim(naa)[[3]],
                              function(x) caa.est.mat(naa[,i,x],saa.tmp[,x],
                                                      waa.catch[,i,x],M[,i,x],tmpcatch,Pope=Pope))
#                                                      waa[,i,x],M[,i,x],tmpcatch,Pope=Pope))   # ここがwaa.catchだとwaa.fun=TRUE&MSEの場合に不都合が生じる？？ちょっと応急処置                                                   
                faa.new <- sweep(saa.tmp,2,sapply(tmp,function(x) x$x),FUN="*")
                caa[,i,] <- sapply(tmp,function(x) x$caa)
                faa[,i,] <- faa.new
            }
            else{
                faa.new <- NULL
            }
            
            ## HCRを使う場合(当年の資源量から当年のFを変更する)
            if(!is.null(HCR) && fyears[i]>=ABC.year
               && is.null(faa.new)) # <- pre.catchで漁獲量をセットしていない
            {
                if(!isTRUE(use.MSE)){
                    tmp <- i+HCR$year.lag
                    if(tmp>0){
                        ssb.tmp <- colSums(naa[,tmp,]*waa[,tmp,]*maa[,tmp,],na.rm=T)*
                            res0$input$unit.waa/res0$input$unit.biom
                    }
                    else{
                        vpayear <- fyears[i]+HCR$year.lag
                        ssb.tmp <- sum(res0$ssb[as.character(vpayear)])
                    }
                    alpha[i,] <- ifelse(ssb.tmp<HCR$Blim,HCR$beta*(ssb.tmp-HCR$Bban)/(HCR$Blim-HCR$Bban),HCR$beta)
                    alpha[i,] <- ifelse(alpha[i,]<0,0,alpha[i,])
                    faa[,i,] <- sweep(faa[,i,],2,alpha[i,],FUN="*")
                    faa[,i,] <- faa_all[,i,] <- ifelse(faa[,i,]<0,0,faa[,i,])

                }
                else{ # when using MSE option
                    ABC.tmp <- get_ABC_inMSE(naa_all,waa_all,maa_all,faa_all,M[,(i-2):(i),],res0,
                                             start_year=i_all-2,nyear=2,
                                             N=MSE.options$N,
                                             recfunc=MSE.options$recfunc,
                                             rec.arg=MSE.options$rec.arg,
                                             Pope=Pope,HCR=HCR,plus.group=plus.group,lag=min.age)
                    y <- colSums(naa[,i,] * waa[,i,])
                    ABC.tmp <- ifelse(ABC.tmp>y*MSE.options$max.ER,y*MSE.options$max.ER,ABC.tmp)
                    ABC.mat[i,] <- ABC.tmp
                    ####
                    saa.tmp <- sweep(faa[,i,],2,apply(faa[,i,],2,max),FUN="/")
                    est.result <- lapply(1:dim(naa)[[3]],
                                  function(x) caa.est.mat(naa[,i,x],saa.tmp[,x],
                                                          waa.catch[,i,x],M[,i,x],ABC.tmp[x],Pope=Pope))
                    fmulti_to_saa <- sapply(est.result,function(x) x$x)
                    faa.new2 <- sweep(saa.tmp,2,fmulti_to_saa,FUN="*")
                    caa[,i,] <- sapply(est.result,function(x) x$caa)
                    faa[,i,] <- faa_all[,i_all,] <- faa.new2
                    ####                    
                    }
            }
            
            ## 漁獲して１年分前進（加入はまだいれていない）
            tmp <- forward.calc.mat2(faa[,i,],naa[,i,],M[,i,],plus.group=plus.group)
            # 既に値が入っているところ（１年目の加入量）は除いて翌年のNAAを入れる
            naa.tmp <- naa[,i+1,]
            naa.tmp[is.na(naa.tmp)] <- tmp[is.na(naa.tmp)]          
            naa[,i+1, ] <- naa_all[,i_all+1,] <- naa.tmp

            # naaを更新するタイミングですかさずwaaを更新するようにしないといけない
            if(waa.fun){
                # 動的なwaaは対応する年のwaaを書き換えた上で使う？
                waa[2:nage,i+1,] <- waa.catch[2:nage,i+1,] <- 
                    waa_all[2:nage,i_all+1,] <-
                    t(sapply(2:nage, function(ii) as.numeric(exp(WAA.b0[ii]+WAA.b1[ii]*log(naa[ii,i+1,])+waa.rand[ii,i+1,]))))
            }            
            
            ## 当年の加入の計算
            if(fyears[i+1]-min.age < start.year){
                # 参照する親魚資源量がVPA期間である場合、VPA期間のSSBをとってくる
                thisyear.ssb[i+1,] <- sum(res0$ssb[,as.character(fyears[i+1]-min.age)],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom
                #              thisyear.ssb <- rep(thisyear.ssb,N)              
                if(!is.null(ssb0)) thisyear.ssb[i+1,] <- colSums(ssb0)
            }
            else{
                # そうでない場合
                thisyear.ssb[i+1,] <- colSums(naa[,i+1-min.age,]*waa[,i+1-min.age,]*maa[,i+1-min.age,],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom            
            }

            thisyear.ssb[i+1,] <- thisyear.ssb[i+1,]+(1e-10)
            rec.tmp <- recfunc(thisyear.ssb[i+1,],res0,
                               rec.resample=rec.tmp$rec.resample,
                               rec.arg=rec.arg)
            if(is.na(naa[1,i+1,1]))  naa[1,i+1,] <- naa_all[1,i_all+1,] <- rec.tmp$rec          
            #          if(!is.null(rec.tmp$rec.arg)) rec.arg <- rec.tmp$rec.arg
            if (waa.fun) {
                waa[1,i+1,] <- waa.catch[1,i+1,] <- waa_all[1,i_all+1,] <- 
                    as.numeric(exp(WAA.b0[1]+WAA.b1[1]*log(naa[1,i+1,])+waa.rand[1,i+1,])) 
            }                        
            rps.mat[i+1,] <- naa[1,i+1,]/thisyear.ssb[i+1,]
            eaa[i+1,] <- rec.tmp$rec.resample[1:N]
            rec.arg$resid <- rec.tmp$rec.resample # ARオプションに対応
            
            i_all <- i_all+1
        }
        
        if (!is.null(rec.arg$rho)) rec.tmp$rec.resample <- NULL

        if(Pope){
            caa[] <- naa*(1-exp(-faa))*exp(-M/2)
        }
        else{
            caa[] <- naa*(1-exp(-faa-M))*faa/(faa+M)
        }
        
        caa <- caa[,-ntime,,drop=F]
        if(isTRUE(waa.fun)){ ## アドホックな対応！ waa.fun=TRUEかつwaa.catchが与えられているとき動かない。また、pre.catchが与えられていてwaa.fun=TRUEの場合も不具合おこる！
            waa.catch <- waa[,-ntime,,drop=F]            
        }
        else{
            waa.catch <- waa.catch[,-ntime,,drop=F]           
            }
        thisyear.ssb <- thisyear.ssb[-ntime,,drop=F]      
        waa <- waa[,-ntime,,drop=F]
        maa <- maa[,-ntime,,drop=F]                
        naa <- naa[,-ntime,,drop=F]
        faa <- faa[,-ntime,,drop=F]
        alpha <- alpha[-ntime,,drop=F]      
        M <- M[,-ntime,,drop=F]
        fyears <- fyears[-ntime]
        
        biom <- naa*waa*res0$input$unit.waa/res0$input$unit.biom
        ssb <- naa*waa*maa*res0$input$unit.waa/res0$input$unit.biom
        
        wcaa <- caa*waa.catch*res0$input$unit.waa/res0$input$unit.biom
        vwcaa <- apply(wcaa,c(2,3),sum,na.rm=T)
        
        ABC <- apply(as.matrix(vwcaa[fyears%in%ABC.year,,drop=F]),2,sum)

        if(!is.null(rec.arg$resample)) if(rec.arg$resample==TRUE) eaa[] <- NA # resamplingする場合にはeaaにはなにも入れない
        
        fres <- list(faa=faa,naa=naa,biom=biom,baa=biom,ssb=ssb,wcaa=wcaa,caa=caa,M=M,rps=rps.mat,recruit=naa[1,,],
                     maa=maa,vbiom=apply(biom,c(2,3),sum,na.rm=T),
                     eaa=eaa,alpha=alpha,thisyear.ssb=thisyear.ssb,
                     waa=waa,waa.catch=waa.catch,currentF=currentF,
                     futureF=futureF,
                     vssb=apply(ssb,c(2,3),sum,na.rm=T),vwcaa=vwcaa,naa_all=naa_all,
                     years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,rec.arg=rec.arg,
                     waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
                     Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)

        if(is.plot){
            par(mfrow=c(2,2))
            plot.future(fres)
        }
        if(waa.fun) fres$waa.reg <- WAA.res

        
        if(outtype=="Det"){
            fres <- list(faa=faa[,,1],M=M[,,1],recruit=naa[1,,],eaa=eaa,baa=biom,
                         maa=maa[,,1],vbiom=apply(biom,c(2,3),sum,na.rm=T),
                         waa=waa[,,1],waa.catch=waa.catch[,,1],currentF=currentF,
                         vssb=apply(ssb,c(2,3),sum,na.rm=T),vwcaa=vwcaa,alpha=alpha,
                         years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,
                         futureF=futureF,                                                  
                         waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
                         Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
        }

        if(outtype=="short"){
            fres <- list(recruit=naa[1,,],eaa=eaa,alpha=alpha,
                         Fsakugen=-(1-faa[1,,]/currentF[1]),ABC.mat=ABC.mat,
                         vbiom=apply(biom,c(2,3),sum,na.rm=T),
                         currentF=currentF,alpha=alpha,
                         futureF=futureF,  
                         vssb=apply(ssb,c(2,3),sum,na.rm=T),vwcaa=vwcaa,
                         years=fyears,fyear.year=fyear.year,ABC=ABC,
                         waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
                         Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
        }      

        ## if(non.det==TRUE){
        ##     fres <- list(faa=faa[,,-1,drop=F],naa=naa[,,-1,drop=F],biom=biom[,,-1,drop=F],
        ##                  ssb=ssb[,,-1,drop=F],wcaa=wcaa[,,-1,drop=F],caa=caa[,,-1,drop=F],
        ##                  M=M[,,-1,drop=F],rps=rps.mat[,-1,drop=F],
        ##                  maa=maa[,,-1,drop=F],vbiom=apply(biom[,,-1,drop=F],c(2,3),sum,na.rm=T),
        ##                  eaa=eaa[,-1,drop=F],
        ##                  waa=waa[,,-1,drop=F],waa.catch=waa.catch[,,-1,drop=F],currentF=currentF,
        ##                  vssb=apply(ssb[,,-1,drop=F],c(2,3),sum,na.rm=T),vwcaa=vwcaa[,-1,drop=F],
        ##                  years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,rec.arg=rec.arg,
        ##                  waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
        ##                  Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
        ## }
        
        class(fres) <- "future"

        invisible(fres)
    }


get_ABC_inMSE <- function(naa_all,waa_all,maa_all,faa_all,M,res0,start_year,nyear,N,recfunc,rec.arg,Pope,HCR,
                          plus.group=plus.group,lag=0){
    ABC.all <- numeric()
#    N <- dim(naa_all)[[3]]
#    naa_dummy <- naa_all
    #    naa_dummy[] <- NA
    naa_dummy <- array(NA,c(dim(naa_all)[[1]],dim(naa_all)[[2]],N))
    faa_dummy <- array(NA,c(dim(faa_all)[[1]],dim(faa_all)[[2]],N))
    M_dummy <- array(M,c(dim(naa_all)[[1]],dim(naa_all)[[2]],N))

#    faa_all_dummy <- faa
#    waa_dummy <- waa_all    
    rec.tmp <- list(rec.resample=NULL,tmparg=NULL)
    
    for(s in 1:dim(naa_all)[[3]]){
        naa_dummy[,1:start_year,] <- naa_all[,1:start_year,s]
        faa_dummy[] <- faa_all[,,s]
        #        waa_dummy[,,s] <- waa_all[,,s]
        #        waa_dummy[] <- waa_all[]        
        for(j in 1:nyear){
            sj <- start_year+j
            naa_dummy[,sj,] <- forward.calc.mat2(faa_dummy[,sj-1,],naa_dummy[,sj-1,],M_dummy[,sj-1,],
                                                 plus.group=plus.group)
            thisyear.ssb <- colSums(naa_dummy[,sj-lag,] * waa_all[,sj-lag,s] *
                                    maa_all[,sj-lag,s],na.rm=T)
            naa_dummy[1,sj,] <- recfunc(thisyear.ssb,res0,
                                        rec.resample=rec.tmp$rec.resample,
                                        rec.arg=rec.arg)$rec
        }

        lastyear <- start_year+nyear
        ssb.tmp <-  colSums(naa_dummy[,lastyear,]*
                           waa_all[,lastyear,s]*
                           maa_all[,lastyear,s],na.rm=T)*
            res0$input$unit.waa/res0$input$unit.biom    
        alpha <- ifelse(ssb.tmp<HCR$Blim,HCR$beta*(ssb.tmp-HCR$Bban)/(HCR$Blim-HCR$Bban),HCR$beta)
        alpha <- ifelse(alpha<0,0,alpha)        
        #faa_dummy[,lastyear,] <- sweep(faa_all[,lastyear,],2,alpha,FUN="*")
        faa_dummy[,lastyear,] <- sweep(faa_dummy[,lastyear,],2,alpha,FUN="*")

    
        if(Pope){
            ABC <- naa_dummy[,lastyear,]*(1-exp(-faa_dummy[,lastyear,]))*exp(-M_dummy[,lastyear,]/2)*waa_all[,lastyear,s]
        }
        else{
            ABC <- naa_dummy[,lastyear,]*(1-exp(-faa_dummy[,lastyear,]-M_dummy[,lastyear,]))*faa_dummy[,lastyear,]/(faa_dummy[,lastyear,]+M_dummy[,lastyear,])*waa_all[,lastyear,s] 
        }
        
        ABC.all[s] <- mean(colSums(ABC))

        if(0){
            boxplot(t(apply(naa_dummy*waa_all*maa_all,c(2,3),sum)),ylim=c(0,200000),col=2)
            locator(1)
#            if(Pope){
#                ABC <- naa_dummy*(1-exp(-faa_dummy))*exp(-M[,nyear,]/2)*waa_all
#            }
#            else{
#                ABC <- naa_dummy*(1-exp(-faa_dummy-M[,nyear,]))*faa_dummy/(faa_dummy+M[,nyear,])*waa_all 
#            }
#            boxplot(t(apply(ABC,c(2,3),sum)),col=2)            
        }
    }
    
    return(ABC.all)

}


#' 再生産optiond

set_SR_options <- function(rec.arg,N=100, silent=TRUE,eaa0=NULL, det.run=TRUE){
        ##---- set S-R functin option -----
        ## 使う関数によっては必要ないオプションもあるが、使わないオプションを入れてもエラーは出ないので、
        # rec.arg$resampleがNULLかどうかで、パラメトリックな誤差分布かそうでないか（残差リサンプリング）を判別する
        if(is.null(rec.arg$rho)){
            rec.arg$rho <- 0
            if(!silent) cat("rec.arg$rho is assumed to be 0...\n")
        }
        if(is.null(rec.arg$sd2)) rec.arg$sd2 <- sqrt(rec.arg$sd^2/(1-rec.arg$rho^2)) #rho込み平均補正用SD # HS.recAR

        ## resampling optionを使わない場合
        if(is.null(rec.arg$resample)|!isTRUE(rec.arg$resample)){
            if(is.null(rec.arg$bias.correction)) rec.arg$bias.correction <- TRUE # HS.recAR, HS.rec0
            if(is.null(rec.arg$rho)){
                rec.arg$rho <- 0 # HS.recAR, HS.rec0
                rec.arg$resid <- 0
            }
            if(!is.null(rec.arg$rho)){
                if(rec.arg$rho>0){
                    if(is.null(eaa0)){
                        if(is.null(rec.arg$resid.year)) rec.arg$resid <- rep(rev(rec.arg$resid)[1],N)
                        else rec.arg$resid <- rep(mean(rev(rec.arg$resid)[1:rec.arg$resid.year]),N)
                    }
                    else{
                        rec.arg$resid <- eaa0
                    }
                }
                else{
                    rec.arg$resid <- rep(0,N)
                }
            }
        }
        else{
            if(rec.arg$rho>0) stop("You set rho is >0. You cannot use resample=TRUE option when rho>0") # resamplingの場合に自己相関は考慮できないのでrhoは強制的にゼロ
        }
        
    if(!is.null(rec.arg$sd) & isTRUE(det.run)) rec.arg$sd <- c(0,rep(rec.arg$sd,N-1))
    if(!is.null(rec.arg$sd2)  & isTRUE(det.run)) rec.arg$sd2 <- c(0,rep(rec.arg$sd2,N-1))    
    return(rec.arg)
}

# リサンプリングする残差の年数をどんどん増やしていく
resample_backward.rec <- function(ssb,vpares,#deterministic=FALSE,
                   rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1, 
                                resid=0,duration=5,
                                resid.list=list(),
                                SR="HS",# or "BH","RI"
                                bias.correction=TRUE)){

    rec.arg$resid[1] <- rec.arg$resid[1]+1
    bias.factor <- mean(exp(unlist(rec.arg$resid.list)))
    
    if(rec.arg$SR=="HS") rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)
    if(rec.arg$SR=="BH") rec0 <- rec.arg$a*ssb/(1+rec.arg$b*ssb)
    if(rec.arg$SR=="RI") rec0 <- rec.arg$a*ssb*exp(-rec.arg$b*ssb)

    if(rec.arg$resid[1]%%5==1){
        max.sample <- min(ceiling(rec.arg$resid[1]/5),6)
        rec.arg$resid[-1] <- sample(1:max.sample,length(rec.arg$resid)-1,replace=TRUE)
    }

    resampled.resid <- sapply(rec.arg$resid[-1],function(x) sample(rec.arg$resid.list[[x]],1))
    
    if(isTRUE(rec.arg$bias.correction)){
        rec <- c(rec0[1],exp(log(rec0[-1])+resampled.resid)/bias.factor)
    }
    else{
        rec <- c(rec0[1],exp(log(rec0[-1])+resampled.resid))
    }
  return(list(rec=rec,rec.resample=rec.arg$resid))
}


forward.calc.mat2 <- function(fav,nav,Mv,plus.group=TRUE){
  nage <- max(which(!is.na(nav[,1])))#length(fav)
  na.age <- which(is.na(nav[-1,1]))
#  naa <- matrix(NA,nage,dim(nav)[[2]])
  naa <- matrix(NA,dim(nav)[[1]],dim(nav)[[2]])  
#  for(a in 2:(nage-1)){
  naa[-c(1,nage,na.age),] <- nav[-c(nage,nage-1,na.age),]*
      exp(-fav[-c(nage,nage-1,na.age),]-Mv[-c(nage,nage-1,na.age),])
#  }
  naa[nage,] <- nav[nage-1,]*exp(-fav[nage-1,]-Mv[nage-1,]) 
  pg <- nav[nage,]*exp(-fav[nage,]-Mv[nage,])
  if(plus.group) naa[nage,] <- naa[nage,] + pg
  return(naa)
}

caa.est.mat <- function(naa,saa,waa,M,catch.obs,Pope){
  saa <- saa/max(saa)
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


# しきい値を設定し、その境界前後でリサンプリングする残差を変える
resample_2block.rec <- function(ssb,vpares,#deterministic=FALSE,
                   rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1, 
                                resid.lower=0, # しきい値よりも小さいときの残差のセット
                                resid.higher=0, # しきい値よりも大きいときの残差のセット
                                ssb.threshold=0, # しきい値
                                SR="HS", # or "BH", "RI"
                                bias.correction=TRUE)){

    if(rec.arg$SR=="HS") rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)
    if(rec.arg$SR=="BH") rec0 <- rec.arg$a*ssb/(1+rec.arg$b*ssb)
    if(rec.arg$SR=="RI") rec0 <- rec.arg$a*ssb*exp(-rec.arg$b*ssb) 
    
    mean.bias.lower <- 1/mean(exp(rec.arg$resid.lower))
    mean.bias.higher <- 1/mean(exp(rec.arg$resid.higher))
    
    resid.set <- rep(0,length(ssb[-1]))
    tmp <- ssb[-1]<rec.arg$ssb.threshold
    if(isTRUE(rec.arg$bias.correction)){
        
        mean.bias <- rep(1,length(ssb[-1]))
        mean.bias[tmp] <- mean.bias.lower
        mean.bias[!tmp] <- mean.bias.higher
        
        resid.set[tmp] <- sample(rec.arg$resid.lower,sum(tmp),replace=TRUE)
        resid.set[!tmp] <- sample(rec.arg$resid.higher,sum(!tmp),replace=TRUE)
        
        det.bias <- ifelse(ssb[1]<rec.arg$ssb.threshold,
                           mean.bias.lower,mean.bias.higher)
#        cat(det.bias,"\n")
#        cat(mean(mean.bias),"\n")
        rec <- c(rec0[1],#*det.bias,
                 exp(log(rec0[-1])+resid.set)*mean.bias)
#        if(det.bias<1) browser()
    }
    else{
        resid.set[tmp] <- sample(rec.arg$resid.lower,sum(tmp),replace=TRUE)
        resid.set[!tmp] <- sample(rec.arg$resid.higher,sum(!tmp),replace=TRUE)
        rec <- c(rec0[1],exp(log(rec0[-1])+resid.set))
    }
    
#    if(isTRUE(rec.arg$bias.correction)){
#        mean.bias <- mean(exp(c(rec.arg$resid.lower,rec.arg$resid.higher)))
#        rec <- c(rec0[1],exp(log(rec0[-1])+resid.set))/mean.bias
#    }
#    else{
#        rec <- c(rec0[1],exp(log(rec0[-1])+resid.set))
#    }
    return(list(rec=rec,rec.resample=rec.arg$resid)) # 暫定的変更
}

#' HSを仮定したときの加入関数
#'
#' @param ssb 親魚資源量
#' @param vpares VPAの出力結果
#' 
#'
#' @export

HS.recAR <- function(ssb,vpares,#deterministic=FALSE,
                      rec.resample=NULL,
                      rec.arg=list(a=1000,b=1000,#gamma=0.01,
                                   sd=0.1, rho=0,
                                   resid=0)#, bias.correction=TRUE)
                      ){
    ## 再生産関係からの予測値
#    rec0 <- rec.arg$a*(ssb+sqrt(rec.arg$b^2+(rec.arg$gamma^2)/4)-sqrt((ssb-rec.arg$b)^2+(rec.arg$gamma^2)/4))
    rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)     
    rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値

    rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
    new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
    return(list(rec=rec,rec.resample=new.resid))
}

#' BHを仮定したときの加入関数
#'
#' @param ssb 親魚資源量
#' @param vpares VPAの出力結果
#' 
#'
#' @export

# Beverton-Holt
BH.recAR <- function(ssb,vpares,deterministic=FALSE,rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1,bias.correction=TRUE)){
  rec0 <- rec.arg$a*ssb/(1+rec.arg$b*ssb)
  rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値
  rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
  new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
  return(list(rec=rec,rec.resample=new.resid))
}

#' RIを仮定したときの加入関数
#'
#' @param ssb 親魚資源量
#' @param vpares VPAの出力結果
#'
#' 
#'
#' @export

RI.recAR <- function(ssb,vpares,deterministic=FALSE,rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1,bias.correction=TRUE)){                   
    rec0 <- rec.arg$a*ssb*exp(-rec.arg$b*ssb) 
    rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値
    rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
    new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
    return(list(rec=rec,rec.resample=new.resid))
}

## Allee effect (depensation) ありの再生産関係
# Hockey-stick
HS.recAR2 <- function(ssb,vpares,#deterministic=FALSE,
                     rec.resample=NULL,
                     rec.arg=list(a=1000,b=1000,
                                  sd=0.1, rho=0, c=1,
                                  resid=0)#, bias.correction=TRUE)
){
  ## 再生産関係からの予測値
  #    rec0 <- rec.arg$a*(ssb+sqrt(rec.arg$b^2+(rec.arg$gamma^2)/4)-sqrt((ssb-rec.arg$b)^2+(rec.arg$gamma^2)/4))
  a <- rec.arg$a
  b <- rec.arg$b
  c <- rec.arg$c
  rec0 <- ifelse(ssb>b,b*a,a*b*(ssb/b)^c)
  # rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)     
  rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値
  
  rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
  new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
  return(list(rec=rec,rec.resample=new.resid))
}


# Beverton-Holt
BH.recAR2 <- function(ssb,vpares,deterministic=FALSE,rec.resample=NULL,
                     rec.arg=list(a=1000,b=1000,sd=0.1,rho=0,c=1,bias.correction=TRUE)){
  a <- rec.arg$a
  b <- rec.arg$b
  c <- rec.arg$c
  rec0 <- (a/b)/(1+1/(b*ssb)^c)
  # rec0 <- rec.arg$a*ssb/(1+rec.arg$b*ssb)
  rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値
  rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
  new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
  return(list(rec=rec,rec.resample=new.resid))
}

# Ricker 
RI.recAR2 <- function(ssb,vpares,deterministic=FALSE,rec.resample=NULL,
                     rec.arg=list(a=1000,b=1000,sd=0.1,rho=0,c=1,bias.correction=TRUE)){                   
  a <- rec.arg$a
  b <- rec.arg$b
  c <- rec.arg$c
  rec0 <- a/(b*exp(1))*(b*ssb)^c*exp(c*(1-b*ssb))
  # rec0 <- rec.arg$a*ssb*exp(-rec.arg$b*ssb) 
  rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値
  rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
  new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
  return(list(rec=rec,rec.resample=new.resid))
}


#' 複数の将来予測の結果をプロットする（ggplotは使わず）
#'
#' @param fres.list future.vpaからの出力結果をリストで並べたもの
#' 
#' 
#'
#' @export

plot.futures <- function(fres.list,conf=c(0.1,0.5,0.9),target="SSB",legend.text="",xlim.tmp=NULL,y.scale=1){
    if(target=="SSB")  aa <- lapply(fres.list,function(x) apply(x$vssb[,-1],1,quantile,probs=conf))
    if(target=="Biomass") aa <- lapply(fres.list,function(x) apply(x$vbiom[,-1],1,quantile,probs=conf))
    if(target=="Catch") aa <- lapply(fres.list,function(x) apply(x$vwcaa[,-1],1,quantile,probs=conf))
    if(target=="Recruit"){
        if(is.null(x$recruit)) x$recruit <- x$naa
        aa <- lapply(fres.list,function(x) apply(x$recruit[,-1],1,quantile,probs=conf))
    }

    if(is.null(xlim.tmp)) xlim.tmp <- as.numeric(range(unlist(sapply(aa,function(x) colnames(x)))))
    plot(0,max(unlist(aa)),type="n",xlim=xlim.tmp,
         ylim=y.scale*c(0,max(unlist(aa))),xlab="Year",ylab=target)
    lapply(1:length(aa),function(i) matpoints(colnames(aa[[i]]),t(aa[[i]]),col=i,type="l",lty=c(2,1,2)))
    legend("bottomright",col=1:length(aa),legend=legend.text,lty=1)
    invisible(aa)
}

#' 一つの将来予測の結果をプロットする（ggplotは使わず）
#'
#' @param fres0 future.vpaからの出力結果
#' 
#' 
#'
#' @export

plot.future <- function(fres0,ylim.tmp=NULL,xlim.tmp=NULL,vpares=NULL,what=c(TRUE,TRUE,TRUE),conf=0.1,N.line=0,
                        label=c("Biomass","SSB","Catch"),is.legend=TRUE,add=FALSE,col=NULL,...){
    ## 暗黙に、vssbなどのmatrixの1列目は決定論的なランの結果と仮定している 
    if(is.null(col)) col <- 1                        
    matplot2 <- function(x,add=FALSE,...){
        if(add==FALSE) matplot(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)
        if(add==TRUE) matpoints(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)    
    }

    if(is.null(xlim.tmp)) xlim.tmp <- range(as.numeric(rownames(fres0$vssb)))
    
    if(what[1]){
        matplot2(x <- t(apply(fres0$vbiom[,-1],1,quantile,probs=c(conf,0.5,1-conf))),
                 ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[1])),
                 xlim=xlim.tmp,
                 ylab=label[1],main=label[1],add=add,...)
        points(rownames(fres0$vbiom),apply(fres0$vbiom[,-1],1,mean),type="b",pch=1)
        points(rownames(fres0$vbiom),as.numeric(fres0$vbiom[,1]),type="b",pch=3)
        if(!is.null(vpares)){
            points(colnames(vpares$baa),colSums(vpares$baa),type="o",pch=20)
        }
        if(N.line>0) matpoints(rownames(fres0$vbiom),fres0$vbiom[,2:(N.line+1)],col="gray",type="l",lty=1)
    }

  if(what[2]){
    matplot2(x <- t(apply(fres0$vssb[,-1],1,quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[2])),
             xlim=xlim.tmp,           
             ylab=label[2],main=label[2],add=add,...)
    points(rownames(fres0$vssb),apply(fres0$vssb[,-1],1,mean),type="b",pch=1)    
    points(rownames(fres0$vssb),as.numeric(fres0$vssb[,1]),type="b",pch=3)
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
    matplot2(x <- t(apply(fres0$vwcaa[,-1],1,quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[3])),
             xlim=xlim.tmp,           
             ylab=label[3],main=label[3],add=add,...)
    points(rownames(fres0$vwcaa),apply(fres0$vwcaa[,-1],1,mean),type="b",pch=1)        
    points(rownames(fres0$vwcaa),as.numeric(fres0$vwcaa[,1]),type="b",pch=3)
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


ref.F2 <- function(res0,target.year=c(2018,2023),current.year=2011,Blim,
                   interval=c(0,3),...){
  ssb <- apply(res0$ssb,2,sum)
  Frec <- numeric()
  Frec[1] <- ssb[current.year]/Blim

  for(i in 1:length(target.year)){
    tmpfunc <- function(x,res0,Blim,...){
      fres <- future.vpa(res0=res0,multi=x,...)
      cat(x," ")    
      return((fres$vssb[rownames(fres$vssb)==target.year[i]]-Blim)^2)
    }
    Frec[i+1] <- optimize(tmpfunc,interval=interval,res0=res0,Blim=Blim,...)$minimum
  }
  return(Frec)
}

# 2012. 8. 3 -- 管理基準値計算は外に出す
getABC <- function(res.vpa, # VPAの結果
                   res.ref, # 管理基準値計算の結果
                   res.future, # 将来予測計算の結果
                   ref.case="all",
                   multi=NULL,
                   N=NULL,                   
                   SSBcur=1000,
                   Blim=1000,Bban=0,                   
                   target.year=NULL, # NULLの場合，ABC.year+4
                   catch.year=NULL, # 2013:2017など、漁獲量の平均を出したい期間、NULLの場合、ABC.year:ABC.year+4
                   is.plot=TRUE){
  if(all(ref.case=="all")) ref.case <- names(res.ref$summary)
  if(all(is.null(multi))) multi <- rep(1,length(ref.case))
                                       
  nref <- length(ref.case)

  ABC.year <- res.future$input$ABC.year
  if(is.null(target.year)) target.year <- ABC.year+4
  ABC <- wariai <- aveF <- catch5u <- catch5l <- upperSSBlim <- upperSSBcur <- SSBlim <- SSBcur.tmp <- rep(NA,nref)
  names(ABC) <- names(wariai) <- names(aveF) <- paste(ref.case,"x",round(multi,3))
  wcatch <- matrix(NA,5,nref,dimnames=list(((min(ABC.year)):(min(ABC.year)+4)),names(aveF)))

  fres <- list()
  i.tmp <- match(ref.case,names(res.ref$summary))
  
  if(any(is.na(i.tmp)))
    stop(paste("ref.case specification of is wrong!"))

  years <- res.future$year
  currentF <- res.ref$Fcurrent["max"] * res.ref$sel
  N <- ifelse(is.null(N),dim(res.future$naa)[[3]],N)

  for(i in 1:nref){
    tmp <- res.ref$summary[i.tmp[i]][1,1] * res.ref$sel
    tmp <- max(tmp,na.rm=T)/max(currentF,na.rm=T)*multi[i]
    tmpF <- tmp * currentF
    aveF[i] <- mean(tmpF,na.rm=T)
    input.tmp <- res.future$input        
    input.tmp$multi <- tmp
    input.tmp$is.plot <- FALSE
    input.tmp$N <- N

    # Frecで使われたシードはとっておかないといけない=> seedはFrecの引数の外に出すこと！
    input.tmp$Frec <- NULL
    
    fres[[i]] <- do.call(future.vpa, input.tmp)
    ABC[i] <- fres[[i]]$ABC[1]
#    browser()    
    if(res.future$input$ts>1){ # ts>2の場合、漁獲量などの計算は暦年を使う
      input.tmp <- res.future$input
      input.tmp$multi <- tmp      
      input.tmp$ts <- 1
      input.tmp$is.plot <- FALSE      
      input.tmp$ABC.year <- ABC.year <- floor(min(input.tmp$ABC.year))
      input.tmp$waa <- input.tmp$maa <- input.tmp$M <- input.tmp$waa.catch <- NULL
      input.tmp$N <- N
      fres[[i]] <- do.call(future.vpa, input.tmp)
      years <- fres[[i]]$year
    }
    wariai[i] <- sum(fres[[i]]$wcaa[,years==ABC.year,1],na.rm=T)/
            sum(fres[[i]]$biom[,years==ABC.year,1],na.rm=T)
    catch.year <- (ABC.year):(ABC.year+4)
    wcatch[,i] <- apply(fres[[i]]$vwcaa[years %in% (catch.year),-1],1,mean,na.rm=T)
    catch5u[i] <- quantile(fres[[i]]$vwcaa[years==max(catch.year),-1],probs=0.9) # catchは2017年
    catch5l[i] <- quantile(fres[[i]]$vwcaa[years==max(catch.year),-1],probs=0.1) 

    tmp.year <- years %in% target.year
    if(is.null(SSBcur)) SSBcur <- fres[[i]]$vssb[years==(ABC.year),1]    
      
    SSBcur.tmp[i] <- SSBcur
    upperSSBlim[i] <- sum(fres[[i]]$vssb[tmp.year,-1]>Blim)/N*100 # SSBは2018年当初まで
    upperSSBcur[i] <- sum(fres[[i]]$vssb[tmp.year,-1]>SSBcur)/N*100
    SSBlim[i] <- Blim
  }

  if(is.plot){
    par(mfrow=c(1,2),mar=c(4,4,2,1))
    vssb <- apply(res.vpa$ssb,2,sum,na.rm=T)/1000
    x <- sapply(fres,function(x) x$vssb[,1])/1000
    plot(range(c(as.numeric(names(vssb)),years)),
         c(0,max(x)*1.1),type="n",xlab="Year",ylab="SSB (x1000)")
    matpoints(years,x,col=1:nref,type="l",lty=1,
            ylim=c(0,max(x)))
    points(as.numeric(names(vssb)),vssb,type="b")
    abline(h=c(SSBlim/1000,SSBcur/1000),col="gray")
    title("SSB in deterministic runs")
    plot(0,axes=F,xlab="",ylab="")
    legend("topleft",col=1:nref,lty=1,legend=names(ABC))
  }
  average <- apply(wcatch,2,mean)
  res.ref$ABC <- rbind(aveF,wariai,catch5l,catch5u,average,
                         upperSSBcur,SSBcur.tmp,upperSSBlim,SSBlim,ABC)
  rownames(res.ref$ABC)[3] <- paste("catch5l during ",min(catch.year),"-",max(catch.year),sep="")
  rownames(res.ref$ABC)[4] <- paste("catch5u during ",min(catch.year),"-",max(catch.year),sep="")  
  rownames(res.ref$ABC)[5] <- paste("average catch during ",min(catch.year),"-",max(catch.year),sep="")    
  rownames(res.ref$ABC)[6] <- paste("upperSSBcur at",target.year)
  rownames(res.ref$ABC)[8] <- paste("upperSSBlim at",target.year)  
  fres0 <- fres
  write.table(round(res.ref$ABC,2),sep="\t")
  save(fres0,file="fres0.R") # 将来予測の全結果はfres0.Rにてセーブされている

  # Kobe chartの作成
  kobe.array <- array(NA,dim=c(length(fres),nrow(fres[[1]]$vssb),5))
  dimnames(kobe.array) <- list(names(ABC),rownames(fres[[1]]$vssb),
                               c("catch","Biomass","SSB","upperBlimit","upperBban"))
  for(i in 1:length(fres)){
      kobe.array[i,,] <- as.matrix(get.kobematrix(fres[[i]],
                                   Blim=Blim,Bban=Bban,ssb=TRUE))
  }
  return(list(ABC=res.ref$ABC,kobe.matrix=kobe.array))
}  

#----------------------------------------------------------------------
#----------   加入に関する関数。魚種specific        -------------------
#----------------------------------------------------------------------

#-------------- VPA mode 用関数 -------------------
caa.est <- function(naa,saa,waa,M,catch.obs,Pope){
  saa <- saa/max(saa)
  tmpfunc <- function(x,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,out=FALSE,Pope=Pope){
    if(isTRUE(Pope)){
      caa <- naa*(1-exp(-saa*x))*exp(-M/2)
    }
    else{
      caa <- naa*(1-exp(-saa*x-M))*saa*x/(saa*x+M)
    }
    wcaa <- caa*waa
    if(out==FALSE){
      return((sum(wcaa,na.rm=T)-catch.obs)^2)
    }
    else{
      return(caa)
    }
  }
  tmp <- optimize(tmpfunc,c(0,5),catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=FALSE)
  tmp2 <- tmpfunc(x=tmp$minimum,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=TRUE)
  return(list(x=tmp$minimum,caa=tmp2))
}


#---------------- 結果の確かめ用関数 ---------------------
# --------USAGE-------
# tdata <- get.tdata("vpa_results.csv")
# check.res(res.pms,list(fres,fres),tdata,digits=2,type="%")

get.data <- function(tfile){
  tmpdata <- read.csv(tfile,header=F,as.is=F,colClasses="character")
  flags <- which(substr(tmpdata[,1],1,1)=="#")
  tlist <- list()
  for(i in 1:(length(flags)-1)){
      tmp <- tmpdata[(flags[i]+1):(flags[i+1]-1),]
      if(dim(tmp)[[1]]>1){
        dimnames(tmp) <- list(tmp[,1],tmp[1,])
        tmp <- tmp[,!apply(tmp=="",2,all)]
        tlist[[i]] <- sapply((tmp[-1,-1]),as.numeric)
      }
     else{
        tlist[[i]] <- as.numeric(tmp[tmp!=""])
      }
  }
  names(tlist)[1:4] <- c("naa","faa","Biomass","Fc.at.age")
  dimnames(tlist[[3]])[[1]] <- c("SSB","Biomass")
  for(i in 1:tlist[[5]]){
    names(tlist)[(4+(i-1)*4+1):(4+(i*4))] <- c("fnaa","ffaa","fwcaa","ABC")
  }
  return(tlist)
}


#' VPA結果をcsvファイルに出力する
#' @param res  VPAの結果
#' @param srres fit.SRの結果
#' @param msyres est.MSYの結果
#' @param fres_current future.vpaの結果(Fcurrent)
#' @param fres_HCR future.vpaの結果(F with HCR)
#' @param kobeII kobeII.matrixの結果
#'
#' @export

out.vpa <- function(res=NULL,    # VPA result
                    srres=NULL,  # fit.SR result
                    msyres=NULL, # est.MSY result
                    fres_current=NULL,   # future projection result
                    fres_HCR=NULL,   # future projection result                    
                    kobeII=NULL, # kobeII result
                    filename="vpa" # filename without extension
                    ){
  old.par <- par()  
  exit.func <- function(){
    dev.off()
    options(warn=0)      
  }
  on.exit(exit.func())

  csvname <- paste(filename,".csv",sep="")
  pdfname <- paste(filename,".pdf",sep="")
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
    x <- rbind(colSums(res$ssb),colSums(res$baa),colSums(res$wcaa))
    rownames(x) <- c("Spawning biomass","Total biomass","Catch biomass")
    write.table2(x,title.tmp="Total and spawning biomass")
  }

  if(!is.null(srres)){
      sr_summary <- 
          as_tibble(srres$pars) %>% mutate(AIC   =srres$AIC,
                                           method=srres$input$method,
                                           type  =srres$input$SR)      
      write("\n# SR fit resutls",file=csvname,append=T)
      write_csv(res_summary,file=csvname,append=T)
  }  
  
  if(!is.null(msyres)){
    write("\n# MSY Reference points",file=csvname,append=T)
    write_csv(msyres$summary,file=csvname,append=T)
  }

  
  if(!is.null(fres)){
    write("\n# future projection under F current (average)",file=csvname,append=T)  
    write("\n# future F at age",file=csvname,append=T)
    write.table2(apply(res_future_current$faa,c(1,2),mean),title.tmp="Future F at age")
    
    write("\n# future numbers at age",file=csvname,append=T)
    write.table2(apply(res_future_current$naa,c(1,2),mean),title.tmp="Future numbers at age")

    write("\n# future total and spawning biomass",file=csvname,append=T)
    x <- rbind(apply(fres$vssb, 1,mean),
               apply(fres$vbiom,1,mean),
               apply(fres$vwcaa,1,mean))
    rownames(x) <- c("Spawning biomass","Total biomass","Catch biomass")
    write.table2(x,title.tmp="Future total, spawning and catch biomass")    
  }

  if(!is.null(fres_HCR)){
    write("\n# future projection under F current (average)",file=csvname,append=T)  
    write("\n# future F at age",file=csvname,append=T)
    write.table2(apply(res_future_current$faa,c(1,2),mean),title.tmp="Future F at age")
    
    write("\n# future numbers at age",file=csvname,append=T)
    write.table2(apply(res_future_current$naa,c(1,2),mean),title.tmp="Future numbers at age")

    write("\n# future total and spawning biomass",file=csvname,append=T)
    x <- rbind(apply(fres_HCR$vssb, 1,mean),
               apply(fres_HCR$vbiom,1,mean),
               apply(fres_HCR$vwcaa,1,mean))
    rownames(x) <- c("Spawning biomass","Total biomass","Catch biomass")
    write.table2(x,title.tmp="Future total, spawning and catch biomass")    
  }  
  
  ## if(!is.null(ABC)){
  ##   write("\n# ABC summary",file=csvname,append=T)
  ##   write.table2(ABC$ABC,title.tmp="Future F at age",is.plot=F)
  ##   write("\n# Kobe matrix",file=csvname,append=T)
  ##   for(i in 1:dim(ABC$kobe.matrix)[[3]]){
  ##       write(paste("\n# ",dimnames(ABC$kobe.matrix)[[3]][i]),
  ##             file=csvname,append=T)        
  ##       write.table2(ABC$kobe.matrix[,,i],
  ##                    title.tmp=dimnames(ABC$kobe.matrix)[[3]][i],is.plot=T)        
  ##   }
  ## }
  
}

#' csvファイルとしてまとめられた資源計算結果を読み込んでRのオブジェクトにする
#' @param tfile 資源計算結果がまとめられたcsvファイルの名前
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
    dres$Fc.at.age <- dres$Fc.at.age[!is.na(dres$Fc.at.age)]
    if(length(dres$Fc.at.age)!=nrow(dres$naa)) stop("Dimension of Fc.at.age and numbers at age is differerent.")
    
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
        diff.pope <- mean(unlist(dres$input$dat$caa/caa.pope))
        
        faa <- dres$faa
        M <- dres$input$dat$M
        caa.bara <- dres$naa*faa/(faa+M)*(1-exp(-faa-M))
        diff.bara <- mean(unlist(dres$input$dat$caa/caa.bara))

        if(abs(1-mean(diff.bara))>abs(1-mean(diff.pope))){
            dres$input$Pope <- TRUE
            cat("Pope is TRUE... OK?\n")
        }
        else{
            dres$input$Pope <- FALSE
            cat("Pope is FALSE... OK?\n")            
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


#type="TorF" # true or false
#type="diff" # excel-RVPA
#type="%" # (excel-RVPA)/excel
check.res <- function(res,fres,tdata,digits=3,type="%"){
    
  check.twomats <- function(mat1,mat2,digits=3,type="%"){
    if(!is.null(colnames(mat1))){
      tmp1 <- mat1[,colnames(mat1)%in%colnames(mat2)]
      tmp2 <- mat2[,colnames(mat2)%in%colnames(mat1)]
    }
    else{
      tmp1 <- mat1
      tmp2 <- mat2
    }
    if(type=="TorF"){
      tmp <- round(tmp1,digits) == round(tmp2,digits)
    }
    if(type=="diff"){
      tmp <- round(tmp1-tmp2,digits)
    }
    if(type=="%"){
      tmp <- round((tmp1-tmp2)/tmp1*100,digits)
    }
    return(tmp)
  }
  
  naa.res <- check.twomats(tdata$naa,res$naa,digits=digits,type=type)
  faa.res <- check.twomats(tdata$faa,res$faa,digits=digits,type=type)
  fcaa.res <- check.twomats(tdata$Fc.at.age,res$Fc.at.age,digits=digits,type=type)
    
  tmp.list <- list(naa=naa.res,faa=faa.res,Fc.at.age=fcaa.res)
  return(tmp.list)       
}


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

forward.calc.simple <- function(fav,nav,Mv,plus.group=TRUE){
    nage <- length(nav)#length(fav)
    naa <- rep(NA,nage)
    naa[c(-1,-nage)] <- nav[c(-nage,-(nage-1))]*exp(-fav[c(-nage,-(nage-1))]-Mv[c(-nage,-(nage-1))])

    naa[nage] <- nav[nage-1]*exp(-fav[nage-1]-Mv[nage-1]) 
    pg <- nav[nage]*exp(-fav[nage]-Mv[nage])
    if(plus.group) naa[nage] <- naa[nage] + pg

    return(naa)
}

forward.calc.mat <- function(fav,nav,Mv,plus.group=TRUE){
  nage <- max(which(!is.na(nav[,1])))#length(fav)
  na.age <- which(is.na(nav[,1]))
#  naa <- matrix(NA,nage,dim(nav)[[2]])
  naa <- matrix(NA,dim(nav)[[1]],dim(nav)[[2]])  
#  for(a in 2:(nage-1)){
    naa[c(-1,-nage,-na.age),] <- nav[c(-nage,-(nage-1),-na.age),]*
        exp(-fav[c(-nage,-(nage-1),-na.age),]-Mv[c(-nage,-(nage-1),-na.age),])
#  }
  naa[nage,] <- nav[nage-1,]*exp(-fav[nage-1,]-Mv[nage-1,]) 
  pg <- nav[nage,]*exp(-fav[nage,]-Mv[nage,])
  if(plus.group) naa[nage,] <- naa[nage,] + pg
  return(naa)
}


pred.RI <- function(SSB,a,b) a*SSB*exp(-b*SSB)
pred.BH <- function(SSB,a,b) a*SSB/(1+b*SSB)
pred.HS <- function(SSB,a,b,gamma) a*(SSB+sqrt(b^2+gamma^2/4)-sqrt((SSB-b)^2+gamma^2/4))
pred.SL <- function(SSB,a) a*SSB

##
get.stat <- function(fout,eyear=0,hsp=NULL,tmp.year=NULL){
    col.target <- ifelse(fout$input$N==0,1,-1) 
    tmp <- as.numeric(fout$vssb[(nrow(fout$vssb)-eyear):nrow(fout$vssb),col.target])
    lhs <- sum(tmp<hsp)/length(tmp)
    if(is.null(tmp.year)) tmp.year <- (nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa)
    
    a <- data.frame("catch.mean"=mean(fout$vwcaa[tmp.year,col.target]),
                    "catch.sd"=sd(fout$vwcaa[tmp.year,col.target]),
                    "catch.geomean"=geomean(fout$vwcaa[tmp.year,col.target]),
                    "catch.median"=median(fout$vwcaa[tmp.year,col.target],na.rm=T),
                    "catch.det"=mean(fout$vwcaa[tmp.year,1],na.rm=T),
                    "catch.L10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.1),
                    "catch.H10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.9),
                    "ssb.mean"=mean(fout$vssb[tmp.year,col.target]),
                    "ssb.sd"=sd(fout$vssb[tmp.year,col.target]),                        
                        "ssb.geomean"=geomean(fout$vssb[tmp.year,col.target]),
                        "ssb.median"=median(fout$vssb[tmp.year,col.target],na.rm=T),
                        "ssb.det"=mean(fout$vssb[tmp.year,1],na.rm=T),
                        "ssb.L10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.1),
                        "ssb.H10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.9),

                        "biom.mean"=mean(fout$vbiom[tmp.year,col.target]),
                        "biom.sd"=sd(fout$vbiom[tmp.year,col.target]),                        
                        "biom.geomean"=geomean(fout$vbiom[tmp.year,col.target]),
                        "biom.median"=median(fout$vbiom[tmp.year,col.target],na.rm=T),
                        "biom.det"=mean(fout$vbiom[tmp.year,1],na.rm=T),
                        "biom.L10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.1),
                        "biom.H10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.9),
                        "lower.HSpoint"=lhs,
                        "Fref2Fcurrent"=fout$multi
                        )
        a$U.mean <- a$catch.mean/a$biom.mean
        a$U.median <- a$catch.median/a$biom.median
        a$U.geomean <- a$catch.geomean/a$biom.geomean
        a$U.det <- a$catch.det/a$biom.det

        a$catch.CV <- a$catch.sd/a$catch.mean
        a$ssb.CV <- a$ssb.sd/a$ssb.mean
        a$biom.CV <- a$biom.sd/a$biom.mean

    #        Faa <- as.data.frame(t(fout$multi * fout$input$res0$Fc.at.age))
        Faa <- as.data.frame(t(fout$multi * fout$currentF))    
        colnames(Faa) <- paste("F",dimnames(fout$naa)[[1]],sep="")
        a <- cbind(a,Faa)
        return(a)
    }

get.stat2 <- function(fout,unit.waa=1,eyear=2,hsp=NULL,tmp.year=NULL){
    col.target <- ifelse(fout$input$N==0,1,-1)     
    if(is.null(tmp.year)) tmp.year <- (nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa)
        nage <- dim(fout$naa)[[1]]
        tb <- fout$naa * fout$waa * unit.waa
        if(is.null(fout$waa.catch)) fout$waa.catch <- fout$waa
        tc <- fout$caa * fout$waa.catch * unit.waa
        ssb <- fout$naa * fout$waa *fout$maa  * unit.waa
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

        # MA; mean, ME; median, GM; geometric mean
        names(tc.mat) <- c(paste("TC-MA-A",1:nage,sep=""),paste("TC-ME-A",1:nage,sep=""),
                           paste("TC-GM-A",1:nage,sep=""),paste("TC-DE-A",1:nage,sep=""),
                           paste("TC-L10-A",1:nage,sep=""),paste("TC-H10-A",1:nage,sep=""))
        names(tb.mat) <- c(paste("TB-MA-A",1:nage,sep=""),paste("TB-ME-A",1:nage,sep=""),
                           paste("TB-GM-A",1:nage,sep=""),paste("TB-DE-A",1:nage,sep=""),
                           paste("TB-L10-A",1:nage,sep=""),paste("TB-H10-A",1:nage,sep=""))
        names(ssb.mat) <- c(paste("SSB-GA-A",1:nage,sep=""),paste("SSB-ME-A",1:nage,sep=""),
                            paste("SSB-GM-A",1:nage,sep=""),paste("SSB-DE-A",1:nage,sep=""),
                            paste("SSB-L10-A",1:nage,sep=""),paste("SSB-H10-A",1:nage,sep=""))        
            
        return(as.data.frame(t(c(tb.mat,tc.mat,ssb.mat))))
    }    


get.stat3 <- function(fout,eyear=0,hsp=NULL,tmp.year=NULL,unit.waa=1){
    col.target <- ifelse(fout$input$N==0,1,-1)
    tmp <- as.numeric(fout$vssb[(nrow(fout$vssb)-eyear):nrow(fout$vssb),col.target])
    lhs <- sum(tmp<hsp)/length(tmp)
    if(is.null(tmp.year)) tmp.year <- (nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa)
    
    a <- data.frame("catch.mean"=mean(fout$vwcaa[tmp.year,col.target]),
                    "catch.sd"=sd(fout$vwcaa[tmp.year,col.target]),
                    "catch.geomean"=geomean(fout$vwcaa[tmp.year,col.target]),
                    "catch.median"=median(fout$vwcaa[tmp.year,col.target],na.rm=T),
                    "catch.det"=mean(fout$vwcaa[tmp.year,1],na.rm=T),
                    "catch.L10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.1),
                    "catch.H10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.9),
                    "ssb.mean"=mean(fout$vssb[tmp.year,col.target]),
                    "ssb.sd"=sd(fout$vssb[tmp.year,col.target]),                        
                    "ssb.geomean"=geomean(fout$vssb[tmp.year,col.target]),
                    "ssb.median"=median(fout$vssb[tmp.year,col.target],na.rm=T),
                    "ssb.det"=mean(fout$vssb[tmp.year,1],na.rm=T),
                    "ssb.L10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.1),
                    "ssb.H10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.9),

                    "biom.mean"=mean(fout$vbiom[tmp.year,col.target]),
                    "biom.sd"=sd(fout$vbiom[tmp.year,col.target]),                        
                    "biom.geomean"=geomean(fout$vbiom[tmp.year,col.target]),
                    "biom.median"=median(fout$vbiom[tmp.year,col.target],na.rm=T),
                    "biom.det"=mean(fout$vbiom[tmp.year,1],na.rm=T),
                    "biom.L10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.1),
                    "biom.H10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.9),
                    
                    "rec.mean"=mean(unlist(fout$naa[1,,])[tmp.year,col.target]),
                    "rec.sd"=sd(unlist(fout$naa[1,,])[tmp.year,col.target]),
                    "rec.geomean"=geomean(unlist(fout$naa[1,,])[tmp.year,col.target]),
                    "rec.median"=median(unlist(fout$naa[1,,])[tmp.year,col.target],na.rm=T),
                    "rec.det"=mean(unlist(fout$naa[1,,])[tmp.year,1],na.rm=T),
                    "rec.L10"=quantile(unlist(fout$naa[1,,])[tmp.year,col.target],na.rm=T,probs=0.1),
                    "rec.H10"=quantile(unlist(fout$naa[1,,])[tmp.year,col.target],na.rm=T,probs=0.9),
                    
                    "lower.HSpoint"=lhs,
                    "Fref2Fcurrent"=fout$multi
                    )
    a$U.mean <- a$catch.mean/a$biom.mean
    a$U.median <- a$catch.median/a$biom.median
    a$U.geomean <- a$catch.geomean/a$biom.geomean
    a$U.det <- a$catch.det/a$biom.det

    a$catch.CV <- a$catch.sd/a$catch.mean
    a$ssb.CV <- a$ssb.sd/a$ssb.mean
    a$biom.CV <- a$biom.sd/a$biom.mean
    a$rec.CV <- a$rec.sd/a$rec.mean

    #    Faa <- as.data.frame(t(fout$multi * fout$input$res0$Fc.at.age))
    Faa <- as.data.frame(t(fout$multi * fout$currentF))    
    colnames(Faa) <- paste("F",dimnames(fout$naa)[[1]],sep="")
    res.stat1 <- cbind(a,Faa) # ここまで、get.stat

    agename <- dimnames(fout$naa)[[1]]
    nage <- dim(fout$naa)[[1]]    
    tb <- fout$naa * fout$waa * unit.waa
    if(is.null(fout$waa.catch)) fout$waa.catch <- fout$waa
    tc <- fout$caa * fout$waa.catch * unit.waa
    ssb <- fout$naa * fout$waa *fout$maa  * unit.waa
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
    res.stat <- cbind(res.stat1,res.stat2)
    return(res.stat)    
}    


geomean <- function(x)
{
  ifelse(all(x > 0), exp(mean(log(x))), NA)
}


show.likeprof <- function(res){
    x <- tapply(res$hs$surface$obj,list(res$hs$surface$b,res$hs$surface$a),function(x) x)
    image(as.numeric(rownames(x)),as.numeric(colnames(x)),log(x/min(x)),col=rev(heat.colors(12)),ylab="a",xlab="b")
    contour(as.numeric(rownames(x)),as.numeric(colnames(x)),log(x/min(x)),add=T,nlevels=10,zlim=c(0,0.3))
    points(res$hs$b,res$hs$a)
    title("Diagnostics")    
}


#'
#' VPA計算結果を入れると、毎年のF at ageがどのくらいのSPR, YPRに相当するかを返す
#'
#' @param dres vpa関数の返り値
#' @param target.SPR 目標とするSPR。この値を入れると、結果の$ysdataでその年のFが目標とするSPRを達成するためのFの何倍になっているかを返す。デフォルトは30が入っている
#' @param Fmax 上の計算をするときに探索するFの乗数の最大値
#' @param max.age SPRやYPRの計算をするときに最大何歳まで考慮するか（デフォルトは無限大) 
#'
#' @export 
get.SPR <- function(dres,target.SPR=30,Fmax=10,max.age=Inf){
    # Fの歴史的な%SPRを見てみる                                                                             
    # 毎年異なるFや生物パラメータに対して、YPR,SPR、SPR0がどのくらい変わっているのか見る(Rコード例2)
    # target.SPRが与えられると、target.SPR（％）として与えた数字に対応するSPR値に対するFの乗数も出力する(与えない場合は30%とする)

    dres$ysdata <- matrix(0,ncol(dres$faa),5)
    dimnames(dres$ysdata) <- list(colnames(dres$faa),c("perSPR","YPR","SPR","SPR0","F/Ftarget"))
    for(i in 1:ncol(dres$faa)){
	dres$Fc.at.age <- dres$faa[,i] # Fc.at.ageに対象年のFAAを入れる
        if(all(dres$Fc.at.age>0)){
            byear <- colnames(dres$faa)[i] # 何年の生物パラメータを使うか                                       
            # RVPAのref.F関数でYPRなどを計算。                                                                  
            # 配布している1.3から1.4にアップデートしているので、新しいほうの関数を使うこと(返り値がちょっと違う)
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
#' 
#'
#' @export
#'
#' 

plot_SRdata <- function(SRdata){
    plot(SRdata$SSB,SRdata$R,xlab="SSB",ylab="R",xlim=c(0,max(SRdata$SSB)),ylim=c(0,max(SRdata$R)))
}

#' VPAの結果からMSY推定を行う
#'
#' @param vpares VPAの結果のオブジェクト
#' @param frag   MSY計算時の将来予測で用いる引数のリスト
#' 
#'
#' @export

est.MSY <- function(vpares,
                    farg,
                   seed=farg$seed,
                   nyear=NULL,
                   eyear=0, # 将来予測の最後のeyear+1年分を平衡状態とする
#                   FUN=median, # 漁獲量の何を最大化するか？
                   FUN=mean, # 漁獲量の何を最大化するか？                   
                   N=1000, # stochastic計算するときの繰り返し回数
                   onlylower.pgy=FALSE,# PGY計算するとき下限のみ計算する（計算時間省略のため）
                   optim.method="optimize",
                   max.target="catch.mean", # method="optimize"以外を使うとき、どの指標を最大化するか。他のオプションとしては"catch.median" (漁獲量のmedianの最大化)
                   calc.yieldcurve=TRUE, # yield curveを正確に計算するかどうか。TRUEだと計算時間が余計にかかる。FALSEだと、yield curveは正確ではない
                   Blimit=0, 
                   trace.multi=c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100), # Fmsyを探索したり、Yield curveを書くときにグリッドサーチをするときのFの刻み。Fcurrentに対する乗数。Fが異常に大きい場合、親魚=0になって加入＝NA
                   is.plot=TRUE,
                   PGY=NULL, # PGY管理基準値を計算するかどうか。計算しない場合はNULLを、計算する場合はc(0.8,0.9,0.95)のように割合を入れる
                   B0percent=NULL, # B0_XX%の管理基準値を計算するかどうか
                   Bempirical=NULL, # 特定の親魚量をターゲットにする場合
                   long.term=20, # 世代時間の何倍年後の状態を平衡状態と仮定するか
                   GT=NULL, # 世代時間を外から与える場合(世代時間の計算は将来予測で使われる年齢別成熟率・自然死亡係数を使っているが、別のパラメータを与えたい場合など、外で計算してここに入れる)
                   mY=5, # 自己相関を考慮して管理基準値を計算する場合、平衡状態から何年進めるか

                   estAR.RP=FALSE, # 平衡状態から近年の残差を考慮した将来予測をおこなったときの管理基準値を計算するか
                   resid.year=0,   # ARありの場合、最近年何年分の残差を平均するか
                   current.resid=NULL # 残差の値を直接入れる場合。上の年数が設定されていてもこちらが設定されたらこの値を使う
                   ){

#    require(tidyverse)
  
    farg$seed <- seed

### 内部で使うための関数定義
    ## 最小化のための関数
    ## シミュレーション回数ぶんの漁獲量のFUN（mean, geomean, median）を最大化するFを選ぶ
    msy.objfun <- function(x,f.arg,FUN=FUN,eyear=eyear){
      f.arg$multi <- x
      fout <- do.call(future.vpa,f.arg)
      return(-FUN(fout$vwcaa[(nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa),-1]))
    }

    trace.func <- function(farg,eyear,hsp=0,trace.N=farg$N,
                           fmulti=c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100)){
        trace.res <- NULL
#        ssb.array <- array(0,dim=c(farg$nyear,farg$N+1,length(fmulti)))
        farg$outtype <- "FULL"
        farg$N <- trace.N
        for(i in 1:length(fmulti)){
            farg$multi <- fmulti[i]
            tmp <- do.call(future.vpa,farg)
#            ssb.array[,,i] <- tmp$vssb
            tmp2 <- get.stat3(tmp,eyear=eyear,hsp=hsp)
            trace.res <- rbind(trace.res,tmp2)
            if(tmp2$"ssb.mean"<trace.res$"ssb.mean"[1]/1000){
                fmulti <- fmulti[1:i]
                break()
            }
          }
        trace.res <- as.data.frame(trace.res)
        trace.res$fmulti <- fmulti
        return(list(table=trace.res))
    }

    which.min2 <- function(x){
        max(which(min(x)==x))
    }

    target.func <- function(fout,faa0=NULL,mY=5,N=2,seed=1,eyear=4,p=1,beta=NULL,delta=NULL,Blim=0,Bban=0,sd0=NULL,current.resid=NULL){
        
        farg <- fout$input
        last.year <- dim(fout$naa)[[2]]

        lag <- as.numeric(dimnames(fout$naa)[[1]])[1]        
        # if(lag==0) SSB.m <- NULL else SSB.m <- fout$ssb[,last.year-lag,]
        SSB.m <- fout$ssb[,last.year-lag,]
        ssb0 <- SSB.m
        
        farg$seed <- seed
        farg$N <- N
        farg$nyear <- mY
        farg$naa0 <- p*fout$naa[,last.year,]
        farg$eaa0 <- fout$eaa[last.year,]+current.resid
        farg$ssb0 <- p*ssb0
        farg$faa0 <- faa0
        farg$beta <- beta
        farg$delta <- delta
        farg$Blim <- Blim
        farg$Bban <- Bban
        farg$start.year <- max(as.numeric(colnames(farg$res0$naa)))+1
        farg$ABC.year <- farg$start.year
        if(!is.null(sd0)) farg$rec.arg$sd <- sd0
        farg$Frec <- NULL
        fout <- do.call(future.vpa,farg)
        out <- get.stat3(fout,eyear=0,hsp=Blimit)
#        out <- cbind(out,get.stat2(fout,eyear=0,hsp=Blimit))
        return(list(out,fout))
    }    

### 関数定義おわり
    ## 世代時間を計算
    if(is.null(GT)){
        GT <- Generation.Time(vpares,maa.year=farg$maa.year,
                              M.year=farg$M.year)  # Generation Time
    }
    if(is.null(nyear)){
        nyear <- round(GT*long.term)
    }
    trace.N <- N        
    years <- sort(as.numeric(rev(names(vpares$naa))[1:5]))
    nY <- nyear+1    # これ必要？？

    ## 引数の調整
    b0 <- numeric() # B0
    fout <- fout0 <- trace <- Fhist <- fout.HS.5par <- list()

    farg.org <- farg.tmp <- farg
    farg.tmp$outtype <- "FULL"
    farg.tmp$nyear <- nyear
    farg.tmp$N <- N
    farg.tmp$silent <- TRUE
    farg.tmp$is.plot <- FALSE
    farg.tmp$ABC.year <- max(years)+1
    farg.tmp$add.year <- 1
    farg.tmp$det.run <- FALSE

    if(!is.null(farg.tmp$pre.catch)){
        farg.tmp$pre.catch <- NULL # pre.catchオプションがあるとうまくいかないのでなかったことにする
        cat("notice: option \"pre.catch\" is turned off in estimating MSY.\n")
    }
    if(!is.null(farg.tmp$rec.new)){
        farg.tmp$rec.new <- NULL # rec.newプションがあるとうまくいかないのでなかったことにする
        cat("notice: option \"rec.new\" is turned off in estimating MSY.\n")            
    }

    # B0の計算
    farg.tmp$multi <- 0
    fout0 <- do.call(future.vpa,farg.tmp)
    B0 <- get.stat3(fout0,eyear=eyear,hsp=Blimit)
#    B0 <- cbind(B0,get.stat2(fout0,eyear=eyear,hsp=Blimit))
    rownames(B0) <- "B0"    
    
    trace <- trace.func(farg.tmp,eyear,hsp=Blimit,fmulti=trace.multi,trace.N=trace.N)

    xx <- which.max(trace$table$catch.mean)+c(-1,1)
    range.tmp <- trace$table$fmulti[xx]
    if(xx[1]==0) range.tmp <- c(0,range.tmp)
    if(is.na(range.tmp[2])) range.tmp[2] <- max(trace$table$fmulti)*10

    farg.tmp$multi <- 1
    cat("Estimating MSY\n")
    if(optim.method=="optimize"){
        tmp <- optimize(msy.objfun,range.tmp,f.arg=farg.tmp,eyear=eyear,FUN=FUN)
        # 壁にあたっている限り続ける
        while(sum(round(tmp$minimum,3)==range.tmp)>0){
            tmp0 <- round(tmp$minimum,3)==range.tmp
            range.tmp <- sort(c(range.tmp[tmp0],
                                range.tmp[tmp0] -2*(mean(range.tmp) - range.tmp[tmp0])))
            range.tmp <- ifelse(range.tmp<0,0,range.tmp)
            tmp <- optimize(msy.objfun,range.tmp,f.arg=farg.tmp,eyear=eyear,FUN=FUN)
        }
        farg.msy <- farg.tmp
        farg.msy$multi <- tmp$minimum # Fc.at.a * multiがFmsy
        cat("F multiplier=",tmp$minimum,"\n")
        fout.msy <- do.call(future.vpa,farg.msy)
        fout.msy$input$multi <- fout.msy$multi
        if(calc.yieldcurve){
            trace$table <- rbind(trace$table,trace.func(farg.msy,eyear,hsp=Blimit,trace.N=trace.N,
                                                    fmulti=tmp$minimum+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075))$table)
            trace$table <- trace$table[order(trace$table$fmulti),]
        }
        F.msy <- fout.msy$input$multi*fout.msy$currentF
    }
    # optimizeでなくgridでやる場合
    else{
        Fmulti <- seq(from=min(range.tmp),to=max(range.tmp),by=0.01)
        trace.tmp <- trace.func(farg.tmp,eyear,hsp=Blimit,fmulti=Fmulti,trace.N=trace.N)
        farg.msy <- farg.tmp        
        farg.msy$multi <- trace.tmp$table$fmulti[which.max(unlist(trace.tmp$table[max.target]))]
        cat("F multiplier=",farg.msy$multi,"\n")        
        fout.msy <- do.call(future.vpa,farg.msy)
        trace$table <- rbind(trace$table,trace.tmp$table)
        trace$table <- trace$table[order(trace$table$fmulti),]        
    }

    MSY <- get.stat3(fout.msy,eyear=eyear)
#    MSY <- cbind(MSY,get.stat2(fout.msy,eyear=eyear))
    rownames(MSY) <- "MSY"
#    cat(" SSB=",MSY$"ssb.mean","\n")    
	
    gc(); gc();
	
    ## PGYの計算
    fout.PGY <- list()
    PGYstat <- NULL
    if(!is.null(PGY)){
        s <- 1
        for(j in 1:length(PGY)){
            cat("Estimating PGY ",PGY[j]*100,"%\n")                        
            ttmp <- trace$table$catch.mean-PGY[j]*MSY$catch.mean
            ttmp <- which(diff(sign(ttmp))!=0)
            frange.list <- list(trace$table$fmulti[ttmp[1]+0:1],
                                trace$table$fmulti[ttmp[2]+0:1])
            if(isTRUE(onlylower.pgy)) i.tmp <- 2  else i.tmp <- 1:2
            for(i in i.tmp){
                farg.pgy <- farg.tmp
                if(sum(is.na(frange.list[[i]]))>0) frange.list[[i]] <- c(0,300)
                farg.pgy$Frec <- list(stochastic=TRUE,
                                      future.year=rev(rownames(fout0$vssb))[1],
                                      Blimit=PGY[j]*MSY$catch.mean,
                                      scenario="catch.mean",Frange=frange.list[[i]])
                fout.PGY[[s]] <- do.call(future.vpa,farg.pgy)
                fout.PGY[[s]]$input$multi <- fout.PGY[[s]]$multi
                PGYstat <- rbind(PGYstat,get.stat3(fout.PGY[[s]]))

                if(calc.yieldcurve){
                    trace$table <- rbind(trace$table,trace.func(farg.msy,eyear,hsp=Blimit,trace.N=trace.N,
                                                                fmulti=fout.PGY[[s]]$multi+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075))$table)
                    trace$table <- trace$table[order(trace$table$fmulti),]
                }
                fout.PGY[[s]][names(fout.PGY[[s]])!="input"] <- NULL  
                s <- s+1
		gc(); gc();
            }
        }
#        PGYstat <- as.data.frame(t(sapply(fout.PGY,get.stat3,eyear=eyear,hsp=Blimit)))
#        PGYstat <- cbind(PGYstat,as.data.frame(t(sapply(fout.PGY,get.stat2,eyear=eyear,hsp=Blimit))))
        rownames(PGYstat) <- names(fout.PGY) <- paste("PGY",rep(PGY,each=length(i.tmp)),
                                                      rep(c("upper","lower")[i.tmp],length(PGY)),sep="_")
    }
    else{
        PGYstat <-  NULL
        }
    ###

    ## B0_%の計算
    fout.B0percent <- list()
    B0stat <- NULL
    if(!is.null(B0percent)){
        for(j in 1:length(B0percent)){
            cat("Estimating B0 ",B0percent[j]*100,"%\n")            
            ttmp <- trace$table$ssb.mean-B0percent[j]*B0$ssb.mean
            ttmp <- which(diff(sign(ttmp))!=0)
            frange.list <- trace$table$fmulti[ttmp[1]+0:1]
            farg.b0 <- farg.tmp
            farg.b0$Frec <- list(stochastic=TRUE,
                                 future.year=rev(rownames(fout0$vssb))[1],
                                 Blimit=B0percent[j]*B0$ssb.mean,
                                 scenario="ssb.mean",Frange=frange.list)
            fout.B0percent[[j]] <- do.call(future.vpa,farg.b0)
            fout.B0percent[[j]]$input$multi <- fout.B0percent[[j]]$multi
            B0stat <- rbind(B0stat,get.stat3(fout.B0percent[[j]]))
            if(calc.yieldcurve){
                trace$table <- rbind(trace$table,trace.func(farg.msy,eyear,hsp=Blimit,trace.N=trace.N,
                                                            fmulti=fout.B0percent[[j]]$multi+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075))$table)
                    trace$table <- trace$table[order(trace$table$fmulti),]
            }
            fout.B0percent[[j]][names(fout.B0percent[[j]])!="input"] <- NULL            
            gc(); gc();
        }
        rownames(B0stat) <- names(fout.B0percent) <- paste("B0-",B0percent*100,"%",sep="")
    }
    else{
        B0stat <-  NULL
        }
###

    ## 特定のSSBを目指す場合
    fout.Bempirical <- list()
    Bempirical.stat <- NULL
    if(!is.null(Bempirical)){
        for(j in 1:length(Bempirical)){
            cat("Estimating B empirical ",Bempirical[j],"\n")            
            ttmp <- trace$table$ssb.mean-Bempirical[j]
            ttmp <- which(diff(sign(ttmp))!=0)
            frange.list <- trace$table$fmulti[ttmp[1]+0:1]
            farg.ben <- farg.tmp
            farg.ben$Frec <- list(stochastic=TRUE,
                                 future.year=rev(rownames(fout0$vssb))[1],
                                 Blimit=Bempirical[j],
                                 scenario="ssb.mean",Frange=frange.list)
            fout.Bempirical[[j]] <- do.call(future.vpa,farg.ben)
            fout.Bempirical[[j]]$input$multi <- fout.Bempirical[[j]]$multi
            Bempirical.stat <- rbind(Bempirical.stat,get.stat3(fout.Bempirical[[j]]))

            if(calc.yieldcurve){
                trace$table <- rbind(trace$table,trace.func(farg.msy,eyear,hsp=Blimit,trace.N=trace.N,
                                                            fmulti=fout.Bempirical[[j]]$multi+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075))$table)
                    trace$table <- trace$table[order(trace$table$fmulti),]
            } 
            fout.Bempirical[[j]][names(fout.Bempirical[[j]])!="input"] <- NULL
            gc(); gc();            
        }
        rownames(Bempirical.stat) <- names(fout.Bempirical) <- paste("Ben-",round(Bempirical),"",sep="")
    }
    else{
        Bempirical.stat <-  NULL
        }
###

    refvalue <- bind_rows(MSY,B0,PGYstat,B0stat,Bempirical.stat) %>% as_tibble %>%
        mutate(RP_name=c("MSY","B0",rownames(PGYstat),rownames(B0stat),rownames(Bempirical.stat)),
               AR=FALSE)
    refvalue <- refvalue %>%
                   mutate(SSB2SSB0=refvalue$ssb.mean/refvalue$ssb.mean[2])
    sumvalue <- refvalue %>% select(RP_name,AR,ssb.mean,SSB2SSB0,biom.mean,U.mean,catch.mean,catch.CV,Fref2Fcurrent)
    colnames(sumvalue) <- c("RP_name","AR","SSB","SSB2SSB0","B","U","Catch","Catch.CV","Fref/Fcur")
    sumvalue <- bind_cols(sumvalue,refvalue[,substr(colnames(refvalue),1,1)=="F"])
    

### ARありの場合の管理基準値の計算（平衡状態から5年分進めたときの値）
    if(isTRUE(estAR.RP)){
        if(resid.year > 0 && is.null(current.resid)){
            current.resid <- mean(rev(fout.msy$input$rec.arg$resid)[1:resid.year]) 
            cat("Residuals of ",resid.year," years are averaged as, ",current.resid,"\n")
        }
        else{
            if(resid.year==0){
                current.resid <- 0
            }
        }

        lag <- as.numeric(rownames(fout.msy$naa))[1]            
        eyear <- mY+(lag > 0)*(lag-1)
        
        MSY2 <- target.func(fout.msy,mY=mY,seed=seed,N=N,eyear=mY,current.resid=current.resid)
        B02 <- target.func(fout0,mY=mY,seed=seed,N=N,eyear=mY,current.resid=current.resid)
        if(!is.null(PGY)){
            PGYstat2 <- lapply(1:length(fout.PGY),
                               function(x) target.func(fout.PGY[[x]],mY=mY,seed=seed,N=N,eyear=mY,current.resid=current.resid))
        }
        else{
            PGYstat2 <- NULL
        }

        if(!is.null(B0percent)){
            B0stat2 <- lapply(1:length(fout.B0percent),
                              function(x) target.func(fout.B0percent[[x]],mY=mY,seed=seed,N=N,eyear=mY,current.resid=current.resid)
                              )
        }
        else{
            B0stat2 <- NULL
        }

        if(!is.null(Bempirical)){
            Bempirical.stat2 <- lapply(1:length(fout.Bempirical),
                                       function(x) target.func(fout.Bempirical[[x]],mY=mY,seed=seed,N=N,eyear=mY,current.resid=current.resid)
                                       )
        }
        else{
            Bempirical.stat2 <- NULL
        }    

        refvalue2 <- bind_rows(MSY2[[1]],B02[[1]],
                               purrr::map_dfr(PGYstat2,function(x) x[[1]]),
                               purrr::map_dfr(B0stat2,function(x) x[[1]]),
                               purrr::map_dfr(Bempirical.stat2,function(x) x[[1]])) %>% as_tibble() %>%
            mutate(RP_name=refvalue$RP_name,AR=TRUE)

        refvalue2 <-  refvalue2 %>%
            mutate(SSB2SSB0=refvalue$ssb.mean/refvalue$ssb.mean[2])
        
        sumvalue2 <- refvalue2 %>% select(RP_name,AR,ssb.mean,SSB2SSB0,biom.mean,U.mean,catch.mean,catch.CV,Fref2Fcurrent)
        colnames(sumvalue2) <- c("RP_name","AR","SSB","SSB2SSB0","B","U","Catch","Catch.CV","Fref/Fcur")
        sumvalue2 <- bind_cols(sumvalue2,refvalue2[,substr(colnames(refvalue2),1,1)=="F"])


        ssb.ar.mean <- cbind(apply(MSY2[[2]]$vssb,1,mean),
                             apply(B02[[2]]$vssb,1,mean),
                             sapply(PGYstat2,function(x) apply(x[[2]]$vssb,1,mean)),
                             sapply(B0stat2,function(x) apply(x[[2]]$vssb,1,mean)),
                             sapply(Bempirical.stat2,function(x) apply(x[[2]]$vssb,1,mean)))
        ssb.ar.mean <- sweep(matrix(as.numeric(ssb.ar.mean),nrow(ssb.ar.mean),ncol(ssb.ar.mean)),
                             2,unlist(sumvalue$SSB),FUN="/")
        colnames(ssb.ar.mean) <- rownames(sumvalue$SSB)
    }
    else{ # estAR.RP==FALSEのとき
        ssb.ar.mean <- NULL
        sumvalue2 <- NULL
        refvalue2 <- NULL
    }
    
    ### 結果のプロットなど

    trace$table <- subset(trace$table,fmulti>0)
    
    if(isTRUE(is.plot)){
        # plot of yield curve
        par(mfrow=c(1,3),mar=c(4,4,2,1))
        plot(trace$table$fmulti,trace$table$"ssb.mean"*1.2,type="n",xlab="Fref/Fcurrent",ylab="SSB")
        abline(v=sumvalue$Fref2Fcurrent,col="gray")
        text(sumvalue$Fref2Fcurrent,max(trace$table$"ssb.mean")*seq(from=1.1,to=0.8,length=nrow(sumvalue)),rownames(sumvalue))
        menplot(trace$table$fmulti,cbind(0,trace$table$"ssb.mean"),col="skyblue",line.col="darkblue")
        title("Equiribrium SSB")
        
        plot(trace$table$fmulti,trace$table$"catch.mean",type="n",xlab="Fref/Fcurrent",ylab="Catch")
        abline(v=sumvalue$Fref2Fcurrent,col="gray")        
        menplot(trace$table$fmulti,cbind(0,trace$table$"catch.mean"),col="lightgreen",line.col="darkgreen")
        title("Equiribrium Catch (Yield curve)")        

        # plot of the effect of AR
        if(isTRUE(estAR.RP)){
            matplot(ssb.ar.mean,type="b",ylab="SSB_MSY_AR/SSB_MSY",xlab="Years from Equiribrium")
            legend("topright",col=1:ncol(ssb.ar.mean),legend=rownames(sumvalue),lty=1:ncol(ssb.ar.mean))
            title("plot of the effect of AR")
        }
    }

    ## kobe II matrix
    #kobe2 <- array(0,dim=c(dim(trace$array)[[1]],dim(trace$array)[[2]],length(sumvalue$SSb)))
    #for(i in 1:length(sumvalue$SSB)){
    #tmp <- trace$array > sumvalue$SSB[i]
    #kobe2[,,i] <- cbind(kobe2,apply(tmp,c(1,2),mean))
    #  }
    #dimnames(kobe2) <- list()

    input.list <- list(B0=fout0$input,
                       msy=fout.msy$input,
                       pgy=lapply(fout.PGY,function(x) x$input),
                       B0percent=lapply(fout.B0percent,function(x) x$input))

    allsum <- bind_rows(sumvalue,sumvalue2)
    allsum$RP.definition <- NA
    allsum$RP.definition[allsum$AR==FALSE&allsum$RP_name=="MSY"] <- "Btarget0"
#    allsum$RP.definition[allsum$AR==FALSE&allsum$RP_name=="PGY_0.9_lower"] <- "Blow0"
    allsum$RP.definition[allsum$AR==FALSE&allsum$RP_name=="PGY_0.6_lower"] <- "Blimit0"    
    allsum$RP.definition[allsum$AR==FALSE&allsum$RP_name=="PGY_0.1_lower"] <- "Bban0"
    allsum <- allsum %>% select(1,ncol(allsum),2:(ncol(allsum)-1))    

    Fvector <- select(allsum,num_range("F",0:40))

    output <- list(summary =allsum,#as.data.frame(as.matrix(sumvalue)),
                   summary_tb=allsum,
#                   summary_tb=allsum,
#                   all.stat=as.data.frame(as.matrix(refvalue)),
                   all.stat=bind_rows(refvalue,refvalue2),                   
                   trace   =trace$table,
                   input.list=input.list,
                   Fvector =Fvector,
                   F.msy   =F.msy)
    
    if(isTRUE(estAR.RP)){
        output$summaryAR   <- as.data.frame(as.matrix(sumvalue2))
        output$all.statAR  <- as.data.frame(as.matrix(refvalue2))
        output$ssb.ar.mean <- ssb.ar.mean
        }
    output$SPR.msy <- calc_MSY_spr(output)
    
    invisible(output)    
}


#### function definition
get.perform <- function(fout0,Blimit=0,longyear=50,smallcatch=0.5,N=NULL,
                        shortyear=c(3,5,10),tmp.year=NULL){
    stat1 <- get.stat(fout0,eyear=0,hsp=Blimit,tmp.year=tmp.year)[c("catch.mean","catch.CV","biom.mean","biom.CV","ssb.mean","lower.HSpoint")]
    stat2 <- get.stat2(fout0,eyear=0,tmp.year=tmp.year)
    stat2 <- data.frame(t(as.data.frame(strsplit(colnames(stat2),"-"))),value=as.numeric(stat2))
    rownames(stat2) <- NULL

    # waaによる加重平均年齢&組成
    xx <- subset(stat2,X1=="TB" & X2=="MA")
    nage <- sum(!is.na(xx$value))
    tmp <- c(rep(2,ceiling(nage/3)),rep(3,ceiling(nage/3)))
    tmp <- c(rep(1,nage-length(tmp)),tmp)
    if(sum(tmp==1)==0 & sum(tmp==2)>1) tmp[1] <- 1

    xx$bvalue <- xx$value * fout0$waa[,1,1]
    xx$waa <- fout0$waa[,1,1]
    large.portion1 <- tapply(xx$bvalue[!is.na(xx$bvalue)],tmp,sum,na.rm=T)
    stat1$largefish.nature <- large.portion1[names(large.portion1)==3]/sum(large.portion1)
    aage.biom <- sum(xx$bvalue * 0:(length(xx$bvalue)-1))/sum(xx$bvalue)
    
    xx <- subset(stat2,X1=="TC" & X2=="MA")
    xx$bvalue <- xx$value * fout0$waa[,1,1]    
    aage.catch <- sum(xx$bvalue * 0:(length(xx$bvalue)-1))/sum(xx$bvalue)
    large.portion2 <- tapply(xx$bvalue[!is.na(xx$bvalue)],tmp,sum,na.rm=T)
    stat1$largefish.catch <- large.portion2[names(large.portion2)==3]/sum(large.portion2)    

    # 漁獲量<0.5平均漁獲量の頻度
    if(is.null(tmp.year)) tmp.year <- nrow(fout0$vwcaa)
    stat1$catch.safe <- 1/mean(fout0$vwcaa[tmp.year,]<smallcatch*mean(fout0$vwcaa[tmp.year,]))
    stat1$catch.safe <- ifelse(stat1$catch.safe>longyear,longyear,stat1$catch.safe)
    
    # 親魚量<Blimitの頻度　→　確率の逆数
    stat1$ssb.safe <- 1/stat1$"lower.HSpoint"
    stat1$ssb.safe <- ifelse(stat1$ssb.safe>longyear,longyear,stat1$ssb.safe)

    # ABC.yearから5年目までの平均累積漁獲量
    short.catch <- numeric()
    for(i in 1:length(shortyear)){
        years <- fout0$input$ABC.year:(fout0$input$ABC.year+shortyear[i])
        short.catch[i] <- mean(apply(fout0$vwcaa[rownames(fout0$vwcaa)%in%years,-1],2,sum))
    }
    names(short.catch) <- paste("short.catch",shortyear,sep="")
    short.catch <- as.data.frame(t(short.catch))

    # 平衡状態になった年
    years <- names(fout0$vssb[,1])[-1]
    heikou.diff <- which(diff(fout0$vssb[,1])/fout0$vssb[-1,1]<0.01)
    if(length(heikou.diff)>0) stat1$eq.year <- years[min(heikou.diff)] else stat1$eq.year <- Inf 
    
    dat <- data.frame(stat1,short.catch,aage.biom=aage.biom,aage.catch=aage.catch,effort=fout0$multi,
                      waa=as.data.frame(t(fout0$waa[,1,1])),meigara=as.data.frame(t(tmp)))
    return(dat)
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

## kobe.matrixの計算
# Pr(B<Btarget)のみを返す単純なやつ
get.kobemat <- function(fout,N=fout$input$N,nyear=fout$input$nyear,Btarget=0,
                      fmulti=seq(from=0.3,to=1,by=0.1)){
    multi.org <- 1
    fres.short <- list()
    farg <- fout$input
    farg$Frec <- NULL
    farg$N <- N
    farg$nyear <- nyear
    for(i in 1:length(fmulti)){
        farg$multi <- multi.org * fmulti[i]
        fres.short[[i]] <- do.call(future.vpa,farg)
    }
  	prob.btarget <- sapply(fres.short,function(x) apply(x$vssb>Btarget,1,mean))
	colnames(prob.btarget) <- fmulti
    
    invisible(prob.btarget)
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



## ちょっと複雑なkobe.plot
# foutsが複数の将来予測の結果。brefsは複数の管理基準値
get.kobemat2 <- function(fouts,brefs,xlim=NULL,target.prob=0.5){
#    brefs <- sort(brefs)
    years <- as.numeric(rownames(fouts[[1]]$vssb))        
    probs <- matrix(0,length(years),length(fouts))
    for(j in 1:ncol(brefs)){
        probs <- probs + foreach(i=1:length(fouts),.combine=cbind) %do%
            as.numeric(rowMeans(fouts[[i]]$vssb > brefs[i,j])>target.prob)
        }
    if(is.null(xlim)) xlim <- range(years)
    plot(range(years),
         range(0.5,nrow(brefs)+1.5),
         type="n",xlab="Years",cex=3,xlim=xlim,
         ylab="Strategies",yaxt="n")
    abline(h=1:ncol(brefs),v=years,col="gray")
    axis(side=2,at=1:nrow(brefs),label=rownames(brefs))

#    require()
    cols <- RColorBrewer::brewer.pal(ncol(brefs), "Paired")

    for(i in 1:length(fouts)){
        points(years,rep(i,length(years)),
               col=cols[probs[,i]],pch=20,cex=3)
    }
    legend("topright",pch=20,cex=1,col=cols,ncol=ceiling(ncol(brefs)/2),
           legend=paste("Prob(B>",colnames(brefs),")>",round(target.prob*100),"%"))
}



Generation.Time <- function(vpares,
  maa.year=2014:2015,
  M.year=2014:2015,
  Plus = 19
){

  maa <- vpares$input$dat$maa
  maa <- rowMeans(maa[,colnames(maa) %in% maa.year,drop=F],na.rm=T)
  maa <- maa[!is.na(maa)]    
  M <- vpares$input$dat$M
  M <- rowMeans(M[,colnames(M) %in% M.year,drop=F],na.rm=T)
  M <- M[!is.na(M)]
    
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

#' 期間内のRPSをサンプリングするときの加入関数
#'
#' @export

RPS.simple.rec <- function(ssb,vpares,
                           rec.arg=list(rps.year=NULL, # 点推定値のrpsを計算する期間
                             upper.ssb=Inf, # 親魚資源量の上限（単位はトン？）
                             upper.recruit=Inf,
			     sample.year = NULL, # リサンプリング期間。rps.yearと異なる範囲を使う場合、設定する
                             bias.correction=TRUE, # stochasticのときに平均値と中央値の比率を使うもの)
                             rpsmean=FALSE),# 決定論的予測をRPSの平均でおこなうか、中央値でおこなうか。スケトウでは平均、その他の魚種は中央値。
                           deterministic=FALSE,rec.resample=NULL # ここは外から指定する必要ない
                           ){ 
#  argname <- c("rps.year","upper.ssb","upper.recruit","sample.year","bias.corrected","rpsmean")
#  tmp <- !(names(rec.arg) %in% argname)
#  if(sum(tmp)>0) stop(paste(names(rec.arg)[tmp]," no such arguments in RPS.simple.rec"))
  
  if(is.null(rec.arg$bias.correction)) rec.arg$bias.correction <- TRUE
  if(is.null(rec.arg$rpsmean)) rec.arg$rpsmean <- FALSE
  if(is.null(rec.arg$rps.year)) rec.arg$rps.year <- as.numeric(colnames(vpares$naa))
  if(is.null(rec.arg$sample.year)) rec.arg$sample.year <- rec.arg$rps.year

#  browser()
  names(rec.arg)

  if(is.null(rec.resample)){
    min.age <- min(as.numeric(rownames(vpares$ssb)))
    if(min.age==0) slide.tmp <- TRUE else slide.tmp <- -1:-min.age    
    rps.data <- data.frame(year=years <- as.numeric(names(colSums(vpares$ssb,na.rm=T))),
                           ssb=SSB <- as.numeric(colSums(vpares$ssb,na.rm=T)),
                           recruit=rec <- as.numeric(c(vpares$naa[1,slide.tmp],
                             rep(NA,min.age))))
    rps.data$rps <- rps <- rps.data$recruit/rps.data$ssb
    
    rps.range <- as.numeric(rps[years %in% rec.arg$rps.year])
    rps.med <- median(rps.range) # 点推定のためのrps
    rps.mean <- mean(rps.range) # 点推定のためのrps

    rec.resample <- as.numeric(rps[years %in% rec.arg$sample.year]) # リサンプリングのためのrps
#    sample.mean <- mean(sample.range) # リサンプリングのためのrps
#    sample.median <- median(sample.range) # リサンプリングのためのrps

#    if(rec.arg$bias.correction==TRUE){
##      rec.resample <- sample.range/rps.mean*rps.med
##      rec.resample <- sample.range/sample.mean*rps.med # ここは本当にsample.meanで良いのか？(サンプル期間が同じ場合、sample.mean=rps=meanなので問題ない。期間が異なる場合、rps.meanを使うとsample.range/rps.meanの平均が1にならないため、やはりsample.meanを使うのが適切) => 分母はrps.medでなくsample.meanでは？
#        rec.resample <- sample.range/sample.mean*sample.median 
#    }
#    else{
#      rec.resample <- sample.range
#    }
  }
    
    ssb.tmp <- ifelse(ssb>rec.arg$upper.ssb,rec.arg$upper.ssb,ssb)
    rec_determine <- ifelse(rec.arg$rpsmean,ssb.tmp * mean(rec.resample),ssb.tmp * median(rec.resample))
    rec_random <- sample(rec.resample,length(ssb.tmp)-1,replace=TRUE) * ssb.tmp[-1]
    if(rec.arg$bias.correction==TRUE){
        rec_random <- rec_random * mean(rec.resample)/median(rec.resample)
    }
    rec <- c(rec_determine,rec_random)
    rec2 <- ifelse(rec>rec.arg$upper.recruit,rec.arg$upper.recruit,rec)
#    rec2 <- min(rec,rec.arg$upper.recruit)
    return(list(rec=rec2,rec.resample=rec.resample))
}


menplot <- function(x,y,line.col=1,...){
    polygon(c(x,rev(x)),c(y[,1],rev(y[,2])),...)
    if(dim(y)[[2]]>2) points(x,y[,3],type="l",lwd=2,col=line.col)
}

menplot2 <- function(xy,probs=c(0.1,0.9),new=FALSE,xlab=NULL,ylab=NULL,...){
    xx <- rownames(xy)
    yy <- t(apply(xy,1,quantile,probs=c(0.1,0.9)))
    if(isTRUE(new)) matplot(xx,yy,type="n",xlab=xlab,ylab=ylab)
    menplot(xx,yy,...)
}


#' 期間内の残差をサンプリングするときの加入関数
#'
#' @export

resample.rec <- function(ssb,vpares,#deterministic=FALSE,
                   rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1, 
                                resid=0, # 残差リサンプリングする場合、resample=TRUEにして、residにリサンプリングする残差（対数）を入れる
                                SR="HS",# or "BH","RI"
                                bias.correction=TRUE)){

    if(rec.arg$SR=="HS") rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)
    if(rec.arg$SR=="BH") rec0 <- rec.arg$a*ssb/(1+rec.arg$b*ssb)
    if(rec.arg$SR=="RI") rec0 <- rec.arg$a*ssb*exp(-rec.arg$b*ssb)
    
    if(!isTRUE(rec.arg$resample)){
        if(isTRUE(rec.arg$bias.correction)){
            rec <- rec0*exp(rnorm(length(ssb),-0.5*(rec.arg$sd)^2,rec.arg$sd))
        }
        else{
            rec <- rec0*exp(rnorm(length(ssb),0,rec.arg$sd))
        }
    }
    else{
        if(isTRUE(rec.arg$bias.correction)){
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE))/mean(exp(rec.arg$resid)))
        }
        else{
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE)))
        }
    }
  return(list(rec=rec,rec.resample=rec.arg$resid)) # 暫定的変更
}

