#'
#' @import magrittr          
#' @import dplyr             
#' @import tidyr             
#' @import tibble            
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' 
NULL

#' VPAの結果から再生産関係推定用のデータを作成する
#'
#' @param vpares VPAの結果のオブジェクト
#' @encoding UTF-8
#' @export
get.SRdata <- function(vpares,R.dat = NULL,
                       SSB.dat = NULL,
                       years = as.numeric(colnames(vpares$naa)),
                       return.df = FALSE){
    # R.datとSSB.datだけが与えられた場合、それを使ってシンプルにフィットする
    if (!is.null(R.dat) & !is.null(SSB.dat)) {
        dat <- data.frame(R = R.dat,SSB = SSB.dat,year = 1:length(R.dat))
    } else {
    # データの整形
        n <- ncol(vpares$naa)
        L <- as.numeric(rownames(vpares$naa)[1])

        dat      <- list()
        dat$R    <- as.numeric(vpares$naa[1,])
        dat$SSB  <- as.numeric(colSums(vpares$ssb,na.rm = TRUE))
        dat$year <- as.numeric(colnames(vpares$ssb))
    # 加入年齢分だけずらす
        dat$R    <- dat$R[(L+1):n]
        dat$SSB  <- dat$SSB[1:(n-L)]
        dat$year <- dat$year[(L+1):n]

    # データの抽出
        dat <- as.data.frame(dat)
        dat <- dat[dat$year%in%years,]
    }

    class(dat) <- "SRdata"
    if (return.df == TRUE) return(data.frame(year = dat$year,
                                             SSB  = dat$SSB,
                                             R    = dat$R))
    return(dat[c("year","SSB","R")])
}


#' 再生産関係の推定
#'
#' 3種類の再生産関係の推定を、最小二乗法か最小絶対値法で、さらに加入の残差の自己相関を考慮して行うことができる
#' @param SRdata \code{get.SRdata}で作成した再生産関係データ
#' @param SR 再生産関係 (\code{"HS"}: Hockey-stick, \code{"BH"}: Beverton-Holt, \code{"RI"}: Ricker)
#' @param method 最適化法（\code{"L2"}: 最小二乗法, \code{"L1"}: 最小絶対値法）
#' @param AR 自己相関を推定するか(1), しないか(0)
#' @param out.AR 自己相関係数を一度再生産関係を推定したのちに、外部から推定するか（1), 内部で推定するか(0)
#' @param length 初期値を決める際のgridの長さ
#' @param rep.opt \code{TRUE}で\code{optim}による最適化を収束するまで繰り返す
#' @param p0 \code{optim}で設定する初期値
#' @encoding UTF-8
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1], 
#'                 method = c("L1","L2")[1], AR = 1, 
#'                 out.AR = TRUE, rep.opt = TRUE)
#' resSR$pars
#' }
#' @return 以下の要素からなるリスト
#' \describe{
#' \item{\code{input}}{使用した引数のリスト}
#' \item{\code{pars}}{推定されたパラメータ}
#' \item{\code{opt}}{\code{optim}の結果オブジェクト}
#' \item{\code{resid}}{再生産関係から予測値からの加入量の残差}
#' \item{\code{resid2}}{自己相関のを推定したうえでの加入の残差（自己相関なしの時\code{resid}と等しくなる)}
#' \item{\code{loglik}}{対数尤度}
#' \item{\code{k}}{推定したパラメータ数}
#' \item{\code{AIC}}{AIC (\code{out.AR=TRUE}のときは自己相関推定前の結果)}
#' \item{\code{AICc}}{AICc (\code{out.AR=TRUE}のときは自己相関推定前の結果)}
#' \item{\code{BIC}}{BIC (\code{out.AR=TRUE}のときは自己相関推定前の結果)}
#' \item{\code{AIC.ar}}{\code{out.AR=TRUE}のときに\code{acf}関数で得られた自己相関を推定しない場合(0)とする場合(1)のAICの差}
#' \item{\code{pred}}{予測された再生産関係}
#' }
#'
#' @export

fit.SR <- function(SRdata,
                   SR="HS",
                   method="L2",
                   AR=1,TMB=FALSE,
                   hessian=FALSE,w=rep(1,length(SRdata$year)),
                   length=20,
                   max.ssb.pred=1.3, # 予測値を計算するSSBの最大値（観測された最大値への乗数）
                   p0=NULL,
                   out.AR = TRUE, #自己相関係数rhoを外で推定するか
                   rep.opt = FALSE
){ 
  
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  
  if (AR==0) out.AR <- FALSE
  rec <- SRdata$R
  ssb <- SRdata$SSB
  
  N <- length(rec)
  NN <- sum(w) #likelihoodを計算するサンプル数
  
  #  if (SR=="HS") SRF <- function(x,a,b) a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  
  obj.f <- function(a,b,rho){
    resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
    resid2 <- NULL
    for (i in 1:N) {
      resid2[i] <- ifelse(i==1, resid[i], resid[i]-rho*resid[i-1])
    }
    
    if (method == "L2") {
      rss <- w[1]*resid2[1]^2*(1-rho^2)
      for(i in 2:N) rss <- rss + w[i]*resid2[i]^2
      sd <- sqrt(rss/NN)
      sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))
      obj <- -sum(w*dnorm(resid2,0,sd2,log=TRUE))
    } else {
      rss <- w[1]*abs(resid2[1])*sqrt(1-rho^2)
      for(i in 2:N) rss <- rss + w[i]*abs(resid2[i])
      sd <- sum(abs(w*resid2))/NN
      sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))
      obj <- -sum(w*sapply(1:N, function(i){-log(2*sd2[i])-abs(resid2[i]/sd2[i])}))
    }
    return(obj)
  }
  
  if (is.null(p0)) {
    a.range <- range(rec/ssb)
    b.range <- range(1/ssb)
    if (SR == "HS") b.range <- range(ssb)
    grids <- as.matrix(expand.grid(
      seq(a.range[1],a.range[2],len=length),
      seq(b.range[1],b.range[2],len=length)
    ))
    init <- as.numeric(grids[which.min(sapply(1:nrow(grids),function(i) obj.f(grids[i,1],grids[i,2],0))),])
    init[1] <- log(init[1])
    init[2] <- ifelse (SR == "HS",-log(max(0.000001,(max(ssb)-min(ssb))/max(init[2]-min(ssb),0.000001)-1)),log(init[2]))
    if (AR != 0 && !isTRUE(out.AR)) init[3] <- 0
  } else {
    init = p0
  }
  
  if (SR == "HS") { 
    if (AR == 0 || out.AR) {
      obj.f2 <- function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),0)
    } else {
      obj.f2 <-  function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),1/(1+exp(-x[3])))
    }
  } else {
    if (AR == 0 || out.AR) {
      obj.f2 <- function(x) obj.f(exp(x[1]),exp(x[2]),0)
    } else {
      obj.f2 <-  function(x) obj.f(exp(x[1]),exp(x[2]),1/(1+exp(-x[3])))
    }
  }
  
  opt <- optim(init,obj.f2)
  if (rep.opt) {
    for (i in 1:100) {
      opt2 <- optim(opt$par,obj.f2)
      if (abs(opt$value-opt2$value)<1e-6) break
      opt <- opt2
    }
  }
  opt <- optim(opt$par,obj.f2,method="BFGS",hessian=hessian)
    
  Res <- list()
  Res$input <- arglist
  Res$opt <- opt
  
  a <- exp(opt$par[1])
  b <- ifelse(SR=="HS",min(ssb)+(max(ssb)-min(ssb))/(1+exp(-opt$par[2])),exp(opt$par[2]))
  rho <- ifelse(AR==0,0,ifelse(out.AR,0,1/(1+exp(-opt$par[3]))))
  resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
  resid2 <- NULL
  for (i in 1:N) {
    resid2[i] <- ifelse(i == 1,resid[i], resid[i]-rho*resid[i-1])
  }
  
  if (method=="L2") {
    rss <- w[1]*resid2[1]^2*(1-rho^2)
    for(i in 2:N) rss <- rss + w[i]*resid2[i]^2
    sd <- sqrt(rss/NN)
  } else {
    rss <- w[1]*abs(resid2[1])*sqrt(1-rho^2)
    for(i in 2:N) rss <- rss + w[i]*abs(resid2[i])
    sd <- sum(abs(w*resid2))/NN
    sd <- sqrt(2)*sd
  }
  # sd <- ifelse(method=="L2",sqrt(sum(w*resid2^2)/(NN-rho^2)),sqrt(2)*sum(abs(w*resid2))/(NN-rho^2))
  
  Res$resid <- resid
  Res$resid2 <- resid2
  
  Res$pars <- c(a,b,sd,rho)
  
  if (method!="L2") {
    if (AR!=0) {
      message("L1 & out.AR=FALSE is NOT recommended")
      arres <- ar(resid,aic=FALSE,order.max=1,demean=FALSE,method="mle")
      Res$pars[3] <- ifelse(arres$ar<0,sd,sqrt(arres$var.pred))
      Res$pars[4] <- ifelse(arres$ar<0,0,arres$ar)
    }
  }
  
  if (AR==1 && out.AR) {
    arres <- ar(resid,aic=FALSE,order.max=1,demean=FALSE,method="mle")
    Res$pars[3] <- sqrt(arres$var.pred)
    Res$pars[4] <- as.numeric(arres$ar)
    Res$resid2[2:length(Res$resid2)] <- arres$resid[-1]
    Res$AIC.ar  <- ar(resid,order.max=1,demean=FALSE,method="mle")$aic
  }
  
  Res$loglik <- loglik <- -opt$value
  
  names(Res$pars) <- c("a","b","sd","rho")
  Res$pars <- data.frame(t(Res$pars))
  #  Res$gamma <- gamma
  
  ssb.tmp <- seq(from=0,to=max(ssb)*max.ssb.pred,length=100)
  R.tmp <- sapply(1:length(ssb.tmp), function(i) SRF(ssb.tmp[i],a,b))
  pred.data <- data.frame(SSB=ssb.tmp,R=R.tmp)
  Res$pred <- pred.data
  
  Res$k <- k <- length(opt$par)+1
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(NN-k-1)
  Res$BIC <- -2*loglik+k*log(NN)
  return(Res)
}


### 西嶋加筆
# Allee effect (depensation)ありの再生産関係の推定用関数 (c.est=FALSEとすればfit.SRと同じ)
# 修正が必要
fit.SR2 <- function(SRdata,
                    SR="HS",
                    method="L2",
                    AR=1,
                    hessian=FALSE,
                    w=rep(1,length(SRdata$year)),
                    length=20, #parameter (a,b) の初期値を決めるときにgrid searchする数
                    c.est = TRUE #Allee effectを推定するかどうか(c>1でdepensation (Allee-like), c<1でcompensation)
){
  
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  
  rec <- SRdata$R
  ssb <- SRdata$SSB
  
  N <- length(rec)
  NN <- sum(w) #sample size for likelihood calculation
  
  if (SR=="HS") SRF <- function(x,a,b,c) ifelse(x>b,b*a,a*b*(x/b)^c)
  if (SR=="BH") SRF <- function(x,a,b,c) (a/b)/(1+1/(b*x)^c)
  if (SR=="RI") SRF <- function(x,a,b,c) a/(b*exp(1))*(b*x)^c*exp(c*(1-b*x))
  
  obj.f <- function(a,b,rho,c){
    resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b,c)))
    resid2 <- NULL
    for (i in 1:N) {
      resid2[i] <- ifelse(i==1,resid[i], resid[i]-rho*resid[i-1])
    }
    
    if (method == "L2") {
      sd <- sqrt(sum(w*resid2^2)/(NN-rho^2))
      sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))
      obj <- -sum(w*dnorm(resid2,0,sd2,log=TRUE))
    } else {
      sd <- sum(abs(w*resid2))/(NN-rho^2)
      sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))
      obj <- -sum(w*sapply(1:N, function(i){-log(2*sd2[i])-abs(resid2[i]/sd2[i])}))
    }
    return(obj)
  }
  
  a.range <- range(rec/ssb)
  b.range <- range(1/ssb)
  if (SR == "HS") b.range <- range(ssb)
  grids <- as.matrix(expand.grid(
    seq(a.range[1],a.range[2],len=length),
    seq(b.range[1],b.range[2],len=length)
  ))
  init <- as.numeric(grids[which.min(sapply(1:nrow(grids),function(i) obj.f(grids[i,1],grids[i,2],0,1))),])
  init[1] <- log(init[1])
  init[2] <- ifelse (SR == "HS",-log(max(0.000001,(max(ssb)-min(ssb))/max(init[2]-min(ssb),0.000001)-1)),log(init[2]))
  if (AR != 0 || isTRUE(c.est)) init[3] <- 0
  if (AR != 0 && isTRUE(c.est)) init[4] <- 0
  
  if (SR == "HS") { 
    if (AR == 0) {
      if (c.est) {
        obj.f2 <- function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),0,exp(x[3]))
      } else {
        obj.f2 <- function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),0,1)
      }
    } else {
      if (c.est) {
        obj.f2 <-  function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),1/(1+exp(-x[3])),exp(x[4]))
      } else {
        obj.f2 <-  function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),1/(1+exp(-x[3])),1)
      }
    }
  } else {
    if (AR == 0) {
      if (c.est) {
        obj.f2 <- function(x) obj.f(exp(x[1]),exp(x[2]),0,exp(x[3]))
      } else {
        obj.f2 <- function(x) obj.f(exp(x[1]),exp(x[2]),0,1)
      }
    } else {
      if (c.est) {
        obj.f2 <-  function(x) obj.f(exp(x[1]),exp(x[2]),1/(1+exp(-x[3])),exp(x[4]))
      } else {
        obj.f2 <-  function(x) obj.f(exp(x[1]),exp(x[2]),1/(1+exp(-x[3])),1)
      }
    }
  }
  
  opt <- optim(init,obj.f2)
  opt <- optim(opt$par,obj.f2,method="BFGS",hessian=hessian)
  
  Res <- list()
  Res$input <- arglist
  Res$opt <- opt
  
  a <- exp(opt$par[1])
  b <- ifelse(SR=="HS",min(ssb)+(max(ssb)-min(ssb))/(1+exp(-opt$par[2])),exp(opt$par[2]))
  rho <- ifelse(AR==0,0,1/(1+exp(-opt$par[3])))
  c <- ifelse(c.est, exp(rev(opt$par)[1]),1)
  resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b,c)))
  resid2 <- NULL
  for (i in 1:N) {
    resid2[i] <- ifelse(i == 1,resid[i], resid[i]-rho*resid[i-1])
  }
  sd <- ifelse(method=="L2",sqrt(sum(w*resid2^2)/(NN-rho^2)),sqrt(2)*sum(abs(w*resid2))/(NN-rho^2))
  
  Res$resid <- resid
  Res$resid2 <- resid2
  
  Res$pars <- c(a,b,sd,rho,c)
  
  if (method!="L2") {
    if (AR!=0) {
      arres <- ar(resid,aic=FALSE,order.max=1)
      Res$pars[3] <- sqrt(arres$var.pred)
      Res$pars[4] <- arres$ar
    }
  }
  
  Res$loglik <- loglik <- -opt$value
  
  names(Res$pars) <- c("a","b","sd","rho","c")
  Res$pars <- data.frame(t(Res$pars))
  
  ssb.tmp <- seq(from=0,to=max(ssb)*1.3,length=100)
  R.tmp <- sapply(1:length(ssb.tmp), function(i) SRF(ssb.tmp[i],a,b,c))
  pred.data <- data.frame(SSB=ssb.tmp,R=R.tmp)
  Res$pred <- pred.data
  
  Res$k <- k <- length(opt$par)+1
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(NN-k-1)
  Res$BIC <- -2*loglik+k*log(NN)
  return(Res)
}

#' parametric bootstrap usnig fit.SR
#' 
#' @encoding UTF-8
#' @export
#' 

boot.SR <- function(Res,n=100,seed=1){
  N <- length(Res$input$SRdata$year)
  
#  if (Res$input$SR=="HS") SRF <- function(x,a,b,gamma=Res$gamma) a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (Res$input$SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a) 
  if (Res$input$SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (Res$input$SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  
  sd <- sapply(1:N, function(i) ifelse(i==1,Res$pars$sd/sqrt(1-Res$pars$rho^2),Res$pars$sd))
  
  set.seed(seed)
  lapply(1:n, function(j){
    N <- length(Res$input$SRdata$SSB)
    resids <- rnorm(N,0,sd)
    pred <- obs <- resid0 <- numeric(N)
    ssb <- Res$input$SRdata$SSB
    
    for(i in 1:N){
      pred[i] <- SRF(ssb[i],Res$pars$a,Res$pars$b)
      if (i==1) {
        obs[i] <- pred[i]*exp(resids[i])
      } else {
        obs[i] <- pred[i]*exp(Res$pars$rho*resid0[i-1])*exp(resids[i])
      }
      resid0[i] <- log(obs[i]/pred[i])
    }
    res.b <- Res
    res.b$input$SRdata$R <- obs
    res.b <- do.call(fit.SR, res.b$input)
    return(res.b)
  })
}

#'  profile likelihood
#' @param Res 再生産関係（\code{fit.SR()}）の結果オブジェクト
#' @encoding UTF-8
#' @export
#' 

prof.lik <- function(Res,a=Res$pars$a,b=Res$pars$b,sd=Res$pars$sd,rho=ifelse(Res$input$out.AR,0,Res$pars$rho)) {
  SRdata <- Res$input$SRdata
  rec <- SRdata$R
  ssb <- SRdata$SSB
  N <- length(rec)
  SR <- Res$input$SR
  gamma <- Res$gamma
  method <- Res$input$method
  w <- Res$input$w
  
#  if (SR=="HS") SRF <- function(x,a,b) a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)   
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  
  resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
  resid2 <- NULL
  for (i in 1:N) {
    resid2[i] <- ifelse(i==1,resid[i], resid[i]-rho*resid[i-1])
  }
  
  obj <- NULL
  if (method == "L2") {
    for (i in 1:N) {
      if (i==1) {
        obj <- c(obj,-0.5*log(2*pi)-log(sd^2/(1-rho^2))-resid2[i]^2/(2*sd^2/(1-rho^2)))
      } else {
        obj <- c(obj, -0.5*log(2*pi)-0.5*log(sd^2)-resid2[i]^2/(2*sd^2))
      }
    }
  } else {
    for (i in 1:N) {
      if (i==1) {
        obj <- c(obj,-log(2*sqrt(sd^2/(1-rho^2)))-abs(resid2[i])/sqrt(sd^2/(1-rho^2)))
      } else {
        obj <- c(obj, -log(2*sd)-abs(resid2[i])/sd)
      }
    }
  }
  obj <- sum(w*obj) # exact likelihood
  return(exp(obj))
}

#' レジーム分けを考慮した再生産関係の推定
#' 
#' レジームシフトが生じた年やレジームであるパラメータが共通する場合やレジームのパターンがA->B->CなのかA->B->Aなのか等が検討できる
#' @param SRdata \code{get.SRdata}で作成した再生産関係データ
#' @param SR 再生産関係 (\code{"HS"}: Hockey-stick, \code{"BH"}: Beverton-Holt, \code{"RI"}: Ricker)
#' @param method 最適化法（\code{"L2"}: 最小二乗法, \code{"L1"}: 最小絶対値法）
#' @param regime.year レジームが変わる年
#' @param regime.key レジームのパターンを表す(\code{0:2}だとA->B->Cで、\code{c(0,1,0)}だとA->B->Aのようなパターンとなる)
#' @param regime.par レジームによって変化するパラメータ(\code{c("a","b","sd")}の中から選ぶ)
#' @param length 初期値を決める際のgridの長さ
#' @param p0 \code{optim}で設定する初期値
#' @inheritParams fit.SR
#' @encoding UTF-8
#' @examples 
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSRregime <- fit.SRregime(SRdata, SR="HS", method="L2", 
#'                             regime.year=c(1977,1989), regime.key=c(0,1,0),
#'                             regime.par = c("a","b","sd")[2:3])
#' resSRregime$regime_pars
#' }
#' @return 以下の要素からなるリスト
#' \describe{
#' \item{\code{input}}{使用した引数のリスト}
#' \item{\code{pars}}{推定されたパラメータ}
#' \item{\code{opt}}{\code{optim}の結果オブジェクト}
#' \item{\code{resid}}{再生産関係から予測値からの加入量の残差}
#' \item{\code{loglik}}{対数尤度}
#' \item{\code{k}}{推定したパラメータ数}
#' \item{\code{AIC}}{AIC}
#' \item{\code{AICc}}{AICc}
#' \item{\code{BIC}}{BIC}
#' \item{\code{regime_pars}}{レジームごとの推定パラメータ}
#' \item{\code{regime_resid}}{レジームごとの残差}
#' \item{\code{pred}}{レジームごとの各親魚量に対する加入量の予測値}
#' \item{\code{pred_to_obs}}{観測値に対する予測値}
#' \item{\code{summary_tbl}}{観測値と予測値を合わせた表}
#' }
#'
#' @export


fit.SRregime <- function(
  SRdata,
  SR = "HS",
  method = "L2",
  regime.year = NULL,
  regime.key = 0:length(regime.year),
  regime.par = c("a","b","sd"),
  use.fit.SR = TRUE,
  length=10,  # parameter (a,b) の初期値を決めるときにgrid searchする数
  p0 = NULL,  # 初期値
  w = rep(1,length(SRdata$year)),
  max.ssb.pred = 1.3
) {
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  
  rec <- SRdata$R
  ssb <- SRdata$SSB
  
  N <- length(rec)
  
  regime.key <- regime.key-min(regime.key)+1
  regime <- a_key <- b_key <- sd_key <- rep(1,N) 
  if (!is.null(regime.year)) {
    for(i in 1:length(regime.year)) {
      regime[SRdata$year>=regime.year[i]] <- regime.key[i+1]
    }
  }
  if ("a" %in% regime.par) a_key <- regime
  if ("b" %in% regime.par) b_key <- regime
  if ("sd" %in% regime.par) sd_key <- regime
  
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  
  obj.f <- function(a,b,out="nll"){ #a,bはベクトル
    resid <- NULL
    for(i in 1:N) {
      pred_rec <- SRF(ssb[i],a[a_key[i]],b[b_key[i]])
      resid <- c(resid, log(rec[i]/pred_rec))
    }
    if (method == "L2") {
      tbl = tibble(resid=resid,sd_key=sd_key,w=w) %>%
        mutate(resid2 = resid^2) %>% 
        group_by(sd_key) %>%
        summarise(rss = sum(w*resid2), n = sum(w)) %>%
        mutate(sd = sqrt(rss/n))
    } else {
      tbl = tibble(resid=resid,sd_key=sd_key,w=w) %>%
        # mutate(resid2 = resid^2) %>% 
        group_by(sd_key) %>%
        summarise(rss = sum(w*abs(resid)), n = sum(w)) %>%
        mutate(sd = rss/n)
    }
    SD = tbl$sd %>% as.numeric()
    nll <- 
      ifelse(method=="L2",
             -sum(sapply(1:N, function(i) dnorm(resid[i],0,SD[sd_key[i]], log=TRUE))),
             -sum(sapply(1:N, function(i) -log(2*SD[sd_key[i]])-abs(resid[i]/SD[sd_key[i]])))
      )
    if (out=="nll") return(nll)
    if (out=="resid") return(resid)
    if (out=="sd") return(SD)
  }
  
  a_grid <- NULL
  for(i in unique(a_key)){
    a_range<-range(rec[a_key==i]/ssb[a_key==i])
    a_grid <- cbind(a_grid,seq(a_range[1],a_range[2],length=length))
  }
  b_grid <- NULL
  for(i in unique(b_key)){
    if (SR=="HS") {
      b_range <- range(ssb[b_key==i])
    } else {b_range <- range(1/ssb[b_key==i])}
    b_grid <- cbind(b_grid,seq(b_range[1],b_range[2],length=length))
  }
  ab_grid <- expand.grid(data.frame(a_grid,b_grid)) %>% as.matrix()
  b_range <- apply(b_grid,2,range)
  
  if (is.null(p0)) {
    if (use.fit.SR) {
      fit_SR_res = fit.SR(SRdata, SR = SR, method = method, w = w, AR = 0)
      init <- c(rep(fit_SR_res$opt$par[1],max(a_key)),rep(fit_SR_res$opt$par[2],max(b_key)))
    } else {
      init_list <- sapply(1:nrow(ab_grid), function(j) {
        obj.f(a=ab_grid[j,1:max(a_key)],b=ab_grid[j,(1+max(a_key)):(max(a_key)+max(b_key))])
      })
      
      ab_init <- as.numeric(ab_grid[which.min(init_list),])
      init <- log(ab_init[1:max(a_key)])
      if (SR=="HS") {
        for(i in unique(b_key)) {
          init <- c(init,-log(max(0.000001,(max(ssb[b_key==i])-min(ssb[b_key==i]))/max(0.000001,(ab_init[max(a_key)+i]-min(ssb[b_key==i])))-1)))
        }
      } else {
        init <- c(init, log(ab_init[(max(a_key)+1):(max(a_key)+max(b_key))]))
      }
    }
  } else {
    init <- p0
  }
  
  if (SR=="HS") {
    obj.f2 <- function(x) {
      a <- exp(x[1:max(a_key)])
      b <- b_range[1,]+(b_range[2,]-b_range[1,])/(1+exp(-x[(max(a_key)+1):(max(a_key)+max(b_key))]))
      return(obj.f(a,b))
    }
  } else {
    obj.f2 <- function(x) obj.f(a=exp(x[1:max(a_key)]),b=exp(x[(1+max(a_key)):(max(a_key)+max(b_key))]))
  }
  
  opt <- optim(init,obj.f2)
  for (i in 1:100) {
    opt2 <- optim(opt$par,obj.f2)
    if (abs(opt$value-opt2$value)<1e-6) break
    opt <- opt2
  }
  opt <- optim(opt$par,obj.f2,method="BFGS")
  
  Res <- list()
  Res$input <- arglist
  Res$opt <- opt
  
  a <- exp(opt$par[1:max(a_key)])
  if (SR=="HS") {
    b <- b_range[1,]+(b_range[2,]-b_range[1,])/(1+exp(-opt$par[(1+max(a_key)):(max(a_key)+max(b_key))]))
  } else {
    b <- exp(opt$par[(1+max(a_key)):(max(a_key)+max(b_key))])
  }
  sd <- obj.f(a,b,out="sd")
  if (method=="L1") sd <- sqrt(2)*sd
  resid <- obj.f(a,b,out="resid")
  
  Res$resid <- resid
  Res$pars$a <- a
  Res$pars$b <- b
  Res$pars$sd <- sd
  
  Res$loglik <- loglik <- -opt$value
  Res$k <- k <- length(opt$par)+max(sd_key)
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(sum(w>0)-k-1)
  Res$BIC <- -2*loglik+k*log(sum(w>0))
  
  Res$regime_pars <- tibble(regime=regime,a=a[a_key],b=b[b_key],sd=sd[sd_key]) %>% distinct()
  Res$regime_resid <- tibble(regime=regime,resid = resid)
  
  ssb.tmp <- seq(from=0,to=max(ssb)*max.ssb.pred,length=100)
  ab_unique <- unique(cbind(a_key,b_key))
  summary_tbl = tibble(Year = SRdata$year,SSB=ssb, R = rec, Regime=regime, Category = "Obs")
  for (i in 1:nrow(ab_unique)) {
    R.tmp <- sapply(1:length(ssb.tmp), function(j) SRF(ssb.tmp[j],a[ab_unique[i,1]],b[ab_unique[i,2]]))
    summary_tbl = bind_rows(summary_tbl,tibble(Year=NA,SSB=ssb.tmp, R=R.tmp, Regime=i, Category="Pred"))
  }
  summary_tbl = summary_tbl %>% mutate(Regime = factor(Regime))
  pred = dplyr::filter(summary_tbl, Category == "Pred") %>%
    dplyr::select(-Year, -Category) %>%
    dplyr::select(Regime,SSB,R)
  Res$pred <- pred
  pred_to_obs = dplyr::filter(summary_tbl, Category == "Obs") %>%
    dplyr::select(-Category) %>%
    mutate(resid = resid) %>%
    mutate(Pred = exp(log(R)-resid)) %>%
    dplyr::select(Year,SSB,R,Regime,Pred,resid)
  Res$pred_to_obs <- pred_to_obs
  Res$summary_tbl
  
  return(Res)
}