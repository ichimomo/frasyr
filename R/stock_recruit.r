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
    assertthat::assert_that(all(dat[["R"]] > 0))
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
#' @param p0 \code{optim}で設定する初期値
#' @encoding UTF-8
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                 out.AR = TRUE)
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
#' 

fit.SR <- function(SRdata,
                   SR="HS",
                   method="L2",
                   AR=1,
                   # TMB=FALSE,
                   hessian=FALSE,w=rep(1,length(SRdata$R)),
                   length=20,
                   max.ssb.pred=1.3, # 予測値を計算するSSBの最大値（観測された最大値への乗数）
                   p0=NULL,
                   out.AR = TRUE #自己相関係数rhoを外で推定するか
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

  if (length(SRdata$R) != length(w)) stop("The length of 'w' is not appropriate!")

  one_max = max(SRdata$year[w>0])
  zero_min =ifelse(sum(w==0)>0, min(SRdata$year[w==0]),one_max)

  if (method == "L2" && AR==1 && out.AR==FALSE && zero_min<one_max) { # For Jackknife

    obj.f = function(a,b,rho,out = "nll") {
      w2 = as.numeric(w>0)
      before_zero = rep(0,N)
      for (i in 1:N) {
        if (w2[i]>0) {
          if (i == 1) {
            w2[i] <- w[i]*(1-rho^2)
          } else {
            for (j in 1:(i-1)) {
              if (rev(w[1:(i-1)])[j] == 0) {
                before_zero[i] = before_zero[i]+1
                w2[i] <- w2[i]+rho^(2*j)
              } else break
            }
            if (before_zero[i] == i-1) {
              w2[i] <- w[i]*(1-rho^2)
            } else {
              w2[i] <- w[i]/w2[i]
            }
          }
        }
      }

      resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
      pred_resid <- NULL
      for (i in 1:N) {
        if (i==1 || before_zero[i] == i-1) {
          pred_resid <- c(pred_resid,0)
        } else {
          pred_resid <- c(pred_resid,rho^(1+before_zero[i])*resid[i-1-before_zero[i]])
        }
      }
      sd2 = sum(w2*(resid-pred_resid)^2)/sum(w2) #variance
      sd = sqrt(sd2) #SD

      sigma = NULL
      for (i in 1:N) {
        if (w2[i]==0) {
          sigma = c(sigma,1000)
        } else {
          sigma = c(sigma,sqrt(sd2/w2[i]))
        }
      }

      nll <- -sum(w*dnorm(resid,pred_resid,sigma,log=TRUE))
      if (out=="nll") return(nll)
      if (out=="resid") return(resid)
      if (out=="resid2") return(resid-pred_resid)
      if (out=="sd") return(sd)
    }
  } else {
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
  #if (rep.opt) {
    for (i in 1:100) {
      opt2 <- optim(opt$par,obj.f2)
      if (abs(opt$value-opt2$value)<1e-6) break
      opt <- opt2
    }
  #}
  opt <- optim(opt$par,obj.f2,method="BFGS",hessian=hessian)

  Res <- list()
  Res$input <- arglist
  Res$obj.f <- obj.f
  Res$obj.f2 <- obj.f2
  Res$opt <- opt

  a <- exp(opt$par[1])
  b <- ifelse(SR=="HS",min(ssb)+(max(ssb)-min(ssb))/(1+exp(-opt$par[2])),exp(opt$par[2]))
  rho <- ifelse(AR==0,0,ifelse(out.AR,0,1/(1+exp(-opt$par[3]))))

  if (method == "L2" && AR==1 && out.AR==FALSE && zero_min<one_max) {
    resid = obj.f(a=a,b=b,rho=rho,out="resid")
    resid2 = obj.f(a=a,b=b,rho=rho,out="resid2")
    sd <- sd.pred <- obj.f(a=a,b=b,rho=rho,out="sd")
  } else {
    resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
    resid2 <- NULL
    for (i in 1:N) {
      resid2[i] <- ifelse(i == 1,resid[i], resid[i]-rho*resid[i-1])
    }

    # if (method=="L2") {
      rss <- w[1]*resid2[1]^2*(1-rho^2)
      for(i in 2:N) rss <- rss + w[i]*resid2[i]^2
      sd <- sd.pred <- sqrt(rss/NN)
    # } else {
      if (method=="L1") {
        rss <- w[1]*abs(resid2[1])*sqrt(1-rho^2)
        for(i in 2:N) rss <- rss + w[i]*abs(resid2[i])
        sd.pred <- sum(abs(w*resid2))/NN
        sd.pred <- sqrt(2)*sd.pred
    }
    # sd <- ifelse(method=="L2",sqrt(sum(w*resid2^2)/(NN-rho^2)),sqrt(2)*sum(abs(w*resid2))/(NN-rho^2))
  }

  Res$resid <- resid
  Res$resid2 <- resid2
  Res$sd.pred = sd.pred
  # if (sd.obs) {
    # if (AR==0 && method != "L2") sd = sqrt(sum(resid[w==1]^2)/sum(w))
  # }
  Res$pars <- c(a,b,sd,rho)

  if (method!="L2") {
    if (AR!=0) {
      if (!isTRUE(out.AR)) {
        message("L1 & out.AR=FALSE is NOT recommended")
      }
      arres <- ar(resid,aic=FALSE,order.max=1,demean=FALSE,method="mle")
      Res$pars[3] <- ifelse(arres$ar<0,sd,sqrt(arres$var.pred))
      Res$pars[4] <- ifelse(arres$ar<0,0,arres$ar)
    }
  }

  if (AR==1 && out.AR) {
    arres <- ar(resid,aic=FALSE,order.max=1,demean=FALSE,method="mle")
    Res$pars[3] <- Res$sd.pred <-sqrt(arres$var.pred)
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

  class(Res) <- "fit.SR"
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

#' 再生産関係のブートストラップ
#'
#' ①残差のパラメトリックブートストラップ、②残差のノンパラメトリックブートストラップ（リサンプリング）、③データのブートストラップ（リサンプリング）が行える
#' @import purrr
#' @param Res \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param method パラメトリック ("p") かノンパラメトリック ("n")
#' @encoding UTF-8
#' @export
#'
boot.SR <- function(Res,method="p",n=100,seed=1){

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  N <- length(Res$input$SRdata$SSB)

  if (Res$input$SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (Res$input$SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (Res$input$SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)

  set.seed(seed)
  if (class(Res) == "fit.SR") { #fit.SR
    N <- length(Res$input$SRdata$SSB)
    RES = lapply(1:n, function(j){
        sd <- sapply(1:N, function(i) ifelse(i==1,Res$pars$sd/sqrt(1-Res$pars$rho^2),Res$pars$sd))
        for (j in 1:100) {
          if (method=="d") { #data bootstrap
            Year0 = Res$input$SRdata$year[Res$input$w==1]
            Year_r = sample(Year0,replace=TRUE)
            w_r = purrr::map_dbl(Year0, function(x) sum(Year_r==x))
            res.b <- Res
            res.b$input$w <- w_r
          } else {
            if (method=="p") { # parametric bootstrap assuming a normal distribution
              resids <- rnorm(N,0,sd)
            } else {# non-parametric bootstrap for residuals
              std.resid = calc.StdResid(Res)$std.resid
              std.resids = sample(std.resid,replace=TRUE)
              resids = std.resids*sd
            }
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
          }
          res.b$input$p0 = Res$opt$par
          res.b <- try(do.call(fit.SR, res.b$input))
          if (class(res.b) != "try-error") break
        }

      return(res.b)
    })
  } else {
    # fit.SRregime
    N <- length(Res$input$SRdata$SSB)
    tmp = calc.StdResid(Res)
    sd = tmp$sigma
    std.resid = tmp$std.resid
    merged = full_join(Res$regime_pars,mutate(Res$regime_resid,Year=Res$input$SRdata$year),by="regime") %>%
      arrange(Year)
    RES = lapply(1:n, function(j){
      for (k in 1:100) {
        if (method=="d") { #data bootstrap
          Year0 = Res$input$SRdata$year[Res$input$w==1]
          Year_r = sample(Year0,replace=TRUE)
          w_r = purrr::map_dbl(Year0, function(x) sum(Year_r==x))
          res.b <- Res
          res.b$input$w <- w_r
        } else {
          if (method=="p") { # parametric bootstrap assuming a normal distribution
            std.resids = rnorm(N,0,1)
          } else {# non-parametric bootstrap for residuals
            std.resids = sample(std.resid,replace=TRUE)
          }
          resids <- std.resids*sd
          pred <- obs <- resid0 <- numeric(N)
          ssb <- Res$input$SRdata$SSB
          for(i in 1:N){
            pred[i] <- SRF(ssb[i],merged$a[i],merged$b[i])
            if (i==1) {
              obs[i] <- pred[i]*exp(resids[i])
            } else {
              obs[i] <- pred[i]*exp(0*resid0[i-1])*exp(resids[i])
            }
            resid0[i] <- log(obs[i]/pred[i])
          }
          res.b <- Res
          res.b$input$SRdata$R <- obs
        }
        res.b$input$p0 = Res$opt$par
        res.b <- try(do.call(fit.SRregime, res.b$input))
        if (class(res.b) != "try-error") break
      }

      return(res.b)
    })
  }
  RES$input = arglist
  return(RES)
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
#'                             regime.year=c(1995,2005), regime.key=c(0,1,0),
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
  # sd.obs = TRUE,
  use.fit.SR = TRUE,
  length=10,  # parameter (a,b) の初期値を決めるときにgrid searchする数
  p0 = NULL,  # 初期値
  w = rep(1,length(SRdata$R)),
  max.ssb.pred = 1.3,
  hessian = FALSE
) {
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  rec <- SRdata$R
  ssb <- SRdata$SSB
  N <- length(rec)

  regime.key0 = regime.key
  unique.key = unique(regime.key)
  regime.key = sapply(1:length(regime.key0), function(i) which(unique.key == regime.key0[i]))
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
             -sum(sapply(1:N, function(i) w[i]*dnorm(resid[i],0,SD[sd_key[i]], log=TRUE))),
             -sum(sapply(1:N, function(i) w[i]*(-log(2*SD[sd_key[i]])-abs(resid[i]/SD[sd_key[i]]))))
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
  opt <- optim(opt$par,obj.f2,method="BFGS",hessian=hessian)

  Res <- list()
  Res$input <- arglist
  Res$opt <- opt

  a <- exp(opt$par[1:max(a_key)])
  if (SR=="HS") {
    b <- b_range[1,]+(b_range[2,]-b_range[1,])/(1+exp(-opt$par[(1+max(a_key)):(max(a_key)+max(b_key))]))
  } else {
    b <- exp(opt$par[(1+max(a_key)):(max(a_key)+max(b_key))])
  }
  sd.pred <- sd <- obj.f(a,b,out="sd")
  if (method=="L1") {
    # L1の場合sd.predとsdは定義が異なる. sdの値は825-834行目は上書きされることに注意
    sd.pred <- sd <- sqrt(2)*sd.pred
  }
  # sd <- sqrt(sum(w*obj.f(a,b,out="resid")^2)/sum(w))
  resid <- obj.f(a,b,out="resid")

  Res$obj.f <- obj.f
  Res$obj.f2 <- obj.f2
  Res$resid <- resid
  Res$pars$a <- a
  Res$pars$b <- b
  Res$pars$sd <- sd

  Res$loglik <- loglik <- -opt$value
  Res$k <- k <- length(opt$par)+max(sd_key)
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(sum(w>0)-k-1)
  Res$BIC <- -2*loglik+k*log(sum(w>0))

  Res$regime_pars <- tibble(regime=regime.key0[regime],a=a[a_key],b=b[b_key],sd=sd[sd_key]) %>% distinct()
  Res$regime_resid <- tibble(regime=regime.key0[regime],resid = resid)
  # if (sd.obs) {
  Res$sd.pred = sd.pred
  if (method=="L1") {
      if ("sd" %in% regime.par) {
        tmp = Res$regime_resid %>% mutate(w = w, squared_resid = w*resid^2) %>%
          group_by(regime) %>% summarise(n = sum(w), RSS = sum(squared_resid)) %>%
          mutate(RMSE = sqrt(RSS/n))
        sd = as.numeric(tmp$RMSE)
      } else {
        sd = sqrt(sum(resid[w==1]^2)/sum(w))
      }
      Res$pars$sd = sd
      Res$regime_pars$sd = sd
    }
  # }
  ssb.tmp <- seq(from=0,to=max(ssb)*max.ssb.pred,length=100)
  ab_unique <- unique(cbind(a_key,b_key))
  summary_tbl = tibble(Year = SRdata$year,SSB=ssb, R = rec, Regime=regime.key0[regime], Category = "Obs")
  for (i in 1:nrow(ab_unique)) {
    R.tmp <- sapply(1:length(ssb.tmp), function(j) SRF(ssb.tmp[j],a[ab_unique[i,1]],b[ab_unique[i,2]]))
    summary_tbl = bind_rows(summary_tbl,tibble(Year=NA,SSB=ssb.tmp, R=R.tmp, Regime=unique.key[i], Category="Pred"))
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
  class(Res) <- "fit.SRregime"

  return(Res)
}

#' 再生産関係の推定における標準化残差を計算する関数
#' @import rmutil
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @encoding UTF-8
#' @export
calc.StdResid = function(resSR) {
  if(class(resSR) == "fit.SR") { #fit.SR
    if (resSR$input$method == "L2") {
      sigma = rep(sqrt(sum(resSR$resid^2)/length(resSR$resid)),length(resSR$resid))
    } else {
      sigma = rep(sqrt(2)*sum(abs(resSR$resid))/length(resSR$resid),length(resSR$resid))
    }
    sigma2 = c(sqrt(resSR$pars$sd^2/(1-resSR$pars$rho^2)), rep(resSR$pars$sd,length(resSR$resid)-1))
    std.resid = resSR$resid/sigma
    std.resid2 = resSR$resid2/sigma2
    if (resSR$input$method == "L2") {
      cumulative.prob = pnorm(std.resid,0,1)
    } else {
      cumulative.prob = rmutil::plaplace(std.resid,0,s=1/sqrt(2))
    }
    if (resSR$input$AR) {
      cumulative.prob2 = pnorm(std.resid2,0,1)
    } else {
      cumulative.prob2 = cumulative.prob
    }
    RES = tibble(sigma,sigma2,std.resid,std.resid2,cumulative.prob,cumulative.prob2)
  } else{ #fit.SRregime
    RES = dplyr::full_join(resSR$regime_pars,
                           resSR$regime_resid %>% mutate(Year = resSR$input$SRdata$year),by="regime") %>%
      dplyr::arrange(Year) %>%
      dplyr::mutate(std.resid = resid/sd) %>%
      dplyr::select(sd,std.resid) %>%
      rename(sigma=sd)
    if (resSR$input$method == "L2") {
      cumulative.prob = pnorm(RES$std.resid,0,1)
    } else {
      cumulative.prob = rmutil::plaplace(RES$std.resid,0,s=1/sqrt(2))
    }
    RES = RES %>% mutate(cumulative.prob=cumulative.prob)
  }
  return(RES)
}


#' 再生産関係の残差の確率分布に関するチェック図を出力する関数
#'
#' 1) 正規性をチェックするための検定結果・ヒストグラム、2) 累積確率分布のヒストグラム、3) 一様分布を使用した累積確率分布のQQプロttの3つが出力される
#' 標準化した残差を使用する
#' 累積確率分布はL2の場合は正規分布、L1の場合はラプラス分布を使用する
#' 自己相関の外側推定の場合は、2段階で推定しているため、ファイルが2つ出力される
#' @import EnvStats
#' @import rmutil
#' @inheritParams calc.StdResid
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param output pngファイルに出力するか否か
#' @encoding UTF-8
#' @export
check.SRdist = function(resSR,test.ks=TRUE,output=FALSE,filename = "SR_error_dist") {
  std.resid.res = calc.StdResid(resSR)
  main_name=paste0(resSR$input$SR," ",resSR$input$method," ")

  if (class(resSR)=="fit.SR" && resSR$input$AR && isTRUE(resSR$input$out.AR)) {
    for(i in 1:2) {
      if (i==1) {
        std.resid = std.resid.res$std.resid
        cumulative.prob = std.resid.res$cumulative.prob
        main_name2 = "Std. Deviance to SR"
      } else {
        std.resid = std.resid.res$std.resid2
        cumulative.prob = std.resid.res$cumulative.prob2
        main_name2 = "Standard. Resid."
      }
      # pdf(file = paste0(filename,"(",main_name2,").pdf"), width=15,height=5)
      if (output) png(file = paste0(filename,"(",main_name2,").png"), width=15, height=5, res=432, units='in')
      par(pch=1,lwd = 2, mfrow=c(1,3),cex=1)
      check1 <- shapiro.test(std.resid)
      hist(std.resid,main = paste0(main_name,main_name2),
           xlab = "Standardized Deviance",freq=FALSE,col="gray")
      X <- seq(min(std.resid)*1.3,max(std.resid)*1.3,length=200)
      points(X,dnorm(X,0,1),col=2,lwd=3,type="l")
      if (resSR$input$method=="L1" && i==1) {
        points(X,rmutil::dlaplace(X,0,1/sqrt(2)),col="blue",lwd=2,type="l",lty=2)
      }
      mtext(text=" P value",adj=1,line=-1,lwd=2,font=2)
      mtext(text=sprintf(" SW: %1.3f",check1$p.value),adj=1,line=-2)
      if (test.ks) {
        if (i==1 && resSR$input$method=="L1") {
          check2 <- ks.test(std.resid,y=function(x) rmutil::plaplace(x,m=0,s=1/sqrt(2)))
        } else {
          check2 <- ks.test(std.resid,y="pnorm",sd=1)
        }
        mtext(text=sprintf(" KS: %1.3f",check2$p.value),adj=1,line=-3)
      }
      hist(cumulative.prob, xlab="Cumulative probability",
           main = paste0(main_name,"Cumlative Prob. Dist."),freq=FALSE,breaks=seq(0,1,by=0.1),col="gray")
      # if (test.ks) {
      #   check2 <- ks.test(cumulative.prob,y="punif",min=0,max=1)
      #   mtext(text=" P value",adj=1,line=-1,lwd=2,font=2)
      #   mtext(text=sprintf(" KS: %1.3f",check2$p.value),adj=1,line=-2)
      # }
      EnvStats::qqPlot(x=cumulative.prob,distribution="unif", param.list = list(min = 0, max = 1),
                       add.line = TRUE,qq.line.type = "0-1",line.lwd=1.5,cex=1.2,line.col="red",
                       main = paste0(main_name,"Uniform QQ plot"),ylab="Quatiles of Cumlative.Prob.",
                       xlab = "Quantiles of Uniform Dist.")
      if (output) dev.off()
    }
  } else {
    if (is.null(std.resid.res$std.resid2)) {
      std.resid = std.resid.res$std.resid
      cumulative.prob = std.resid.res$cumulative.prob
    } else {
      std.resid = std.resid.res$std.resid2
      cumulative.prob = std.resid.res$cumulative.prob2
    }
    # pdf(file = paste0(filename,".pdf"), width=15,height=5)
    if (output) png(file = paste0(filename,".png"), width=15, height=5, res=432, units='in')
    par(pch=1,lwd = 2, mfrow=c(1,3),cex=1)
    check1 <- shapiro.test(std.resid)
    hist(std.resid,main = paste0(main_name,"Standard. Resid."),
         xlab = "Standardized Residual",freq=FALSE,col="gray")
    X <- seq(min(std.resid)*1.3,max(std.resid)*1.3,length=200)
    points(X,dnorm(X,0,1),col=2,lwd=3,type="l")
    if (resSR$input$method=="L1") {
      points(X,rmutil::dlaplace(X,0,1/sqrt(2)),col="blue",lwd=2,type="l",lty=2)
    }
    mtext(text=" P value",adj=1,line=-1,lwd=2,font=2)
    mtext(text=sprintf(" SW: %1.3f",check1$p.value),adj=1,line=-2)
    if (test.ks) {
      if (resSR$input$method=="L1"){
        check2 <- ks.test(std.resid,y=function(x) rmutil::plaplace(x,m=0,s=1/sqrt(2)))
      } else {
        check2 <- ks.test(std.resid,y="pnorm",mean=0,sd=1)
      }
      mtext(text=sprintf(" KS: %1.3f",check2$p.value),adj=1,line=-3)
    }
    hist(cumulative.prob, xlab="Cumulative probability",
         main = paste0(main_name,"Cumlative Prob. Dist."),freq=FALSE,breaks=seq(0,1,by=0.1),col="gray")
    # if (test.ks) {
    #   check2 <- ks.test(cumulative.prob,y="punif",min=0,max=1)
    #   mtext(text=" P value",adj=1,line=-1,lwd=2,font=2)
    #   mtext(text=sprintf(" KS: %1.3f",check2$p.value),adj=1,line=-2)
    # }
    EnvStats::qqPlot(x=cumulative.prob,distribution="unif", param.list = list(min = 0, max = 1),
                     add.line = TRUE,qq.line.type = "0-1",line.lwd=1.5,cex=1.2,line.col="red",
                     main = paste0(main_name,"Uniform QQ plot"),ylab="Quatiles of Cumlative.Prob.",
                     xlab = "Quantiles of Uniform Dist.")
    if (output) dev.off()
  }
}

#' 再生産関係の残差から事後的に1次の自己相関係数を推定する関数
#'
#' 計算は\code{fit.SR}の\code{out.AR}と同じであるが、AICcも見れるほか、\code{fit.SRregime}にも対応している
#' \code{fit.SRregime}の場合、各レジームの最初の年は初期値となる（つまり前のレジームの最後の年からの残差を引きずらない）
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param per_regime 自己相関係数をレジームごとに推定するか (\code{TRUE}) 否か
#' @encoding UTF-8
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[1], AR = 0,
#'                 out.AR = FALSE)
#' resSR_post = calc.residAR(resSR)
#' resSR_post$AICc
#' resSR_post$pars
#' resSRregime <- fit.SRregime(SRdata, SR="HS", method="L2",
#'                             regime.year=c(1995,2005), regime.key=c(0,1,0),
#'                             regime.par = c("a","b","sd")[2:3])
#' resSRregime_post = calc.residAR(resSRregime, per_regime = TRUE)
#' resSRregime_post$AICc
#' resSRregime_post$regime$pars
#' }
#'
#' @export
calc.residAR = function(resSR, per_regime=TRUE, output=TRUE, filename="residARouter") {
  RES = list()
  if (class(resSR) == "fit.SR") { #fit.SR
    if (resSR$input$AR && !isTRUE(resSR$input$out.AR)) {
      warning("This function is meaningless when AR=TRUE & out.AR=FALSE")
    }
    deviance = resSR$resid
    arimares0 = arima(deviance,method="ML",include.mean=FALSE,order=c(0,0,0))
    arimares1 = arima(deviance,method="ML",include.mean=FALSE,order=c(1,0,0))
    nn = length(deviance)
    k = c("AR(0)"=1,"AR(1)"=2)
    loglik = c("AR(0)"=arimares0$loglik,"AR(1)"=arimares1$loglik)
    AIC = c("AR(0)"=arimares0$aic,"AR(1)"=arimares1$aic)
    AICc = AIC + 2*k*(k+1)/(nn-k-1)
    BIC = AIC + k*log(nn)
    rho = as.numeric(arimares1$coef)
    sd = sqrt(arimares1$sigma2)
    resid2 = arimares1$residuals
    RES$pars = resSR$pars
    RES$pars$rho = rho
    RES$pars$sd = sd
    RES$arima0 = arimares0
    RES$arima1 = arimares1
    RES$resid = resSR$resid
    RES$resid2 = resid2
    RES$loglik = loglik
    RES$k = k
    RES$AIC = AIC
    RES$AICc = AICc
    RES$BIC = BIC
  } else {  #fit.SRregime
    tbl = resSR$regime_pars %>% full_join(resSR$regime_resid %>% mutate(Year = resSR$input$SRdata$year)) %>%
      arrange(Year) %>% mutate(Shift = if_else(Year %in% c(min(Year), resSR$input$regime.year),1,0))
    tbl = tbl %>%
      mutate(regime_id = purrr::map_dbl(regime, function(x) which(x==unique(resSR$regime_pars$regime))))
    deviance = tbl$resid; shift = tbl$Shift; regime_id = tbl$regime_id;
    w = resSR$input$w;nn = nrow(tbl);
    obj.f = function(x,out="nll") {
      rho = (exp(x)-1.0)/(exp(x)+1.0);resid2 = c();
      if (per_regime) {
        for (i in 1:nn) {
          if (shift[i]==1) {
            resid2 = c(resid2,sqrt(1-rho[regime_id[i]]^2)*deviance[i])
          } else {
            resid2 = c(resid2,deviance[i]-rho[regime_id[i-1]]*deviance[i-1])
          }
        }
      } else {
        for (i in 1:nn) {
          if (shift[i]==1) {
            resid2 = c(resid2,sqrt(1-rho^2)*deviance[i])
          } else {
            resid2 = c(resid2,deviance[i]-rho*deviance[i-1])
          }
        }
      }
      RSS = tbl %>% mutate(resid2 = resid2,w=w) %>% group_by(sd) %>%
        summarise(var=sum(w*resid2^2),n=sum(w)) %>%
        mutate(sigma = sqrt(var/n))
      tbl2 = left_join(tbl,RSS,by="sd")
      sigma = tbl2$sigma
      if (out=="nll") {
        return(-sum(dnorm(0,resid2,sigma,log=TRUE))) } else {
          if (out=="resid") {
            return(resid2) } else {
              return(sigma)
            }
        }
    }
    if (per_regime) {
      opt = optim(rep(0,max(tbl$regime_id)),obj.f)
      for (i in 1:100) {
        opt2 <- optim(opt$par,obj.f)
        if (abs(opt$value-opt2$value)<1e-6) break
        opt <- opt2
      }
      opt <- optim(opt$par,obj.f,method="BFGS")
    } else {
      opt = optim(0,obj.f,method="Brent",lower=-20,upper=20)
    }
    rho = (exp(opt$par)-1)/(exp(opt$par)+1)
    loglik1 = -opt$value
    resid2 = obj.f(x=opt$par,out="resid")
    sigma = obj.f(x=opt$par,out="sigma")
    tbl3 = tbl %>% mutate(resid2 = resid2,sigma=sigma)
    RES$pars = resSR$pars
    RES$pars$rho = rho
    if (per_regime) rho = rho[tbl3$regime_id]
    tbl3 = tbl3 %>% mutate(rho = rho)
    RES$regime_pars = tbl3 %>% dplyr::select(regime,a,b,sigma,rho) %>% distinct() %>% rename(sd = sigma)
    RES$regime_resid = tbl3 %>% dplyr::select(regime,resid,resid2)
    RES$pars$sd = unique(RES$regime_pars$sd)
    RES$opt = opt
    if (per_regime) {
      loglik0 = -obj.f(rep(0,max(tbl$regime_id)))
    } else {
      loglik0 = -obj.f(0)
    }
    k = c("AR(0)"=length(resSR$pars$sd),"AR(1)"=length(resSR$pars$sd)+length(opt$par))
    loglik = c("AR(0)"=loglik0,"AR(1)"=loglik1)
    AIC = -2*loglik+2*k
    AICc = AIC + 2*k*(k+1)/(nn-k-1)
    BIC = AIC + k*log(nn)
    RES$loglik = loglik
    RES$k = k
    RES$AIC = AIC
    RES$AICc = AICc
    RES$BIC = BIC
  }
  if (output) {
    capture.output(RES, file = paste0(filename,".txt"))
  }
  return(RES)
}

#' 再生産関係における残差の時系列の自己相関等についてのプロット
#'
#' 1) 残差のトレンド、2) \code{acf}関数による自己相関係数のプロット、3) Ljung-Box検定におけるP値の3つの図を出力
#' @inheritParams calc.StdResid
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param use.resid 再生産関係との残差 (deviance \code{resid}) を使うか (1:default)、自己相関を除いた残差 (\code{resid2}) を使うか (2)
#' @param output pngファイルに出力するか否か
#' @encoding UTF-8
#' @export
autocor.plot = function(resSR,use.resid=1,lag.max=NULL,output = FALSE,filename = "Residual_trend",pch=16,lwd=2,cex=1.2,cex.main=1.2,cex.lab=1.2,...){
  Year = resSR$input$SRdata$year
  if (output) png(file = paste0(filename,".png"), width=15, height=5, res=432, units='in')
  par(pch=pch,lwd = lwd, mfrow=c(1,3),cex=cex)
  if (class(resSR) == "fit.SR") { #fit.SR
    if (use.resid==1) {
      Resid = resSR$resid
      plot(Year,Resid,pch=pch,main="",xlab="Year",ylab="Deviance",cex.lab=cex.lab,...)
      title("Time series of deviance to SR",cex.main=cex.main)
    } else {
      Resid = resSR$resid2
      Resid[1] = sqrt(1-resSR$pars$rho^2)*Resid[1]
      plot(Year,Resid,pch=pch,main="",xlab="Year",ylab="Residual",cex.lab=cex.lab,...)
      title("Time series of Residuals",cex.main=cex.main)
    }
  } else { #fit.SRregime
    message("Standardized residuals are used for 'fit.SRregime'")
    table = calc.StdResid(resSR)
    Resid = table$std.resid
    plot(Year,Resid,pch=pch,main="",xlab="Year",ylab="Std.Residual",cex.lab=cex.lab,...)
    title("Time series of Standardized Residuals",cex.main=cex.main)
    abline(v=c(resSR$input$regime.year)-0.5,lty=3,col="blue")
  }
  abline(0,0,lty=2)
  par(new=T)
  scatter.smooth(Year, Resid, lpars=list(col="red",lwd=lwd),ann=F,axes=FALSE)

  if (is.null(lag.max)) lag.max = 10*log10(length(Resid))
  ac.res <- acf(Resid,plot=FALSE,lag.max=lag.max)
  plot(ac.res,main="",lwd=lwd,cex=cex,cex.lab=cex.lab,...)
  title("Autocorrelation (rho vs. lag)",cex.main=cex.main)

  p = c()
  col = c()
  for(i in 1:lag.max){
    LBtest = Box.test(Resid,lag = i,type="Ljung")
    p = c(p,LBtest$p.value)
    col = c(col,isTRUE(LBtest$p.value<0.05)+1)
  }
  plot(p,pch=pch,...,ylim=c(0,max(c(p,0.2))),col=col,xlab="Lag",ylab="P value",cex.lab=cex.lab,...)
  abline(0.05,0.,lty=2,col="blue")
  title("Ljung-Box test",cex.main=cex.main)
  if (output) dev.off()
}

#' 再生産関係の残差ブートストラップをプロットする関数
#'
#' @param boot.res \code{boot.SR}のオブジェクト
#' @param CI プロットする信頼区間
#' @param output pngファイルに出力するか否か
#' @param filename ファイル名
#' @encoding UTF-8
#' @export
bootSR.plot = function(boot.res, CI = 0.8,output = FALSE,filename = "boot",lwd=1.2,pch=1,...) {
  res_base = boot.res$input$Res
  if (class(boot.res$input$Res)=="fit.SR") {
    # for fit.SR
    if (output) png(file = paste0(filename,"_pars.png"), width=10, height=10, res=432, units='in')
    par(pch=pch,lwd = lwd, mfrow=c(2,2))
    jmax = ifelse(boot.res$input$Res$pars$rho==0,3,4)
    for (j in 1:jmax) {
      par0 = c("a","b","sd","rho")[j]

      hist(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$pars[,par0]),xlab=par0,ylab="Frequency",main="",col="gray")
      abline(v=boot.res$input$Res$pars[,par0],col=2,lwd=3)
      abline(v=median(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$pars[,par0])),col=3,lwd=3,lty=2)
      arrows(quantile(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$pars[,par0]),0.5*(1-CI)),0,
             quantile(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$pars[,par0]),0.5*(1-CI)+CI),0,
             col=4,lwd=3,code=3)
      legend("topright",
             legend=c("Estimate","Median","CI(0.8)"),lty=1:2,col=2:4,lwd=2,ncol=1,cex=1)
      if (boot.res$input$method=="d") {
        title(paste0(par0," in Data Bootstrap"))
      } else {
        if (boot.res$input$method=="p") {
          title(paste0(par0," in Parametric Bootstrap"))
        } else {
          title(paste0(par0," in Non-Parametric Bootstrap"))
        }
      }
    }
    if (output) dev.off()

    par(mfrow=c(1,1),pch=pch,lwd = lwd)
    if (output) png(file = paste0(filename,"_SRcurve.png"), width=10, height=7.5, res=432, units='in')
    data_SR = boot.res$input$Res$input$SRdata
    plot(data_SR$R ~ data_SR$SSB, cex=2, type = "p",xlab="SSB",ylab="R",pch=1,
         main="",ylim=c(0,max(data_SR$R)*1.3),xlim=c(0,max(data_SR$SSB)*1.3))
    if (boot.res$input$method=="d") {
      title("Data Bootstrap")
    } else {
      if (boot.res$input$method=="p") {
        title("Parametric Bootstrap for Residuals")
      } else {
        title("Non-Parametric Bootstrap for Residuals")
      }
    }
    points(rev(data_SR$SSB)[1],rev(data_SR$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
    for (i in 1:boot.res$input$n) {
      points(boot.res[[i]]$pred$SSB,boot.res[[i]]$pred$R,type="l",lwd=2,col=rgb(0,0,1,alpha=0.1))
    }
    points(res_base$pred$SSB,res_base$pred$R,col=2,type="l",lwd=3)
    if (output) dev.off()
  } else {
    # fit.SRregime
    regime_unique = boot.res$input$Res$regime_pars$regime
    obs_data = boot.res$input$Res$pred_to_obs
    if (output) png(file = paste0(filename,"_pars.png"), width=15, height=5*nrow(boot.res$input$Res$regime_pars), res=432, units='in')
    par(lwd = lwd, mfrow=c(nrow(boot.res$input$Res$regime_pars),3))
    for (ii in 1:nrow(boot.res$input$Res$regime_pars)) {
      regime = boot.res$input$Res$regime_pars$regime[ii]
      jmax = 3
      for (j in 1:jmax) {
        par0 = c("a","b","sd","rho")[j]
        boot_pars = sapply(1:boot.res$input$n, function(i) as.numeric(boot.res[[i]]$regime_pars[ii,par0]))
        hist(boot_pars,xlab=par0,ylab="Frequency",main="",col="gray")
        abline(v=as.numeric(boot.res$input$Res$regime_pars[ii,par0]),col=2,lwd=3)
        abline(v=median(boot_pars),col=3,lwd=3,lty=2)
        arrows(quantile(boot_pars,0.5*(1-CI)),0,
               quantile(boot_pars,0.5*(1-CI)+CI),0,
               col=4,lwd=3,code=3)
        legend("topright",
               legend=c("Estimate","Median","CI(0.8)"),lty=1:2,col=2:4,lwd=2,ncol=1,cex=1)
        if (boot.res$input$method=="d") {
          title(paste0(par0," of Regime",regime," in Data Bootstrap"))
        } else {
          if (boot.res$input$method=="p") {
            title(paste0(par0," of Regime",regime," in Para Bootstrap"))
          } else {
            title(paste0(par0," of Regime",regime," in Non-Para Bootstrap"))
          }
        }

      }
    }
    if (output) dev.off()

    if (output) png(file = paste0(filename,"_SRcurve.png"), width=7.5*length(regime_unique), height=7.5, res=432, units='in')
    par(mfrow=c(1,length(regime_unique)))
    for(i in 1:length(regime_unique)) {
      use_data = dplyr::filter(obs_data,Regime == regime_unique[i])
      plot(use_data$R ~ use_data$SSB, cex=2, type = "p",pch=1,xlab="SSB",ylab="R",
           main="",ylim=c(0,max(use_data$R)*1.3),xlim=c(0,max(obs_data$SSB)*1.3))
      if (boot.res$input$method=="d") {
        title(paste0("Data Bootstrap for Regime ",regime_unique[i]))
      } else {
        if (boot.res$input$method=="p") {
          title(paste0("Para Bootstrap for Regime ",regime_unique[i]))
        } else {
          title(paste0("Non-para Bootstrap for Regime ",regime_unique[i]))
        }
      }
      for (j in 1:boot.res$input$n) {
        pred_data = boot.res[[j]]$pred %>% filter(Regime == regime_unique[i])
        points(pred_data$SSB,pred_data$R,type="l",lwd=2,col=rgb(0,0,1,alpha=0.1))
      }
      pred_data = res_base$pred %>% filter(Regime == regime_unique[i])
      points(pred_data$SSB,pred_data$R,col=2,type="l",lwd=3)
    }
    if (output) dev.off()
  }
}

#' 再生産関係のジャックナイフ解析
#'
#' 結果のプロットもこの関数で行う
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param is.plot プロットするかどうか
#' @param output pngファイルに出力するか否か
#' @param filename ファイル名
#' @encoding UTF-8
#' @export
jackknife.SR = function(resSR,is.plot=TRUE,output=FALSE,filename = "jackknife",ylim.range = c(0.5,1.5),pch=19,cex=1.1,...) {
  RES = lapply(1:length(resSR$input$SRdata$SSB), function(i){
    jack <- resSR
    jack$input$w[i] <- 0
    jack$input$p0 <- resSR$opt$par
    if (class(resSR)=="fit.SR") {
      do.call(fit.SR,jack$input)
    } else {
      do.call(fit.SRregime,jack$input)
    }
  })
  if (is.plot) {
    jack.res <- RES
    data_SR = resSR$input$SRdata
    if (class(resSR)=="fit.SR") {
      if (output) png(file = paste0(filename,"_pars.png"), width=10, height=10, res=432, units='in')
      par(mfrow=c(2,2),mar=c(3,3,2,2),oma=c(3,3,2,2),pch=pch,cex=cex)
      plot(data_SR$year,sapply(1:length(data_SR$year), function(i) jack.res[[i]]$pars$a),type="b",
           xlab="Year removed",ylab="",main="a in jackknife",ylim=resSR$pars$a*ylim.range)
      abline(resSR$pars$a,0,lwd=2,col=2,lty=2)

      plot(data_SR$year,sapply(1:length(data_SR$year), function(i) jack.res[[i]]$pars$b),type="b",
           xlab="Year removed",ylab="",main="b in jackknife",ylim=resSR$pars$b*ylim.range)
      abline(resSR$pars$b,0,lwd=2,col=2,lty=2)

      plot(data_SR$year,sapply(1:length(data_SR$year), function(i) jack.res[[i]]$pars$sd),type="b",
           xlab="Year removed",ylab="",main="sd in jackknife",ylim=resSR$pars$sd*ylim.range)
      abline(resSR$pars$sd,0,lwd=2,col=2,lty=2)

      if (resSR$input$AR==1) {
        plot(data_SR$year,sapply(1:length(data_SR$year), function(i) jack.res[[i]]$pars$rho),type="b",
             xlab="Year removed",ylab="",main="rho in jackknife",ylim=resSR$pars$rho*ylim.range)
        abline(resSR$pars$rho,0,lwd=2,col=2,lty=2)
      }
      if (output) dev.off()

      if (output) png(file = paste0(filename,"_SRcurve.png"), width=8, height=5, res=432, units='in')
      par(mar=c(3,3,2,2),oma=c(3,3,2,2),mfrow=c(1,1))
      plot(data_SR$R ~ data_SR$SSB, cex=2, type = "p",xlab="SSB",ylab="R",pch=1,
           main="jackknife SR functions",ylim=c(0,max(data_SR$R)*1.3),xlim=c(0,max(data_SR$SSB)*1.3))
      points(rev(data_SR$SSB)[1],rev(data_SR$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
      for (i in 1:length(data_SR$R)) {
        points(jack.res[[i]]$pred$SSB,jack.res[[i]]$pred$R,type="l",lwd=2,col=rgb(0,0,1,alpha=0.1))
      }
      points(resSR$pred$SSB,resSR$pred$R,col=2,type="l",lwd=3,lty=2)
      if (output) dev.off()

    } else { #fit.SRregime
      regime_unique = resSR$regime_resid$regime %>% unique()
      if (output) png(file = paste0(filename,"_pars.png"), width=15, height=7.5, res=432, units='in')
      par(mar=c(3,3,2,2),oma=c(3,3,2,2),mfrow=c(length(regime_unique),3),pch=pch,cex=cex)
      for(j in 1:length(regime_unique)) {
        plot(data_SR$year,sapply(1:length(data_SR$year), function(i) jack.res[[i]]$regime_pars$a[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("a in jackknife (regime ",regime_unique[j],")"),ylim=resSR$regime_pars$a[j]*ylim.range)
        abline(resSR$regime_pars$a[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
        plot(data_SR$year,sapply(1:length(data_SR$year), function(i) jack.res[[i]]$regime_pars$b[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("b in jackknife (regime ",regime_unique[j],")"),ylim=resSR$regime_pars$b[j]*ylim.range)
        abline(resSR$regime_pars$b[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
        plot(data_SR$year,sapply(1:length(data_SR$year), function(i) jack.res[[i]]$regime_pars$sd[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("sd in jackknife (regime ",regime_unique[j],")"),ylim=resSR$regime_pars$sd[j]*ylim.range)
        abline(resSR$regime_pars$sd[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
      }
      if (output) dev.off()
      obs_data = resSR$pred_to_obs
      if (output) png(file = paste0(filename,"_SRcurve.png"), width=12, height=6, res=432, units='in')
      par(mar=c(3,3,2,2),oma=c(3,3,2,2),mfrow=c(1,length(regime_unique)))
      for(i in 1:length(regime_unique)) {
        use_data = dplyr::filter(obs_data,Regime == regime_unique[i])
        plot(use_data$R ~ use_data$SSB, cex=2, type = "p",pch=1,xlab="SSB",ylab="R",
             main=paste0("Jackknife SR for Regime ",regime_unique[i]),ylim=c(0,max(use_data$R)*1.3),xlim=c(0,max(obs_data$SSB)*1.3))
        for (j in 1:nrow(obs_data)) {
          pred_data = jack.res[[j]]$pred %>% filter(Regime == regime_unique[i])
          points(pred_data$SSB,pred_data$R,type="l",lwd=2,col=rgb(0,0,1,alpha=0.1))
        }
        pred_data = resSR$pred %>% filter(Regime == regime_unique[i])
        points(pred_data$SSB,pred_data$R,col=2,type="l",lwd=3,lty=2)
      }
      if (output) dev.off()
    }
  }
  return(invisible(RES))
}

#' プロファイル尤度を計算する関数
#'
#' 結果のプロットもこの関数で行う
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param output pngファイルに出力するか否か
#' @param filename ファイル名
#' @param a_range 推定値に\code{a_range}を乗じた範囲で尤度計算を行う (\code{NULL}のとき\code{c(0.5,2)})
#' @param b_range 推定値に\code{a_range}を乗じた範囲で尤度計算を行う (\code{NULL}のとき\code{c(0.5,2)})
#' @param HS_b_restrict Hockey-StickのときにbをSSBの観測値の範囲にするか否か (\code{b_range}より優先される)
#' @param length 範囲を区切る数
#' @encoding UTF-8
#' @export
prof.likSR = function(resSR,output=FALSE,filename="Profile_Likelihood",a_range = NULL,b_range = NULL,HS_b_restrict = TRUE,length=50) {
  RES = list()
  if (is.null(a_range)) a_range = c(0.5,2)
  if (is.null(b_range)) b_range = c(0.5,2)
  if (class(resSR) == "fit.SR") {
    a.grid <- seq(resSR$pars$a*a_range[1],resSR$pars$a*a_range[2],length=length)
    if (resSR$input$SR!="HS" || !isTRUE(HS_b_restrict)) {
      b.grid <- seq(resSR$pars$b*b_range[1],resSR$pars$b*b_range[2],length=length)
    } else {
      b.grid <- seq(min(resSR$input$SRdata$SSB),max(resSR$input$SRdata$SSB),length=length)
    }
    ba.grid = expand.grid(b=b.grid,a=a.grid)

    if (resSR$pars$rho==0 || resSR$input$out.AR==TRUE) {
      obj.f = function (a,b) resSR$obj.f(a=a,b=b,rho=0)
      prof.lik.res <- exp(-sapply(1:nrow(ba.grid), function(i) obj.f(a=ba.grid[i,2],b=ba.grid[i,1])))
    } else {
      obj.f = function (a,b,x) resSR$obj.f(a=a,b=b,rho=x)
      prof.lik.res <- exp(-sapply(1:nrow(ba.grid), function(i) {
        optimize(function(x) obj.f(a=ba.grid[i,2],b=ba.grid[i,1],x),interval=c(-0.999,0.999))$objective
      }))
    }
    if (output) png(file = paste0(filename,".png"), width=7.5, height=5, res=432, units='in')
    image(b.grid,a.grid,matrix(prof.lik.res,nrow=length),ann=F,col=cm.colors(12),
          ylim=range(a.grid),xlim=range(b.grid))
    par(new=T, xaxs="i",yaxs="i")
    contour(b.grid,a.grid,matrix(prof.lik.res,nrow=length),
            ylim=range(a.grid),xlim=range(b.grid),
            xlab="b",ylab="a",main="Profile Likelihood")
    points(resSR$pars$b,resSR$pars$a,cex=2,pch=4,col=2,lwd=3)
    if (output) dev.off()
    ba.grid.res = ba.grid
  } else { #fit.SRregime
    if (output) png(file = paste0(filename,".png"), width=15, height=5, res=432, units='in')
    par(mfrow=c(1,nrow(resSR$regime_pars)),mar=c(4,4,2,2))
    # par(mfrow=c(1,nrow(resSR$regime_pars)))
    prof.lik.res = NULL
    ba.grid.res = list()
    for (j in 1:nrow(resSR$regime_pars)) {
      a.grid <- seq(resSR$regime_pars$a[j]*a_range[1],resSR$regime_pars$a[j]*a_range[2],length=length)
      if (resSR$input$SR!="HS"|| !isTRUE(HS_b_restrict)) {
        b.grid <- seq(resSR$regime_pars$b[j]*b_range[1],resSR$regime_pars$b[j]*b_range[2],length=length)
      } else {
        if ("b" %in% resSR$input$regime.par) {
          ssb_range=range(resSR$pred_to_obs %>% dplyr::filter(Regime==resSR$regime_pars$regime[j]) %>%
                            dplyr::select(SSB))
          b.grid <- seq(ssb_range[1],ssb_range[2],length=length)
        } else {
          b.grid <- seq(min(resSR$input$SRdata$SSB),max(resSR$input$SRdata$SSB),length=length)
        }
      }
      ba.grid = expand.grid(b=b.grid,a=a.grid)
      ab_order = c("a"=which(resSR$regime_pars$a[j]==resSR$pars$a),"b"=which(resSR$regime_pars$b[j]==resSR$pars$b))
      a_fix = resSR$regime_pars$a[j]
      b_fix = resSR$regime_pars$b[j]
      ab = c(resSR$pars$a,resSR$pars$b)
      x = ab[!(ab %in% c(a_fix,b_fix))]
      obj.f = function(par_a,par_b,x) {
        a = resSR$pars$a
        xa_length=length(a[-ab_order[1]])
        if (xa_length>0) a[-ab_order[1]] <- x[1:xa_length]
        b = resSR$pars$b
        xb_length=length(b[-ab_order[2]])
        if (xb_length>0) b[-ab_order[2]] <- x[xa_length+(1:xb_length)]
        a[ab_order[1]] <- par_a
        b[ab_order[2]] <- par_b
        resSR$obj.f(a=a,b=b)
      }
      if (length(x)==1) {
        prof.lik.res <- cbind(prof.lik.res,exp(-sapply(1:nrow(ba.grid), function(i) {
          # opt = optim(x,obj.f,par_a=ba.grid[i,2],par_b=ba.grid[i,1],lower=x*1.0e-3,upper=x*1.0e+3,method="Brent")
          opt = optim(x,obj.f,par_a=ba.grid[i,2],par_b=ba.grid[i,1],method="BFGS")
          opt$value
        })))
      } else {
        prof.lik.res <- cbind(prof.lik.res,exp(-sapply(1:nrow(ba.grid), function(i) {
          # opt = optim(x,obj.f,par_a=ba.grid[i,2],par_b=ba.grid[i,1],lower=x*0.001,
          #             upper=x*1000,method="L-BFGS-B")
          opt = optim(x,obj.f,par_a=ba.grid[i,2],par_b=ba.grid[i,1],method="BFGS")
          opt$value
        })))
      }
      ba.grid.res[[j]] <- ba.grid

      image(b.grid,a.grid,matrix(prof.lik.res[,j],nrow=length),ann=F,col=cm.colors(12),
            ylim=range(a.grid),xlim=range(b.grid))
      par(new=T, xaxs="i",yaxs="i")
      contour(b.grid,a.grid,matrix(prof.lik.res[,j],nrow=length),
              ylim=range(a.grid),xlim=range(b.grid),
              xlab="b",ylab="a",main=paste0("Profile Likelihood for Regime ",resSR$regime_pars$regime[j]))
      points(resSR$regime_pars$b[j],resSR$regime_pars$a[j],cex=2,pch=4,col=2,lwd=3)
    }
    if (output) dev.off()
  }
  RES$prof.lik <- prof.lik.res
  RES$ba.grid <- ba.grid.res
  return(invisible(RES))
}

#' 再生産関係の推定結果をtxtファイルに出力する関数
#'
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param filename ファイル名('.txt')がつく
#' @encoding UTF-8
#' @export
out.SR = function(resSR,filename = "resSR") {
  RES = list()
  RES$SR = resSR$input$SR
  RES$method = resSR$input$method
  if (class(resSR) == "fit.SR") {
    RES$AR = resSR$input$AR
    RES$out.AR = resSR$input$out.AR
    RES$pars = resSR$pars
  } else {
    RES$regime.year = resSR$input$regime.year
    RES$regime.key = resSR$input$regime.key
    RES$regime.par = resSR$input$regime.par
    RES$pars = resSR$regime_pars
  }
  RES$n = sum(resSR$input$w)
  RES$k = resSR$k
  RES$loglik = resSR$loglik
  RES$AIC = resSR$AIC
  if (!is.null(resSR$AIC.ar)) RES$AIC.ar = resSR$AIC.ar
  RES$AICc = resSR$AICc
  RES$BIC = resSR$BIC
  RES$opt = resSR$opt
  if (class(resSR) == "fit.SR") {
    RES$pred_to_obs = as_tibble(resSR$input$SRdata) %>%
      dplyr::rename(Year = year) %>%
      dplyr::mutate(resid = resSR$resid,resid2 = resSR$resid2) %>%
      dplyr::mutate(Pred_from_SR = R/exp(resid),Pred_from_AR=R/exp(resid2)) %>%
      dplyr::select(Year,SSB,R,Pred_from_SR,resid,Pred_from_AR,resid2)
  } else {
    RES$pred_to_obs = resSR$pred_to_obs
  }
  RES$pred_to_obs =  as.data.frame(RES$pred_to_obs)

  capture.output(RES, file = paste0(filename,".txt"))
}

#' 再生産関係推定が収束しているかや最適解を得られているかを診断する関数
#'
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param n 初期値を変えてパラメータ推定する回数
#' @param sigma 初期値を変えるときの生起乱数の標準偏差
#' @param seed \code{set.seed}で使用するseed
#' @param output テキストファイルに結果を出力するか
#' @param filename ファイル名('.txt')がつく
#' @encoding UTF-8
#' @export
check.SRfit = function(resSR,n=100,sigma=5,seed = 1,output=FALSE,filename="checkSRfit") {
  opt = resSR$opt
  SRdata = resSR$input$SRdata

  RES = list()
  # check convergence
  if (opt$convergence==0) {
    cat(RES$convergence <- "Successful convergence","\n")
  } else {
    message(RES$convergence <- "False convergencen")
  }
  # hessian
  resSR2 = resSR
  if (is.null(resSR$opt$hessian)) {
    input = resSR$input
    input$p0 = resSR$opt$par
    input$hessian = TRUE
    if (class(resSR) == "fit.SR") {
      #input$rep.opt = TRUE
      resSR2 = do.call(fit.SR, input)
    } else {
      resSR2 = do.call(fit.SRregime, input)
    }
  }
  if (all(diag(resSR2$opt$hessian) > 0)) {
    cat(RES$hessian <- "Hessian successfully having positive definite","\n")
  } else {
    message(RES$hessian <- "Hessian NOT having positive definite")
  }

  # check boundary
  RES$boundary <- NULL
  if (class(resSR) == "fit.SR") {
    if (resSR$input$SR == "HS") {
      if (resSR$pars$b > max(SRdata$SSB)*0.99) message(RES$boundary <- "Parameter b reaching the maximum SSB")
      if (resSR$pars$b < min(SRdata$SSB)*1.01) message(RES$boundary <- "Parameter b reaching the minimum SSB")
    } else {
      if (1/resSR$pars$b > 10*max(SRdata$SSB)) message(RES$boundary <- "Proportional recruitment to SSB (no density-dependence)")
      if (1/resSR$pars$b < 0.1*min(SRdata$SSB)) message(RES$boundary <- "Extremely strong density-dependence of recruitment against SSB")
    }
  } else {
    for (i in 1:nrow(resSR$regime_pars)) {
      if ("b" %in% resSR$input$regime.par) {
        SRdata_r = dplyr::filter(resSR$pred_to_obs,Regime==resSR$regime_pars$regime[i])
      } else {
        SRdata_r = SRdata
      }
      if (resSR$input$SR == "HS") {
        if (resSR$regime_pars$b[i] > max(SRdata_r$SSB)*0.99) RES$boundary <- c(RES$boundary,paste0("Parameter b of regime ",resSR$regime_pars$regime[i], " reaching the maximum SSB"))
        if (resSR$regime_pars$b[i] < min(SRdata_r$SSB)*1.01) RES$boundary <- c(RES$boundary,paste0("Parameter b of regime ",resSR$regime_pars$regime[i], " reaching the minimum SSB"))
      } else {
        if (1/resSR$regime_pars$b[i] > 10*max(SRdata_r$SSB)) RES$boundary <- c(RES$boundary,paste0("Proportional recruitment of regime ",resSR$regime_pars$regime[i], " to SSB (no density-dependence)"))
        if (1/resSR$regime_pars$b[i] < 0.1*min(SRdata_r$SSB)) RES$boundary <- c(RES$boundary,paste0("Extremely strong density-dependence of recruitment against SSB in ",resSR$regime_pars$regime[i]))
      }
    }
  }
  if (is.null(RES$boundary)) {
    cat(RES$boundary <- "Parameters not reaching boundaries (successful)","\n")
  } else {
    for (i in 1:length(RES$boundary)) {
      message(RES$boundary[i])
    }
  }

  set.seed(seed)
  loglik = NULL
  pars = NULL
  resSR_list = list()
  for (i in 1:n) {
    input = resSR$input
    for (j in 1:100) {
      input$p0 <- opt$par + rnorm(length(opt$par),0,sigma)
      if (class(resSR) == "fit.SR") {
        #input$rep.opt = TRUE
        resSR2 = try(do.call(fit.SR, input),silent=TRUE)
      } else {
        resSR2 = try(do.call(fit.SRregime, input),silent=TRUE)
      }
      if (class(resSR2) != "try-error") break
    }
    resSR_list[[i]] = resSR2
    loglik = c(loglik, resSR2$loglik)
    pars = rbind(pars,resSR2$opt$par)
  }
  max_loglik = max(loglik)
  optimal = NULL
  if (resSR$loglik-max_loglik < -0.001) {
    message(RES$optim <- "NOT achieving the global optimum")
    diff_loglik = abs(resSR$loglik-max_loglik)
    message(paste0("Maximum difference of log-likelihood is ",round(diff_loglik,6)))
    optimal = resSR_list[[which.max(loglik)]]
    RES$loglik_diff = diff_loglik
  } else {
    cat(RES$optim <- "Successfully achieving the global optimum","\n")
    # global optimumに達している場合のみ
    loglik_diff = purrr::map_dbl(loglik, function(x) abs(diff(c(x,max(loglik)))))
    problem = NULL
    diff_threshold = 1.0e-6
    # a_diff = NULL; b_diff = NULL; sd_diff = NULL; rho_diff = NULL
    for (i in 1:n) {
      if (loglik_diff[i] < diff_threshold) {
        if (all(abs(pars[i,] - resSR$opt$par) < 0.001)) {
          problem = c(problem,FALSE)
          # a_diff = c(a_diff,0); b_diff = c(b_diff,0); sd_diff = c(sd_diff,0)
          # if (class(resSR)=="fit.SR" && resSR$pars$rho != 0) rho_diff = c(rho_diff,NULL)
        } else {
          problem = c(problem,TRUE)
          # a_diff = c(a_diff,max(abs(resSR_list[[i]]$pars$a/resSR$pars$a-1))*100)
          # b_diff = c(b_diff,max(abs(resSR_list[[i]]$pars$b/resSR$pars$b-1))*100)
          # sd_diff = c(sd_diff,max(abs(resSR_list[[i]]$pars$sd/resSR$pars$sd-1))*100)
          # if (class(resSR)=="fit.SR" && resSR$pars$rho != 0) {
          #   rho_diff = c(rho_diff,max(abs(resSR_list[[i]]$pars$rho/resSR$pars$rho-1))*100)
          # }
        }
      } else {
        problem = c(problem,FALSE)
      }
    }
    if (sum(problem)>0) {
      message(RES$pars <- "Different parameter values achieving the global optimum")
      # RES$percent_bias = c("a"=max(a_diff),"b"=max(b_diff),"sd" = max(sd_diff))
      # message("Maximum percent bias of 'a' is ", round(as.numeric(RES$percent_bias["a"]),6),"%")
      # message("Maximum percent bias of 'b' is ", round(as.numeric(RES$percent_bias["b"]),6),"%")
      # message("Maximum percent bias of 'sd' is ", round(as.numeric(RES$percent_bias["sd"]),6),"%")
      # if (class(resSR)=="fit.SR" && resSR$pars$rho != 0) {
      #   RES$percent_bias = c(RES$percent_bias,"rho" = max(rho_diff))
      #   message("Maximum percent bias of 'rho' is ", round(as.numeric(RES$percent_bias["rho"]),6),"%")
      # }
      par_list = t(sapply(1:n, function(i) unlist(resSR_list[[i]]$pars)[unlist(resSR$pars) != 0]))
      par_list = par_list[loglik_diff<diff_threshold,]
      bias_list = t(sapply(1:n, function(i) 100*(unlist(resSR_list[[i]]$pars)[unlist(resSR$pars) != 0]/unlist(resSR$pars)[unlist(resSR$pars)!=0]-1)))
      bias_list = bias_list[loglik_diff<diff_threshold,]
      par_summary = apply(par_list,2,summary)
      percent_bias_summary = apply(bias_list,2,summary)
      RES$par_summary <- par_summary
      RES$percent_bias_summary <- percent_bias_summary
    } else {
      cat(RES$pars <- "Parameters successfully achieving the single solution","\n")
    }
  }
  if (output) {
    capture.output(RES,file=paste0(filename,".txt"))
  }
  if (!is.null(optimal)) RES$optimum = optimal
  # RES$loglik = loglik
  # RES$par_list = par_list
  # RES$percent_bias_list = bias_list
  return(RES)
}

#' 再生産関係のブートストラップの結果をggplotで生成する関数する関数
#'
#' @param boot.res \code{boot.SR}のオブジェクト
#' @param CI プロットする信頼区間
#' @encoding UTF-8
#' @export

bootSR.ggplot = function(boot.res, CI=0.80) {
  if (class(boot.res$input$Res) == "fit.SRregime") {
    stop("This function cannot handle 'fit.SRregime' at present")
  }
  estimate_data = boot.res$input$Res$pred
  obs_data = as_tibble(boot.res$input$Res$input$SRdata)

  for (i in 1:boot.res$input$n) {
    tmp_tbl = boot.res[[i]]$pred %>% mutate(simID = i)
    if (i == 1) {
      boot_pred = tmp_tbl
    } else {
      boot_pred = bind_rows(boot_pred,tmp_tbl)
    }
  }

  summary_boot = boot_pred %>%
    group_by(SSB) %>%
    summarise(median = median(R), lower=quantile(R, probs=0.5*(1-CI))[1], upper=quantile(R, probs=1-0.5*(1-CI))[1])

  g1 = ggplot(data = NULL)+
    geom_ribbon(data = summary_boot, aes(x=SSB,ymin=lower,ymax=upper),alpha=0.3,fill="red")+
    geom_path(data = estimate_data,aes(x=SSB,y=R),size=2) +
    geom_point(data = obs_data,aes(x=SSB,y=R),size=2)+
    geom_path(data = summary_boot, aes(x=SSB,y=median),colour="red",size=2,linetype="dashed")+
    theme_SH()
  g1
}
