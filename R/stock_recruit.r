#' VPAの結果から再生産関係推定用のデータを作成する
#'
#' @param vpares VPAの結果のオブジェクト
#' 
#'
#' @export

get.SRdata <- function(vpares,R.dat=NULL,SSB.dat=NULL,years=as.numeric(colnames(vpares$naa))){
    # R.datとSSB.datだけが与えられた場合、それを使ってシンプルにフィットする
    if(!is.null(R.dat) & !is.null(SSB.dat)){
        dat <- data.frame(R=R.dat,SSB=SSB.dat,year=1:length(R.dat))
    }
    else{
    # データの整形
        n <- ncol(vpares$naa)
        L <- as.numeric(rownames(vpares$naa)[1])

        dat <- list()
        dat$R <- as.numeric(vpares$naa[1,])
        dat$SSB <- as.numeric(colSums(vpares$ssb,na.rm=TRUE))
        dat$year <- as.numeric(colnames(vpares$ssb))
    # 加入年齢分だけずらす
        dat$R <- dat$R[(L+1):n]
        dat$SSB <- dat$SSB[1:(n-L)]
        dat$year <- dat$year[(L+1):n]

                                        # データの抽出
        dat <- as.data.frame(dat)
        dat <- dat[dat$year%in%years,]
    }

    class(dat) <- "SRdata"
    return(dat[c("year","SSB","R")])
}


#' 加入の残差の自己相関を考慮した再生産関係の推定
#' L1ノルム（最小絶対値）も推定できる (sigmaはSD)
#' TMB = TRUEでmarginal likelihood (.cppファイルが必要)
#'
#' @param SRdata get.SRdataで作成した再生産関係データ
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
                   out.AR = FALSE #自己相関係数rhoを外で推定するか
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
      arres <- ar(resid,aic=FALSE,order.max=1)
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
#' 
#' @export
#' 

prof.lik <- function(Res,a=Res$pars$a,b=Res$pars$b,sd=Res$pars$sd,rho=Res$pars$rho) {
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

