#'
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr select

#'
NULL

#' VPAの結果から再生産関係推定用のデータを作成する
#'
#' VPA結果に放流尾数が含まれている場合、天然加入尾数＝加入尾数ー放流尾数として計算する
#'
#' @param vpares VPAの結果のオブジェクト
#' @param R.dat VPAの結果オブジェクトでなく直接データを与えたい場合の加入の値
#' @param SSB.dat VPAの結果オブジェクトでなく直接データを与えたい場合の親魚量の値
#' @param release.dat VPAの結果オブジェクトでなく直接データを与えたい場合の放流尾数の値
#' @param return.df データフレームとして結果を返すか。このオプション関係なくデータフレームとして返すように変更したので、近いうちに廃止予定
#' @param weight.year fit.SRに渡すとき、フィットの対象とする年を指定する。特例として０を与えると、全部の年のデータが使われるようになる。返り値にweightという列が加わる。
#'
#' @encoding UTF-8
#'
#' @export
get.SRdata <- function(vpares=NULL,
                       R.dat = NULL,
                       SSB.dat = NULL,
                       release.dat=NULL,
                       years = NULL,
                       weight.year = NULL,
                       weight.data = NULL,
                       return.df = TRUE){

    is.release <- !is.null(release.dat) | !is.null(vpares$input$dat$release.alive)
    if(is.null(years) && !is.null(vpares)) years <- as.numeric(colnames(vpares$naa))

    # R.datとSSB.datだけが与えられた場合、それを使ってシンプルにフィットする
    if (!is.null(R.dat) & !is.null(SSB.dat)) {
        assertthat::assert_that(length(R.dat)==length(SSB.dat))
        if(!is.null(release.dat)){
            assertthat::assert_that(length(R.dat)==length(release.dat))
            R.dat <- R.dat - release.dat
        }
        if(is.null(years)){
            years <- 1:length(R.dat)
        }
        else{
            assertthat::assert_that(length(R.dat)==length(years))
        }
        dat <- tibble(R = R.dat,SSB = SSB.dat,year = years)
        if(!is.null(release.dat)) R.dat$release <- release.dat
    }

    # vpa dataが与えられた場合にはデータの整形をする
    if(!is.null(vpares)) {
        n <- ncol(vpares$naa)
        L <- as.numeric(rownames(vpares$naa)[1])

        dat      <- list()
        dat$SSB  <- as.numeric(colSums(vpares$ssb,na.rm = TRUE))
        dat$year <- as.numeric(colnames(vpares$ssb))
        dat$R   <-  as.numeric(vpares$naa[1,])

        # 加入年齢分だけずらす
        dat$R    <- dat$R[(L+1):n]
        dat$SSB  <- dat$SSB[1:(n-L)]
        dat$year <- dat$year[(L+1):n]

        dat <- as_tibble(dat)

        if(!is.null(vpares$input$dat$release.dat)) vpares$input$dat$release.alive <- vpares$input$dat$release.dat
        if(!is.null(vpares$input$dat$release.alive)){

          tmpfunc <- function(rdat, data_name){
            if(!is.null(rdat)){
              dat <- rdat[rownames(rdat) %in% as.character(L), ]
              rownames(dat) <- data_name
              dat %>% t() %>% as.data.frame() %>% rownames_to_column(var="year") %>%
                mutate(year=as.numeric(year)) %>% as_tibble()
            }
          }

          dat <- left_join(dat, tmpfunc(vpares$input$dat$release.alive, "release_alive"))
          dat <- dat %>% mutate(R=R-release_alive)
          if(!is.null(vpares$input$dat$release.all)){
            dat <- left_join(dat, tmpfunc(vpares$input$dat$release.all, "release_all"))
            dat <- left_join(dat, tmpfunc(vpares$input$dat$release.ratealive, "release_ratealive"))
          }
        }

        dat$R[dat$R<0] <- 0.001

        # データの抽出
        dat <- dat %>% dplyr::filter(year%in%years)
    }

    assertthat::assert_that(all(dat[["R"]] > 0))

    #if (return.df == TRUE) return(data.frame(year = dat$year,
    #SSB  = dat$SSB,
    #                                         R    = dat$R,
    #                                             release=dat$release))
    #return(dat[c("year","SSB","R","release")])

#    dat.df <- data.frame(year = dat$year,
#                         SSB  = dat$SSB,
#                         R    = dat$R)
#    if(is.release){
#        dat.df$release_alive <- dat$release_alive
#        if(!is.null(dat$release_all)) dat.df$release_all <- dat$release_all
#        if(!is.null(dat$release_ratealive))  dat.df$release_ratealive <- dat$
#    }

  if(!is.null(weight.year)){
    dat <- dat %>% mutate(weight=0)
    dat$weight[dat$year %in% weight.year]  <- 1

    if(length(weight.year)==1 && weight.year==0){
      dat <- dat %>% mutate(weight=1)
    }
  }
  if(is.null(weight.year)){
    dat <- dat %>% mutate(weight=1)
  }

  if(!is.null(weight.data)){
    assertthat::assert_that(all(dat$year==weight.data$year))
    if("weight" %in% names(dat)) dat <- dat %>% select(-weight)
    dat <- dat %>% left_join(weight.data)
  }

  dat <- dat %>% select(year, SSB, R, weight, everything())
  return(dat)

}

#' 再生産関係のパラメタが想定した形になっているかを確認
#'
#' @inheritParams fit.SR
#' @param res_sr fit.SR の結果
validate_sr <- function(SR = NULL, method = NULL, AR = NULL, out.AR = NULL, res_sr = NULL) {
  if (!is.null(SR)) {
    assertthat::assert_that(
      length(SR) == 1,
      SR %in% c("HS", "BH", "RI","Mesnil", "Shepherd", "Cushing","BHS")
    )
  }
  if (!is.null(method)) {
    method = as.character(method)
    assertthat::assert_that(method %in% c("L1", "L2"))
  }
  if (!is.null(AR)) {
    assertthat::assert_that(
      is.numeric(AR),
      AR %in% c(0, 1)
    )
  }
  if (!is.null(out.AR)) {
    assertthat::assert_that(is.logical(out.AR))
  }
  if (!is.null(res_sr)) {
    validate_sr(SR     = res_sr$input$SR,
                method = res_sr$input$method,
                AR     = res_sr$input$AR,
                out.AR = res_sr$input$out.AR)
  }
}

#' 再生産関係の推定
#'
#' 3種類の再生産関係の推定を、最小二乗法か最小絶対値法で、さらに加入の残差の自己相関を考慮して行うことができる
#' @param SRdata \code{get.SRdata}で作成した再生産関係データ
#' @param SR 再生産関係 (\code{"HS"}: Hockey-stick, \code{"BH"}: Beverton-Holt, \code{"RI"}: Ricker, \code{"Mesnil"}: Continuous HS)
#' @param method 最適化法（\code{"L2"}: 最小二乗法, \code{"L1"}: 最小絶対値法）
#' @param AR 自己相関を推定するか(1), しないか(0)
#' @param out.AR 自己相関係数を一度再生産関係を推定したのちに、外部から推定するか（1), 内部で推定するか(0)
#' @param length 初期値を決める際のgridの長さ
#' @param p0 \code{optim}で設定する初期値
#' @param bio_par data.frame(waa=c(100,200,300),maa=c(0,1,1),M=c(0.3,0.3,0.3)) のような生物パラメータをあらわすデータフレーム。waaは資源量を計算するときのweight at age, maaはmaturity at age, Mは自然死亡係数。これを与えると、steepnessやR0も計算して返す
#' @param plus_group hなどを計算するときに、プラスグループを考慮するか
#' @param HS_fix_b SR=HSのとき、折れ点bを固定して計算したい時に値を代入。デフォルットはNULL。ただし、min(SSB)/100より小さい値は入れないこと。
#' @param gamma \code{SR="Mesnil"}のときに使用するsmoothing parameter
#'
#' @encoding UTF-8
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                 out.AR = TRUE)
#' resSR$pars
#'
#' # When outputting steepness, create a bio_par object using derive_biopar function with the res_vpa object and the corresponding year as its argument, and put the bio_par object in the argument of fit.SR.
#' bio_par <- derive_biopar(res_obj=res_vpa,derive_year = 2010)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                 out.AR = TRUE,bio_par=bio_par)
#' resSR$pars
#' resSR$steepness
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
#' \item{\code{steepness}}{bio_parを与えたときに、steepness (h) やR0(漁獲がない場合の平均加入尾数), SB0(漁獲がない場合の平均親魚量)なども追加的に返す}
#' }
#'
#' @export
#'

fit.SR <- function(SRdata,
                   SR="HS",
                   method="L2",
                   AR=1,
                   # TMB=FALSE,
                   hessian=FALSE,w=NULL,
                   length=20,
                   max.ssb.pred=1.3, # 予測値を計算するSSBの最大値（観測された最大値への乗数）
                   p0=NULL,
                   out.AR = TRUE, #自己相関係数rhoを外で推定するか
                   bio_par = NULL,
                   plus_group = TRUE,
                   is_jitter = FALSE,
                   HS_fix_b = NULL,
                   gamma=0.01,
                   bias_correct=FALSE # only for test and L2 option
){
  validate_sr(SR = SR, method = method, AR = AR, out.AR = out.AR)

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  tmp <- check_consistent_w(w, SRdata)
  SRdata <- tmp$SRdata
  w <- arglist$w <- tmp$w


  if (AR==0) out.AR <- FALSE
  rec <- SRdata$R
  ssb <- SRdata$SSB

  N <- length(rec)
  NN <- sum(w) #likelihoodを計算するサンプル数

  #  if (SR=="HS") SRF <- function(x,a,b) a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  if (SR=="Mesnil") SRF <- function(x,a,b) 0.5*a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="Shepherd") SRF <- function(x,a,b) a*x/(1+(b*x)^gamma)
  if (SR=="Cushing") SRF <- function(x,a,b) a*x^b
  if (SR=="BHS") SRF <- function(x,a,b) ifelse(x<b, a*b*(x/b)^{1-(x/b)^gamma}, a*b)

  if (length(SRdata$R) != length(w)) stop("The length of 'w' is not appropriate!")

  if (!is.null(HS_fix_b) && SR!="HS" ) stop("Set SR type as HS if you use HS_fix_b option!")
  if (!is.null(HS_fix_b) && HS_fix_b < (min(ssb)/100) ) stop("Set larger HS_fix_b than min(SSB)/100 if you use HS_fix_b option!")

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
        if(bias_correct==FALSE){
          rss <- w[1]*resid2[1]^2*(1-rho^2)
          for(i in 2:N) rss <- rss + w[i]*resid2[i]^2
          sd <- sqrt(rss/NN)
          sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))
          obj <- -sum(w*dnorm(resid2,0,sd2,log=TRUE))
        }
        else{
#          resid2 <- resid-mean(resid)
          rss <- w[1]*resid2[1]^2*(1-rho^2)
          for(i in 2:N) rss <- rss + w[i]*resid2[i]^2
          sd <-  sqrt(rss/NN)
          sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))          
          obj <- -sum(w*dnorm(resid2,-0.5*sd2^2,sd2,log=TRUE))
        }
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

  if(is.null(HS_fix_b)){
      if (is.null(p0)) {
    a.range <- range(rec/ssb)
    b.range <- range(1/ssb)
    if (SR == "HS" | SR=="Mesnil" | SR=="BHS") b.range <- range(ssb)
    grids <- as.matrix(expand.grid(
      seq(a.range[1],a.range[2],len=length),
      seq(b.range[1],b.range[2],len=length)
    ))
    init <- as.numeric(grids[which.min(sapply(1:nrow(grids),function(i) obj.f(grids[i,1],grids[i,2],0))),])
    init[1] <- log(init[1])
    init[2] <- ifelse (SR == "HS" | SR =="Mesnil" | SR=="BHS",-log(max(0.000001,(max(ssb)-min(ssb))/max(init[2]-min(ssb),0.000001)-1)),log(init[2]))
    if (AR != 0 && !isTRUE(out.AR)) init[3] <- 0
  } else {
    init = p0
  }

  if (SR == "HS" | SR == "Mesnil" | SR=="BHS") {
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
      if(is_jitter==FALSE) par2 <- opt$par
      if(is_jitter==TRUE)  par2 <- opt$par + rnorm(length(opt$par),0,5)
      opt2 <- optim(par2,obj.f2)
      if (abs(opt$value-opt2$value)<1e-6) break
      opt <- opt2
    }
    #}
    opt <- optim(opt$par,obj.f2,method="BFGS",hessian=hessian)

  }
  else{
    if (is.null(p0)) {
      a.range <- c(min(rec/ssb),max(rec)/(min(ssb)/100))
      grids <- as.matrix(seq(a.range[1],a.range[2],len=length))
      init <- as.numeric(grids[which.min(sapply(1:length(grids),function(i) obj.f(grids[i,1],HS_fix_b,0))),])
      init[1] <- log(init[1])
      if (AR != 0 && !isTRUE(out.AR)) init[2] <- 0
    } else {
      init = p0
    }

    if (AR == 0 || out.AR) {
      obj.f2 <- function(x) obj.f(exp(x[1]),HS_fix_b,0)
      opt <- optim(init,obj.f2,method ="Brent",upper = log(max(a.range)), lower = log(min(a.range)))
      #if (rep.opt) {
      #}

    } else {
      obj.f2 <-  function(x) obj.f(exp(x[1]),HS_fix_b,1/(1+exp(-x[2])))

      opt <- optim(init,obj.f2)
      #if (rep.opt) {
      for (i in 1:100) {
        if(is_jitter==FALSE) par2 <- opt$par
        if(is_jitter==TRUE)  par2 <- opt$par + rnorm(length(opt$par),0,5)
        opt2 <- optim(par2,obj.f2)
        if (abs(opt$value-opt2$value)<1e-6) break
        opt <- opt2
      }
      #}
      opt <- optim(opt$par,obj.f2,method="BFGS",hessian=hessian)

    }

  }

  Res <- list()
  Res$input <- arglist
  Res$obj.f <- obj.f
  Res$obj.f2 <- obj.f2
  Res$opt <- opt

  if(is.null(HS_fix_b)){
    a <- exp(opt$par[1])
    b <- ifelse(SR=="HS"|SR=="Mesnil",min(ssb)+(max(ssb)-min(ssb))/(1+exp(-opt$par[2])),exp(opt$par[2]))
    rho <- ifelse(AR==0,0,ifelse(out.AR,0,1/(1+exp(-opt$par[3]))))
  } else{
    a <- exp(opt$par[1])
    b <- HS_fix_b
    rho <- ifelse(AR==0,0,ifelse(out.AR,0,1/(1+exp(-opt$par[2]))))
  }

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

    if(bias_correct==TRUE){
      resid <- resid-mean(resid)
      resid2 <- resid2-mean(resid2)
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
  #Res$gamma <- gamma

  ssb.tmp <- seq(from=0,to=max(ssb)*max.ssb.pred,length=100)
  R.tmp <- sapply(1:length(ssb.tmp), function(i) SRF(ssb.tmp[i],a,b))
  pred.data <- data.frame(SSB=ssb.tmp,R=R.tmp)
  Res$pred <- pred.data

  Res$k <- k <- length(opt$par)+1
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(NN-k-1)
  Res$BIC <- -2*loglik+k*log(NN)

  if(!is.null(bio_par)){
    if(SR!="Mesnil") Res$steepness <- calc_steepness(SR=SR,Res$pars,bio_par$M,bio_par$waa,bio_par$maa,plus_group=plus_group) else{ #add gamma to Res$pars if SR=Mesnil
      pars_and_gamma <- data.frame(Res$pars,gamma)
      Res$steepness <- calc_steepness(SR=SR,pars_and_gamma,bio_par$M,bio_par$waa,bio_par$maa,plus_group=plus_group)
    }
  }

  class(Res) <- "fit.SR"

  return(Res)
}



#' L1とL2のmixtureによる再生産関係の推定（beta版なのでexportしない）
#'
#' 3種類の再生産関係の推定を、最小二乗法か最小絶対値法で、さらに加入の残差の自己相関を考慮して行うことができる
#' @param SRdata \code{get.SRdata}で作成した再生産関係データ
#' @param SR 再生産関係 (\code{"HS"}: Hockey-stick, \code{"BH"}: Beverton-Holt, \code{"RI"}: Ricker)
#' @param alpha alpha:(1-alpha)の比でL1:L2を混ぜる
#' @param length 初期値を決める際のgridの長さ
#' @param rep.opt \code{TRUE}で\code{optim}による最適化を収束するまで繰り返す
#' @param p0 \code{optim}で設定する初期値
#' @encoding UTF-8
# #' @export
#' @noRd
#'
fit.SRalpha <- function(SRdata,
                   SR="HS",
                   alpha=0,
                   hessian=FALSE,
                   w=rep(1,length(SRdata$R)),
                   length=20,
                   max.ssb.pred=1.3, # 予測値を計算するSSBの最大値（観測された最大値への乗数）
                   p0=NULL,
                   rep.opt = TRUE
){

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  # if (AR==0) out.AR <- FALSE
  rec <- SRdata$R
  ssb <- SRdata$SSB

  N <- length(rec)
  NN <- sum(w) #likelihoodを計算するサンプル数

  validate_sr(SR = SR)
  #  if (SR=="HS") SRF <- function(x,a,b) a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)

  obj.f <- function(a,b,out="nll"){ #rhoは無し
    resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))

    # L2 part
    rss_L2 = sum(w*resid^2)
    sd = sqrt(rss_L2/NN)
    obj_L2 = -sum(w*dnorm(resid,0,sd,log=TRUE))

    # L1 part
    rsa_L1 = sum(w*abs(resid))
    phi = rsa_L1/NN
    obj_L1 = -sum(w*(-log(2*phi)-abs(resid/phi)))

    obj = alpha*obj_L2 + (1-alpha)*obj_L1
    SD = sqrt(alpha*sd^2 + (1-alpha)*2*phi^2)
    if (out=="nll") return(obj)
    if (out=="sd") return(SD)
    if (out=="resid") return(resid)
  }

  if (is.null(p0)) {
    a.range <- range(rec/ssb)
    b.range <- range(1/ssb)
    if (SR == "HS") b.range <- range(ssb)
    grids <- as.matrix(expand.grid(
      seq(a.range[1],a.range[2],len=length),
      seq(b.range[1],b.range[2],len=length)
    ))
    init <- as.numeric(grids[which.min(sapply(1:nrow(grids),function(i) obj.f(grids[i,1],grids[i,2]))),])
    init[1] <- log(init[1])
    init[2] <- ifelse (SR == "HS",-log(max(0.000001,(max(ssb)-min(ssb))/max(init[2]-min(ssb),0.000001)-1)),log(init[2]))
    # if (AR != 0 && !isTRUE(out.AR)) init[3] <- 0
  } else {
    init = p0
  }

  if (SR == "HS") {
    obj.f2 <- function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])))
  } else {
    obj.f2 <- function(x) obj.f(exp(x[1]),exp(x[2]))
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
  Res$obj.f <- obj.f
  Res$obj.f2 <- obj.f2
  Res$opt <- opt

  a <- exp(opt$par[1])
  b <- ifelse(SR=="HS",min(ssb)+(max(ssb)-min(ssb))/(1+exp(-opt$par[2])),exp(opt$par[2]))
  # rho <- ifelse(AR==0,0,ifelse(out.AR,0,1/(1+exp(-opt$par[3]))))
  rho <- 0

  sd = obj.f(a=a,b=b,out="sd")
  resid = obj.f(a=a,b=b,out="resid")

  Res$resid <- resid
  Res$resid2 <- resid
  Res$pars <- c(a,b,sd,rho)

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

  class(Res) <- "fit.SRalpha"
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

  validate_sr(SR = SR, method = method, AR = AR)
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
#'
#' 残差のパラメトリックブートストラップはmethod="p"で実行可能で、自己相関なしの場合、fit.SRのpars$sdを標準偏差とする正規分布からランダムに乱数を発生させ、予測値からのずれを加えて加入量のブートストラップデータを生成し、再推定している
#'
#' 自己相関ありの場合、fit.SRのpars$sdを1から自己相関係数rho^2を引き2乗根をとって除したもの(pars$sd/sqrt(1-pars$rho^2))を標準偏差とする正規分布からランダムに乱数（epsilon_t）を発生させ、 毎年の残差(resid_t) を resid_(t+1) = rho × resid_t + epsilon_t とし、resid_tを予測値に加えたブートストラップデータを生成し、再推定している
#'
#' 自己相関ありの場合はノンパラメトリックブートストラップは使わずにパラメトリックブートストラップを用いること
#'
#' 残差のノンパラメトリックブートストラップはmethod="n"で実行可能で、残差の確率分布を仮定せず、残差を重複ありでリサンプリングして、加入量のブートストラップデータを生成する
#'
#' データのブートストラップはmethod="d"で実行可能で、データを重複ありでリサンプリングしたデータを使用して、再生産関係の再推定を行う
#'
#' 親魚量データもリサンプリングにより変化するため、親魚量の不確実性も考慮されることになることから、親魚量データに偏りがあったり、データ数が少なかったり、あるデータ点に推定値が大きく依存している場合はバイアスや不確実性が大きくなりやすい
#'
#' 図のプロットにはbootSR.plotを使用する
#'
#' 自己相関を推定していない場合は最後のrhoの図は表示されない
#'
#' fit.SRの引数にbio_parを入れてスティープネスを計算した場合、SB0、R0、B0、hの図もされる
#'
#' 図の出力bootSR.plotのオプションoutput=Tとすると、各パラメータのヒストグラムが出力される
#'
#' @import purrr
#' @param Res \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @param method 残差ブートストラップ(パラメトリック ("p") かノンパラメトリック ("n")) もしくはデータブートストラップ("d")
#' @param n ブートストラップの回数（例では100回だが500回あるいは1000回を推奨）
#' @encoding UTF-8
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' bio_par <- derive_biopar(res_obj=res_vpa,derive_year = 2010)
#' resL1outer = fit.SR(SRdata, SR = "HS", method = "L1", out.AR = TRUE, AR = 1,bio_par=bio_par)
#'
#' # example if parametric bootstrap
#' boot.res1 = boot.SR(resL1outer, n = 100, method = "p")
#' bootSR.plot(boot.res1)
#'
#' # example if non-parametric bootstrap
#' boot.res2 = boot.SR(resL1outer, n = 100, method = "n")
#' bootSR.plot(boot.res2)
#'
#' #' # example if data bootstrap (optional method)
#' boot.res3 = boot.SR(resL1outer, n = 100, method = "d")
#' bootSR.plot(boot.res3)

#' }
#' @export
#'
boot.SR <- function(Res,method="p",n=100,seed=1){

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  N <- length(Res$input$SRdata$SSB)

  validate_sr(res_sr = Res)
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
            res.b$input$w <- res.b$input$SRdata$weight <-w_r
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
          res.b$input$w <- res.b$input$SRdata$weight <- w_r
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
  validate_sr(SR = SR)
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
#' @param bio_par data.frame(waa=c(100,200,300),maa=c(0,1,1),M=c(0.3,0.3,0.3)) のような生物パラメータをあらわすデータフレーム。waaは資源量を計算するときのweight at age, maaはmaturity at age, Mは自然死亡係数。これを与えると、steepnessやR0も計算して返す
#' @param plus_group hなどを計算するときに、プラスグループを考慮するか
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
#' \item{\code{steepness}}{bio_parを与えたときに、steepness (h) やR0(漁獲がない場合の平均加入尾数), SB0(漁獲がない場合の平均親魚量)なども追加的に返す}
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
  w = NULL,
  max.ssb.pred = 1.3,
  hessian = FALSE,
  bio_par = NULL,
  plus_group = TRUE,
  gamma = 0.001
) {
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  #  if(!is.null(SRdata$weight)) w <- SRdata$w
  #  if(is.null(w)) w <- rep(1,length(SRdata$R))
  tmp <- check_consistent_w(w, SRdata)
  SRdata <- tmp$SRdata
  w <- arglist$w <- tmp$w

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

  validate_sr(SR = SR, method = method)
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  if (SR=="Mesnil") SRF <- function(x,a,b) 0.5*a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="Shepherd") SRF <- function(x,a,b) a*x/(1+(b*x)^gamma)
  if (SR=="Cushing") SRF <- function(x,a,b) a*x^b

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
        summarise(rss = sum(w*resid2), n = sum(w),.groups="drop") %>%
        mutate(sd = sqrt(rss/n))
    } else {
      tbl = tibble(resid=resid,sd_key=sd_key,w=w) %>%
        # mutate(resid2 = resid^2) %>%
        group_by(sd_key) %>%
        summarise(rss = sum(w*abs(resid)), n = sum(w),.groups="drop") %>%
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
    if (SR=="HS" | SR=="Mesnil") {
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
      if (SR=="HS" | SR=="Mesnil") {
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

  if (SR=="HS" | SR=="Mesnil") {
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
  if (SR=="HS" | SR=="Mesnil") {
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

  if(!is.null(bio_par)){
    if(SR!="Mesnil") par.matrix <- Res$regime_pars[c("a","b")]
    else par.matrix <- tibble(Res$regime_pars[c("a","b")],gamma=gamma)
      Res$steepness <- purrr::map_dfr(seq_len(nrow(par.matrix)),
                                function(i){
                                    calc_steepness(SR=SR,rec_pars=par.matrix[i,],M=bio_par$M,waa=bio_par$waa,maa=bio_par$maa,
                                                   plus_group=plus_group) %>%
                                        mutate(regime=Res$regime_pars$regime[i])
                                })
  }

  Res$pars <- NULL
  return(Res)
}

#' 再生産関係の推定における標準化残差を計算する関数
#' @import rmutil
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクト
#' @encoding UTF-8
#' @export
calc.StdResid = function(resSR) {
  validate_sr(res_sr = resSR)
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
    validate_sr(res_sr = resSR)
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
    validate_sr(res_sr = resSR)
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
    validate_sr(res_sr = resSR)
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
#' ブートストラップによる各パラメータのヒストグラムと、それらパラメータをつかった再生産関係が描画される
#'
#' ヒストグラムの描画ではとその中央値、上側・下側CIが（デフォルトのggplt=Tでは加えて平均値が）推定値とともに表示される
#'
#' ヒストグラムのプロットにはggplotをつかうオプションggplt=TRUEをデフォルトにしているが、ggplotを使わない作図もggplt=Fで実行可能
#'　
#' ggplotで描画する場合でRstudio利用時にはplotをZoomで表示にしないと描画されないことがあるので注意
#' @param boot.res \code{boot.SR}のオブジェクト
#' @param CI プロットする信頼区間
#' @param output pngファイルに出力するか否か
#' @param filename ファイル名
#' @encoding UTF-8
#' @export
bootSR.plot = function(boot.res, CI = 0.8,output = FALSE,filename = "boot",lwd=1.2,pch=1,ggplt=TRUE,...) {
  res_base = boot.res$input$Res
  if (class(boot.res$input$Res)=="fit.SR") {
    validate_sr(res_sr = boot.res$input$Res)
    # for fit.SR ----
    # parameter histogram ----
    if(ggplt){ # ggplot (plot histograms of a,b,B0,h,R0,rho,SB0,sd)
      if(is.null(convert_SR_tibble(boot.res[[1]]))) print("Do fit.SR with argument(bio.par) to calculate steepness.")
      N <- boot.res$input$n
      res_boot_par_tibble <- purrr::map_dfr(boot.res[1:N], convert_SR_tibble) %>% dplyr::filter(type=="parameter"&name!="SPR0")
      # AR=0のときrhoのプロットを除く
      if(boot.res$input$Res$input$AR==0) res_boot_par_tibble <- res_boot_par_tibble %>% dplyr::filter(name!="rho")

      res_boot_par_tibble_summary <- res_boot_par_tibble %>%  group_by(name) %>% summarise(median=median(value),mean=mean(value),CI10=quantile(value,0.1),CI90=quantile(value,0.9)) %>% mutate(name_with_CI = str_c(name," (",round(CI10,2),"-",round(CI90,2),")")) %>% pivot_longer(cols=-c(name,name_with_CI), names_to="stats")

      res_boot_par_base <- res_boot_par_tibble_summary %>% filter(stats=="median") #推定結果のために中央値から形だけ持ってくる
      res_fitSRtibble <- convert_SR_tibble(boot.res[["input"]]$Res) %>% dplyr::filter(type=="parameter"&name!="SPR0")


      for(i in 1:length(res_boot_par_base)){
        for(j in 1:length(res_boot_par_base)){
          if(res_boot_par_base$name[i] ==res_fitSRtibble$name[j])
            res_boot_par_base$value[i] <- res_fitSRtibble$value[j]
        }
      }
      res_boot_par_base$stats <- "estimated"

      res_boot_par_tibble <- res_boot_par_tibble_summary %>% select(name, name_with_CI) %>% right_join(res_boot_par_tibble)

      plot_col <- c("blue","blue","red","darkgreen","green")
      base_linetype <- c("solid")
      bootsr_linetype <- c("dashed","dashed","longdashed","solid")
      plot_bootsr_linetype <- rep(bootsr_linetype,length(levels(as.factor(res_boot_par_tibble$name))))
      #legend.labels <- c("estimated", "CI10", "CI90","mean","median")

      if (boot.res$input$method=="d") {
        plot_bootsr_title<- paste0("Data Bootstrap")
      } else {
        if (boot.res$input$method=="p") {
          plot_bootsr_title<- paste0("Parametric Bootstrap")
        } else {
          plot_bootsr_title<- paste0("Non-Parametric Bootstrap")
        }
      }

      boot_par_hist <-ggplot(res_boot_par_tibble) + geom_histogram(aes(x=value)) + facet_wrap(.~fct_inorder(name_with_CI), scale="free")+theme_SH(legend.position="bottom")  + labs(title=plot_bootsr_title)+
        geom_vline(data=res_boot_par_base, mapping = aes(xintercept=value,color=stats),linetype=base_linetype) +
        geom_vline(data=res_boot_par_tibble_summary, mapping = aes(xintercept=value,color=stats),linetype="dashed") +
        scale_color_manual(name="stats",values = plot_col) #+
      #scale_linetype_manual(name="",values = plot_bootsr_linetype) #+ #scale_color_discrete(name="stats",breaks=legend.labels)

      print(boot_par_hist)

      if (output) ggsave(file = paste0(filename,"_pars.png"), plot=boot_par_hist, width=10, height=10,  units='in')

    }
    else { #plot not using ggplot
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

      # steepness histogram (if available) ----
      if(!is.null(boot.res[[1]]$steepness)){
        if (output) png(file = paste0(filename,"_pars_steepness.png"), width=10, height=10, res=432, units='in')
        par(pch=pch,lwd = lwd, mfrow=c(2,2))
        for (j in 1:4) {
          par0 = c("SB0","R0","B0","h")[j]

          hist(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$steepness[,par0]),xlab=par0,ylab="Frequency",main="",col="gray")
          abline(v=boot.res$input$Res$steepness[,par0],col=2,lwd=3)
          abline(v=median(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$steepness[,par0])),col=3,lwd=3,lty=2)
          arrows(quantile(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$steepness[,par0]),0.5*(1-CI)),0,
                 quantile(sapply(1:boot.res$input$n, function(i) boot.res[[i]]$steepness[,par0]),0.5*(1-CI)+CI),0,
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
      }

    }


    # SR curve ----
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
    # fit.SRregime ----

    regime_unique = boot.res$input$Res$regime_pars$regime
    obs_data = boot.res$input$Res$pred_to_obs

    if(ggplt){ # plot using ggplot
      res_base = boot.res$input$Res
      N <- boot.res$input$n
      res_boot_par_tibble <- purrr::map_dfr(boot.res[1:N], convert_SR_tibble) %>% dplyr::filter(type=="parameter"&name!="SPR0")

      legend.labels <- c("estimated", "CI10", "CI90","mean","median")
      plot_col <- c("blue","blue","red","darkgreen","green")
      base_linetype <- c("solid")
      bootsr_linetype <- c("dashed","dashed","longdashed","solid")
      plot_bootsr_linetype <- rep(bootsr_linetype,length(levels(as.factor(res_boot_par_tibble$name))))

      regime.num <- length(levels(as.factor(res_boot_par_tibble$regime)))-1

      for(k in 0:regime.num){
        if (boot.res$input$method=="d") {
          plot_bootsr_title<- paste0("Data Bootstrap regime",k)
        } else {
          if (boot.res$input$method=="p") {
            plot_bootsr_title<- paste0("Parametric Bootstrap regime",k)
          } else {
            plot_bootsr_title<- paste0("Non-Parametric Bootstrap regime",k)
          }
        }


        res_boot_par_tibble_regime<- res_boot_par_tibble %>% filter(regime==k)

        res_boot_par_tibble_summary <- res_boot_par_tibble_regime %>%  group_by(name) %>% summarise(median=median(value),mean=mean(value),CI10=quantile(value,0.1),CI90=quantile(value,0.9)) %>% mutate(name_with_CI = str_c(name," (",round(CI10,2),"-",round(CI90,2),")")) %>% pivot_longer(cols=-c(name,name_with_CI), names_to="stats")

        res_boot_par_base <- res_boot_par_tibble_summary %>% filter(stats=="median")
        res_fitSRtibble <- convert_SR_tibble(res_base) %>% dplyr::filter(type=="parameter"&name!="SPR0") %>% filter(regime==k)
        for(i in 1:nrow(res_boot_par_base)){
          for(j in 1:nrow(res_boot_par_base)){
            if(res_boot_par_base$name[i] == res_fitSRtibble$name[j])
              res_boot_par_base$value[i] <- res_fitSRtibble$value[j]
          }
        }
        res_boot_par_base$stats <- "estimated"

        res_boot_par_tibble_regime <- res_boot_par_tibble_summary %>% select(name, name_with_CI) %>% right_join(res_boot_par_tibble_regime)

        boot_par_hist <- ggplot(res_boot_par_tibble_regime) + geom_histogram(aes(x=value)) + facet_wrap(.~fct_inorder(name_with_CI), scale="free")+theme_SH(legend.position="bottom")  + labs(title=plot_bootsr_title)+
          geom_vline(data=res_boot_par_base, mapping = aes(xintercept=value,color=stats),linetype=base_linetype) +
          geom_vline(data=res_boot_par_tibble_summary, mapping = aes(xintercept=value,color=stats),linetype="dashed") +
          scale_color_manual(name="stats",values = plot_col)

        print(boot_par_hist)

        if (output) ggsave(file = paste0(filename,"_regime",k,"_pars.png"), plot=boot_par_hist, width=10, height=10,  units='in')

      }

    }
    else{ # plot not using ggplot

    # histogram ----
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

    # histogram (steepness) ----
    if(!is.null(boot.res[[1]]$steepness)){
      if (output) png(file = paste0(filename,"_pars_steepness.png"), width=15, height=5*nrow(boot.res$input$Res$regime_pars), res=432, units='in')
      par(lwd = lwd, mfrow=c(nrow(boot.res$input$Res$regime_pars),4))
      for (ii in 1:nrow(boot.res$input$Res$regime_pars)) {
        regime = boot.res$input$Res$regime_pars$regime[ii]
        jmax = 4
        for (j in 1:jmax) {
          par0 = c("SB0","R0","B0","h")[j]
          boot_pars = sapply(1:boot.res$input$n, function(i) as.numeric(boot.res[[i]]$steepness[ii,par0]))
          hist(boot_pars,xlab=par0,ylab="Frequency",main="",col="gray")
          abline(v=as.numeric(boot.res$input$Res$steepness[ii,par0]),col=2,lwd=3)
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
    }
    }

    # SR curve -----
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
#' @param use.p0 初期値に\code{resSR}の結果を使うか否か
#' @param is.plot プロットするかどうか
#' @param output pngファイルに出力するか否か
#' @param filename ファイル名
#' @encoding UTF-8
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                 out.AR = TRUE)
#' res_jackSR <- jackknife.SR(resSR,output = T)
#'
#' # if calculate steepness
#' bio_par <- derive_biopar(res_obj=res_vpa,derive_year = 2010)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                 out.AR = TRUE,bio_par=bio_par)
#' res_jackSR <- jackknife.SR(resSR,output = T)
#' }
#' @export
#'

jackknife.SR = function(resSR,is.plot=TRUE,use.p0 = TRUE, output=FALSE,filename = "jackknife",ylim.range = c(0,1.5),pch=19,cex=1.1,...) {

  used_data <- which(resSR$input$w==1)
  if(is.null(resSR$input$SRdata$weight)) resSR$input$SRdata$weight <- resSR$input$w
  RES = lapply(used_data, function(i){
    jack <- resSR
    jack$input$w[i] <- jack$input$SRdata$weight[i] <- 0
    if (use.p0) jack$input$p0 <- resSR$opt$par
    if (class(resSR)=="fit.SR") {
      validate_sr(res_sr = resSR)
      do.call(fit.SR,jack$input)
    } else {
      if (class(resSR)=="fit.SRregime") {
        do.call(fit.SRregime,jack$input)
      } else {
        do.call(fit.SRalpha,jack$input)
      }
    }
  })
  if (is.plot) {
    jack.res <- RES
    data_SR = resSR$input$SRdata
    years <- data_SR$year[used_data]
    # no regime ----
    if (class(resSR)=="fit.SR" || class(resSR)=="fit.SRalpha") {         # plot parameters ----
      if (output) png(file = paste0(filename,"_pars.png"), width=10, height=10, res=432, units='in')
      par(mfrow=c(2,2),mar=c(3,3,2,2),oma=c(3,3,2,2),pch=pch,cex=cex)
      plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$pars$a),type="b",
           xlab="Year removed",ylab="",main="a in jackknife",ylim=resSR$pars$a*ylim.range)
      abline(resSR$pars$a,0,lwd=2,col=2,lty=2)

      plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$pars$b),type="b",
           xlab="Year removed",ylab="",main="b in jackknife",ylim=resSR$pars$b*ylim.range)
      abline(resSR$pars$b,0,lwd=2,col=2,lty=2)

      plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$pars$sd),type="b",
           xlab="Year removed",ylab="",main="sd in jackknife",ylim=resSR$pars$sd*ylim.range)
      abline(resSR$pars$sd,0,lwd=2,col=2,lty=2)

      if (class(resSR)=="fit.SR") {
        if (resSR$input$AR==1) {
          plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$pars$rho),type="b",
               xlab="Year removed",ylab="",main="rho in jackknife",ylim=resSR$pars$rho*ylim.range)
          abline(resSR$pars$rho,0,lwd=2,col=2,lty=2)
        }
      }
      if (output) dev.off()
      if(!is.null(resSR$steepness)){ #steepness----
        if (output) png(file = paste0(filename,"_steepness.png"), width=10, height=10, res=432, units='in')
        par(mfrow=c(2,2),mar=c(3,3,2,2),oma=c(3,3,2,2),pch=pch,cex=cex)
        plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$steepness$h),type="b",
             xlab="Year removed",ylab="",main="h in jackknife",ylim=resSR$steepness$h*ylim.range)
        abline(resSR$steepness$h,0,lwd=2,col=2,lty=2)

        plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$steepness$B0),type="b",
           xlab="Year removed",ylab="",main="B0 in jackknife",ylim=resSR$steepness$B0*ylim.range)
        abline(resSR$steepness$B0,0,lwd=2,col=2,lty=2)

        plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$steepness$R0),type="b",
         xlab="Year removed",ylab="",main="R0 in jackknife",ylim=resSR$steepness$R0*ylim.range)
        abline(resSR$steepness$R0,0,lwd=2,col=2,lty=2)

        plot(years,sapply(1:length(used_data), function(i) jack.res[[i]]$steepness$SB0),type="b",
           xlab="Year removed",ylab="",main="SB0 in jackknife",ylim=resSR$steepness$SB0*ylim.range)
        abline(resSR$steepness$SB0,0,lwd=2,col=2,lty=2)
        if (output) dev.off()
      }

    # plot SR curve ----
      if (output) png(file = paste0(filename,"_SRcurve.png"), width=8, height=5, res=432, units='in')
      par(mar=c(3,3,2,2),oma=c(3,3,2,2),mfrow=c(1,1))
      plot(data_SR$R ~ data_SR$SSB, cex=2, type = "p",xlab="SSB",ylab="R",pch=1,
           main="jackknife SR functions",ylim=c(0,max(data_SR$R)*1.3),xlim=c(0,max(data_SR$SSB)*1.3))
      points(rev(data_SR$SSB)[1],rev(data_SR$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
      for (i in 1:length(years)) {
        points(jack.res[[i]]$pred$SSB,jack.res[[i]]$pred$R,type="l",lwd=2,col=rgb(0,0,1,alpha=0.1))
      }
      points(resSR$pred$SSB,resSR$pred$R,col=2,type="l",lwd=3,lty=2)
      legend("topright", legend=c("Estimated SR function","A jackknife SR function"),col=c("red",rgb(0,0,1,alpha=0.1)),lty=c(2,1), cex=0.8)
      if (output) dev.off()
    }
    #fit.SRregime ----
    else {
      regime_unique = resSR$regime_resid$regime %>% unique()
      # plot parameters ----
      if (output) png(file = paste0(filename,"_pars.png"), width=15, height=7.5, res=432, units='in')
      par(mar=c(3,3,2,2),oma=c(3,3,2,2),mfrow=c(length(regime_unique),3),pch=pch,cex=cex)
      for(j in 1:length(regime_unique)) {
        plot(years,sapply(1:length(years), function(i) jack.res[[i]]$regime_pars$a[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("a in jackknife (regime ",regime_unique[j],")"),ylim=resSR$regime_pars$a[j]*ylim.range)
        abline(resSR$regime_pars$a[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
        plot(years,sapply(1:length(years), function(i) jack.res[[i]]$regime_pars$b[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("b in jackknife (regime ",regime_unique[j],")"),ylim=resSR$regime_pars$b[j]*ylim.range)
        abline(resSR$regime_pars$b[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
        plot(years,sapply(1:length(years), function(i) jack.res[[i]]$regime_pars$sd[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("sd in jackknife (regime ",regime_unique[j],")"),ylim=resSR$regime_pars$sd[j]*ylim.range)
        abline(resSR$regime_pars$sd[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
      }
      if (output) dev.off()
      if(!is.null(resSR$steepness)){
        if (output) png(file = paste0(filename,"_steepness.png"), width=15, height=7.5, res=432, units='in')
        par(mar=c(3,3,2,2),oma=c(3,3,2,2),mfrow=c(length(regime_unique),4),pch=pch,cex=cex)
        for(j in 1:length(regime_unique)) {
        plot(years,sapply(1:length(years), function(i) jack.res[[i]]$steepness$h[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("h in jackknife (regime ",regime_unique[j],")"),ylim=resSR$steepness$h[j]*ylim.range)
        abline(resSR$steepness$h[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")

        plot(years,sapply(1:length(years), function(i) jack.res[[i]]$steepness$B0[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("B0 in jackknife (regime ",regime_unique[j],")"),ylim=resSR$steepness$B0[j]*ylim.range)
        abline(resSR$steepness$B0[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
        plot(years,sapply(1:length(years), function(i) jack.res[[i]]$steepness$R0[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("R0 in jackknife (regime ",regime_unique[j],")"),ylim=resSR$steepness$R0[j]*ylim.range)
        abline(resSR$steepness$R0[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
        plot(years,sapply(1:length(years), function(i) jack.res[[i]]$steepness$SB0[j]),type="b",
             xlab="Year removed",ylab="",
             main=paste0("SB0 in jackknife (regime ",regime_unique[j],")"),ylim=resSR$steepness$SB0[j]*ylim.range)
        abline(resSR$steepness$SB0[j],0,lwd=2,col=2,lty=2)
        abline(v=resSR$input$regime.year-0.5,lwd=1,lty=3,col="blue")
        }
        if (output) dev.off()
      }

      # plot SRcurve ----
      obs_data = resSR$pred_to_obs
      if (output) png(file = paste0(filename,"_SRcurve.png"), width=12, height=6, res=432, units='in')
      par(mar=c(3,3,2,2),oma=c(3,3,2,2),mfrow=c(1,length(regime_unique)))
      for(i in 1:length(regime_unique)) {
        use_data = dplyr::filter(obs_data,Regime == regime_unique[i])
        plot(use_data$R ~ use_data$SSB, cex=2, type = "p",pch=1,xlab="SSB",ylab="R",
             main=paste0("Jackknife SR for Regime ",regime_unique[i]),ylim=c(0,max(use_data$R)*1.3),xlim=c(0,max(obs_data$SSB)*1.3))
        for (j in 1:length(used_data)) {
          pred_data = jack.res[[j]]$pred %>% filter(Regime == regime_unique[i])
          points(pred_data$SSB,pred_data$R,type="l",lwd=2,col=rgb(0,0,1,alpha=0.1))
        }
        pred_data = resSR$pred %>% filter(Regime == regime_unique[i])
        points(pred_data$SSB,pred_data$R,col=2,type="l",lwd=3,lty=2)
       if(output) legend("topright", legend=c("Estimated SR function","A jackknife SR function"),col=c("red",rgb(0,0,1,alpha=0.1)),lty=c(2,1), cex=0.8)
       else if(i==1) legend("topleft", legend=c("Estimated SR function","A jackknife SR function"),col=c("red",rgb(0,0,1,alpha=0.1)),lty=c(2,1), cex=0.8)
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
    validate_sr(res_sr = resSR)
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
      ab_order = c("a"=which(resSR$regime_pars$a[j]==resSR$regime_pars$a),"b"=which(resSR$regime_pars$b[j]==resSR$regime_pars$b))
      a_fix = resSR$regime_pars$a[j]
      b_fix = resSR$regime_pars$b[j]
      ab = c(resSR$regime_pars$a,resSR$regime_pars$b)
      x = ab[!(ab %in% c(a_fix,b_fix))]
      obj.f = function(par_a,par_b,x) {
        a = resSR$regime_pars$a
        xa_length=length(a[-ab_order[1]])
        if (xa_length>0) a[-ab_order[1]] <- x[1:xa_length]
        b = resSR$regime_pars$b
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
    validate_sr(res_sr = resSR)
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
  if(!is.null(resSR$input$bio_par)) RES$steepness=resSR$steepness
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
#' @param fun_when_check5_replace check5にひっかかったとき、候補となるパラメータの中からどのパラメータを代表としてとってくるか。median, max, minなどの関数を入れる。デフォルト（今までの実装）はmedian。保守的な結果を出すならmaxのほうがよいかも（エクセルではmaxに近い値が出力されるみたい）
#' @encoding UTF-8
#' @export

check.SRfit = function(resSR,n=100,sigma=5,seed = 1,output=FALSE,filename="checkSRfit",
                       fun_when_check5_replace=median) {
  opt = resSR$opt
  SRdata = resSR$input$SRdata
  flag = rep(0,5)

  RES = list()
  # check convergence
  if (opt$convergence==0) {
    cat(RES$convergence <- "1. 収束しています (OK)","\n")
  } else {
    flag[1] <- 1
    stop(RES$convergence <- "1. 収束していないので初期値(引数p0)を変えて計算しなおしてください")
  }
  # hessian
  resSR2 = resSR
  if (is.null(resSR$opt$hessian)) {
    input = resSR$input
    input$p0 = resSR$opt$par
    input$hessian = TRUE
    if (class(resSR) == "fit.SR") {
      validate_sr(res_sr = resSR)
      #input$rep.opt = TRUE
      resSR2 = do.call(fit.SR, input)
    } else {
      if (class(resSR) == "fit.SRregime") {
        resSR2 = do.call(fit.SRregime, input)
      } else {
        input$rep.opt = TRUE
        resSR2 = do.call(fit.SRalpha, input)
      }
    }
  }
  if (all(diag(resSR2$opt$hessian) > 0)) {
    cat(RES$hessian <- "2. Hessian行列の対角成分が正定値になっています (OK)","\n")
  } else {
    flag[2] <- 1
    message(RES$hessian <- "2. Hessian行列の対角成分が正定値になっていません (HSかつ3以降のチェックがOKであれば問題ありません)")
  }

  # check boundary
  RES$boundary <- NULL
  if (class(resSR) == "fit.SR" || class(resSR) == "fit.SRalpha") {
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
    cat(RES$boundary <- "3. どの推定パラメータも壁(boundaries)にあたっていないのでOKです (OK)","\n")
  } else {
    flag[3] <- 1
    for (i in 1:length(RES$boundary)) {
      message("3. パラメータは壁にあたっています(HSで折れ点が過去最小・最大親魚量になっているときにそうなります。HSでない場合は推定の不安定性を示唆します。)", RES$boundary[i])
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
        if (class(resSR) == "fit.SRregime") {
          resSR2 = try(do.call(fit.SRregime, input),silent=TRUE)
        } else {
          input$rep.opt = TRUE
          resSR2 = try(do.call(fit.SRalpha, input),silent=TRUE)
        }
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
    flag[4] <- 1
    diff_loglik = abs(resSR$loglik-max_loglik)
    message(paste0("Maximum difference of log-likelihood is ",round(diff_loglik,6)))
    optimal = resSR_list[[which.max(loglik)]]
    message(RES$optim <- str_c("4. パラメータが大域解に達していません (fit.SRtolで再計算をおこなうか、本関数の返り値のoptimumに結果を置き換えてください。大域解を得るための初期値は,",str_c(optimal$input$p0,collapse=", "),"です)"))
    RES$loglik_diff = diff_loglik
  } else {
    cat(RES$optim <- "4. パラメータが大域解に達しているのでOKです (OK)","\n")
    # global optimumに達している場合のみ
    loglik_diff = purrr::map_dbl(loglik, function(x) abs(diff(c(x,max_loglik))))
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
      flag[5] <- 1
      RES$loglik_diff <- loglik_diff
      message(RES$pars <- str_c("5. 同じ最大尤度(",diff_threshold,"よりも小さい違い)を持つ複数のパラメータが見つかりました（L1かつHSでよく見られます）。"))
      # RES$percent_bias = c("a"=max(a_diff),"b"=max(b_diff),"sd" = max(sd_diff))
      # message("Maximum percent bias of 'a' is ", round(as.numeric(RES$percent_bias["a"]),6),"%")
      # message("Maximum percent bias of 'b' is ", round(as.numeric(RES$percent_bias["b"]),6),"%")
      # message("Maximum percent bias of 'sd' is ", round(as.numeric(RES$percent_bias["sd"]),6),"%")
      # if (class(resSR)=="fit.SR" && resSR$pars$rho != 0) {
      #   RES$percent_bias = c(RES$percent_bias,"rho" = max(rho_diff))
      #   message("Maximum percent bias of 'rho' is ", round(as.numeric(RES$percent_bias["rho"]),6),"%")
      # }
      selected_set <- which(loglik_diff<diff_threshold)
      if(class(resSR)!="fit.SRregime"){
        par_list = t(sapply(1:n, function(i) unlist(resSR_list[[i]]$pars)[unlist(resSR$pars) != 0]))
        par_list0 <- par_list
        par_list = par_list[selected_set,]
        bias_list = t(sapply(1:n, function(i) 100*(unlist(resSR_list[[i]]$pars)[unlist(resSR$pars) != 0]/unlist(resSR$pars)[unlist(resSR$pars)!=0]-1)))
        bias_list = bias_list[selected_set,]
      }else{
        par_list = t(sapply(1:n, function(i) unlist(resSR_list[[i]]$regime_pars)[unlist(resSR$regime_pars) != 0]))
        par_list0 <- par_list
        par_list = par_list[selected_set,]
        bias_list = t(sapply(1:n, function(i) 100*(unlist(resSR_list[[i]]$regime_pars)[unlist(resSR$regime_pars) != 0]/unlist(resSR$regime_pars)[unlist(resSR$regime_pars)!=0]-1)))
        bias_list = bias_list[selected_set,]
      }
      par_summary = apply(par_list,2,summary)
      percent_bias_summary = apply(bias_list,2,summary)
      RES$par_summary <- par_summary
      RES$percent_bias_summary <- percent_bias_summary
      RES$par_list <- par_list
      RES$par_list0 <- par_list0
      # すべてのパラメータのfun_when_check5_replaceに最も近いパラメータセットを持つものを選んでoptimalに入れちゃう
      #      x <- sweep(par_list,2,apply(par_list,2,median),FUN="/") %>% apply(1,mean)
     if(class(resSR)!="fit.SRregime"){
        x <- sweep(par_list[,c("a","b")],2,apply(par_list[,c("a","b")],2,fun_when_check5_replace),FUN="/") %>% apply(1,mean)
     }else{
        tmp <- 2:(1+2*length(unique(resSR$input$regime.key)))
        x <- sweep(par_list[,tmp],2,apply(par_list[,tmp],2,fun_when_check5_replace),FUN="/") %>% apply(1,mean)
      }
      selected <- which.min(abs(x-1))
      cat("ほとんど同じ尤度を持つパラメータの範囲 (",n,"回試行のうち",nrow(par_list),"回分),\n")
      print(apply(par_list,2,summary))
      optimal <- resSR_list[selected_set][[selected]]
      cat(round(as.numeric(unlist(optimal$pars)),4), "をoptimumに出力します(そのときの初期値は",str_c(optimal$input$p0,collapse="-"),"です)\n")
    } else {
      cat(RES$pars <- "5. パラメータが唯一の解として推定されているのでOKです (OK)","\n")
    }
  }
  if (output) {
    capture.output(RES,file=paste0(filename,".txt"))
  }
  if (!is.null(optimal)) RES$optimum = optimal
  # RES$loglik = loglik
  # RES$par_list = par_list
  # RES$percent_bias_list = bias_list
  RES$flag <- flag
  RES$loglik <- loglik
  return(RES)
}

#' 再生産関係推定でSR＝HSにおいて同一尤度を持つ複数のパラメータセットがえられたとき、SR＝Mesnilとしてcheck.SRfitを呼び出し、gammaを変えながら尤度が一意に定まる値を探す関数
#'
#' @param resSR \code{fit.SR}か\code{fit.SRregime}のオブジェクトで、かつSR=Mesnilのもの
#' @param n 初期値を変えてパラメータ推定する回数 (check.SRfitの引数)
#' @param sigma 初期値を変えるときの生起乱数の標準偏差(check.SRfitの引数)
#' @param seed \code{set.seed}で使用するseed(check.SRfitの引数)
#' @param gamma_ini SR=Mesnilの曲がり具合を決めるパラメータ デフォルトは10で、9,8,…,2,1,0.9,0.8,…と複数のパラメータが得られるまで探していく。
#' @encoding UTF-8
#' @export
specify.Mesnil.gamma <- function(resSR,n=100,seed = 1,sigma=5,gamma_ini=10){
  opt = resSR$opt
  SRdata = resSR$input$SRdata
  resSR2 = resSR
  RES = list()

  if (resSR$input$SR!="Mesnil") stop("This function checks convergence of the fitting SR function if SR=Mesnil.")

  identical.likely<-TRUE
  gamma_post <- gamma_ini
  while(identical.likely){
    set.seed(seed)
    loglik = NULL
    pars = NULL
    resSR_list = list()
    for (i in 1:n) {
      input = resSR$input
      input$gamma = gamma_post

      for (j in 1:100) {
        input$p0 <- opt$par + rnorm(length(opt$par),0,sigma)
        if (class(resSR) == "fit.SR") {
          #input$rep.opt = TRUE
          resSR2 = try(do.call(fit.SR, input),silent=TRUE)
        } else {
          if (class(resSR) == "fit.SRregime") {
            resSR2 = try(do.call(fit.SRregime, input),silent=TRUE)
          } else {
            input$rep.opt = TRUE
            resSR2 = try(do.call(fit.SRalpha, input),silent=TRUE)
          }
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
      diff_loglik = abs(resSR$loglik-max_loglik)
      message(paste0("Maximum difference of log-likelihood is ",round(diff_loglik,6)))
      optimal = resSR_list[[which.max(loglik)]]
      RES$loglik_diff = diff_loglik
    } else {
      # global optimumに達している場合のみ
      loglik_diff = purrr::map_dbl(loglik, function(x) abs(diff(c(x,max_loglik))))
      problem = NULL
      diff_threshold = 1.0e-6
      # a_diff = NULL; b_diff = NULL; sd_diff = NULL; rho_diff = NULL
      for (i in 1:n) {
        if (loglik_diff[i] < diff_threshold) {
          if (all(abs(pars[i,] - resSR$opt$par) < 0.001)) {
            problem = c(problem,FALSE)
          } else {
            problem = c(problem,TRUE)
          }
        } else {
          problem = c(problem,FALSE)
        }
      }
      if (sum(problem)>0) {
        RES$loglik_diff <- loglik_diff
        message(RES$pars <- str_c("同じ最大尤度(",diff_threshold,"よりも小さい違い)を持つ複数のパラメータが見つかりました（L1かつHSでよく見られます）。"))
        selected_set <- which(loglik_diff<diff_threshold)
        if(class(resSR)!="fit.SRregime"){
          par_list = t(sapply(1:n, function(i) unlist(resSR_list[[i]]$pars)[unlist(resSR$pars) != 0]))
          par_list0 <- par_list
          par_list = par_list[selected_set,]
          bias_list = t(sapply(1:n, function(i) 100*(unlist(resSR_list[[i]]$pars)[unlist(resSR$pars) != 0]/unlist(resSR$pars)[unlist(resSR$pars)!=0]-1)))
          bias_list = bias_list[selected_set,]
        }else{
          par_list = t(sapply(1:n, function(i) unlist(resSR_list[[i]]$regime_pars)[unlist(resSR$regime_pars) != 0]))
          par_list0 <- par_list
          par_list = par_list[selected_set,]
          bias_list = t(sapply(1:n, function(i) 100*(unlist(resSR_list[[i]]$regime_pars)[unlist(resSR$regime_pars) != 0]/unlist(resSR$regime_pars)[unlist(resSR$regime_pars)!=0]-1)))
          bias_list = bias_list[selected_set,]
        }
        par_summary = apply(par_list,2,summary)
        percent_bias_summary = apply(bias_list,2,summary)
        RES$par_summary <- par_summary
        RES$percent_bias_summary <- percent_bias_summary
        RES$par_list <- par_list
        RES$par_list0 <- par_list0
        # すべてのパラメータのmedianに最も近いパラメータセットを持つものを選んでoptimalに入れちゃう
        #      x <- sweep(par_list,2,apply(par_list,2,median),FUN="/") %>% apply(1,mean)
        if(class(resSR)!="fit.SRregime"){
          x <- sweep(par_list[,c("a","b")],2,apply(par_list[,c("a","b")],2,median),FUN="/") %>% apply(1,mean)
        }else{
          tmp <- 2:(1+2*length(unique(resSR$input$regime.key)))
          x <- sweep(par_list[,tmp],2,apply(par_list[,tmp],2,median),FUN="/") %>% apply(1,mean)
        }
        selected <- which.min(abs(x-1))
        optimal <- resSR_list[selected_set][[selected]]

        identical.likely<-FALSE
      }else{
        gamma_pre <- input$gamma
        degit_n <- floor(log10(gamma_pre))
        ssd <- gamma_pre %/% 10^(degit_n)
        if(ssd!=1){
          gamma_post <- gamma_pre - 10^(degit_n)
        }else{
          gamma_post <- gamma_pre - 10^(degit_n-1)
        }
        input$gamma <- gamma_post
      }
    }

  }

  if(gamma_pre==gamma_ini)stop("Mesnil with gamma_ini has multi-identical likelyhood. Set gamma_ini larger than current value.")
  input$gamma = gamma_pre
  if (class(resSR) == "fit.SR") {
    #input$rep.opt = TRUE
    optimal = try(do.call(fit.SR, input),silent=TRUE)
  } else {
    if (class(resSR) == "fit.SRregime") {
      optimal = try(do.call(fit.SRregime, input),silent=TRUE)
    } else {
      input$rep.opt = TRUE
      optimal = try(do.call(fit.SRalpha, input),silent=TRUE)
    }
  }
  print(paste0("unique SR parameter was given where gamma = ",gamma_pre,"."))
  RES$optimum = optimal
  RES$loglik <- loglik
  RES$gamma <- gamma_pre
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

#' 再生産関係の推定パラメータの相関を出力する関数
#'
#' @inheritParams fit.SR
#' @inheritParams fit.SRregime
#' @param resSR \code{fit.SR}または\code{fit.SRregime}のオブジェクト
#' @return 以下の要素からなるリスト
#' \describe{
#' \item{\code{hessian}}{ヘッセ行列}
#' \item{\code{cov}}{推定されたパラメータの分散共分散行列}
#' \item{\code{cor}}{推定されたパラメータの相関行列}
#' }
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                 out.AR = TRUE)
#' corRes = corSR(resSR)
#' corRes$cor
#' }
#' @encoding UTF-8
#' @export
corSR = function(resSR) {
  if (!resSR$input$hessian) {
    resSR$input$hessian <- TRUE
    resSR$input$p0 = resSR$opt$par
    if (class(resSR) == "fit.SR") resSR = do.call(fit.SR, resSR$input)
    if (class(resSR) == "fit.SRregime") resSR = do.call(fit.SRregime, resSR$input)
  }
  hessian = resSR$opt$hessian
  cov = solve(hessian)
  cor = stats::cov2cor(cov)
  return (list(hessian=hessian,cov=cov,cor=cor))
}

#' Steepness (h) と関連するパラメータ (SB0,R0,B0)を計算する関数
#'
#' @param SR "HS", "BH", "RI"のいずれか
#' @param rec_pars 再生産関係のパラメータで\code{rec_pars$a},\code{rec_pars$b}で使用する
#' @param M 年齢別自然死亡係数 (ベクトルで与えるか、年齢共通の場合\code{M=0.4}のようにしてもよい)
#' @param waa （親魚量の）年齢別体重
#' @param maa 年齢別親魚量
#' @param plus_group 最高齢がプラスグループかどうか
#' @return 以下の要素からなるデータフレーム
#' \describe{
#' \item{\code{SPR0}}{F=0のときのSPR(この逆数がreplacement lineの傾き)}
#' \item{\code{SB0}}{F=0のときの親魚量}
#' \item{\code{R0}}{F=0のときの加入量}
#' \item{\code{B0}}{F=0のときの資源量}
#' \item{\code{h}}{steepness: BHかRIのときは0.2×SB0のときの加入量がh×R0, HSのときはh=1-b/SB0}
#' }
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                out.AR = TRUE)
#' rec_pars = resSR$pars
#' year <- "2017"
#' M = res_vpa$input$dat$M[,year]
#' waa = res_vpa$input$dat$waa[,year]
#' maa = res_vpa$input$dat$maa[,year]
#' Res_h = calc_steepness(SR="HS",rec_pars=rec_pars,M=M,waa=waa,maa=maa,plus_group=TRUE)
#' Res_h
#' }
#' @encoding UTF-8
#' @export
calc_steepness = function(SR="HS",rec_pars,M,waa,maa,plus_group=TRUE,faa = NULL, Pope=TRUE) {
  if (length(M)==1) {
    M = rep(M,length(waa))
  }
  if (length(waa) != length(maa) || length(M) != length(maa)) {
    stop("The lengths of 'waa' and 'maa' must be equal")
  }
  is_MSY <- ifelse(is.null(faa),0,1)
  if (is.null(faa)) {
    faa <- rep(0,length(waa))
  }
  est_detMSY = function(x) {
    NAA0 = 1
    for (i in 1:(length(waa)-1)) {
      NAA0 = c(NAA0,rev(NAA0)[1]*exp(-M[i]-x*faa[i]))
    }
    if (plus_group) {
      NAA0[length(NAA0)] = rev(NAA0)[1]/(1-exp(-1*rev(M)[1]-x*rev(faa)[1]))
    }
    BAA0 = NAA0*waa
    SSB0 = BAA0*maa
    SPR0 = sum(SSB0) #get.SRRと一致 (testに使える)

    # 再生産関係とy=(1/SPR0)*xの交点を求める
    rec_a = rec_pars$a
    rec_b = rec_pars$b
    validate_sr(SR = SR)
    if (SR == "HS") {
      R0 = rec_pars$a*rec_pars$b
      SB0 = R0*SPR0
      if (SB0<rec_b) {
        # warning("Virgin equilibrium does not exist!")
        R0 <- 0
        SB0 <- 0
      }
      h = (SB0-rec_b)/SB0
    }
    if (SR == "Mesnil") {
      gamma <- rec_pars$gamma
      K = sqrt(rec_b^2+gamma^2/4)
      SB0 = (2*K/(SPR0*rec_a/2)-2*rec_b-2*K)/(1/(SPR0*rec_a/2)^2-2/(SPR0*rec_a/2)) #Mesnil and Rochet 2010 ICES JMSより
      R0 = SB0/SPR0
      if (1/SPR0>rec_a) {
        # warning("Virgin equilibrium does not exist!")
        R0 <- 0
        SB0 <- 0
      }
      h = (SB0-rec_b)/SB0 # steepnessの定義はHSと同じでよい？
    }
    if (SR == "BH") {
      SB0 = (rec_a*SPR0-1)/rec_b
      R0 = SB0/SPR0
      h = (rec_a*0.2*SB0/(1+rec_b*0.2*SB0))/R0
    }
    if (SR == "RI") {
      SB0 = (1/rec_b)*log(rec_a*SPR0)
      R0 = SB0/SPR0
      h = (rec_a*0.2*SB0*exp(-rec_b*0.2*SB0))/R0
    }
    B0 = sum(R0*BAA0)
    Res = data.frame(SPR0 = SPR0, SB0 = SB0, R0 = R0, B0 = B0, h = h)

    if(is_MSY==1){
      ypr.spr = ref.F(Fcurrent=x*faa,M=M,waa=waa,waa.catch = waa,maa =maa,
                      Pope=Pope,pSPR=NULL,F.range=NULL,plot=FALSE)
      ypr.spr <- ypr.spr$ypr.spr[1,]
      Yield <- as.numeric(R0*ypr.spr["ypr"])
      Res = cbind(Res,data.frame(Yield=Yield))
    }
    return(Res)
  }

  RES = est_detMSY(0)
  if (is_MSY==1) RES = RES %>% dplyr::select(-Yield)
  if(is_MSY==1){
    obj_fun = function(x) {
      RES = est_detMSY(x)
      return(-RES$Yield)
    }
    x_grid = seq(0,10,length=101)
    tmp = x_grid %>% purrr::map_dbl(.,obj_fun)
    Opt = optimize(obj_fun,x_grid[c(max(1,which.min(tmp)-1),min(which.min(tmp)+1,length(x_grid)))])
    Fmsy2F=Opt$minimum
    RES2 = est_detMSY(Fmsy2F)
    RES2 = RES2 %>% dplyr::select(-h) %>% dplyr::mutate(Fmsy2F=Fmsy2F)
    colnames(RES2) <- c("SPRmsy","SBmsy","Rmsy","Bmsy","MSY","Fmsy2F")
    RES = cbind(RES,RES2)
  }
  return(RES)
}


#' 再生産関係の推定結果からbootstrapを実行し、hの分布を計算する関数
#'
#' @param res_SR fit.SRまたはfit.SRregime、モデル平均の結果
#' @param M 年齢別自然死亡係数 (ベクトルで与えるか、年齢共通の場合\code{M=0.4}のようにしてもよい)
#' @param waa （親魚量の）年齢別体重
#' @param maa 年齢別親魚量
#' @param plus_group 最高齢がプラスグループかどうか
#'
#' @export
#'


boot_steepness <- function(res_SR, M, waa, maa, n=100, plus_group=TRUE){

    is.model.average <- class(res_SR)!="fit.SRregime" && class(res_SR)!="fit.SR"

    if(is.model.average){ # model average case (calculate recursively)
      res_steepness <- purrr::map_dfr(res_SR, boot_steepness, n=n, M=M, waa=waa, maa=maa, plus_group=plus_group, .id="id")
    }else{
      res_boot <- boot.SR(res_SR, n=n)
      if(class(res_boot[[1]])=="fit.SRregime"){ # regime shift
          res_steepness <- purrr::map_dfr(res_boot[1:n], function(x){
              par.matrix <- x$regime_pars[c("a","b")]
              tmplist <- purrr::map_dfr(seq_len(nrow(par.matrix)),
                                    function(i){
                                        calc_steepness(SR=res_SR$input$SR,rec_pars=par.matrix[i,],M=M,waa=waa,maa=maa,
                                                       plus_group=plus_group)
                                    },.id="id")
              tmplist <- bind_cols(tmplist,par.matrix)
          })
      }
      if(class(res_boot[[1]])=="fit.SR"){ # normal
          res_steepness <- purrr::map_dfr(res_boot[1:n], function(x){
              calc_steepness(SR=res_SR$input$SR,rec_pars=x$pars,M=M,waa=waa,maa=maa,plus_group=plus_group)
          })
          res_steepness$id <- 1
          res_steepness <- res_steepness %>%
              bind_cols(purrr::map_dfr(res_boot[1:n], function(x) x$pars))
      }}

    res_steepness
}

#' fit.SRとfit.SRregimeの重み付け（どのデータを使うか）の設定方法は、引数wで与える場合とSRdata$weightで与える場合の２パターンある。両者が矛盾なく設定されているかを確認するための関数

check_consistent_w <- function(w, SRdata){
    if(is.null(SRdata$weight)  && is.null(w)) SRdata$weight <- w <- rep(1,length(SRdata$R))
    if(is.null(SRdata$weight)  && !is.null(w)) SRdata$weight <- w
    if(!is.null(SRdata$weight) &&  is.null(w)) w <- SRdata$weight
    if(!is.null(SRdata$weight) && !is.null(w) && sum(unlist(SRdata$weight)!=unlist(w))>0) stop("SRdata$weight と引数wの両方に重み付けの指定がなされていて、両者が違います")
    lst(w,SRdata)
}


#'
#' 再生産関係を網羅的にフィットする
#'
#' regimeなしのものだけに対応
#'
#' @export
#'

tryall_SR <- function(data_SR, plus_group=TRUE, bio_par=NULL, tol=FALSE, detail=FALSE){

  SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), L.type = c("L1", "L2")) %>%
    as_tibble()

  if(tol==TRUE) fit.SR <- fit.SR_tol
  allres <- purrr::map2(SRmodel.list$SR.rel, SRmodel.list$L.type,
                                   function(x,y){
                                     tmp0 <- fit.SR(data_SR, SR = x, method = y,
                                                    AR = 0, hessian = FALSE, out.AR=TRUE,
                                                    bio_par = bio_par, plus_group = plus_group)
                                     res1 <- unlist(tmp0[c("pars","AICc","steepness")])
                                     tmp <- fit.SR(data_SR, SR = x, method = y,
                                                   AR = 1, hessian = FALSE, out.AR=TRUE,
                                                   bio_par = bio_par, plus_group = plus_group)
                                     res2 <- unlist(tmp[c("pars","AICc","steepness")])
                                     res2 <- c(res2,"deltaAIC(AIC_AR-AIC_noAR)"=tmp$AIC.ar[2]-tmp$AIC.ar[1])
                                     tmp2 <- fit.SR(data_SR, SR = x, method = y,
                                                    AR = 1, hessian = FALSE, out.AR=FALSE,
                                                    bio_par = bio_par, plus_group = plus_group)
                                     res3 <- unlist(tmp2[c("pars","AICc","steepness")])
                                     list(pars=bind_rows(res1,res2,res3,.id="id"),
                                          model=list(tmp0, tmp, tmp2))
                                   })
  SRmodel.list$pars <- purrr::map(allres, function(x) x$pars)
  result.list <- SRmodel.list %>%
    unnest(col=pars) %>%
    left_join(tibble(id=as.character(1:3),AR.type=c("non","outer","inner"))) 
  
  if(detail==TRUE){
   SRmodel.list <- SRmodel.list %>% mutate(model=purrr::map(allres, function(x) x$model))  %>%
    unnest(col=model) 
   result.list$model <- SRmodel.list$model
  }
  
  return(arrange(result.list, AICc, AR.type))

}

#'
#' 再生産関係のフィット（頑健版）
#'
#' fit.SRのwrapper。与えられたinputはそのままfit.SRに渡すが、計算したあとにcheck.SRfitを回し、大域解に達していない場合や同じ尤度を持つ複数の大域解がある場合、その中央値をとるres_SRを返す
#'
#' @param ... fit.SRへの引数
#' @param n_check check.SRfitに渡す引数n（初期値を変えて計算を繰り返す回数）
#' @param seed check.SRfitに渡す引数seed。初期値を変えるときに使う乱数。toolのスクリプトではずっと12345を利用していた関係で、デフォルトの引数は12345にしている。
#'
#' @details fit.SR_tolの動作は以下の通りです。1) fit.SRを実施 2) fit.SRの結果(便宜的にres0とする)にcheck.SRfitをあてはめる。check.SRfitでは以下のことが実施されます。2-1) 初期値を変えてn回パラメータ推定を繰り返す。2-2) res0から得られている対数尤度と、n回パラメータ推定を繰り返したときに得られた対数尤度を比較し、res0の対数尤度よりも大きいものが見つかったらres0は「4. 大域解ではない」という判定をする。また、帰り値$optimumの中に、最大の対数尤度を得たときの推定結果を返す。2-3) 4.がOKの場合、n回分の対数尤度の中で最大値との差が1-06以下だけど、パラメータの推定値が異なる（0.001以上）場合には、「5. 同じ対数尤度を持つ複数のパラメータが見つかりました」というメッセージが出る。この場合、対数尤度の差が1-06以下を示すパラメータの範囲（最小、最大、平均など）が出力される。また、中央値に最も近いパラメータのセットを$optimumに返す。3) check.SRfitの返り値を見て、4のフラグが立っている場合、res0をcheck.SRfitの返り値$optimumに置き換える。そうなった場合、もう一度check.SRfitを実行する（これによって4のフラグは立たなくなることが想定）。4) check.SRfitの返り値を見て、5のフラグが立っている場合、res0をcheck.SRfitの返り値$optimumに置き換える。5) res0を返す。
#'
#' @export
#'

fit.SR_tol <- function(...,n_check=100,is_regime=FALSE,seed=12345,fun_when_check5_replace=median){

  if(is_regime==FALSE)  res_SR <- fit.SR(...)
  if(is_regime==TRUE )  res_SR <- fit.SRregime(...)

  cat("...check.SRfit 1 回目....\n")
  check <- check.SRfit(res_SR, output=FALSE, n=n_check, fun_when_check5_replace=fun_when_check5_replace)
  # まず大域解に達しているかどうかのflagを確認し、flagが立っている場合には$optimと置き換える
  if(check$flag[4]==1) {
    res_SR <- check$optimum
    cat("大域解を得るための初期値に変えたres_SRの結果に置き換えます\n")

    # 再度check.SRfit
    cat("...check.SRfit 2回目....\n")
    check <- check.SRfit(res_SR,output=FALSE,n=n_check, fun_when_check5_replace=fun_when_check5_replace)
  }

  # 5番目のflag(尤度が同じパラメータの範囲が存在する)が立っているかどうかも確認し、flagが立っていればoptimumと置き換える
  if(check$flag[5]==1){
    res_SR <- check$optimum
    cat("尤度が同じパラメータの範囲で中央値をとる結果に置き換えます\n")
  }
  return(res_SR)
}


#'
#' fit.SRのwrapper関数でhに制約を置いたもの
#'
#' 現状は推定されたパラメータのみが出力される
#'
#' @export
#'
#'

fit.SR_pen <- function(bio_par, h_upper=Inf, h_lower=0.2, plus_group=TRUE, ...){
    res1 = fit.SR(...)
#    ri1$pars
#    ri1$steepness
#    h_upper = 1
    init = res1$opt$par

    obj_pen = function(x,out=FALSE) {
        a = exp(x[1]); b= exp(x[2]);
        h0 = calc_steepness(SR="RI",rec_pars=data.frame("a"=a,"b"=b),waa=bio_par$waa,maa=bio_par$maa,M=bio_par$M,plus_group=plus_group)
        h <- h0[1,"h"]
        if(out==FALSE) {
            res1$obj.f2(x)+1000*max(0,h-h_upper) + 1000*max(0,h_lower-h)
        } else{
            bind_cols(data.frame("a"=a,"b"=b),h0)
        }
    }

    opt = optim(init,obj_pen)
    obj_pen(opt$par,out=TRUE)

}


#'
#' 再生産関係のレジームシフトを隠れマルコフモデルで推定する
#'
#' Tang et al. (2021) ICESJMSの論文を改変
#'
#' @inheritParams fit.SRregime
#' @params k_regime 推定するレジームの数（2以上の整数）
#' @params b_range パラメータbの範囲
#' @params p0 パラメータの初期値（リスト）
#' @params overwrite cppファイルのコンパイルを上書きして行うか
#' @encoding UTF-8
#' @export
#'
#
hmm_SR = function(SRdata,SR="BH",k_regime=2,gamma=0.01,b_range=NULL,p0=NULL,overwrite=FALSE,max.ssb.pred = 1.3) {
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  k_regime<-round(k_regime)

  if(k_regime<2) stop("Incorrect 'k_regime'")

  data = SRdata
  st = data$SSB
  # st = st/(10^4)
  yt = log(data$R/data$SSB)

  SRcode = case_when(SR=="RI" ~ 1,SR=="BH" ~ 2,SR=="HS"~3,SR=="Mesnil"~4,TRUE~5)
  if (SRcode==5) stop("SR not recognized")
  if (is.null(b_range)) {
    if (SRcode<3) { # Ricker or BH
      b_range = range(1/st)
    } else {
      b_range = range(st)
    }
  }

  tmb_data = list(st=st,yt=yt,
                  alpha_u=max(yt),alpha_l=min(yt),
                  beta_u=max(b_range),beta_l=min(b_range),
                  sigma_u=sd(yt),SRcode=SRcode,gamma=gamma)

  if (!is.null(p0)) {
    parameters = p0
  } else {
    parameters = list(
      lalpha = -log(k_regime+1-1:k_regime),
      lbeta = rep(0,k_regime),
      lsigma = rep(0,k_regime),
      pi1_tran = rep(0,k_regime-1),
      qij_tran = matrix(0,nrow=k_regime,ncol=k_regime-1)
    )
  }

  if (overwrite==TRUE){
    use_rvpa_tmb("HMM_SR",overwrite=TRUE)
  }
  if (!file.exists("HMM_SR.dll")) {
    use_rvpa_tmb("HMM_SR")
  }

  obj = TMB::MakeADFun(tmb_data,parameters,DLL="HMM_SR",inner.control=list(maxit=50000,trace=F),silent=TRUE)
  if (length(obj$par)>length(st)) {
    stop("NOT estimable because k > n (k: parameter number, n: sample size")
  }
  opt = nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,iter.max=10000,eval.max=10000,sing.tol=1e-20))

  Res = list()

  Res$tmb_data = tmb_data
  Res$p0 = parameters
  Res$obj = obj
  Res$opt = opt

  rep = obj$report(opt$par)
  Res$rep = rep

  alpha = rep[["alpha"]]
  a = exp(alpha)
  b = rep[["beta"]]
  sigma = rep[["sigma"]]
  pars=data.frame(regime=1:length(a),a=a,b=b,sd=sigma)
  Res$pars=pars
  regime_prob = (rep[["r_pred"]])
  colnames(regime_prob) = data$year
  Res$regime_prob=round(t(regime_prob),3)
  regime=Rfast::colMaxs(regime_prob)
  names(regime) = data$year
  Res$regime = regime
  Res$trans_prob = rep[["qij"]]

  Res$loglik <- loglik <- -opt$objective
  Res$k <- k <- length(opt$par)
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(length(data$year)-k-1)
  Res$BIC <- -2*loglik+k*log(length(data$year))
  Res$input <- arglist

  # calculate residuals and predicted values
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)
  SRFV = Vectorize(SRF,vectorize.args = "x")

  pred = sapply(1:nrow(pars), function(i) {
    regime_prob[i,] *log(SRFV(SRdata$SSB,pars$a[i],pars$b[i]))
  }) %>% rowSums %>% exp()

  pred_to_obs = SRdata %>%
    dplyr::rename(Year=year) %>%
    dplyr::mutate(Regime=factor(regime),
                  Pred=pred) %>%
    dplyr::mutate(resid=log(R/pred))

  Res$pred_to_obs <- pred_to_obs
  Res$resid <- pred_to_obs$resid

  ssb.tmp <- seq(from=0,to=max(SRdata$SSB)*max.ssb.pred,length=100)
  pred2 = sapply(1:nrow(pars), function(i) SRFV(ssb.tmp,pars$a[i],pars$b[i])) %>% as.data.frame
  colnames(pred2) <- pars$regime
  pred2 = pred2 %>% mutate(SSB = ssb.tmp) %>%
    pivot_longer(., cols=-SSB,values_to="R") %>%
    mutate(Regime=factor(name)) %>%
    arrange(Regime,SSB) %>%
    dplyr::select(Regime,SSB,R)

  Res$pred <- pred2

  class(Res) <- "hmm_SR"

  return( Res )
}
