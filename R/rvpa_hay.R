#' csvデータを読み込んでvpa.hay用のデータを作成する関数
#'
#' @details `data.handler`の`vpa.hay`関数バージョン。
#'
#' @encoding UTF-8
#'
#' @export

data.handler.hay <- function(
    caa1,
    caa2,
    waa1,
    waa2,
    maa,
    maa_ssb=NULL,
    waa_ssb=NULL,
    M = 0.295,
    omega1=NULL,
    omega2=NULL,
    index=NULL
)
{
  years <- as.numeric(sapply(strsplit(names(caa1[1,]),"X"), function(x) x[2]))

  if (is.null(dim(waa1)) | dim(waa1)[2]==1) waa1 <- as.data.frame(matrix(unlist(waa1), nrow=nrow(caa1), ncol=ncol(caa1)))
  if (is.null(dim(waa2)) | dim(waa2)[2]==1) waa2 <- as.data.frame(matrix(unlist(waa2), nrow=nrow(caa2), ncol=ncol(caa2)))
  if (is.null(dim(maa)) | dim(maa)[2]==1) maa <- as.data.frame(matrix(unlist(maa), nrow=nrow(caa1), ncol=ncol(caa1)))

  if (is.null(dim(M))) M <- as.data.frame(matrix(M, nrow=nrow(caa1), ncol=ncol(caa1)))

  if (is.null(maa_ssb)) maa_ssb <- maa else {
    maa_ssb <- as.data.frame(matrix(unlist(maa_ssb), nrow=nrow(caa1), ncol=ncol(caa1)))
  }

  if (is.null(waa_ssb)) waa_ssb <- waa2 else {
    waa_ssb <- as.data.frame(matrix(unlist(waa_ssb), nrow=nrow(caa1), ncol=ncol(caa1)))
  }

  colnames(caa1) <- colnames(waa1) <- colnames(caa2) <- colnames(waa2) <- colnames(maa) <- colnames(maa_ssb) <- colnames(waa_ssb) <- colnames(M) <- years
  rownames(maa_ssb) <- rownames(waa_ssb) <- rownames(M) <- rownames(caa1)

  if (!is.null(index)) colnames(index) <- years
  if (!is.null(omega1)) {rownames(omega1) <- rownames(caa1); colnames(omega1) <- years}
  if (!is.null(omega2)) {rownames(omega2) <- rownames(caa1); colnames(omega2) <- years}

  res <- list(caa1=caa1, caa2=caa2, waa1=waa1, waa2=waa2, maa=maa, maa_ssb=maa_ssb, waa_ssb=waa_ssb, M=M, omega1=omega1, omega2=omega2, index=index)

  invisible(res)
}

sel.func <- function(faa, def="maxage") {
  if(def=="maxage") saa <- apply(faa, 2, function(x) x/x[length(x[!is.na(x)])])
  if(def=="max") saa <- apply(faa, 2, function(x) x/max(x,na.rm=TRUE))
  if(def=="mean") saa <- apply(faa, 2, function(x) x/sum(x,na.rm=TRUE))

  return(saa)
}

zc_logit <- function(x) -log(2/(1+x)-1)
zc_expit <- function(x) 2/(1+exp(-x))-1


#' 半期VPAで資源量推定を行う関数
#'
#' @details ホッケ道北系群の資源評価を種に使われているコード。\code{frasyr::vpa}と対応していないコードもあったりするが、利用する人が適宜開発しながら、関数をアップデートしていくといいと思います。
#'
#' @encoding UTF-8
#'
#' @export

vpa.hay <- function(
    dat,  # data for vpa
    con_tf=3,  # constraint of terminal F
    con_fc=3,  # Fcurrent years from the last
    rec=NULL, # 最新年の加入を外から
    alpha=1,   # F_{a-1} = alpha*F_a
    tune = FALSE,  # tuningをするかどうか
    abund = "B",   # tuningの際，何の指標に対応するか
    min.age = 0,  # tuning指標の年齢参照範囲の下限
    max.age = 0,  # tuning指標の年齢参照範囲の上限
    ki = 1,
    link = "id",  # tuningのlink関数
    base = NA,  # link関数が"log"のとき，底を何にするか
    af = NA,  # 資源量指数が年の中央のとき，af=0なら漁期前，af=1なら漁期真ん中，af=2なら漁期後となる
    p.m = 0.5,  # Popeの近似式でどこで漁獲が起こるか（0.5は年の真ん中）
    stat.tf="mean",
    index.w = NULL,  # tuning indexの重み
    use.index = "all",
    scale=1000,
    sel.def = "max",
    b_est = TRUE,
    b_fix = 1,
    lambda=0,
    no.est=FALSE,
    sel.const = FALSE,
    W_s=1000,
    start_row=1,
    ssb_type=1,
    p_ssb=0.25,
    p_slide=1,
    est.method="ml",
    optimizer="nlm",
    Lower=-Inf,
    Upper=Inf,
    hessian=FALSE,
    dd=10^(-10),
    wf=5,
    AR=FALSE,
    Lo_AR=c(-Inf,-Inf,0,-Inf),
    Up_AR=c(Inf,Inf,Inf,Inf),
    do_compile=FALSE,
    rho_init=0,
    rho_est=TRUE,
    waa_new=NULL,
    rho_se=TRUE,
    cov_q=FALSE,
    X_q=NULL,
    rec_cor=NULL,
    eta=1,
    p.init=0.5
){


  # inputデータをリスト化

  argname <- ls()  # 関数が呼び出されたばかりのときのls()は引数のみが入っている
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  #

  if (do_compile){
    require(TMB)
    compile("ar1.cpp")
    dyn.load(dynlib("ar1"))
  }

  #

  caa1 <- dat$caa1
  caa2 <- dat$caa2
  maa <- dat$maa
  waa1 <- dat$waa1
  waa2 <- dat$waa2
  M <- dat$M
  index <- dat$index
  omega1 <- dat$omega1
  omega2 <- dat$omega2
  maa_ssb <- dat$maa_ssb
  waa_ssb <- dat$waa_ssb

  min.age <- min.age + 1
  max.age <- max.age + 1

  caa <- waa <- M2 <- NULL
  for (i in 1:5) {
    caa <- rbind(caa, rbind(caa1[i,], caa2[i,]))
    waa <- rbind(waa, rbind(waa1[i,], waa2[i,]))
    M2 <- rbind(M2, rbind(M[i,]/2, M[i,]/2))
  }

  M <- M2

  if (length(abund) > 1 & length(b_est)==1) {
    b_est <- rep(b_est, length(abund))
    b_fix <- rep(b_fix, length(abund))
  }
  if (length(abund) > 1 & length(AR)==1) {
    AR <- rep(AR, length(abund))
  }
  if (any(AR)){
    if (length(abund) > 1 & length(rho_est)==1) {
      rho_est <- rep(rho_est, length(abund))
    }
    if (length(abund) > 1 & length(rho_init)==1) {
      rho_init <- rep(rho_init, length(abund))
    }
  }
  if (length(abund) > 1 & length(p_slide)==1){
    p_slide <- rep(p_slide, length(abund))
  }
  if (is.null(index.w)) index.w <- rep(1,length(abund))

  nr <- nrow(caa)
  nc <- ncol(caa)

  naa <- baa <- ssb <- faa <- saa <- matrix(NA, nrow=nr, ncol=nc)

  if (tune){
    if (use.index[1] != "all") index <- index[use.index,]
    n_ind <- nrow(index)
    q <- b <- sigma <- rep(NA, n_ind)
    d <- rho <- logit_rho_se <- rep(0, n_ind)
    num_index <- sapply(1:n_ind, function(i) length(which(!is.na(index[i,]))))
    avail <- list()
    for (i in 1:nrow(index)){
      avail[[i]] <- which(!is.na(as.numeric(index[i,])))
    }

    ### new code
    zzz <- matrix(1, nrow=nrow(caa1), ncol=ncol(caa1))
    if (!is.null(rec_cor[1])) zzz[1, as.numeric(colnames(caa1)) %in% rec_cor] <- eta
    colnames(zzz) <- colnames(caa1)
    ###
  } else {
    q <- b <- sigma <- d <- rho <- logit_rho_se <- NULL
    num_index <- NULL
    zzz <- NULL
  }

  if (is.null(waa_new)) waa_new <- rowMeans(waa1[,(nc-(wf-1)):nc])/scale

  if (length(abund) > 1 & length(cov_q)==1) {
    cov_q <- rep(cov_q, length(abund))
  }

  if (!any(cov_q)) X_q <- matrix(0,nrow=nrow(index),ncol=ncol(index))

  tf.year <- 1:con_tf

  objective.func <- function(p,output="p"){

    if (tune) {
      faa[seq(2,(nr-2),by=2), nc] <- exp(p)
      faa[nr, nc] <- alpha*exp(p[nr/2-1])
    } else{
      faa[nr,nc] <- exp(p)
    }

    naa[nr,nc] <- caa[nr,nc]/(1-exp(-faa[nr,nc]))*exp(M[nr,nc]/2)

    naa[nr-1,nc] <- naa[nr,nc]*exp(M[nr-1,nc])+caa[nr-1,nc]*exp(M[nr-1,nc]/2)

    faa[nr-1,nc] <- -log(1-caa[nr-1,nc]*exp(M[nr-1,nc]/2)/naa[nr-1,nc])

    for (i in 1:(nc-1)){
      naa[c(nr-2,nr),nc-i] <- caa[c(nr-2,nr),nc-i]*exp(M[c(nr-2,nr),nc-i]/2)+naa[nr-1,nc-(i-1)]*exp(M[c(nr-2,nr),nc-i])*caa[c(nr-2,nr),nc-i]/sum(caa[c(nr-2,nr),nc-i])
      faa[nr-2,nc-i] <- -log(1-caa[nr-2,nc-i]*exp(M[nr-2,nc-i]/2)/naa[nr-2,nc-i])
      faa[nr,nc-i] <- alpha*faa[nr-2,nc-i]
      naa[nr-1,nc-i] <- naa[nr,nc-i]*exp(M[nr-1,nc-i])+caa[nr-1,nc-i]*exp(M[nr-1,nc-i]/2)
      faa[nr-1,nc-i] <- -log(1-caa[nr-1,nc-i]*exp(M[nr-1,nc-i]/2)/naa[nr-1,nc-i])
    }

    for (k in 1:((nr-2)/2)){
      if (!tune) faa[nr-2*k,nc] <- mean(faa[nr-2*k,nc-tf.year])

      naa[nr-2*k,nc] <- caa[nr-2*k,nc]/(1-exp(-faa[nr-2*k,nc]))*exp(M[nr-2*k,nc]/2)

      for (i in 0:(nc-1)){
        naa[nr-2*k-1,nc-i] <- naa[nr-2*k,nc-i]*exp(M[nr-2*k-1,nc-i])+caa[nr-2*k-1,nc-i]*exp(M[nr-2*k-1,nc-i]/2)
        faa[nr-2*k-1,nc-i] <- -log(1-caa[nr-2*k-1,nc-i]*exp(M[nr-2*k-1,nc-i]/2)/naa[nr-2*k-1,nc-i])
        if (i > 0){
          naa[nr-2*k-2,nc-i] <- naa[nr-2*k-1,nc-(i-1)]*exp(M[nr-2*k-2,nc-i])+caa[nr-2*k-2,nc-i]*exp(M[nr-2*k-2,nc-i]/2)
          faa[nr-2*k-2,nc-i] <- -log(1-caa[nr-2*k-2,nc-i]*exp(M[nr-2*k-2,nc-i]/2)/naa[nr-2*k-2,nc-i])
        }
      }
    }

    if (!is.null(rec)){
      naa[1,nc] <- rec
      faa[1,nc] <- -log(1-caa[1,nc]*exp(M[1,nc]/2)/naa[1,nc])
      naa[2,nc] <- (naa[1,nc]*exp(-M[1,nc]/2)-caa[1,nc])*exp(-M[1,nc]/2)
      faa[2,nc] <- -log(1-caa[2,nc]*exp(M[2,nc]/2)/naa[2,nc])
    }

    naa1 <- naa[seq(1,nr,by=2),]
    naa2 <- naa[seq(2,nr,by=2),]
    baa1 <- waa1*naa1/scale
    baa2 <- waa2*naa2/scale
    faa1 <- faa[seq(1,nr,by=2),]
    faa2 <- faa[seq(2,nr,by=2),]
    saa1 <- sel.func(faa1,def=sel.def)
    saa2 <- sel.func(faa2,def=sel.def)

    if (output=="p" | output=="q"){
      obj <- 0
      if (tune){
        for (i in 1:nrow(index)){
          if (abund[i]=="Bo"){
            if (ki[i]==1){
              saa_o1 <- sel.func(zzz*faa1*omega1, def=sel.def)
              saa_o1 <- sweep(saa_o1,2,colSums(saa_o1),FUN="/")
              pred <- colSums((saa_o1*baa1)[min.age[i]:max.age[i],,drop=FALSE])
            }
            if (ki[i]==2){
              saa_o2 <- sel.func(zzz*faa2*omega2, def=sel.def)
              saa_o2 <- sweep(saa_o2,2,colSums(saa_o2),FUN="/")
              pred <- colSums((saa_o2*baa2)[min.age[i]:max.age[i],,drop=FALSE])
            }
          }
          if (abund[i]=="N"){
            if (ki[i]==1){
              pred <- colSums((naa1*exp(-p_slide[i]*(faa1+0.5*dat$M)))[min.age[i]:max.age[i],,drop=FALSE])
            }
            if (ki[i]==2){
              pred <- colSums((naa2*exp(-p_slide[i]*(faa2+0.5*dat$M)))[min.age[i]:max.age[i],,drop=FALSE])
            }
          }
          if (abund[i]=="B"){
            if (ki[i]==1){
              pred <- colSums((baa1*exp(-p_slide[i]*(faa1+0.5*dat$M)))[min.age[i]:max.age[i],,drop=FALSE])
            }
            if (ki[i]==2){
              pred <- colSums((baa2*exp(-p_slide[i]*(faa2+0.5*dat$M)))[min.age[i]:max.age[i],,drop=FALSE])
            }
          }
          if (abund[i]=="Bf"){
            pred <- colSums((baa1)[min.age[i]:max.age[i],,drop=FALSE])
            pred_new <- (naa2[,ncol(naa2)]*exp(-1*(faa2[,ncol(naa2)]+0.5*dat$M[,ncol(naa2)])))
            pred_new <- sum((waa_new*c(0,pred_new[1:(nrow(naa2)-2)],sum(pred_new[(nrow(naa2)-1):nrow(naa2)])))[min.age[i]:max.age[i]])
            pred <- c(pred[-1],pred_new)
          }
          if (abund[i]=="Nf"){
            pred <- colSums((baa1)[min.age[i]:max.age[i],,drop=FALSE])
            pred_new <- (naa2[,ncol(naa2)]*exp(-1*(faa2[,ncol(naa2)]+0.5*dat$M[,ncol(naa2)])))
            pred_new <- sum((c(0,pred_new[1:(nrow(naa2)-2)],sum(pred_new[(nrow(naa2)-1):nrow(naa2)])))[min.age[i]:max.age[i]])

            pred <- c(pred[-1],pred_new)
          }

          if (b_est[i]) b[i] <- cov(log(as.numeric(index[i,avail[[i]]])),log(as.numeric(pred[avail[[i]]])))/var(log(as.numeric(pred[avail[[i]]]))) else b[i] <- b_fix[i]
          q[i] <- exp(mean(log(as.numeric(index[i,avail[[i]]]))-b[i]*log(as.numeric(pred[avail[[i]]]))))
          sigma[i] <- sqrt(sum((log(as.numeric(index[i,avail[[i]]]))-log(q[i])-b[i]*log(as.numeric(pred[avail[[i]]])))^2)/length(avail[[i]]))
          if(sigma[i]==0 | is.na(sigma[i])) sigma[i] <- dd

          if (AR[i]) {
            data_ar <- list(OBS=log(as.numeric(index[i,avail[[i]]])), PRED=log(as.numeric(pred[avail[[i]]])),X_q=as.numeric(X_q[i,avail[[i]]]))
            params <- list(log_q=log(q[i]), log_b=log(b[i]), logit_rho=zc_logit(rho_init[i]), d=0)
            map <- NULL

            if (!cov_q[i]) map$d <- factor(NA)
            if (!b_est[i]) map$log_b <- factor(NA)
            if (!rho_est[i]) map$logit_rho <- factor(NA)
            obj_tmb <- MakeADFun(data_ar, params, map=map, DLL="ar1", silent=TRUE)
            opt <- nlminb(obj_tmb$par, obj_tmb$fn, obj_tmb$gr, lower=Lo_AR, upper=Up_AR)

            q[i] <- obj_tmb$report()$q
            b[i] <- obj_tmb$report()$b
            sigma[i] <- obj_tmb$report()$sigma
            rho[i] <- obj_tmb$report()$rho
            d[i] <- obj_tmb$report()$d

            if (output=="q" & rho_se & rho_est[i] & rho[i] > 0) logit_rho_se[i] <- summary(sdreport(obj_tmb))["logit_rho",2]

            obj <- obj + index.w[i]*opt$objective
          } else{
            if (est.method=="ml") obj <- obj-index.w[i]*sum(as.numeric(na.omit(dnorm(log(as.numeric(index[i,avail[[i]]])),log(q[i])+b[i]*log(as.numeric(pred[avail[[i]]])),sigma[i],log=TRUE)))) else obj <- obj-index.w[i]*sum(as.numeric(na.omit((log(as.numeric(index[i,avail[[i]]]))-(log(q[i])+b[i]*log(as.numeric(pred[avail[[i]]]))))^2)))
          }
        }

        obj <- (1-lambda)*obj + lambda*sum(exp(p)^2)
        if (sel.const) {
          nr2 <- nrow(faa2)
          obj <- obj+W_s*sum((abs(faa2[start_row:(nr2-2),nc,drop=FALSE]-faa2[nr2,nc]*apply(sweep(faa2[start_row:(nr2-2),nc-1:con_tf,drop=FALSE],2,faa2[nr2,nc-1:con_tf],FUN="/"),1,get(stat.tf))))^2)
        }
        if (output=="q") obj <- list(q=q, b=b, sigma=sigma, d=d, rho=rho, logit_rho_se=logit_rho_se, obj=obj)
      } else obj <- (faa[nr-2,nc]-faa[nr,nc])^2
    }

    if (output=="all"){

      baa <- waa*naa/scale
      saa <- sel.func(faa,def=sel.def)

      uaa1 <- caa1/naa1
      uaa2 <- caa2/naa2

      rownames(naa1) <- rownames(faa1) <- rownames(uaa1) <- rownames(saa1) <- rownames(caa1)
      rownames(naa2) <- rownames(faa2) <- rownames(uaa2) <- rownames(saa2) <- rownames(caa2)
      colnames(naa) <- colnames(baa) <- colnames(ssb) <- colnames(faa) <- colnames(saa) <- colnames(naa1) <- colnames(faa1) <- colnames(uaa1) <-colnames(saa1) <- colnames(caa1)
      colnames(naa2) <- colnames(faa2) <- colnames(uaa2) <- colnames(saa2) <- colnames(caa2)

      if (ssb_type==1) ssb <- naa1*maa*waa2/scale
      if (ssb_type==2) ssb <- naa2*maa_ssb*waa_ssb*exp(-p_ssb*dat$M)/scale
      if (ssb_type==3) ssb <- naa2*maa*waa2*exp(-faa2-0.5*dat$M)/scale

      obj <- list(waa_new=waa_new, naa=naa, naa1=naa1, naa2=naa2, baa=baa, baa1=baa1, baa2=baa2, ssb=ssb, faa=faa, faa1=faa1, faa2=faa2, uaa1=uaa1, uaa2=uaa2, saa=saa, saa1=saa1, saa2=saa2)
    }

    return(obj)
  }

  if (no.est) {
    res <- NULL
    p <- log(p.init)
  } else {
    if (optimizer=="nlm") res <- nlm(objective.func, log(p.init), hessian=hessian)
    if (optimizer=="nlminb") {
      res <- nlminb(log(p.init), objective.func, hessian=hessian, lower=Lower, upper=Upper)
      res$estimate <- res$par
      res$minimum <- res$objective
      res$gradient <- NA
      res$code <- res$convergence
    }
    p <- res$estimate
  }

  res$opts <-  objective.func(p,output="q")

  if (any(AR)) nq <- sum(res$opts$rho > 0) else nq <- 0

  res$aic <- res$minimum*2+2*(length(res$estimate)+sum(b_est)+nq)

  res$outputs <- objective.func(p,output="all")

  res$input <- arglist

  res$AR <- AR
  res$X_q <- X_q
  res$num_index <- num_index
  res$wcaa <- as.data.frame(caa*waa)
  res$naa <- as.data.frame(res$outputs$naa)
  res$faa <- as.data.frame(res$outputs$faa)
  res$baa <- as.data.frame(res$outputs$baa)
  res$ssb <- as.data.frame(res$outputs$ssb)
  res$saa <- as.data.frame(res$outputs$saa)
  res$Fc.at.age <- apply(res$faa[,ncol(res$faa)-0:(con_fc-1),drop=FALSE],1,mean)
  res$zzz <- zzz

  class(res) <- "hvpa"
  return(res)
} #vpa.hay


zenki <- function(res){
  caa1 <- res$input$dat$caa1
  caa2 <- res$input$dat$caa2

  naa1 <- res$outputs$naa1
  naa2 <- res$outputs$naa2
  M <- res$input$dat$M

  faa <- -log(1-(caa1*exp(M/4)+caa2*exp(3/4*M))/naa1)

  colnames(faa) <- colnames(naa1)

  list(naa=naa1,faa=faa)
}

#' 半期VPAのレトロスペクティブ解析用の関数
#'
#' @details `do_retrospective_vpa`関数内でこの関数を動かしてレトロスペクティブ解析の図示もしている。
#'
#' @encoding UTF-8
#'
#' @export

hretro_est <- function(res,n=5,b_fix=TRUE,rho_fix=TRUE){
  res.c <- res

  if (!is.null(res$num_index) & n >= min(res$num_index)) stop("n >= min(num_index)")

  if (b_fix){
    if (class(res$opts)=="list"){
      res.c$input$b_fix <- res$opts$b
      res.c$input$b_est <- rep(FALSE,length(res$opts$b))
    }
  }

  if (rho_fix){
    if (class(res$opts)=="list"){
      if (!is.null(res$opts$rho)){
        res.c$input$rho_init <- res$opts$rho
        res.c$input$rho_est <- rep(FALSE,length(res$opts$rho))
      }
    }
  }

  Res <- list()

  for (i in 1:n){
    nc <- ncol(res.c$input$dat$caa1)

    res.c$input$dat$caa1 <- res.c$input$dat$caa1[,-nc]
    res.c$input$dat$caa2 <- res.c$input$dat$caa2[,-nc]
    res.c$input$dat$maa <- res.c$input$dat$maa[,-nc]
    res.c$input$dat$maa_ssb <- res.c$input$dat$maa_ssb[,-nc]
    res.c$input$dat$waa1 <- res.c$input$dat$waa1[,-nc]
    res.c$input$dat$waa2 <- res.c$input$dat$waa2[,-nc]
    res.c$input$dat$waa_ssb <- res.c$input$dat$waa_ssb[,-nc]
    res.c$input$dat$omega1 <- res.c$input$dat$omega1[,-nc]
    res.c$input$dat$omega2 <- res.c$input$dat$omega2[,-nc]
    res.c$input$dat$M <- res.c$input$dat$M[,-nc]
    res.c$input$dat$index <- res.c$input$dat$index[,-nc,drop=FALSE]

    res11 <- do.call(vpa.hay,res.c$input)

    Res[[i]] <- res11
  }

  retro_res <- NULL

  nc <- ncol(res$input$dat$caa1)

  for (i in 1:n){
    retro_res <- rbind(retro_res, c(last((colSums(Res[[i]]$outputs$naa1)- colSums(res$outputs$naa1)[-((nc-i+1):nc)])/colSums(res$outputs$naa1)[-((nc-i+1):nc)]),last((colSums(Res[[i]]$outputs$naa2)- colSums(res$outputs$naa2)[-((nc-i+1):nc)])/colSums(res$outputs$naa2)[-((nc-i+1):nc)]),last((colSums(Res[[i]]$outputs$baa1)- colSums(res$outputs$baa1)[-((nc-i+1):nc)])/colSums(res$outputs$baa1)[-((nc-i+1):nc)]),last((colSums(Res[[i]]$outputs$baa2)- colSums(res$outputs$baa2)[-((nc-i+1):nc)])/colSums(res$outputs$baa2)[-((nc-i+1):nc)]),last((colSums(Res[[i]]$ssb)- colSums(res$ssb)[-((nc-i+1):nc)])/colSums(res$ssb)[-((nc-i+1):nc)]),last((Res[[i]]$outputs$naa1[1,]- res$outputs$naa1[1,-((nc-i+1):nc)])/res$outputs$naa1[1,-((nc-i+1):nc)]),last((Res[[i]]$outputs$naa2[1,]- res$outputs$naa2[1,-((nc-i+1):nc)])/res$outputs$naa2[1,-((nc-i+1):nc)]),last((colSums(Res[[i]]$outputs$faa1)- colSums(res$outputs$faa1)[-((nc-i+1):nc)])/colSums(res$outputs$faa1)[-((nc-i+1):nc)]),last((colSums(Res[[i]]$outputs$faa2)- colSums(res$outputs$faa2)[-((nc-i+1):nc)])/colSums(res$outputs$faa2)[-((nc-i+1):nc)])))
  }

  colnames(retro_res) <- c("N1","N2","B1","B2","SSB","R1","R2","F1","F2")

  list(Res=Res, retro_res=retro_res, mohn=colMeans(retro_res))
}

sel_const_check <- function(res){
  faa2 <- res$outputs$faa2

  nr <- nrow(faa2)
  nc <- ncol(faa2)

  cf <- res$input$con_tf

  target <- rowMeans(sweep(faa2[,nc-1:cf],2,faa2[nr,nc-1:cf],FUN="/"))*faa2[nr,nc]

  relative_bias <- (faa2[,nc]-target)/abs(target)*100

  list(faa_target=target, faa_terminal=faa2[,nc], percent_relative_bias=relative_bias)
}

#' 半期VPAの作図関数
#'
#' @export
#' @encoding UTF-8

hvpa_plot <- function(res,res0=NULL,years=NULL,age=NULL,rel_heights=c(0.7,0.3),labels=c("Abundance","Fishing Rate")){
  assertthat::assert_that(class(res) == "hvpa")

  if (is.null(age)) age <- 0:4

  dat <- res$input$dat

  caa1 <- dat$caa1
  caa2 <- dat$caa2
  naa1 <- res$outputs$naa1
  naa2 <- res$outputs$naa2

  if (!is.null(age)){
    Age <- age+1
    caa1 <- caa1[Age,,drop=FALSE]
    caa2 <- caa2[Age,,drop=FALSE]
    naa1 <- naa1[Age,,drop=FALSE]
    naa2 <- naa2[Age,,drop=FALSE]
  }

  caa_all <- rbind(cbind(caa1,caa2),rep(1:ncol(caa1),2))

  caa <- caa_all[1:nrow(caa1),order(caa_all[nrow(caa_all),])]

  naa_all <- rbind(cbind(naa1,naa2),rep(1:ncol(naa1),2))

  naa <- naa_all[1:nrow(naa1),order(naa_all[nrow(naa_all),])]

  if (is.null(res0)){

    data_for_plot <- data.frame(year=rep(colnames(caa1),each=nrow(caa1)*2),season=rep(c(rep(1,nrow(caa1)),rep(2,nrow(caa1))),ncol(caa1)),age=rep(age,2*ncol(caa1)),caa=as.numeric(unlist(caa)),naa=as.numeric(unlist(naa)))

    data_for_plot <- data_for_plot %>% mutate(ys=as.Date(if_else(season==1,paste0(year,"-01-01"),paste0(year,"-06-01"))))

    data_for_plot2 <- data_for_plot %>% group_by(year,season) %>% summarize(C=sum(caa),N=sum(naa),ys=mean(ys))
  } else {
    naa01 <- res0$outputs$naa1
    naa02 <- res0$outputs$naa2

    if (!is.null(age)){
      naa01 <- naa01[Age,,drop=FALSE]
      naa02 <- naa02[Age,,drop=FALSE]
    }

    naa0_all <- rbind(cbind(naa01,naa02),rep(1:ncol(naa01),2))

    naa0 <- naa0_all[1:nrow(naa01),order(naa0_all[nrow(naa0_all),])]

    data_for_plot <- data.frame(year=rep(colnames(caa1),each=nrow(caa1)*2),season=rep(c(rep(1,nrow(caa1)),rep(2,nrow(caa1))),ncol(caa1)),age=rep(age,2*ncol(caa1)),caa=as.numeric(unlist(caa)),naa=as.numeric(unlist(naa)),naa0=as.numeric(unlist(naa0)))

    data_for_plot <- data_for_plot %>% mutate(ys=as.Date(if_else(season==1,paste0(year,"-01-01"),paste0(year,"-06-01"))))

    data_for_plot2 <- data_for_plot %>% group_by(year,season) %>% summarize(C=sum(caa),N=sum(naa),N0=sum(naa0),ys=mean(ys))
  }

  if (!is.null(years)) data_for_plot2 <- subset(data_for_plot2, year %in% years)

  p <- ggplot(data_for_plot2, mapping=aes(x=ys, y=N))

  p1 <- p + geom_line(size=2)+labs(x="Year", y=labels[1])+theme_bw()

  if (!is.null(res0)) p1 <- p1 + geom_line(aes(x=ys, y=N0),size=1,color="blue",linetype="twodash")

  pp <- ggplot(data_for_plot2, mapping=aes(x=ys, y=C/N))

  pp1 <- pp + geom_line(size=2, color="red")+labs(x="Year", y=labels[2])+theme_bw()

  if (!is.null(res0)) pp1 <- pp1 + geom_line(aes(x=ys, y=C/N0),size=1, color="green",linetype="twodash")

  cowplot::plot_grid(p1, pp1, nrow=2, rel_heights=rel_heights,align="v")
}

#' 半期VPAの作図関数
#'
#' @export
#' @encoding UTF-8

hvpa_plot2 <- function(res,res0=NULL,years=NULL,age=NULL,rel_heights=c(0.7,0.3),labels=c("Abundance","Fishing Rate")){
  assertthat::assert_that(class(res) == "hvpa")

  if (is.null(age)) age <- 0:4

  dat <- res$input$dat

  caa1 <- dat$caa1
  caa2 <- dat$caa2
  naa1 <- res$outputs$naa1
  naa2 <- res$outputs$naa2

  if (!is.null(age)){
    Age <- age+1
    caa1 <- caa1[Age,,drop=FALSE]
    caa2 <- caa2[Age,,drop=FALSE]
    naa1 <- naa1[Age,,drop=FALSE]
    naa2 <- naa2[Age,,drop=FALSE]
  }

  if (is.null(res0)){

    data_for_plot <- data.frame(year=rep(colnames(caa1),each=nrow(caa1)),season=rep(c(rep(1,nrow(caa1))),ncol(caa1)),age=rep(age,ncol(caa1)),caa=as.numeric(unlist(caa1+caa2)),naa=as.numeric(unlist(naa1)))

    data_for_plot <- data_for_plot %>% mutate(ys=as.Date(paste0(year,"-01-01")))

    data_for_plot2 <- data_for_plot %>% group_by(year,season) %>% summarize(C=sum(caa),N=sum(naa),ys=mean(ys))
  } else {
    naa01 <- res0$naa

    if (!is.null(age)){
      naa01 <- naa01[Age,,drop=FALSE]
    }

    data_for_plot <- data.frame(year=rep(colnames(caa1),each=nrow(caa1)),season=rep(c(rep(1,nrow(caa1))),ncol(caa1)),age=rep(age,ncol(caa1)),caa=as.numeric(unlist(caa1+caa2)),naa=as.numeric(unlist(naa1)),naa0=as.numeric(unlist(naa01)))

    data_for_plot <- data_for_plot %>% mutate(ys=as.Date(paste0(year,"-01-01")))

    data_for_plot2 <- data_for_plot %>% group_by(year,season) %>% summarize(C=sum(caa),N=sum(naa),N0=sum(naa0),ys=mean(ys))
  }

  if (!is.null(years)) data_for_plot2 <- subset(data_for_plot2, year %in% years)

  p <- ggplot(data_for_plot2, mapping=aes(x=ys, y=N))

  p1 <- p + geom_line(size=2)+labs(x="Year", y=labels[1])+theme_bw()

  if (!is.null(res0)) p1 <- p1 + geom_line(aes(x=ys, y=N0),size=1,color="blue",linetype="twodash")

  pp <- ggplot(data_for_plot2, mapping=aes(x=ys, y=C/N))

  pp1 <- pp + geom_line(size=2, color="red")+labs(x="Year", y=labels[2])+theme_bw()

  if (!is.null(res0)) pp1 <- pp1 + geom_line(aes(x=ys, y=C/N0),size=1, color="green",linetype="twodash")

  cowplot::plot_grid(p1, pp1, nrow=2, rel_heights=rel_heights,align="v")
}

#' 半期VPA専用のレトロスペクティブ解析の作図関数
#'
#' @export
#' @encoding UTF-8

hretro_plot <- function(
    res,
    retro,
    start_year=2014,
    target="naa",
    age=0:4,
    season=1,
    add_mohn=TRUE,
    digits=3,
    hanki=TRUE,
    out="plot",
    title_size=NULL,
    ylim=NULL,
    xlab=NULL,
    ylab=NULL,
    size=1.5
){
  assertthat::assert_that(class(res) == "hvpa")
  age <- age + 1

  if (hanki) res00 <- res$outputs else res00 <- res
  dat <- res00[names(res00)==paste0(target,season)]

  years <- as.numeric(dimnames(dat[[1]])[[2]])

  dat1 <- data.frame(type=1,year=years,value=colSums(as.data.frame(dat)[age,]))

  n <- length(retro$Res)
  ny <- length(years)

  mohn_rho <- NULL

  for (k in 1:n){
    res1 <- retro$Res[[k]]
    if (hanki) res00 <- res1$outputs else res00 <- res1
    dat <- res00[names(res00)==paste0(target,season)]
    years <- as.numeric(dimnames(dat[[1]])[[2]])

    dat1 <- rbind(dat1, data.frame(type=k+1,year=years,value=colSums(as.data.frame(dat)[age,])))

    mohn_rho <- c(mohn_rho, last((subset(dat1,type==k+1)$value-subset(dat1,type==1)$value[1:(ny-k)])/subset(dat1,type==1)$value[1:(ny-k)]))
  }

  dat1 <- as_tibble(dat1)

  dat1$type <- factor(dat1$type)

  age <- age - 1

  if (add_mohn) mohn_title <- paste0("\nMohn's rho = ", round(mean(mohn_rho),digits=digits)) else mohn_title <- NULL

  if (is.null(season)) season_title <- NULL else season_title <- paste0(", season = ",season)

  if (length(age)==1) age_title <- paste0(", age = ", age) else age_title <- paste0(", age = ", age[1], "-", last(age))

  p1 <- ggplot(subset(dat1, year >= start_year), aes(x=year, y=value, color=type))+geom_line(size=size)+theme_bw()+
    labs(x = "Year", y=paste0(target,season), title=paste0("target = ",target, season_title, age_title, mohn_title))+theme(plot.title = element_text(size=title_size))+guides(color = "none")

  if (!is.null(ylim)) p1 <- p1 + coord_cartesian(ylim=ylim)
  if (!is.null(xlab)) p1 <- p1 + labs(x=xlab)
  if (!is.null(ylab)) p1 <- p1 + labs(y=ylab)

  print(p1)

  if (out=="dat") return(list(dat=dat1,mohn=mohn_rho)) else return(p1)
}

#' @export
#' @encoding UTF-8

make_dat <- function(res){

  AR <- res$input$AR

  abund <- res$input$abund
  index <- res$input$dat$index
  use.index <- res$input$use.index
  if (use.index[1] != "all") index <- index[use.index,]
  ki <- res$input$ki
  sel.def <- res$input$sel.def
  min.age <- res$input$min.age+1
  max.age <- res$input$max.age+1
  b_est <- res$input$b_est
  if (length(abund) > 1 & length(b_est)==1) b_est <- rep(b_est, length(abund))
  p_slide <- res$input$p_slide
  if (length(abund) > 1 & length(p_slide)==1) p_slide <- rep(p_slide, length(abund))

  waa1 <- res$input$dat$waa1
  waa2 <- res$input$dat$waa2
  omega1 <- res$input$dat$omega1
  omega2 <- res$input$dat$omega2
  M <- res$input$dat$M

  wf <- res$input$wf
  scale <- res$input$scale

  faa1 <- res$outputs$faa1
  faa2 <- res$outputs$faa2
  naa1 <- res$outputs$naa1
  naa2 <- res$outputs$naa2
  baa1 <- res$outputs$baa1
  baa2 <- res$outputs$baa2

  b <- res$opts$b
  q <- res$opts$q
  d <- res$opts$d

  X_q <- res$X_q

  zzz <- res$zzz

  avail <- list()
  for (i in 1:nrow(index)){
    avail[[i]] <- which(!is.na(as.numeric(index[i,])))
  }

  waa_new <- res$outputs$waa_new

  dat_for_plot <- NULL

  for (i in 1:nrow(index)){
    if (abund[i]=="Bo"){
      if (ki[i]==1){
        saa_o1 <- sel.func(zzz*faa1*omega1, def=sel.def)
        saa_o1 <- sweep(saa_o1,2,colSums(saa_o1),FUN="/")
        pred <- colSums((saa_o1*baa1)[min.age[i]:max.age[i],,drop=FALSE])
      }
      if (ki[i]==2){
        saa_o2 <- sel.func(zzz*faa2*omega2, def=sel.def)
        saa_o2 <- sweep(saa_o2,2,colSums(saa_o2),FUN="/")
        pred <- colSums((saa_o2*baa2)[min.age[i]:max.age[i],,drop=FALSE])
      }
    }
    if (abund[i]=="N"){
      if (ki[i]==1){
        pred <- colSums((naa1*exp(-p_slide[i]*(faa1+0.5*M)))[min.age[i]:max.age[i],,drop=FALSE])
      }
      if (ki[i]==2){
        pred <- colSums((naa2*exp(-p_slide[i]*(faa2+0.5*M)))[min.age[i]:max.age[i],,drop=FALSE])
      }
    }
    if (abund[i]=="B"){
      if (ki[i]==1){
        pred <- colSums((baa1*exp(-p_slide[i]*(faa1+0.5*M)))[min.age[i]:max.age[i],,drop=FALSE])
      }
      if (ki[i]==2){
        pred <- colSums((baa2*exp(-p_slide[i]*(faa2+0.5*M)))[min.age[i]:max.age[i],,drop=FALSE])
      }
    }
    if (abund[i]=="Bf"){
      pred <- colSums((baa1)[min.age[i]:max.age[i],,drop=FALSE])
      pred_new <- (naa2[,ncol(naa2)]*exp(-1*(faa2[,ncol(naa2)]+0.5*dat$M[,ncol(naa2)])))
      pred_new <- sum((waa_new*c(0,pred_new[1:(nrow(naa2)-2)],sum(pred_new[(nrow(naa2)-1):nrow(naa2)])))[min.age[i]:max.age[i]])

      pred <- c(pred[-1],pred_new)
    }
    if (abund[i]=="Nf"){
      pred <- colSums((baa1)[min.age[i]:max.age[i],,drop=FALSE])
      pred_new <- (naa2[,ncol(naa2)]*exp(-1*(faa2[,ncol(naa2)]+0.5*dat$M[,ncol(naa2)])))
      pred_new <- sum((c(0,pred_new[1:(nrow(naa2)-2)],sum(pred_new[(nrow(naa2)-1):nrow(naa2)])))[min.age[i]:max.age[i]])

      pred <- c(pred[-1],pred_new)
    }

    pred <- log(q[i])+d[i]*as.numeric(X_q[i,avail[[i]]])+b[i]*log(as.numeric(pred[avail[[i]]]))

    obs <- as.numeric(index[i,avail[[i]]])

    dat_for_plot <- rbind(dat_for_plot, data.frame(year=as.numeric(names(index)[avail[[i]]]),num=i,obs=obs,pred=exp(pred)))
  }

  if (any(AR)) {
    rho <- res$opts$rho
    resid <- NULL
    for (j in 1:nrow(index)){
      dat_for_plot_sub <- subset(dat_for_plot, num==j)
      n <- nrow(dat_for_plot_sub)
      resid <- c(resid,log(dat_for_plot_sub$obs[1])-log(dat_for_plot_sub$pred[1]))
      eps <- last(resid)
      for (i in 2:n){
        resid <- c(resid,log(dat_for_plot_sub$obs[i])-log(dat_for_plot_sub$pred[i])-rho[j]*eps)
        eps <- log(dat_for_plot_sub$obs[i])-log(dat_for_plot_sub$pred[i])
      }
    }
    dat_for_plot$residuals <- resid
  } else dat_for_plot$residuals <- log(dat_for_plot$obs)-log(dat_for_plot$pred)

  return(dat_for_plot)
}

#' 半期VPA専用の残差プロット作図の関数
#'
#' @export
#' @encoding UTF-8

residual_plot <- function(
    res,
    rel_heights=c(0.6,0.4),
    add_log_plot=TRUE,
    out="plot",
    abund_name=NULL,
    stat=TRUE
){
  assertthat::assert_that(class(res) == "hvpa")

  dat_for_plot <- make_dat(res)

  if(is.null(abund_name)) abund_name <- res$input$abund
  min.age <- res$input$min.age+1
  max.age <- res$input$max.age+1
  age_name <- rownames(res$outputs$naa1)

  dat_for_plot$name <- as.factor(paste0("Index ",dat_for_plot$num,": ",abund_name[dat_for_plot$num]," (age ",ifelse(min.age[dat_for_plot$num]==max.age[dat_for_plot$num],age_name[min.age[dat_for_plot$num]],paste0(age_name[min.age[dat_for_plot$num]],"-",age_name[max.age[dat_for_plot$num]])),")"))

  p1 <- ggplot(dat_for_plot,aes(x=year, y=obs))+geom_point()+geom_line(aes(x=year, y=pred),color="red")+facet_wrap(~name, scales = "free")+labs(x="Year",y="Abundance")+theme_bw()

  if (add_log_plot) p2 <- ggplot(dat_for_plot,aes(x=year, y=log(obs)))+geom_point()+geom_line(aes(x=year, y=log(pred)),color="red")+facet_wrap(~name, scales = "free")+labs(x="Year",y="log Abundance")+theme_bw() else p2 <- NULL

  p3 <- ggplot(dat_for_plot,aes(x=year, y=residuals))+geom_point()+geom_hline(yintercept=0,color="red",linetype="dashed")+geom_smooth(method="lm",se = TRUE, level = 0.95)+facet_wrap(~name)+labs(x="Year",y="Residuals")+theme_bw()
  if (stat) p3 <- p3+ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(stat(eq.label), stat(rr.label), stat(p.value.label), sep = "~~~")), parse = TRUE)

  if (out=="dat") return(dat_for_plot) else if (add_log_plot) {
    return(cowplot::plot_grid(p1, p2, p3, nrow=3, rel_heights=rel_heights,align="v"))
  } else return(cowplot::plot_grid(p1, p3, nrow=2, rel_heights=rel_heights,align="v"))
}

#' 半期VPA専用のバブルプロット作図の関数
#'
#' @export
#' @encoding UTF-8

bubble_plot <- function(
    res,
    target="naa",
    season=1,
    hanki=TRUE,
    years=NULL,
    scale=100,
    fix_ratio=2,
    range=seq(0,100,len=5),
    digits=0,
    max_size=10
){
  assertthat::assert_that(class(res) == "hvpa")
  if (hanki) res00 <- res$outputs else res00 <- res
  dat <- res00[names(res00)==paste0(target,season)]

  dat <- dat[[1]] %>% as.data.frame() %>% tibble::rownames_to_column("age") %>% pivot_longer(!age, names_to = "year", values_to = "value") %>% mutate_at(vars(year), as.factor)

  dat <- dat %>% mutate_at(vars(age), as.factor)
  dat$age <- factor(dat$age, levels=rev(levels(dat$age)))

  range <- range[1:(length(range)-1)]

  if(!is.null(years)) dat <- subset(dat, year %in% years)

  ggplot(dat,aes(x=year, y=age)) + geom_point(aes(size=value, fill=age, color=age), alpha = 0.75, shape = 21) + coord_fixed(ratio=fix_ratio) + scale_size_continuous(limits = c(0.0001, 1.0001)*max(dat$value,na.rm=TRUE), range = c(1,max_size), breaks = round(range*max(dat$value,na.rm=TRUE)/scale,digits=digits))+theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + labs(x="Year", y="Age",size=paste0(target,season)) + guides(color="none", fill="none")
}

future_projection <- function(
    res,
    sr,
    Y=20,
    beta=1,
    scale=1000){

  faa1 <- res$outputs$faa1
  faa2 <- res$outputs$faa2

  nc <- ncol(faa1)

  faa1_current <- faa1[,nc]
  faa2_current <- faa2[,nc]

  caa1 <- res$input$dat$caa1
  caa2 <- res$input$dat$caa2
  naa1 <- res$outputs$naa1
  naa2 <- res$outputs$naa2
  maa <- res$input$dat$maa
  waa1 <- res$input$dat$waa1/scale
  waa2 <- res$input$dat$waa2/scale
  SSB <- res$outputs$ssb
  M <- res$input$dat$M
  ssb_type <- res$input$ssb_type
  p_ssb <- res$input$p_ssb
  M <- M[,nc]

  faaz_current <- -log(1-(caa1[,nc]*exp(M/4)+caa2[,nc]*exp(3/4*M))/naa1[,nc])

  a <- SR$pars[1]
  b <- SR$pars[2]
  SRtype <- SR$input$SR

  pred_SR <- function(SSB,SRtype="HS"){
    if (SRtype=="HS") out <- as.numeric(a*min(b,colSums(SSB)))
    if (SRtype=="BH") out <- as.numeric(a*SSB/(1+b*SSB))
    if (SRtype=="RI") out <- as.numeric(a*SSB*exp(-b*SSB))

    return(out)
  }

  nage <- nrow(naa1)

  Naa1 <- c(pred_SR(SSB,SRtype=SRtype),naa2[,nc]*exp(-faa2[,nc]-M/2))
  Naa1[nage] <- sum(Naa1[nage:(nage+1)])
  naa1_f <- as.matrix(Naa1[1:nage],ncol=1)
  caa1_f <- caa2_f <- naa2_f <- NULL

  Naaz <- Naa1
  naaz_f <- naa1_f
  caaz_f <- NULL

  for (i in 1:Y){
    Naa2 <- naa1_f[,i]*exp(-beta*faa1_current-M/2)
    Caa1 <- naa1_f[,i]*exp(-M/4)*(1-exp(-beta*faa1_current))
    naa2_f <- cbind(naa2_f,Naa2)
    caa1_f <- cbind(caa1_f,Caa1)
    if (ssb_type==1) {
      SSB <- Naa1*maa*waa2
      SSBz <- Naaz*maa*waa2
    }
    if (ssb_type==2) {
      SSB <- Naa2*maa_ssb*waa_ssb*exp(-p_ssb*M)
      SSBz <- Naaz*maa_ssb*waa_ssb*exp(-p_ssb*M)
    }
    if (ssb_type==3) {
      SSB <- Naa2*maa*waa2*exp(-beta*faa2_current-0.5*M)
      SSBz <- Naaz*maa*waa2*exp(-beta*faaz_current-M)
    }
    Naa1 <- c(pred_SR(SSB,SRtype=SRtype),naa2_f[,i]*exp(-beta*faa2_current-M/2))
    Naa1[nage] <- sum(Naa1[nage:(nage+1)])
    Naa1 <- Naa1[1:nage]
    Caa2 <- naa2_f[,i]*exp(-M/4)*(1-exp(-beta*faa2_current))
    naa1_f <- cbind(naa1_f,Naa1)
    caa2_f <- cbind(caa2_f,Caa2)

    Naaz <- c(pred_SR(SSBz,SRtype=SRtype),naaz_f[,i]*exp(-beta*faaz_current-M))
    Naaz[nage] <- sum(Naaz[nage:(nage+1)])
    Naaz <- Naaz[1:nage]
    Caaz <- naaz_f[,i]*exp(-M/2)*(1-exp(-beta*faaz_current))

    naaz_f <- cbind(naaz_f, Naaz)
    caaz_f <- cbind(caaz_f,Caaz)
  }

  naa1_f <- naa1_f[,1:Y]
  naaz_f <- naaz_f[,1:Y]

  years <- as.numeric(colnames(naa1)[nc])+1:Y

  colnames(naa1_f) <- colnames(naa2_f) <- colnames(naaz_f) <- colnames(caa1_f) <- colnames(caa2_f) <- colnames(caaz_f) <- years

  list(faa1_current=faa1_current,faa2_current=faa2_current,faaz_current=faaz_current,naa1_f=naa1_f, naa2_f=naa2_f, caa1_f=caa1_f, caa2_f=caa2_f, caa12_f=caa1_f+caa2_f, naaz_f=naaz_f, caaz_f=caaz_f)
}

jackknife_test <- function(res){
  assertthat::assert_that(class(res) == "hvpa")

  index <- res$input$dat$index
  use.index <- res$input$use.index
  if (use.index[1] != "all") index <- index[use.index,]

  res.c <- res

  ResJ <- list()

  k <- 1

  for (i in 1:nrow(index)){
    for (j in 1:ncol(index)){
      if (!is.na(index[i,j])){
        res.c$input$dat$index[use.index[i],j] <- NA
        res11 <- do.call(vpa.hay,res.c$input)

        ResJ[[k]] <- res11
        k <- k + 1
        res.c <- res
      }
    }
  }

  list(res=res,ResJ=ResJ)
}

jackknife_plot <- function(
    jack_res,
    target="naa",
    age=0,
    season=1,
    num_id=1,
    year=2022,
    years=2005:2022,
    scale=100,
    fix_ratio=2,
    range=seq(0,100,len=5),
    digits=1,
    absolute=FALSE,
    max_value=NULL,
    max_size=10
){
  res <- jack_res$res

  index <- res$input$dat$index
  use.index <- res$input$use.index
  if (use.index[1] != "all") index <- index[use.index,]

  impact <- index

  Years <- as.numeric(colnames(index))
  Ages <- 1:nrow(res$outputs$naa1)-1

  k <- 1

  for (i in 1:nrow(index)){
    for (j in 1:ncol(index)){
      if (!is.na(index[i,j])){
        if (target=="b" | target=="q" | target=="sigma"){
          resJ <- jack_res$ResJ[[k]]$opts
          resJ <- resJ[names(resJ)==paste0(target)]
          impact[i,j] <- (resJ[[1]])[num_id]
        } else {
          resJ <- jack_res$ResJ[[k]]$outputs
          resJ <- resJ[names(resJ)==paste0(target,season)]
          impact[i,j] <- (resJ[[1]])[Ages %in% age, Years %in% year]
        }

        k <- k + 1
      }
    }
  }

  dat <- impact %>% as.data.frame() %>% tibble::rownames_to_column("index") %>% pivot_longer(!index, names_to = "year", values_to = "value") %>% mutate_at(vars(year), as.factor)

  dat <- dat %>% mutate_at(vars(index), as.factor)
  dat$index <- factor(dat$index, levels=rev(levels(dat$index)))

  range <- range[1:(length(range)-1)]

  dat <- subset(dat, year %in% years)

  colours <- c("blue","red")

  if(is.null(max_value)) max_value <- max(abs(dat$value),na.rm=TRUE)

  if (target=="b" | target=="q" | target=="sigma") {
    Title <- paste0("target = ",target, ", index = ",num_id)
    season <- NULL
    threshold <- ((res$opts)[names(res$opts)==paste0(target)][[1]])[num_id]
  } else {
    Title <- paste0("target = ",target, ", season = ",season, ", age = ", age, ", year = ", year)
    threshold <- ((res$outputs)[names(res$outputs)==paste0(target,season)][[1]])[Ages %in% age, Years %in% year]
  }

  ggplot(dat,aes(x=year, y=index)) + geom_point(aes(size=abs(value)), fill=colours[(dat$value > threshold)*0.5+1.5], color=colours[(dat$value > threshold)*0.5+1.5], alpha = 0.75, shape = 21) + coord_fixed(ratio=fix_ratio) + scale_size_continuous(limits = c(0.0001, 1.0001)*max_value, range = c(1,max_size), breaks = round(range*max_value/scale,digits=digits))+theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + labs(x="Year", y="Index", size=paste0(target,season),title=Title) + guides(color="none", fill="none")
}

sensitivity_test <- function(res,Delta=2,SD=NULL){
  assertthat::assert_that(class(res) == "hvpa")

  index <- res$input$dat$index
  use.index <- res$input$use.index
  if (use.index[1] != "all") index <- index[use.index,]

  sigma <- res$opts$sigma
  if(!is.null(SD)) sigma <- rep(SD,length(sigma))

  res.c <- res

  ResP <- ResN <- list()

  k <- 1

  for (i in 1:nrow(index)){
    for (j in 1:ncol(index)){
      if (!is.na(index[i,j])){
        res.c$input$dat$index[use.index[i],j] <- exp(log(res$input$dat$index[use.index[i],j])+Delta*sigma[i])
        res11 <- do.call(vpa.hay,res.c$input)
        res.c$input$dat$index[use.index[i],j] <- exp(log(res$input$dat$index[use.index[i],j])-Delta*sigma[i])
        res12 <- do.call(vpa.hay,res.c$input)

        ResP[[k]] <- res11
        ResN[[k]] <- res12
        k <- k + 1
        res.c <- res
      }
    }
  }

  list(Delta=Delta,sigma=sigma,res=res, ResP=ResP, ResN=ResN)
}


sensitivity_plot <- function(
    sens_res,
    target="naa",
    age=0,
    season=1,
    year=2022,
    num_id=1,
    years=2005:2022,
    scale=100,
    fix_ratio=2,
    range=seq(0,100,len=5),
    digits=1,
    absolute=FALSE,
    max_value=NULL,
    max_size=10
){
  res <- sens_res$res
  Delta <- sens_res$Delta

  index <- res$input$dat$index
  use.index <- res$input$use.index
  if (use.index[1] != "all") index <- index[use.index,]

  sigma <- sens_res$sigma

  impact <- index

  Years <- as.numeric(colnames(index))
  Ages <- 1:nrow(res$outputs$naa1)-1

  k <- 1

  for (i in 1:nrow(index)){
    for (j in 1:ncol(index)){
      if (!is.na(index[i,j])){
        if (target=="b" | target=="q" | target=="sigma"){
          resP <- sens_res$ResP[[k]]$opts
          resN <- sens_res$ResN[[k]]$opts
          resP <- resP[names(resP)==paste0(target)]
          resN <- resN[names(resN)==paste0(target)]
          impact[i,j] <- (log((resP[[1]])[num_id])-log((resN[[1]])[num_id]))/(2*Delta*sigma[i])
          if (absolute) impact[i,j] <- impact[i,j]*(res$opts[names(res$opts)==paste0(target)][[1]])[num_id]
        } else {
          resP <- sens_res$ResP[[k]]$outputs
          resN <- sens_res$ResN[[k]]$outputs

          resP <- resP[names(resP)==paste0(target,season)]
          resN <- resN[names(resN)==paste0(target,season)]
          impact[i,j] <- (log((resP[[1]])[Ages %in% age, Years %in% year])-log((resN[[1]])[Ages %in% age, Years %in% year]))/(2*Delta*sigma[i])

          if (absolute) impact[i,j] <- impact[i,j]*(res$outputs[names(res$outputs)==paste0(target,season)][[1]])[Ages %in% age, Years %in% year]
        }

        k <- k + 1
      }
    }
  }

  dat <- impact %>% as.data.frame() %>% tibble::rownames_to_column("index") %>% pivot_longer(!index, names_to = "year", values_to = "value") %>% mutate_at(vars(year), as.factor)

  dat <- dat %>% mutate_at(vars(index), as.factor)
  dat$index <- factor(dat$index, levels=rev(levels(dat$index)))

  range <- range[1:(length(range)-1)]

  dat <- subset(dat, year %in% years)

  colours <- c("blue","red")

  if(is.null(max_value)) max_value <- max(abs(dat$value),na.rm=TRUE)

  if (target=="b" | target=="q" | target=="sigma") {
    Title <- paste0("target = ",target, ", index = ",num_id)
    season <- NULL
    threshold <- ((res$opts)[names(res$opts)==paste0(target)][[1]])[num_id]
  } else {
    Title <- paste0("target = ",target, ", season = ",season, ", age = ", age, ", year = ", year)
    threshold <- ((res$outputs)[names(res$outputs)==paste0(target,season)][[1]])[Ages %in% age, Years %in% year]
  }

  ggplot(dat,aes(x=year, y=index)) + geom_point(aes(size=abs(value)), fill=colours[sign(dat$value)*0.5+1.5], color=colours[sign(dat$value)*0.5+1.5], alpha = 0.75, shape = 21) + coord_fixed(ratio=fix_ratio) + scale_size_continuous(limits = c(0.0001, 1.0001)*max_value, range = c(1,max_size), breaks = round(range*max_value/scale,digits=digits))+theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + labs(x="Year", y="Index", size=paste0(target,season),title=Title) + guides(color="none", fill="none")
}

boot_hvpa <- function(
    res,
    B=100,
    seed0=1,
    method="P",
    p.init_rewrite=FALSE){
  assertthat::assert_that(class(res) == "hvpa")

  set.seed(seed0)

  AR <- res$AR
  rho <- res$opts$rho
  logit_rho_se <- res$opts$logit_rho_se

  dat_for_plot <- make_dat(res)

  res.b <- res

  Res_boot <- list()
  nJ <- max(dat_for_plot$num)

  if(p.init_rewrite) res.b$input$p.init <- res$outputs$faa2[1:(nrow(res$outputs$faa2)-1),ncol(res$outputs$faa2)]

  eps_boot <- list()
  rho_boot <- list()

  for (j in 1:nJ){
    if (method=="P") {
      resid <- dat_for_plot$residuals[dat_for_plot$num==j]
      resid[1] <- resid[1]*sqrt(1-rho[j]^2)
      eps_boot[[j]] <- matrix(rnorm(B*length(resid),0.0,sd(resid)),nrow=B)
      if (AR[j]) {
        eps_boot[[j]][,1] <- eps_boot[[j]][,1]*1/sqrt(1-rho[j]^2)
        rho_boot[[j]] <- zc_expit(rnorm(B,zc_logit(rho[j]),logit_rho_se[j]))
      }
    } else {
      eps_boot[[j]] <- matrix(sample(dat_for_plot$residuals[dat_for_plot$num==j],B*length(dat_for_plot$residuals[dat_for_plot$num==j]),replace=TRUE),nrow=B)
    }
  }

  for (i in 1:B){
    for (j in 1:nJ){
      new_index <- NULL
      if (AR[j]){
        new_index <- exp(log(dat_for_plot$pred[dat_for_plot$num==j][1])+eps_boot[[j]][i,1])
        for (k in 2:ncol(eps_boot[[j]])){
          new_index <- c(new_index, exp(log(dat_for_plot$pred[dat_for_plot$num==j][k])+rho_boot[[j]][i]*(log(new_index[k-1]) - log(dat_for_plot$pred[dat_for_plot$num==j][k-1])+eps_boot[[j]][i,k])))
        }

      } else {
        new_index <- exp(log(dat_for_plot$pred[dat_for_plot$num==j])+eps_boot[[j]][i,])
      }

      avail <- which(!is.na(res.b$input$dat$index[res$input$use.index[j],]))
      res.b$input$dat$index[res$input$use.index[j],avail] <- new_index
    }

    res2 <- do.call(vpa.hay,res.b$input)

    Res_boot[[i]] <- res2
  }

  list(res=res, boot=Res_boot)
}

boot_hvpa_plot <- function(
    res_boot,
    target="naa",
    season=1,
    prob=c(0.05,0.95),
    log_scale=FALSE,
    years=2014:2022,
    size=1.5,
    sum_age=FALSE,
    y_label=NULL,
    Title=NULL,
    out="plot"
){
  assertthat::assert_that(class(res) == "hvpa")

  B <- length(res_boot$boot)

  dat_boot <- NULL
  for (i in 1:B){
    res00 <- res_boot$boot[[i]]$outputs
    dat <- res00[names(res00)==paste0(target,season)]

    dat <- dat[[1]] %>% as.data.frame() %>% tibble::rownames_to_column("age") %>% pivot_longer(!age, names_to = "year", values_to = "value")

    dat <- dat %>% mutate_at(vars(age), as.factor)
    dat$age <- factor(dat$age, levels=rev(levels(dat$age)))

    dat_boot <- rbind(dat_boot, cbind(dat, iteration=i))
  }

  if (sum_age) {
    dat_ci <- dat_boot %>% group_by(year, iteration) %>% summarize(sval=sum(value)) %>% group_by(year) %>% summarize(lo=quantile(sval,probs=prob[1]),med=quantile(sval,probs=0.5),up=quantile(sval,probs=prob[2]),mean=mean(sval),sd=sd(sval),cv=sd(sval)/mean(sval))
    dat_ci$year <- as.numeric(dat_ci$year)

    dat_ci <- cbind(dat_ci, obs=colSums(res_boot$res$outputs[names(res_boot$res$outputs)==paste0(target,season)][[1]]))

    pp1 <- NULL
  } else {
    dat_ci <- dat_boot %>% group_by(year, age) %>% summarize(lo=quantile(value,probs=prob[1]),med=quantile(value,probs=0.5),up=quantile(value,probs=prob[2]),mean=mean(value),sd=sd(value),cv=sd(value)/mean(value))
    dat_ci$year <- as.numeric(dat_ci$year)

    dat_ci$age <- factor(dat_ci$age, levels=rev(levels(dat_ci$age)))
    dat_ci <- dat_ci %>% arrange(year, age)

    dat_ci <- cbind(dat_ci, obs=unlist(res_boot$res$outputs[names(res_boot$res$outputs)==paste0(target,season)]))

    age_name <- paste0("age ", dat_ci$age)

    dat_ci$age_name <- age_name

    pp1 <- facet_wrap(~age_name,scale="free")
  }

  if(!is.null(years)) dat_ci <- subset(dat_ci, year %in% years)

  if (is.null(y_label)) y_label <- paste0(if_else(log_scale,"log-",""),target,season)

  if (log_scale) p1 <- ggplot(dat_ci,aes(x=year,y=log(med)))+geom_ribbon(aes(ymin=log(lo),ymax=log(up)),alpha=0.2)+geom_line(size=size)+geom_point(aes(x=year,y=log(obs)), size=size,color="red")+pp1+labs(x="Year", y=y_label, title=Title)+theme_bw() else p1 <- ggplot(dat_ci,aes(x=year,y=med))+geom_ribbon(aes(ymin=lo,ymax=up),alpha=0.2)+geom_line(size=size)+geom_point(aes(x=year,y=obs), color="red")+pp1+labs(x="Year", y=y_label, title=Title)+theme_bw()

  if (out=="dat") list(dat_boot=dat_boot, dat_ci=dat_ci) else print(p1)
}

#' 半期VPA専用のジャックナイフ法の作図関数
#'
#' @export
#' @encoding UTF-8

out_vpa2 <- function(res=NULL, # VPA result
                     rres=NULL, # reference point
                     fres=NULL, # future projection result (not nessesarily)
                     ABC=NULL,
                     filename="vpa" # filename without extension
){
  assertthat::assert_that(class(res) == "hvpa")

  old.par <- par()
  exit.func <- function(){
    #    par(old.par)
    dev.off()
    options(warn=0)
  }
  on.exit(exit.func())

  csvname <- paste(filename,".csv",sep="")
  pdfname <- paste(filename,".pdf",sep="")
  pdf(pdfname)
  par(mfrow=c(3,2),mar=c(3,3,2,1))
  options(warn=-1)

  write.table2 <- function(x,title.tmp="",is.plot=FALSE,...){
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

  write(paste("# RVPA outputs at ",date()," & ",getwd()),file=csvname)

  if(!is.null(res)){
    write("# VPA results",file=csvname, append=T)

    write("\n# catch at age 1",file=csvname,append=T)
    write.table2(res$input$dat$caa1,title.tmp="Catch at age1")

    write("\n# catch at age 2",file=csvname,append=T)
    write.table2(res$input$dat$caa2,title.tmp="Catch at age2")

    write("\n# maturity at age",file=csvname,append=T)
    write.table2(res$input$dat$maa,title.tmp="Maturity at age")

    write("\n# weight at age",file=csvname,append=T)
    write.table2(res$input$dat$waa,title.tmp="Weight at age")

    write("\n# M at age",file=csvname,append=T)
    write.table2(res$input$dat$M,title.tmp="M at age")

    write("\n# fishing mortality at age 1",file=csvname,append=T)
    write.table2(res$outputs$faa1,title.tmp="F at age1")

    write("\n# fishing mortality at age 2",file=csvname,append=T)
    write.table2(res$outputs$faa2,title.tmp="F at age2")

    #    write("\n# Current F",file=csvname,append=T)
    #    write.table2(res$Fc.at.age,title.tmp="Current F")

    write("\n# numbers at age 1",file=csvname,append=T)
    write.table2(res$outputs$naa1,title.tmp="Numbers at age1")

    write("\n# numbers at age 1",file=csvname,append=T)
    write.table2(res$outputs$naa2,title.tmp="Numbers at age1")

    write("\n# total and spawning biomass ",file=csvname,append=T)
    x <- rbind(colSums(res$ssb),colSums(res$outputs$baa1),colSums(res$outputs$baa2),colSums(res$input$dat$caa1*res$input$dat$waa1),colSums(res$input$dat$caa2*res$input$dat$waa2))
    rownames(x) <- c("Spawning biomass","Total biomass 1","Total biomass 2","Catch biomass 1","Catch biomass 2")
    write.table2(x,title.tmp="Total and spawning biomass")
  }

  if(!is.null(rres)){
    write("\n# Reference points",file=csvname,append=T)
    write.table2(rres$summary,title.tmp="Future F at age",is.plot=F)
  }

  if(!is.null(fres)){
    write("\n# future projection results",file=csvname,append=T)
    write("\n# future F at age",file=csvname,append=T)
    write.table2(fres$faa[,,1],title.tmp="Future F at age")

    write("\n# future numbers at age",file=csvname,append=T)
    write.table2(fres$naa[,,1],title.tmp="Future numbers at age")

    write("\n# future total and spawning biomass",file=csvname,append=T)
    x <- rbind(fres$vssb[,1],fres$vbiom[,1],fres$vwcaa[,1])
    rownames(x) <- c("Spawning biomass","Total biomass","Catch biomass")
    write.table2(x,title.tmp="Future total, spawning and catch biomass")
  }

  if(!is.null(ABC)){
    write("\n# ABC summary",file=csvname,append=T)
    write.table2(ABC$ABC,title.tmp="Future F at age",is.plot=F)
    write("\n# Kobe matrix",file=csvname,append=T)
    for(i in 1:dim(ABC$kobe.matrix)[[3]]){
      write(paste("\n# ",dimnames(ABC$kobe.matrix)[[3]][i]),
            file=csvname,append=T)
      write.table2(ABC$kobe.matrix[,,i],
                   title.tmp=dimnames(ABC$kobe.matrix)[[3]][i],is.plot=T)
    }
  }
}

##

run <- FALSE

if (run){
  ## Initial Settings

  source("rvpa_morita_R.r")
  caa1 <- read.csv("caa1.csv",row.names=1)
  caa2 <- read.csv("caa2.csv",row.names=1)
  waa1 <- read.csv("waa1.csv",row.names=1)
  waa2 <- read.csv("waa2.csv",row.names=1)
  waa_s <- read.csv("waa2.csv",row.names=1)
  maa <- read.csv("maa.csv",row.names=1)
  M <- read.csv("M.csv",row.names=1)
  oki1 <- t(read.csv("oki_ratio2_half1.csv",row.names=1))
  oki2 <- t(read.csv("oki_ratio2_half2.csv",row.names=1))
  index <- read.csv("index2021_2.csv",row.names=1)

  dat <- data.handler.hay(caa1, caa2, waa1, waa2, maa, maa_ssb=maa, waa_ssb=waa_s, M=M, omega1=oki1, omega2=oki2, index=index)

  ## Model Run

  res0 <- vpa.hay(dat)
  res_d <- vpa.hay(dat, rec=100)

  p_init <- res0$outputs$faa2[1:4,"2022"]

  source("rvpa1.9.4.r")

  waa <- read.csv("waa.csv",row.names=1)

  oki.ratio <- read.csv("oki_ratio2r.csv")
  omega <- t(oki.ratio[,2:6])

  dat1 <- data.handler(caa1+caa2, waa, maa, maa.tune=maa, waa.catch=NULL, M=M, index=index)

  res1 <- vpa(dat1,tf.year=2018:2019,tune=TRUE,Pope=TRUE,sel.update=TRUE,use.index=c(1),abund=c("Bo"),min.age=c(0),max.age=c(4),omega=omega,fc.year=2017:2019,sel.def="max",alpha=1,est.method="ls",p.init=c(0.8),b.est=TRUE,plot=TRUE,plot.year=1985:2019)

  res1r <- vpa(dat1,tf.year=2018:2019,tune=TRUE,Pope=TRUE,sel.update=TRUE,use.index=c(1),abund=c("Bo"),min.age=c(0),max.age=c(4),omega=omega,fc.year=2017:2019,sel.def="max",alpha=1,est.method="ls",p.init=c(0.8),b.est=TRUE,plot=TRUE,plot.year=1985:2019,ww=1/CVs^2)

  res2 <- vpa.hay(dat,tune=TRUE,use.index=c(1,40),abund=c("Bo","B"),min.age=c(0,0),max.age=c(4,0),ki=c(1,2),p.ini=p_init,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4))

  res290 <- vpa.hay(dat,tune=TRUE,use.index=c(1,46),abund=c("Bo","Bf"),min.age=c(0,1),max.age=c(4,1),ki=c(1,NA),p.ini=p_init,p_slide=1,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4))

  res290w <- vpa.hay(dat,tune=TRUE,use.index=c(1,46),abund=c("Bo","Bf"),min.age=c(0,1),max.age=c(4,1),ki=c(1,NA),p.ini=p_init,p_slide=1,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),waa_lm=TRUE)

  res2901 <- vpa.hay(dat,tune=TRUE,use.index=c(1,46),abund=c("Bo","Bf"),min.age=c(0,1),max.age=c(4,1),ki=c(1,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=c(FALSE,TRUE))

  res2901r <- vpa.hay(dat,tune=TRUE,use.index=c(1,46),abund=c("Bo","Bf"),min.age=c(0,1),max.age=c(4,1),ki=c(1,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.3,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res2900 <- vpa.hay(dat,tune=TRUE,use.index=c(1,2),abund=c("Bo","B"),min.age=c(0,0),max.age=c(4,0),ki=c(1,1),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res2902 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,54,55),abund=c("Bo","Bo","Bf","B"),min.age=c(0,0,1,1),max.age=c(4,4,1,1),ki=c(1,2,NA,2),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=FALSE)

  res2903 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,54,55),abund=c("Bo","Bo","Bf","B"),min.age=c(0,0,1,1),max.age=c(4,4,1,1),ki=c(1,2,NA,2),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res2904 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,46),abund=c("Bo","Bo","Bf"),min.age=c(0,0,1),max.age=c(4,4,1),ki=c(1,2,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=FALSE)

  res2905 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,46),abund=c("Bo","Bo","Bf"),min.age=c(0,0,1),max.age=c(4,4,1),ki=c(1,2,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res2905r <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,46),abund=c("Bo","Bo","Bf"),min.age=c(0,0,1),max.age=c(4,4,1),ki=c(1,2,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.5,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res2905rw <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,46),abund=c("Bo","Bo","Bf"),min.age=c(0,0,1),max.age=c(4,4,1),ki=c(1,2,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.5,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE,waa_lm=TRUE)

  res2906 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,2),abund=c("Bo","Bo","B"),min.age=c(0,0,0),max.age=c(4,4,0),ki=c(1,2,1),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res2907 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,2),abund=c("Bo","Bo","B"),min.age=c(0,0,0),max.age=c(4,4,0),ki=c(1,2,1),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=FALSE)

  res2908 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,52),abund=c("Bo","Bo","Bf"),min.age=c(0,0,1),max.age=c(4,4,1),ki=c(1,2,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=FALSE)

  res2909 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,52),abund=c("Bo","Bo","Bf"),min.age=c(0,0,1),max.age=c(4,4,1),ki=c(1,2,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res2903 <- vpa.hay(dat,tune=TRUE,use.index=c(3,4,54,55),abund=c("Bo","Bo","Bf","B"),min.age=c(0,0,1,1),max.age=c(4,4,1,1),ki=c(1,2,NA,2),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  res290D <- vpa.hay(dat,tune=TRUE,use.index=c(1,2),abund=c("Bo","B"),min.age=c(0,0),max.age=c(4,0),ki=c(1,1),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4))

  res290C <- vpa.hay(dat,tune=TRUE,use.index=c(1,2,46),abund=c("Bo","B","Bf"),min.age=c(0,0,1),max.age=c(4,0,1),ki=c(1,1,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4))

  res290Ca <- vpa.hay(dat,tune=TRUE,use.index=c(1,2,46),abund=c("Bo","B","Bf"),min.age=c(0,0,1),max.age=c(4,0,1),ki=c(1,1,NA),p.ini=p_init,p_slide=0,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),AR=TRUE)

  ###

  res290z0 <- vpa.hay(dat,tune=TRUE,use.index=c(1,46),abund=c("Bo","Bf"),min.age=c(0,1),max.age=c(4,1),ki=c(2,NA),p.ini=p_init,p_slide=1,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=1000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),eta=1,rec_cor=2016:2020)

  res290z1 <- vpa.hay(dat,tune=TRUE,use.index=c(4,5),abund=c("Bo","Bf"),min.age=c(0,1),index.w=c(1,1),max.age=c(4,1),ki=c(2,NA),p.ini=p_init,p_slide=1,no.est=FALSE,b_est=TRUE,lambda=0.0,sel.const=TRUE,start_row=2,W_s=50000,optimizer="nlminb",Lower=c(log(0.0001),rep(-3,3)),Upper=rep(3,4),eta=10,rec_cor=2016:2020)

  res <- res290

  sel_const_check(res)   # selectivity check

  out_vpa2(res)

  res$outputs$naa1
  res$outputs$naa2
  res$outputs$faa1
  res$outputs$faa2
  colSums(res$outputs$ssb,na.rm=TRUE)

  ## Retrospective Analysis

  Res <- list()
  res_retro_tb <- NULL
  k <- 1
  lambda_range <- seq(0.0,0.9,by=0.1)
  for (lambda in lambda_range){
    dat1 <- res$input
    dat1$lambda <- lambda

    Res[[k]] <- do.call(vpa.hay, dat1)
    res_retro_tb <- rbind(res_retro_tb, hretro_est(Res[[k]], n=4, b_fix=TRUE)$mohn)
    k <- k + 1
  }

  rownames(res_retro_tb) <- lambda_range

  res_retro0 <- hretro_est(res0, n=5, b_fix=TRUE)
  res_retro2 <- hretro_est(res, n=5, b_fix=TRUE)

  retro_res <- rbind(res_retro0, res_retro2)
  rownames(retro_res) <- c("no_tuned","tuned")

  res_z <- zenki(res)
  res_z0 <- zenki(res0)

  ## Comparison Plot

  hvpa_plot(res, res0, years=2014:2022)
  write.csv(res_z$faa,"res_z_faa.csv")
  write.csv(res_z$naa,"res_z_naa.csv")

  ## Retrospective Pattern

  res <- res2905r
  res_retro <- res_retro2905r

  p0 <- hretro_plot(res, res_retro, target="faa",season=2,age=0)
  p1 <- hretro_plot(res, res_retro, target="faa",season=2,age=1)
  p2 <- hretro_plot(res, res_retro, target="faa",season=2,age=2)
  p3 <- hretro_plot(res, res_retro, target="faa",season=2,age=3:4)

  cowplot::plot_grid(p0, p1, p2, p3, nrow=2)

  hretro_plot(res, res_retro, target="ssb",season=NULL,age=0:4)

  p0 <- hretro_plot(res, res_retro, target="naa",season=1,age=0)
  p1 <- hretro_plot(res, res_retro, target="naa",season=1,age=1)
  p2 <- hretro_plot(res, res_retro, target="naa",season=1,age=2)
  p3 <- hretro_plot(res, res_retro, target="naa",season=1,age=3:4)

  cowplot::plot_grid(p0, p1, p2, p3, nrow=2)

  ## Bubble Plot

  q1 <- bubble_plot(res,target="naa",season=1,digits=0,max_size=6,fix_ratio=2)
  q2 <- bubble_plot(res,target="naa",season=2,digits=0,max_size=6,fix_ratio=2)
  q3 <- bubble_plot(res,target="uaa",season=1,digits=1,max_size=6,fix_ratio=2)
  q4 <- bubble_plot(res,target="uaa",season=2,digits=1,max_size=6,fix_ratio=2)
  cowplot::plot_grid(q1,q2,q3,q4,nrow=2)

  ## Residual Plot

  residual_plot(res,add_log_plot=TRUE,rel_heights=c(0.35,0.35,0.3))
  d1 <- residual_plot(res,add_log_plot=TRUE,out="dat")

  colours <- c("blue","red")
  d1$num <- factor(d1$num)
  ggplot(d1,aes(x=year,y=num))+geom_point(aes(size=abs(residuals)),color=colours[sign(d1$residuals)*0.5+1.5], fill=colours[sign(d1$residuals)*0.5+1.5], alpha = 0.75, shape = 21)+coord_fixed(ratio=2) + scale_size_continuous(limits = c(10^(-10), 1.0001)*max(d1$residuals), range = c(1,20), breaks = round(c(0,25,50,75)*max(d1$residuals)/(100),digits=2))+labs(x="Year", y="Index",size="Residuals")+theme_bw()

  ## Future Prediction (Deterministic)

  library(frasyr)
  SR <- fit.SR(SRdata=list(year=colnames(res$outputs$naa1),R=as.numeric(res$outputs$naa1[1,]),SSB=colSums(res$outputs$ssb)),SR = "HS",method = "L1")

  Y <- 20
  dat_for_future <- NULL

  for(beta in seq(0,1,by=0.1)){
    res_f <- future_projection(res,SR,beta=beta,Y=Y)

    dat_for_future <- rbind(dat_for_future, data.frame(beta=beta,type=rep(c("hanki","nenki"),each=Y),naa=c(colSums(res_f$naa1_f),colSums(res_f$naaz_f)),caa=c(colSums(res_f$caa12_f),colSums(res_f$caaz_f)),years=c(colnames(res_f$naa1_f),colnames(res_f$naaz_f))))
  }

  dat_for_future$year <- factor(dat_for_future$year)

  ggplot(dat_for_future,aes(x=years,y=naa,color=type))+geom_point(position = position_dodge(width = .9))+facet_wrap(~beta)+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggplot(dat_for_future,aes(x=years,y=caa,color=type))+geom_point(position = position_dodge(width = .9))+facet_wrap(~beta)+theme_bw()+theme(axis.text.x = element_text(angle=90, hjust=1))

  ## Sensitivity Test

  res <- res2905r

  jack_test1 <- jackknife_test(res)

  jackknife_plot(jack_test1,targe="b",num_id=3)

  sen_test1 <- sensitivity_test(res,SD=0.2)
  sen_test2 <- sensitivity_test(res2,SD=0.2)

  x1 <- sensitivity_plot(sen_test1,year=2020,digits=2,max_value=0.65,range=seq(0,100,len=8),max_size=10)
  x1 <- sensitivity_plot(sen_test1,target="b",year=2020,digits=2,max_value=0.65,range=seq(0,100,len=8),max_size=10)
  x2 <- sensitivity_plot(sen_test2,year=2020,digits=2,max_value=0.65,range=seq(0,100,len=8),max_size=10)
  cowplot::plot_grid(x1,x2,nrow=2)

  ## Bootstrap

  # res_boot2905 <- boot_hvpa(res2905, B=200)
  res_boot2905r <- boot_hvpa(res2905r, B=200)
  res_boot2903 <- boot_hvpa(res2903, B=200)
  res_boot2 <- boot_hvpa(res2, B=20)

  res_boot2904 <- boot_hvpa(res2904, B=200)
  res_boot2905 <- boot_hvpa(res2905, B=200)

  bb1 <- boot_hvpa_plot(res_boot1)
  bb2 <- boot_hvpa_plot(res_boot2)
  cowplot::plot_grid(bb1, bb2)
}


#res_boot2901 <- boot_hvpa(res2901, B=200)
#res_boot290 <- boot_hvpa(res290, B=200)
