#' csvデータを読み込んでvpa用のデータを作成する関数
#'
#' @param caa catch at age
#' @param waa weight at age
#' @param maa maturity at age
#' @param maa.tune チューニングに使うmaturity at ageが異なる場合
#' @param waa.catch abundanceとcatchでwaaが異なる場合
#' @encoding UTF-8
#'
#' @export

data.handler <- function(
  caa,
  waa,
  maa,
  index=NULL,
  M = 0.4,
  maa.tune=NULL,
  waa.catch=NULL,
  catch.prop=NULL,
  release.dat=NULL
)
{
  years <- as.numeric(sapply(strsplit(names(caa[1,]),"X"), function(x) x[2]))

  if (is.null(dim(waa)) | dim(waa)[2]==1) waa <- as.data.frame(matrix(unlist(waa), nrow=nrow(caa), ncol=ncol(caa)))

  if (is.null(dim(maa)) | dim(maa)[2]==1) maa <- as.data.frame(matrix(unlist(maa), nrow=nrow(caa), ncol=ncol(caa)))

  colnames(caa) <- colnames(waa) <- colnames(maa) <- years

  if (!is.null(waa.catch)) {
    if (is.null(dim(waa.catch)) | dim(waa.catch)[2]==1) waa.catch <- as.data.frame(matrix(unlist(waa.catch), nrow=nrow(caa), ncol=ncol(caa)))
    colnames(waa.catch) <- years
    assertthat::assert_that(
      all(rownames(caa) == rownames(waa.catch)),
      all(colnames(caa) == colnames(waa.catch))
    )
  }

  if (!is.null(maa.tune)) {
    if (is.null(dim(maa.tune)) | dim(maa.tune)[2]==1) maa.tune <- as.data.frame(matrix(unlist(maa.tune), nrow=nrow(caa), ncol=ncol(caa)))
    colnames(maa.tune) <- years
    assertthat::assert_that(
      all(rownames(caa) == rownames(maa.tune)),
      all(colnames(caa) == colnames(maa.tune))
    )
  }

  if (!is.null(catch.prop)) colnames(catch.prop) <- years

  if (!is.null(index)) colnames(index) <- years

  if (!is.null(release.dat)) colnames(release.dat) <- years

  if (is.null(dim(M))) M <- as.data.frame(matrix(M, nrow=nrow(caa), ncol=ncol(caa)))

  colnames(M) <- years
  rownames(M) <- rownames(caa)

  assertthat::assert_that(
    all(
      rownames(caa) == purrr::flatten_chr(purrr::map(list(maa,
                                                          waa,
                                                          M),
                                                     rownames))
    ),
    all(
      colnames(caa) == purrr::flatten_chr(purrr::map(list(maa,
                                                          waa,
                                                          M),
                                                     colnames))
    )
  )

  res <- list(caa=caa, maa=maa, waa=waa, index=index, M=M, maa.tune=maa.tune, waa.catch=waa.catch, catch.prop=catch.prop, release.dat=release.dat)

  invisible(res)
}

# miscellaneous functions

#max.age <- function(x) max(which(!is.na(x)))
max.age.func <- function(x) max(which(!is.na(x)))

vpa.core <- function(caa,faa,M,k){
  out <- caa[,k]/(1-exp(-faa[,k]-M[,k]))*(faa[,k]+M[,k])/faa[,k]
  return(out)
}

vpa.core.Pope <- function(caa,faa,M,k,p=0.5){
  out <- caa[, k]*exp(p*M[, k])/(1-exp(-faa[, k]))
  return(out)
}

ik.est <- function(caa,naa,M,i,k,min.caa=0.01,maxit=5,d=0.0001){
  K <- 1
  it <- 0

  f0 <- 1
  f1 <- NA

  if (!is.na(naa[i+1,k+1])){
    while(it < maxit & K > d){
      it <- it + 1
      f1 <- log(1+max(caa[i,k],min.caa)/naa[i+1,k+1]*exp(-M[i,k])*(f0+M[i,k])*(1-exp(-f0))/(f0*(1-exp(-f0-M[i,k]))))
      if(f1==0) stop("Fの推定値が0になり計算不能です。初期値p.initを高めに設定して計算しなおしてみてください")
      K <- sqrt((f1-f0)^2)
      f0 <- f1
    }
  }

  return(f1)
}

hira.est <- function(caa,naa,M,i,k,alpha=1,min.caa=0.01,maxit=5,d=0.0001){
  K <- 1
  it <- 0

  f0 <- 1
  f1 <- NA

  if (!is.na(naa[i+1,k+1])){
    while(it < maxit & K > d){
      it <- it + 1
      f1 <- log(1+(1-exp(-f0))*exp(-M[i,k])/(naa[i+1,k+1]*f0)*(max(caa[i+1,k],min.caa)*(alpha*f0+M[i+1,k])/(alpha*(1-exp(-alpha*f0-M[i+1,k])))*exp((1-alpha)*f0)+max(caa[i,k],min.caa)*(f0+M[i,k])/(1-exp(-f0-M[i,k]))))
	   if(f1==0) stop("Fの推定値が0になり計算不能です。初期値p.initを高めに設定して計算しなおしてみてください")
      K <- sqrt((f1-f0)^2)
      f0 <- f1
    }
  }

  return(f1)
}

f.forward.est <- function(caa,naa,M,i,k,maxit=5,d=0.0001){
  K <- 1
  it <- 0

  f0 <- f1 <- 1

  while(it < maxit & K > d){
    it <- it + 1
    f1 <- caa[i,k]/naa[i,k]*(f0+M[i,k])/f0*1/(1-exp(-f0-M[i,k]))
    K <- sqrt((f1-f0)^2)
    f0 <- f1
  }

  return(f1)
}

fp.forward.est <- function(caa,naa,M,i,k,alpha=1,maxit=5,d=0.0001){
  K <- 1
  it <- 0

  f0 <- f1 <- 1

  while(it < maxit & K > d){
    it <- it + 1
    f1 <- 1/(1+alpha)*(caa[i,k]/naa[i,k]*(f0+M[i,k])*1/(1-exp(-f0-M[i,k]))+caa[i+1,k]/naa[i+1,k]*(alpha*f0+M[i+1,k])*1/(1-exp(-alpha0*f0-M[i+1,k])))
    K <- sqrt((f1-f0)^2)
    f0 <- f1
  }

  return(f1)
}

#' @export

backward.calc <- function(caa,naa,M,na,k,min.caa=0.001,p=0.5,plus.group=TRUE,sel.update,alpha,use.equ){
  out <- rep(NA, na[k])
  if(na[k+1] > na[k]){
    if(isTRUE(sel.update)){stop("Selectivity update method is currently not supported for the plus group change scenario")}
    if (isTRUE(plus.group) & use.equ=="new"){
      for (i in 1:(na[k]-2)){
        out[i] <- naa[i+1,k+1]*exp(M[i,k])+caa[i,k]*exp(p*M[i,k])
      }
      out[na[k]-1]<- (caa[na[k]-1,k]*alpha  * (naa[na[k],k+1]+naa[na[k+1],k+1]) * exp(M[na[k]-1,k]))/(caa[na[k]-1,k]*alpha +caa[na[k],k]) + caa[na[k]-1,k] * exp(p * M[na[k]-1,k])
      out[na[k]]  <- (caa[na[k],k] * (naa[na[k],k+1]+naa[na[k+1],k+1]) * exp(M[na[k],k]))/(caa[na[k]-1,k]*alpha +caa[na[k],k]) + caa[na[k],k] * exp(p * M[na[k],k])
    }
    else if (isTRUE(plus.group) & use.equ=="old"){
      for (i in 1:(na[k])){
        out[i] <- naa[i+1,k+1]*exp(M[i,k])+caa[i,k]*exp(p*M[i,k])
      }
      #out[na[k]-1]<- pmax(caa[na[k]-1,k],min.caa)*alpha/(pmax(caa[na[k]-1,k],min.caa)*alpha +pmax(caa[na[k],k],min.caa))* naa[na[k],k+1] * exp(M[na[k]-1,k])+ caa[na[k]-1,k] * exp(p * M[na[k]-1,k])
      #out[na[k]]<- pmax(caa[na[k],k], min.caa)/(pmax(caa[na[k]-1,k], min.caa)*alpha + pmax(caa[na[k],k], min.caa)) * naa[na[k], k+1] * exp(M[na[k],k]) + caa[na[k], k] * exp(p * M[na[k],k])
    }
    else{
      for (i in 1:(na[k]-2)){
        out[i] <- naa[i+1,k+1]*exp(M[i,k])+caa[i,k]*exp(p*M[i,k])
      }
      out[na[k]-1] <- naa[na[k],k+1]*exp(M[na[k]-1,k])+caa[na[k]-1,k]*exp(p*M[na[k]-1,k])
      out[na[k]] <- out[na[k]-1]*caa[na[k],k]/caa[na[k]-1,k]*exp(p*(M[na[k],k]-M[na[k]-1,k]))
    }
  }
  else{
    for (i in 1:(na[k+1]-2)){
      out[i] <- naa[i+1,k+1]*exp(M[i,k])+caa[i,k]*exp(p*M[i,k])
    }
    if (isTRUE(plus.group)){
      if(na[k]-na[k+1]<=0){
        pp <- c(1, 1/alpha)
        out[(na[k+1]-1):na[k]] <- pp*pmax(caa[(na[k+1]-1):na[k],k],min.caa)/sum(pp*pmax(caa[(na[k+1]-1):na[k],k],min.caa))*naa[na[k+1],k+1]*exp(M[(na[k+1]-1):na[k],k])+caa[(na[k+1]-1):na[k],k]*exp(p*M[(na[k+1]-1):na[k],k])
      }
      else{
        pp <- c(1, 1, 1/alpha)
        out[(na[k+1]-1):na[k]] <- pp*pmax(caa[(na[k+1]-1):na[k],k],min.caa)/sum(pp*pmax(caa[(na[k+1]-1):na[k],k],min.caa))*naa[na[k+1],k+1]*exp(M[(na[k+1]-1):na[k],k])+caa[(na[k+1]-1):na[k],k]*exp(p*M[(na[k+1]-1):na[k],k])
      }
    }
    else{
      out[na[k+1]-1] <- naa[na[k+1],k+1]*exp(M[na[k+1]-1,k])+caa[na[k+1]-1,k]*exp(p*M[na[k+1]-1,k])
      out[na[k]] <- out[na[k+1]-1]*caa[na[k+1],k]/caa[na[k+1]-1,k]*exp(p*(M[na[k+1],k]-M[na[k+1]-1,k]))
    }
  }
  return(out)
}


forward.calc <- function(faa,naa,M,na,k,plus.group=plus.group){
  out <- rep(NA, na[k])
  for (i in 2:(na[k]-1)){
    out[i] <- naa[i-1,k-1]*exp(-faa[i-1,k-1]-M[i-1,k-1]) #最終年の最低年齢より上，最高年齢より下のnaaを計算してる
  }
  if (isTRUE(plus.group)){
    out[na[k]] <- sum(sapply(seq(na[k]-1,max(na[k], na[k-1])), plus.group.eq, naa=naa, faa=faa, M=M, k=k)) #最終年最高年齢のnaa
  }
  else
  {
    out[na[k]] <- sapply(na[k]-1, plus.group.eq, naa=naa, faa=faa, M=M, k=k)
  }


  return(out)
}

plus.group.eq <- function(x, naa, faa, M, k) naa[x,k-1]*exp(-faa[x,k-1]-M[x,k-1])

f.at.age <- function(caa,naa,M,na,k,p=0.5,alpha=1,use.equ) {
 if(na[k+1]>na[k]){
 if (use.equ=="new"){
 out <- -log(1-caa[1:(na[k]-1),k]*exp(p*M[1:(na[k]-1),k])/naa[1:(na[k]-1),k])
  c(out, alpha*out[length(out)])
 }
 else
 {
  out <- -log(1-caa[1:(na[k]),k]*exp(p*M[1:(na[k]),k])/naa[1:(na[k]),k])
   c(out)
 }
 }
 else
 {
  out <- -log(1-caa[1:(na[k]-1),k]*exp(p*M[1:(na[k]-1),k])/naa[1:(na[k]-1),k])
  c(out, alpha*out[length(out)])

  }
}

sel.func <- function(faa, def="maxage") {
  if(def=="maxage") saa <- apply(faa, 2, function(x) x/x[length(x[!is.na(x)])])
  if(def=="max") saa <- apply(faa, 2, function(x) x/max(x,na.rm=TRUE))
  if(def=="mean") saa <- apply(faa, 2, function(x) x/sum(x,na.rm=TRUE))

  return(saa)
}

ff <- function(x, z) get(x)(z)

abund.extractor <- function(
  abund="SSB",
  naa,
  faa,
  dat,
  min.age=0,
  max.age=0,
  link="id",
  base=exp(1),
  af=1,
  catch.prop=NULL,
  sel.def="maxage",
  p.m=0.5,
  omega=NULL,
  scale=1000
){
# abund = "N": abundance
# abund = "Nm": abundance at the middle of the year
# abund = "B": biomass
# abund = "Bm": abundance at the middle of the year
# abund = "SSB": SSB
# # abund = "SSB": SSB at the middle of the year

  naa <- as.data.frame(naa)
  faa <- as.data.frame(faa)

  waa <- dat$waa/scale
  maa <- dat$maa

  min.age <- min.age + 1
  max.age <- max.age + 1

  maa.tune <- dat$maa.tune

 if (abund=="N") res <- colSums(naa[min.age:max.age,], na.rm=TRUE)
 if (abund=="Nm") res <- colSums(naa[min.age:max.age,]*exp(-p.m*dat$M[min.age:max.age,]-p.m*af*faa[min.age:max.age,]), na.rm=TRUE)
 if (abund=="B") res <- colSums((naa*waa)[min.age:max.age,], na.rm=TRUE)
 if (abund=="Bm") res <- colSums((naa*waa)[min.age:max.age,]*exp(-p.m*dat$M[min.age:max.age,]-p.m*af*faa[min.age:max.age,]), na.rm=TRUE)
 if (abund=="SSB"){
   if (is.null(maa.tune)) ssb <- naa*waa*maa else ssb <- naa*waa*maa.tune
   res <- colSums(ssb,na.rm=TRUE)
 }
 if (abund=="Bs"){
       saa <- sel.func(faa, def=sel.def)
       res <- colSums((naa*waa*saa)[min.age:max.age,], na.rm=TRUE)
 }
 if (abund=="Bo"){
        saa <- sel.func(faa*omega, def=sel.def)
        saa <- sweep(saa,2,colSums(saa),FUN="/")
       res <- colSums((naa*waa*saa)[min.age:max.age,], na.rm=TRUE)
 }
 if (abund=="Ns"){
       saa <- sel.func(faa, def=sel.def)
       res <- colSums((naa*saa)[min.age:max.age,], na.rm=TRUE)
 }
 if (abund=="SSBm"){
   if (is.null(maa.tune)) ssb <- naa*waa*maa*exp(-p.m*dat$M-p.m*af*faa) else ssb <- naa*waa*maa.tune*exp(-p.m*dat$M-p.m*af*faa)
   res <- colSums(ssb,na.rm=TRUE)
 }
 if (abund=="SSBmsj"){
   if (is.null(maa.tune)) ssb <- naa*waa*maa*exp(-p.m*dat$M) else ssb <- naa*waa*maa.tune*exp(-p.m*dat$M)
   res <- colSums(ssb,na.rm=TRUE)
 }

 if (abund=="N1sj") res <- colSums(naa[1,]*exp(dat$M[1,]), na.rm=TRUE)
 if (abund=="N0sj") res <- colSums(naa[1,]*exp(dat$M[1,]*2), na.rm=TRUE)

 if (abund=="F") if (is.null(catch.prop)) res <- colMeans(faa[min.age:max.age,], na.rm=TRUE) else res <- colMeans(catch.prop[min.age:max.age, ]*faa[min.age:max.age,], na.rm=TRUE)

 if (link=="log") res <- log(res, base=base)

  return(invisible(res))
}

#

tmpfunc2 <- function(x=1,y=2,z=3){
  argname <- ls()  # 関数が呼び出されたばかりのときのls()は引数のみが入っている
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  value <- x+y+z
  return(list(value=value,args=arglist))
}

#
##


qbs.f <- function(q.const, b.const, sigma.constraint, index, Abund, nindex, index.w, max.dd=0.0001, max.iter=100){

  np.q <- length(unique(q.const[q.const > 0]))
  np.b <- length(unique(b.const[b.const > 0]))
  np.s <- length(unique(sigma.constraint[sigma.constraint > 0]))

  q <- b <- sigma <- numeric(nindex)

  q[1:nindex] <- b[1:nindex] <- sigma[1:nindex] <- 1

  delta <- 1

  obj <- NULL

  NN <- 0

  while(delta > max.dd & NN < max.iter){
    NN <- NN+1
    q0 <- q
    b0 <- b
    sigma0 <- sigma

    if (np.q > 0){
      for(i in 1:np.q){
        id <- which(q.const==i)
        num <- den <- 0
        for (j in id){
          avail <- which(!is.na(as.numeric(index[j,])))
          num <- num+index.w[j]*mean(log(as.numeric(index[j,avail]))-b[j]*log(as.numeric(Abund[j,avail])))/sigma[j]^2
          den <- den+index.w[j]/sigma[j]^2
        }
        q[i] <- num/den
      }
    }
    if (np.b > 0){
      for(i in 1:np.b){
        id <- which(b.const==i)
        num <- den <- 0
        for (j in id){
          avail <- which(!is.na(as.numeric(index[j,])))
          num <- num+index.w[j]*cov(log(as.numeric(index[j,avail])),log(as.numeric(Abund[j,avail])))/var(log(as.numeric(Abund[j,avail])))/sigma[j]^2
          den <- den+index.w[j]/sigma[j]^2
        }
        b[i] <- num/den
      }
    }
    if (np.s > 0){
      for(i in 1:np.s){
        id <- which(sigma.constraint==i)
        num <- den <- 0
        for (j in id){
          avail <- which(!is.na(as.numeric(index[j,])))
          nn <- length(avail)
          num <- num+index.w[j]*sum((log(as.numeric(index[j,avail]))-q[j]-b[j]*log(as.numeric(Abund[j,avail])))^2)
          den <- den+index.w[j]*nn
        }
        sigma[i] <- sqrt(num/den)
      }
    }

    q[which(q.const>0)] <- q[q.const[which(q.const>0)]]
    b[which(b.const>0)] <- b[b.const[which(b.const>0)]]

    sigma[which(sigma.constraint>0)] <- sigma[sigma.constraint[which(sigma.constraint>0)]]

    delta <- max(c(sqrt((q-q0)^2),sqrt((b-b0)^2),sqrt((sigma-sigma0)^2)))
  }

  for (i in 1:nindex){
    avail <- which(!is.na(as.numeric(index[i,])))
    obj <- c(obj, index.w[i]*(-as.numeric(na.omit(dnorm(log(as.numeric(index[i,avail])),log(q[i])+b[i]*log(as.numeric(Abund[i,avail])),sigma[i],log=TRUE)))))
  }

  convergence <- ifelse(delta <= max.dd, 1, 0)

  return(list(q=q, b=b, sigma=sigma, obj=sum(obj), convergence=convergence))
}

qbs.f2 <- function(p0,index, Abund, nindex, index.w, fixed.index.var=NULL){

  if (is.null(fixed.index.var)) fixed.index.var <- matrix(0, nrow=nrow(index), ncol=ncol(index))
  if (class(fixed.index.var)=="numeric") fixed.index.var <- matrix(fixed.index.var, nrow=1)
  if (class(fixed.index.var)=="matrix" | class(fixed.index.var)=="data.frame") fixed.index.var <- array(fixed.index.var, dim=c(dim(fixed.index.var),1))

  np <- min(nindex,sum(index.w>0))

  p <- vector(length=2*np)
  q <- sigma <- vector(length=np)

  obj.f <- function(p){
    obj <- j <- 0

    for (i in 1:nindex){
      if (index.w[i] > 0 ){
      j <- j + 1
      q[j] <- exp(p[2*j-1])
      sigma[j] <- exp(p[2*j])

      avail <- which(!is.na(as.numeric(index[i,])))
      obj <- obj+index.w[i]*(-as.numeric(na.omit(dmvnorm(log(as.numeric(index[i,avail])),log(q[j])+log(as.numeric(Abund[i,avail])),as.matrix(fixed.index.var[avail,avail,j])+sigma[j]^2*diag(length(avail)),log=TRUE))))
      }
    }

    sum(obj)
  }

  res <- nlm(obj.f,p0)

  q <- res$estimate[1:np]
  sigma <- exp(res$estimate[1:np+np])
  obj <- res$minimum
  convergence <- res$code

   return(list(q=q, b=rep(1,np), sigma=sigma, obj=obj, convergence=convergence))
}


#' VPAによる資源計算を実施する
#'
#' @param dat    data for vpa
#' @param sel.f  最終年の選択率
#' @param tf.year terminal Fをどの年の平均にするか
#' @param rec.new 翌年の加入を外から与える
#' @param rec     rec.yearの加入
#' @param rec.year 加入を代入する際の年
#' @param rps.year 翌年のRPSをどの範囲の平均にするか
#' @param fc.year   Fcurrentでどの範囲を参照するか
#' @param last.year vpaを計算する最終年を指定（retrospective analysis）
#' @param last.catch.zero TRUEなら強制的に最終年の漁獲量を0にする
#' @param faa0  sel.update=TRUEのとき，初期値となるfaa
#' @param naa0  sel.update=TRUEのとき，初期値となるnaa
#' @param f.new
#' @param Pope  Popeの近似式を使うかどうか
#' @param p.pope  Popeの式でどこで漁獲するか（0.5 = 真ん中の時期）
#' @param tune    tuningをするかどうか
#' @param abund  tuningの際，何の指標に対応するか. "N"=年の初めの資源尾数，"Nm"=年の中間での資源尾数，"B"=年の初めの資源重量，"Bm"=年の中間での資源重量，"SSB"=産卵親魚量，"Bs"=資源重量×選択率，"Bo"=資源重量×オメガで調節された選択率，"Ns"=資源尾数×選択率，"SSBm"=年の中間での産卵親魚量，"N1sj"=1歳の資源尾数（日本海スケトウダラ用）,"N0sj"=0歳の資源尾数（日本海スケトウダラ用），"F"=漁獲係数の平均
#' @param min.age tuning指標の年齢参照範囲の下限
#' @param max.age tuning指標の年齢参照範囲の上限
#' @param link    tuningのlink関数
#' @param base    link関数が"log"のとき，底を何にするか
#' @param af      資源量指数が年の中央のとき，af=0なら漁期前，af=1なら漁期真ん中，af=2なら漁期後となる
#' @param p.m
#' @param index.w  tuning indexの重み
#' @param use.index   indexのなにを使うか．c(1,3)というような与え方ができる.デフォルトは"all"
#' @param scale  重量のscaling
#' @param hessian
#' @param alpha  最高齢と最高齢-1のFの比 F_a = alpha*F_{a-1}
#' @param maxit  石岡・岸田/平松の方法の最大繰り返し数
#' @param d  石岡・岸田/平松の方法の収束判定基準
#' @param min.caa  caaに0があるとき，0をmin.caaで置き換える
#' @param plot  tuningに使った資源量指数に対するフィットのプロット
#' @param plot.year  上のプロットの参照年
#' @param term.F  terminal Fの何を推定するか: "max" or "all"
#' @param plus.group
#' @param stat.tf  最終年のFを推定する統計量（年齢で同じとする：changed on Nov/2021）
#' @param add.p.est  追加で最高齢以外のfaaを推定する際．年齢を指定する．
#' @param add.p.ini
#' @param sel.update チューニングVPAにおいて，選択率を更新しながら推定
#' @param sel.def   sel.update=TRUEで選択率を更新していく際に，選択率をどのように計算するか．"max"=選択率が一番大きい年齢の選択率を１とする，"maxage"=最高齢の選択率を１とする，"mean"=全体に対する割合として選択率を決める（saa=faa/sum(faa))
#' @param max.dd  sel.updateの際の収束判定基準
#' @param ti.scale  資源量の係数と切片のscaling
#' @param tf.mat terminal Fの平均をとる年の設定．NA行列に平均をとる箇所に1を入れる．
#' @param eq.tf.mean terminal Fの平均値を過去のFの平均値と等しくする
#' @param no.est  パラメータ推定しない．
#' @param est.method 推定方法 （ls = 最小二乗法，ml = 最尤法,  ls_nolog =最小二乗法で実数)
#' @param b.est bを推定するかどうか
#' @param est.constraint  制約付き推定をするかどうか
#' @param q.const  qパラメータの制約（0は推定しないで1にfix）
#' @param b.const  bパラメータの制約（0は推定しないで1にfix）
#' @param q.fix
#' @param b.fix
#' @param fixed.index.var
#' @param max.iter  q,b,sigma計算の際の最大繰り返し数
#' @param optimizer
#' @param Lower
#' @param Upper
#' @param p.fix
#' @param lambda  ridge penaltyの大きさ
#' @param beta  penaltyのexponent  (beta = 1: lasso, 2: ridge)
#' @param penalty  ridgeVPAの際に与えるpenaltyの種類．"p"=指定した年齢範囲の最終年の年齢別Fのbeta乗の和，"f"=｛最終年の年齢aのF－（tf.yearで指定した年のa歳の平均のF)｝のbeta乗の和，"s"=｛最終年の年齢aのs－（tf.yearで指定した年のa歳の平均のs)｝のbeta乗の和
#' @param ssb.def  i: 年はじめ，m: 年中央, l: 年最後
#' @param ssb.lag  0: no lag, 1: lag 1
#' @param TMB  TMBで高速計算する場合TMB=TRUE (事前にuse_rvpa_tmb()を実行)　全Ｆ推定法，POPE=TRUE, alpha=1, 途中でプラスグループが変化しない場合のみのときに使用可能
#' @param sel.rank
#' @param p.init
#' @param sigma.constraint  sigmaパラメータの制約.使い方としては，指標が５つあり，2番目と3番目の指標のsigmaは同じとしたい場合はc(1,2,2,3,4)と指定する
#' @param eta  Fのpenaltyを分けて与えるときにeta.ageで指定した年齢への相対的なpenalty (0~1)
#' @param eta.age  Fのpenaltyを分けるときにetaを与える年齢(0 = 0歳（加入）,0:1 = 0~1歳)
#' @param tmb.file  TMB=TRUEのとき使用するcppファイルの名前
#' @param remove.abund ある値を引いた値に対してチューニングに使用する（トラフグ伊勢・三河湾系群で放流魚を引いて天然魚のみに対してチューニングするため）.
#' @param madara  マダラ太平洋系群で用いているチューニングのやり方
#' @param penalty_age  選択率更新法でridgeVPAをする際のpenalty="p"のときで，etaがNULLで，p_by_age=TRUE (年齢別にペナルテイーを与えたい）ときにペナルテイーを与える年齢範囲．0歳から2歳までなら0：２とする．
#' @param no_eta_age 選択率更新法でridgeVPAの際にpenalty="p"のときで，etaがNULLでなく，p_by_age=TRUE (年齢別にペナルテイーを与えたい）ときに，etaがかからないほうの年齢範囲
#' @param p_by_age 選択率更新法でridgeVPAの際にpenalty="p"のときに年齢別にペナルテイーを与えるか与えないか．与えたい場合はTRUEとして,penalty_age（eta=NULLのとき）もしくはno_eta_age(etaがNULLでないとき）に年齢範囲を指定する．
#' @param sdreport \code{TMB=TRUE}のときに\code{sdreport()}を実行するかどうか（naa, faa, 資源量, 親魚量, Fの平均, 漁獲割合のSDを計算する）
#' @param use.equ plus groupが途中で変わる場合の計算式の選択 （old = 従来の方法（プラスグループが延長している年はプラスグループのＦが一歳若い年齢のＦと等しいという仮定は置かない），new= 新しい方法（プラスグループが延長している年はプラスグループのＦが一歳若い年齢のＦと等しいという仮定を置く)
#' @param　ave_S 選択率更新法において，最終年の選択率の仮定が通常(TRUE)とは違い，ヒラメ瀬戸内のように，最終年の選択率がtf.yearで指定した年の平均のF（つまりSUM（F,a)/SUM(F,maxage))に等しいと仮定する場合はFALSEにする.この場合，sel.def="maxage"必須．
#' @return list object:
#' \describe{
#' \item{\code{input}}{解析に用いたデータや仮定}
#' \item{\code{term.f}}{推定されたターミナルF}
#' \item{\code{np}}{推定されたターミナルFの数}
#' \item{\code{minimum}}{最適解における目的関数の値（合計）}
#' \item{\code{minimum.c}}{最適解における目的関数の値（個々）}
#' \item{\code{logLik}}{対数尤度}
#' \item{\code{gradient}}{最適解での傾き}
#' \item{\code{code}}{最適化法から返されるコード（どのような理由で最適化が停止したのかがわかる）}
#' \item{\code{q}}{推定されたq（資源量指標値の比例定数））}
#' \item{\code{b}}{推定されたb（資源量指標値の非線形性））}
#' \item{\code{sigma}}{資源量指標値の分散の平方根}
#' \item{\code{convergence}}{解が収束していれば1，そうでないと0}
#' \item{\code{message}}{最適化に関する注意（あれば）}
#' \item{\code{hessian}}{ヘッセ行列の値}
#' \item{\code{Ft}}{最終年のFの平均値}
#' \item{\code{Fc.at.age}}{Fcurrentで指定した年における平均のF}
#' \item{\code{Fc.mean}}{Fc.at.ageを平均したもの}
#' \item{\code{Fc.max}}{Fc.at.ageの最大値}
#' \item{\code{last.year}}{vpaを計算する最終年を別に指定した場合}
#' \item{\code{Pope}}{popeの近似式を用いたか否か}
#' \item{\code{ssb.coef}}{産卵親魚量の計算時期（年始めなら0,年中央なら0.5,年最後なら1)}
#' \item{\code{pred.index}}{推定された資源量指標値}
#' \item{\code{wcaa}}{caa*waa.catch}
#' \item{\code{naa}}{推定された年別年齢別資源尾数}
#' \item{\code{faa}}{推定された年別年齢別漁獲係数}
#' \item{\code{baa}}{推定された年別年齢別資源重量}
#' \item{\code{sbb}}{推定された年別年齢別産卵親魚量}
#' \item{\code{saa}}{推定された年別年齢別選択率}
#' }
#' @encoding UTF-8
#'
#' @export
#'

vpa <- function(
  dat,  # data for vpa
  sel.f = NULL,  # 最終年の選択率
  tf.year = 2008:2010, # terminal Fをどの年の平均にするか
  rec.new = NULL, # 翌年の加入を外から与える
  rec=NULL, # rec.yearの加入
  rec.year=2010,  # 加入を代入する際の年
  rps.year = 2001:2010, # 翌年のRPSをどの範囲の平均にするか
  fc.year = 2009:2011, # Fcurrentでどの範囲を参照するか
  last.year = NULL,   # vpaを計算する最終年を指定（retrospective analysis）
  last.catch.zero = FALSE,   # TRUEなら強制的に最終年の漁獲量を0にする
  faa0 = NULL,  # sel.update=TRUEのとき，初期値となるfaa
  naa0 = NULL,    # sel.update=TRUEのとき，初期値となるnaa
  f.new = NULL,
  Pope = TRUE,  # Popeの近似式を使うかどうか
  p.pope = 0.5, # Popeの式でどこで漁獲するか（0.5 = 真ん中の時期）
  tune = FALSE,  # tuningをするかどうか
  abund = "B",   # tuningの際，何の指標に対応するか
  min.age = 0,  # tuning指標の年齢参照範囲の下限
  max.age = 0,  # tuning指標の年齢参照範囲の上限
  link = "id",  # tuningのlink関数
  base = NA,  # link関数が"log"のとき，底を何にするか
  af = NA,  # 資源量指数が年の中央のとき，af=0なら漁期前，af=1なら漁期真ん中，af=2なら漁期後となる
  p.m = 0.5,  # Popeの近似式でどこで漁獲が起こるか（0.5は年の真ん中）
  omega = NULL,  # チューニングの際にselectivityを補正する行列
  index.w = NULL,  # tuning indexの重み
  use.index = "all",
  scale = 1000,  # 重量のscaling
  hessian = TRUE,
  alpha = 1,  # 最高齢と最高齢-1のFの比 F_a = alpha*F_{a-1}
  maxit = 5,  # 石岡・岸田/平松の方法の最大繰り返し数
  d = 0.0001,  # 石岡・岸田/平松の方法の収束判定基準
  min.caa = 0.001,   # caaに0があるとき，0をmin.caaで置き換える
  plot = FALSE,   # tuningに使った資源量指数に対するフィットのプロット
  plot.year = NULL,   # 上のプロットの参照年
  term.F = "max",   # terminal Fの何を推定するか: "max" or "all"
  plus.group = TRUE,
  stat.tf = "mean",  # 最終年のFを推定する統計量（年齢で同じとする）
  add.p.est = NULL,  # 追加で最高齢以外のfaaを推定する際．年齢を指定する．
  add.p.ini = NULL,
  sel.update=FALSE,  # チューニングVPAにおいて，選択率を更新しながら推定
  sel.def = "max",  #  sel.update=TRUEで選択率を更新していく際に，選択率をどのように計算するか．最大値を1とするか，平均値を1にするか...
  max.dd = 0.000001,  # sel.updateの際の収束判定基準
  ti.scale = NULL,   # 資源量の係数と切片のscaling
  tf.mat = NULL,   # terminal Fの平均をとる年の設定．0-1行列．
  eq.tf.mean = FALSE, # terminal Fの平均値を過去のFの平均値と等しくする
  no.est = FALSE,   # パラメータ推定しない．
  est.method = "ls",  # 推定方法 （ls = 最小二乗法，ml = 最尤法, ls_nolog =最小二乗法で実数） #option追加_by miyagawa(2020/10)
  b.est = FALSE,  #  bを推定するかどうか
  est.constraint = FALSE,   # 制約付き推定をするかどうか
  q.const = 1:length(abund),   # qパラメータの制約（0は推定しないで1にfix）
  b.const = 1:length(abund),   # bパラメータの制約（0は推定しないで1にfix）
  q.fix = NULL,
  b.fix = NULL,
  sigma.const = 1:length(abund),
  fixed.index.var = NULL,
  max.iter = 100,    # q,b,sigma計算の際の最大繰り返し数
  optimizer = "nlm",
  Lower = -Inf,
  Upper = Inf,
  p.fix = NULL,
  lambda = 0,   # ridge回帰係数
  beta = 2,   # penaltyのexponent  (beta = 1: lasso, 2: ridge)
  penalty = "p",
  ssb.def = "i",  # i: 年はじめ，m: 年中央, l: 年最後
  ssb.lag = 0,   # 0: no lag, 1: lag 1
  TMB=FALSE,
  # TMB.compile=FALSE,
  sel.rank=NULL,
  p.init = 0.2,   # 推定パラメータの初期値
  sigma.constraint = 1:length(abund),
  eta = NULL,
  eta.age = 0,
  tmb.file = "rvpa_tmb",
  remove.abund = NULL,
  madara=FALSE, #マダラ太平洋系群で用いているチューニングのやり方
  p_by_age = FALSE, #penalty="p"で選択率更新法のときに年齢別にペナルテイーを与えるか与えないか．
  penalty_age=NULL, #penalty="p"でp_by_age = ＴＲＵＥで選択率更新法を採用しているときの年齢参照範囲．0歳から2歳までなら0：２とする．
  no_eta_age = NULL, #etaがNULLでなく，penalty="p"で，選択率更新法を採用していて，年齢別にペナルテイーを与えたいときに，etaがかからないほうの年齢範囲
  sdreport = FALSE,
  use.equ = "new" ,#plus-groupが途中で変わる場合の計算方法の指定．従来の方法でないものを用いる場合は"new"を指定する
  ave_S=TRUE #ヒラメ瀬戸内海のように，選択率更新法において，最終年の選択率がtf.yearで指定した年の平均のF（つまりSUM（F,a)/SUM(F,maxage))に等しいと仮定する場合．注意：tf.yearで指定した年の平均の選択率とは異なる
)
{
  #sigma.constで引数を指定してしまったときは，sigma.constraintで引数を指定しなおしてもらうようにする
  if(length(sigma.const)>length(unique(sigma.const))){print("Try again!: please set sigma.const as sigma.constraint in the argument.");stop()}

  # if (TMB.compile) {
  #   library(TMB)
  #   cpp_name = paste0(tmb.file, ".cpp")
  #   compile(cpp_name)
  #   dyn.load(dynlib(tmb.file))
  # }

  # inputデータをリスト化

  argname <- ls()  # 関数が呼び出されたばかりのときのls()は引数のみが入っている
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  if (isTRUE(TMB) & isTRUE(no.est) ) TMB <- FALSE

  # data handling

  caa <- dat$caa    # catch-at-age
  waa <- dat$waa    # weight-at-age
  maa <- dat$maa    # maturity-at-age
  if (!is.null(dat$maa.tune)) maa.tune <- dat$maa.tune
  if (!is.null(dat$catch.prop)) catch.prop <- dat$catch.prop
  index <- dat$index   # abundance indices
  M <- dat$M    # natural mortality-at-age
  if(is.null(dat$waa.catch)) waa.catch <- waa else waa.catch <- dat$waa.catch

  if (isTRUE(tune) & is.null(index)) {print("Check!: There is no abundance index."); stop()}

  years <- dimnames(caa)[[2]]  # 年
  ages <- dimnames(caa)[[1]]  # 年齢

  if (class(index)=="numeric") index <- t(as.matrix(index))

# tuningの際のパラメータが1個だけ指定されている場合は，nindexの数だけ増やす
  if (isTRUE(tune)){

    nindex <- nrow(index)

    if (nindex > length(abund) & length(abund)==1) abund <- rep(abund, nindex)
    if (nindex > length(min.age) & length(min.age)==1) min.age <- rep(min.age, nindex)
    if (nindex > length(max.age) & length(max.age)==1) max.age <- rep(max.age, nindex)
    if (nindex > length(link) & length(link)==1) link <- rep(link, nindex)
    if (nindex > length(base) & length(base)==1) base <- rep(base, nindex)

    if (is.null(index.w)) index.w <- rep(1, nindex)
    if (!is.na(af[1])) if(nindex > length(af) & length(af)==1) af <- rep(af, nindex)

    q <- rep(NA, nindex)
  }

# nindexの数だけtuningの際のパラメータを増やしたのちに，use.indexを使用する場合は，use.indexで指定した部分だけを抜き出して計算する．
  if (use.index[1]!="all") {
    index <- index[use.index,,drop=FALSE]
	nindex <- nrow(index)
	q <- rep(NA, nindex)

  #以前追加したwarningは削除．何故なら，この前の部分で， nrow(index)と同じ長さだけのベクトルをtuningの際のパラメータに与える仕様にしているため．

    if(length(use.index)!=length(abund)){
      if (length(abund)>1) abund <- abund[use.index]
      if (length(min.age)>1) min.age <- min.age[use.index]
      if (length(max.age)>1) max.age <- max.age[use.index]
      if (length(link)>1) link <- link[use.index]
      if (length(base)>1) base <- base[use.index]
      if (length(af)>1) af <- af[use.index]
      if (length(index.w)>1) index.w <- index.w[use.index]
    }
  }

  if (!is.null(omega)) omega <- omega[,1:ncol(caa)]

  #

  if(!is.null(fixed.index.var)) require(mvtnorm)

  # 最終年 last.yearに値が入っている場合は，それ以降のデータを削除する（retrospective analysis）
  if (!is.null(last.year)) {
    caa <- caa[,years <= last.year]
    waa <- waa[,years <= last.year]
    maa <- maa[,years <= last.year]
    if (!is.null(dat$maa.tune)) maa.tune <- maa.tune[,years <= last.year]
    if (!is.null(dat$cathc.prop)) maa.tune <- catch.prop[,years <= last.year]
    M <- M[,years <= last.year]
    if(!is.null(index)) index <- index[,years <= last.year,drop=FALSE]
    years <- dimnames(caa)[[2]]
    dat <- list(caa=caa, waa=waa, maa=maa, M=M, index=index)
  }

  na <- apply(caa, 2, max.age.func)  # 年ごとの最大年齢（年によって最大年齢が違う場合に対応するため）
  ny <- ncol(caa)  # 年の数

  if (isTRUE(last.catch.zero)) {caa[,ny] <- 0; ny <- ny - 1; n.add <- 1; saa.new <- NULL} else n.add <- 0  # 最終年の漁獲量を0とし，年を1個減らす

  if (term.F=="max") {  # 最高齢のFだけを推定する場合
    p.init <- ifelse(is.na(p.init[1]), M[na[ny],ny], p.init[1])   # 初期値がNAなら，最終年最高齢の自然死亡係数を初期値とする
    if (!is.null(add.p.est) & is.null(add.p.ini)) {add.p.est <- add.p.est + 1; p.init <- rep(p.init, length(add.p.est)+1)}  # add.p.estに数字があれば，その年齢を追加のパラメータとして推定する
    if (!is.null(add.p.est) & !is.null(add.p.ini)) {add.p.est <- add.p.est + 1; p.init <- c(add.p.ini,p.init)}  # add.p.estに数字があれば，その年齢を追加のパラメータとして推定する
  }
  if (term.F=="all"){  # 最終年のすべての年齢のFを推定
    if(length(p.init)==0) p.init <- rep(M[na[ny],ny], na[ny]-1)  # 初期値がNAなら，最終年最高齢の自然死亡係数を0~na-1の初期値とする
    if(length(p.init) < na[ny]-1) p.init <- rep(p.init[1], na[ny]-1)  # 初期値の成分数が年齢数-1より小さい場合は，初期値の最初の値を要素に持つ年齢数-1の大きさのベクトルを初期値とする
    if(length(p.init) >= na[ny]-1) p.init <- p.init[1:na[ny]-1]  # 初期値の成分数が年齢数-1以上であれば，年齢数以上の値は使用しない
  }

  assertthat::assert_that(length(stat.tf) == 1) # stat.tfがベクトルで与えられた場合にエラーを出す

  # selectivityを更新する場合にfaa0，naa0が与えられていれば，それを使う
   if (!isTRUE(sel.update)){
       faa <- naa <- matrix(NA, nrow=max(na), ncol=ny+n.add, dimnames=list(ages, years))
   }else {
     if(is.null(faa0) | is.null(naa0)) faa <- naa <- matrix(1, nrow=max(na), ncol=ny+n.add, dimnames=list(ages, years))
     else {faa <- as.matrix(faa0); naa <- as.matrix(naa0)}
   }

  if (is.null(p.fix)) p.fix <- 1:length(p.init)

  # warnings

  if (!tune & sel.update) print("sel.update = TRUE but tune=FALSE. So, the results are unreliable.")
  if (tune & is.null(sel.f) & (!madara) & (!sel.update & term.F=="max")) print("sel.f=NULL although tune=TRUE & sel.update=FALSE & term.F=max. The results are unreliable.")
  if (tune) if(length(abund)!=nrow(index)) print("Check!: The number of abundance definition is different from the number of indices.")
  if (Pope & alpha!=1) print("Warning! The estimated F for the older ages may not be accurate if C<<N is not satisfied for the older ages.")
#  ssb.def

  if (ssb.def=="i") ssb.coef <- 0
  if (ssb.def=="m") ssb.coef <- 0.5
  if (ssb.def=="l") ssb.coef <- 1



# core function for optimization

  p.est <- function(log.p, out=FALSE){

    p <- exp(log.p)

    # sel.f==NULLで，パラメータpが1個なら，最終年最高齢のfaaとnaaを推定
    if (is.null(sel.f) & length(p) == 1){
      faa[na[ny], ny] <- p
      if (isTRUE(Pope)) naa[na[ny], ny] <- caa[na[ny], ny]*exp(M[na[ny], ny]/2)/(1-exp(-faa[na[ny], ny]))
      else  naa[na[ny], ny] <- caa[na[ny], ny]/(1-exp(-faa[na[ny], ny]-M[na[ny], ny]))*(faa[na[ny],ny]+M[na[ny],ny])/faa[na[ny],ny]
    }

    # sel.f!=NULLで，パラメータが年齢-1より少ない場合，sel.fを使って，最終年/全年齢のfaaとnaaを計算
    if (!is.null(sel.f) & length(p) < na[ny]-1){
      if(length(p)==1) faa[, ny] <- sel.f*p
      if(length(p) > 1) {   # パラメータ数が1より大きい場合，add.p.estの分，推定パラメータ数を増やす
      faa[,ny] <- sel.f*p[length(p)]
        for (i in 1:(length(p)-1)){
          faa[add.p.est[i],ny] <- p[i]*p[length(p)]
        }
      }
      if (isTRUE(Pope)) naa[, ny] <- vpa.core.Pope(caa,faa,M,ny,p=p.pope)
      else  naa[, ny] <- vpa.core(caa,faa,M,ny)
    }

   #
   if (is.null(sel.f) & isTRUE(sel.update)){
      if(length(p)==1) faa[, ny] <- p
      if(length(p) > 1) {   # パラメータ数が1より大きい場合，add.p.estの分，推定パラメータ数を増やす
      faa[,ny] <- p[length(p)]
      }
    }

   # パラメータが年齢-1であれば，それらをパラメータとして最終年/全年齢のfaaとnaaを計算
   if (length(p) == na[ny]-1){
     faa[1:(na[ny]-1), ny] <- p[p.fix]
     faa[na[ny], ny] <- alpha*p[na[ny]-1]
     if (isTRUE(Pope)) naa[, ny] <- vpa.core.Pope(caa,faa,M,ny,p=p.pope)
     else naa[, ny] <- vpa.core(caa,faa,M,ny)
   }

   # selctivityを更新しながら推定する場合
   if (isTRUE(sel.update)){
     dd <- itt <- 1
     while(dd > max.dd & itt < max.iter){
      saa <- sel.func(faa, def=sel.def)   # sel.defに従って選択率を計算
      for (i in (na[ny]-1):1){
	  if(isTRUE(ave_S)){
        saa[i, ny] <- get(stat.tf)(saa[i, years %in% tf.year])
		}
		else  saa[i, ny] <- get(stat.tf)(faa[i, years %in% tf.year])/get(stat.tf)(faa[na[ny], years %in% tf.year])
      }

      if(isTRUE(ave_S)){
      saa[na[ny], ny] <- get(stat.tf)(saa[na[ny], years %in% tf.year])
	  }

	  else  saa[na[ny], ny] <- get(stat.tf)(faa[na[ny], years %in% tf.year])/get(stat.tf)(faa[na[ny], years %in% tf.year])

      if(length(p)==1) faa[1:na[ny], ny] <- p*sel.func(saa, def=sel.def)[1:na[ny],ny] else faa[1:na[ny], ny] <- p[length(p)]*sel.func(saa, def=sel.def)[1:na[ny],ny]

      if (isTRUE(Pope)) naa[ , ny] <- vpa.core.Pope(caa,faa,M,ny,p=p.pope)
      else naa[, ny] <- vpa.core(caa,faa,M,ny)

      if (isTRUE(Pope)){
        for (i in (ny-1):1){
         naa[1:na[i], i] <- backward.calc(caa,naa,M,na,i,min.caa=min.caa,p=p.pope,plus.group=plus.group,sel.update=sel.update,alpha=alpha, use.equ=use.equ)
         if(na[1]>na[ny]){
		 if(is.na(naa[na[ny]+1, i])){
 naa[na[ny]+1, i]<-NA} else if(naa[na[ny]+1, i]==1.00000){
 naa[na[ny]+1, i]<-NA} else{
 naa[na[ny]+1, i]<-naa[na[ny]+1, i]}}#adjustment for sawara

		 faa[1:na[i], i] <- f.at.age(caa,naa,M,na,i,p=p.pope,alpha=alpha, use.equ=use.equ)
         if(na[1]>na[ny]){
		if(is.na(faa[na[ny]+1, i])){
 faa[na[ny]+1, i]<-NA} else if(faa[na[ny]+1, i]==1.00000){
 faa[na[ny]+1, i]<-NA} else{
 faa[na[ny]+1, i]<-faa[na[ny]+1, i]}}#adjustment for sawara
	   }
     }
     else{
       for (i in (ny-1):1){
         for (j in 1:(na[i]-2)){
           faa[j, i] <- ik.est(caa,naa,M,j,i,min.caa=min.caa,maxit=maxit,d=d)
         }
         if (isTRUE(plus.group)){
           faa[na[i]-1, i] <- hira.est(caa,naa,M,na[i]-1,i,alpha=alpha,min.caa=min.caa,maxit=maxit,d=d)
         }
         else faa[na[i]-1, i] <- ik.est(caa,naa,M,na[i]-1,i,min.caa=min.caa,maxit=maxit,d=d)

         faa[na[i], i] <- alpha*faa[na[i]-1, i]

         naa[1:na[i], i] <- vpa.core(caa,faa,M,i)
       }
     }

 if(na[1]>na[ny]){
 if(is.na(faa[na[ny]+1, ny-1]))faa[na[ny]+1, ny]<-NA else faa[na[ny]+1, ny-1]<-faa[na[ny]+1, ny-1]} #adjustment for sawara

	 faa1 <- faa
     saa1 <- sel.func(faa1, def=sel.def)

     for (i in (na[ny]-1):1){
	  if(isTRUE(ave_S)){
       saa1[i, ny] <- get(stat.tf)(saa1[i, years %in% tf.year])
	   }
	   else saa1[i, ny] <- get(stat.tf)(faa1[i, years %in% tf.year])/get(stat.tf)(faa1[na[ny], years %in% tf.year])
     }
	   if(isTRUE(ave_S)){
     saa1[na[ny], ny] <- get(stat.tf)(saa1[na[ny], years %in% tf.year])
	 }

	 else saa1[na[ny], ny] <- get(stat.tf)(faa1[na[ny], years %in% tf.year])/get(stat.tf)(faa1[na[ny], years %in% tf.year])

     if(length(p)==1) faa1[1:na[ny], ny] <- p*sel.func(saa1, def=sel.def)[1:na[ny],ny] else  faa1[1:na[ny], ny] <- p[length(p)]*sel.func(saa1, def=sel.def)[1:na[ny],ny]
     faa1[na[ny], ny] <- alpha*faa1[na[ny]-1, ny]

     dd <- max(sqrt((saa1[,ny] - saa[,ny])^2),na.rm=TRUE)
     itt <- itt + 1

     faa <- faa1
     }

     # トラフグ伊勢三河湾用オプション
     # 特定の年齢だけチューニングせずに最近数年のFの平均
     # 平均する年齢・年に1を入れた行列tf.matを使用（それ以外はNA）
     if (!is.null(tf.mat)) {
       for (i in 1:nrow(faa)) {
         if (sum(!is.na(tf.mat[i,]))) faa[i, ny] <- get(stat.tf)(faa[i, !is.na(tf.mat[i,])])
       }
       naa[, ny] <- vpa.core.Pope(caa,faa,M,ny,p=p.pope)
       for (i in (ny-1):(ny-na[ny]+1)){
         naa[1:na[i], i] <- backward.calc(caa,naa,M,na,i,min.caa=min.caa,p=p.pope,plus.group=plus.group,sel.update=sel.update,alpha=alpha,use.equ=use.equ)
       }
     }

   saa <- sel.func(faa, def=sel.def)

     if(length(p) > 1) {   # パラメータ数が1より大きい場合，add.p.estの分，推定パラメータ数を増やす
      for (i in 1:(length(p)-1)){
        faa[add.p.est[i],ny] <- p[i]
      }
      naa[, ny] <- vpa.core.Pope(caa,faa,M,ny,p=p.pope)
      for (i in (ny-1):(ny-na[ny]+1)){
        naa[1:na[i], i] <- backward.calc(caa,naa,M,na,i,min.caa=min.caa,p=p.pope,plus.group=plus.group,sel.update=sel.update,alpha=alpha,use.equ=use.equ)
        faa[1:na[i], i] <- f.at.age(caa,naa,M,na,i,p=p.pope,alpha=alpha, use.equ=use.equ)
      }
    }
 }

   if (!isTRUE(sel.update)){
   if (isTRUE(Pope)){
     for (i in (ny-1):1){
       naa[1:na[i], i] <- backward.calc(caa,naa,M,na,i,min.caa=min.caa,p=p.pope,plus.group=plus.group,sel.update=sel.update,alpha=alpha, use.equ=use.equ)
       faa[1:na[i], i] <- f.at.age(caa,naa,M,na,i,p=p.pope,alpha=alpha, use.equ=use.equ)
      }
   }
  else{
     for (i in (ny-1):1){
       for (j in 1:(na[i]-2)){
         faa[j, i] <- ik.est(caa,naa,M,j,i,min.caa=min.caa,maxit=maxit,d=d)
       }
       if (isTRUE(plus.group)) faa[na[i]-1, i] <- hira.est(caa,naa,M,na[i]-1,i,alpha=alpha,min.caa=min.caa,maxit=maxit,d=d)
       else faa[na[i]-1, i] <- ik.est(caa,naa,M,na[i]-1,i,min.caa=min.caa,maxit=maxit,d=d)
       faa[na[i], i] <- alpha*faa[na[i]-1, i]
       naa[1:na[i], i] <- vpa.core(caa,faa,M,i)
     }
   }

if (isTRUE(madara)){
     denom <- get(stat.tf)(faa[na[ny], years %in% tf.year])
     for (i in (na[ny]-1):1){
          faa[i,ny] <- get(stat.tf)(faa[i, years %in% tf.year])/denom*faa[na[ny],ny]
          naa[i, ny] <- caa[i, ny]*exp(M[i, ny]/2)/(1-exp(-faa[i, ny]))
          k <- 0
          for (j in (i-1):1){
            k <- k + 1
            if (i-k > 0){
              naa[j,ny-k] <- naa[j+1,ny-k+1]*exp(M[j,ny-k])+caa[j,ny-k]*exp(M[j,ny-k]/2)
              faa[j,ny-k] <- -log(1-caa[j,ny-k]*exp(M[j,ny-k]/2)/naa[j,ny-k])
            }
          }

     }
   }


    if (is.na(naa[na[ny]-1,ny])){
      if(isTRUE(Pope)){
        for (i in (na[ny]-1):1){
          if (is.null(tf.mat)) faa[i, ny] <- get(stat.tf)(faa[i, years %in% tf.year])
          else faa[i, ny] <- get(stat.tf)(faa[i, !is.na(tf.mat[i,])])
          naa[i, ny] <- caa[i, ny]*exp(M[i, ny]/2)/(1-exp(-faa[i, ny]))
          k <- 0
          for (j in (i-1):1){
            k <- k + 1
            if (i-k > 0){
              naa[j,ny-k] <- naa[j+1,ny-k+1]*exp(M[j,ny-k])+caa[j,ny-k]*exp(M[j,ny-k]/2)
              faa[j,ny-k] <- -log(1-caa[j,ny-k]*exp(M[j,ny-k]/2)/naa[j,ny-k])
            }
          }
        }
      }
      else{
        for (i in (na[ny]-1):1){
          faa[i, ny] <- get(stat.tf)(faa[i, years %in% tf.year])
          naa[i, ny] <- caa[i, ny]/(1-exp(-faa[i, ny]-M[i, ny]))*(faa[i, ny]+M[i, ny])/faa[i, ny]
          k <- 0
          for (j in (i-1):1){
            k <- k + 1
            if (i-k > 0){
              faa[j,ny-k] <- ik.est(caa,naa,M,j,ny-k,min.caa=min.caa,maxit=maxit,d=d)
              naa[j,ny-k] <- caa[j, ny-k]/(1-exp(-faa[j, ny-k]-M[j, ny-k]))*(faa[j, ny-k]+M[j, ny-k])/faa[j, ny-k]
            }
          }
        }
      }
    }
   }

   if (!is.null(rec)){
     naa[1, years %in% rec.year] <- rec
     if(isTRUE(Pope)) faa[1, years %in% rec.year] <- -as.numeric(log(1-caa[1, years %in% rec.year]/naa[1, years %in% rec.year]*exp(M[1, years %in% rec.year]/2)))
     else{
       for (j in which(years %in% rec.year)){
         faa[1,j] <- f.forward.est(caa,naa,M,1,j,maxit=maxit,d=d)
       }
     }

     terminal.year <- as.numeric(years[ny])
     for (kk in 1:length(rec.year)){
       for (i in rec.year[kk]:terminal.year){
         if(terminal.year-i > 0 & i-rec.year[kk]+1 <= max(ages)){
           naa[i-rec.year[kk]+2, years %in% (i+1)] <- naa[i-rec.year[kk]+1, years %in% i]*exp(-faa[i-rec.year[kk]+1, years %in% i]-M[i-rec.year[kk]+1, years %in% i])
           if (isTRUE(Pope)) faa[i-rec.year[kk]+2, years %in% (i+1)] <- -log(1-caa[i-rec.year[kk]+2, years %in% (i+1)]/naa[i-rec.year[kk]+2, years %in% (i+1)]*exp(M[i-rec.year[kk]+2, years %in% (i+1)]/2))
           else {
             for (j in which(years %in% (i+1))){
               if(i-rec.year[kk]+2 < na[j]-1) faa[i-rec.year[kk]+2, j] <- f.forward.est(caa,naa,M,i-rec.year[kk]+2,j,maxit=maxit,d=d)
               if (isTRUE(plus.group)){
                 if(i-rec.year[kk]+2 == na[j]-1) faa[i-rec.year[kk]+2, j] <- fp.forward.est(caa,naa,M,i-rec.year[kk]+2,j,alpha,maxit=maxit,d=d)
                 if(i-rec.year[kk]+2 == na[j]) faa[i-rec.year[kk]+2, j] <- alpha*fp.forward.est(caa,naa,M,i-rec.year[kk]+1,j,alpha,maxit=maxit,d=d)
               }
               else{
                 if(i-rec.year[kk]+2 == na[j]-1) faa[i-rec.year[kk]+2, j] <- f.forward.est(caa,naa,M,i-rec.year[kk]+2,j,maxit=maxit,d=d)
                 if(i-rec.year[kk]+2 == na[j]) faa[i-rec.year[kk]+2, j] <- alpha*f.forward.est(caa,naa,M,i-rec.year[kk]+1,j,maxit=maxit,d=d)
               }
             }
           }
         }
       }
     }
   }

  # next year

    if (isTRUE(tune)){
      if (n.add==1 & !is.na(mean(index[,ny+n.add],na.rm=TRUE))){

        new.naa <- forward.calc(faa,naa,M,na,ny+n.add,plus.group=plus.group)

        naa[,ny+n.add] <- new.naa
        baa <- naa*waa
        ssb <- baa*maa*exp(-ssb.coef*(faa+M))

        if (is.null(rec.new) & !isTRUE(last.catch.zero)) {
          new.naa[1] <- NA #median((naa[1,]/colSums(ssb))[years %in% rps.year])*sum(ssb[,ny+n.add],na.rm=TRUE)
        }else if(!is.null(rec.new)) {new.naa[1] <- rec.new
        }else if(is.null(rec.new) & isTRUE(last.catch.zero)){new.naa[1] <-NA}

        naa[1,ny+n.add] <- new.naa[1]
        baa[1,ny+n.add] <- naa[1,ny+n.add]*waa[1,ny+n.add]

        if (!is.null(f.new) & !is.null(saa.new)) faa[,ny+n.add] <- f.new*saa.new else faa[,ny+n.add] <- 0
         if (isTRUE(Pope)) caa[,ny+n.add] <- naa[,ny+n.add]*(1-exp(-faa[,ny+n.add]))*exp(-M[,ny+n.add]/2) else caa[,ny+n.add] <- naa[,ny+n.add]*(1-exp(-faa[,ny+n.add]-M[,ny+n.add]))*faa[,ny+n.add]/(faa[,ny+n.add]+M[,ny+n.add])

        ssb[1,ny+n.add] <- baa[1,ny+n.add]*maa[1,ny+n.add]*exp(-ssb.coef*(faa[1,ny+n.add]+M[1,ny+n.add]))

        if (ssb.lag==1) ssb <- cbind(NA, ssb[,-ncol(ssb)])
      }

  # tuning

    obj <- NULL

   if (tune){
     if (est.constraint | !is.null(fixed.index.var)){

       Abund <- NULL

       for (i in 1:nindex){
         abundance <- abund.extractor(abund=abund[i], naa, faa, dat, min.age=min.age[i], max.age=max.age[i], link=link[i], base=base[i], af=af[i], catch.prop=catch.prop, sel.def=sel.def, p.m=p.m, omega=omega, scale=scale)
         if (!is.null(remove.abund)) {
           if (nrow(remove.abund)!=nindex || ncol(remove.abund)!=ncol(dat$caa)) {
             stop("The dimension of 'remove.abund' must be the same as that of 'caa'")
           }
           abundance <- abundance - as.numeric(remove.abund[i,])
           abundance[abundance<1e-6] <- 1e-6
           }
         Abund <- rbind(Abund, abundance)
       }

       if (is.null(fixed.index.var)) est.qbs <- qbs.f(q.const, b.const, sigma.constraint, index, Abund, nindex, index.w, max.dd, max.iter) else {
       p00 <- c(log(q.const[which(index.w >0)]), log(sigma.constraint[which(index.w >0)]))
       est.qbs <- qbs.f2(p00, index, Abund, nindex, index.w, fixed.index.var)
       }

       q <- exp(est.qbs$q)
       b <- est.qbs$b
       sigma <- est.qbs$sigma
       obj <- est.qbs$obj
       convergence <- est.qbs$convergence
       obj0 <- obj

       rownames(Abund) <- 1:nindex
     }
     else{

    if (est.method=="ls")
    {
        Abund <- nn <- sigma <- b <- NULL
        for (i in 1:nindex)
        {
            abundance <- abund.extractor(abund=abund[i], naa, faa, dat, min.age=min.age[i], max.age=max.age[i], link=link[i], base=base[i], af=af[i], catch.prop=catch.prop, sel.def=sel.def, p.m=p.m, omega=omega, scale=scale)
            if (!is.null(remove.abund)) {
              abundance <- abundance - as.numeric(remove.abund[i,])
              abundance[abundance<1e-6] <- 1e-6
            }
            Abund <- rbind(Abund, abundance)
            avail <- which(!is.na(as.numeric(index[i,])))

            if (b.est)
            {
                if (is.null(b.fix))
                {
                    b[i] <- cov(log(as.numeric(index[i,avail])),log(as.numeric(abundance[avail])))/var(log(as.numeric(abundance[avail])))
                }else
                {
                    if (is.na(b.fix[i])) b[i] <- cov(log(as.numeric(index[i,avail])),log(as.numeric(abundance[avail])))/var(log(as.numeric(abundance[avail]))) else b[i] <- b.fix[i]
                }
            }else
            {
                if (is.null(b.fix)) b[i] <- 1 else b[i] <- b.fix[i]
            }
            if (is.null(q.fix))
            {
                q[i] <- exp(mean(log(as.numeric(index[i,avail]))-b[i]*log(as.numeric(abundance[avail]))))
            }else
            {
                q[i] <- q.fix[i]
            }
            obj <- c(obj,index.w[i]*sum((log(as.numeric(index[i,avail]))-log(q[i])-b[i]*log(as.numeric(abundance[avail])))^2))
        }
    }
	if (est.method=="ls_nolog") #option added for the use of madara (2020/10)
    {
        Abund <- nn <- sigma <- b <- NULL
        for (i in 1:nindex)
        {
            abundance <- abund.extractor(abund=abund[i], naa, faa, dat, min.age=min.age[i], max.age=max.age[i], link=link[i], base=base[i], af=af[i], catch.prop=catch.prop, sel.def=sel.def, p.m=p.m, omega=omega, scale=scale)
            Abund <- rbind(Abund, abundance)
            avail <- which(!is.na(as.numeric(index[i,])))
            #avail <- as.numeric(index[i,])
            if (b.est)
            {stop("no option to estimate b when using real number instead of log !!!!")#実数平方和を用いる場合はｂ推定のオプションはなし。
                if (is.null(b.fix))
                {
                    b[i] <- cov(log(as.numeric(index[i,avail])),log(as.numeric(abundance[avail])))/var(log(as.numeric(abundance[avail])))
                }else
                {
                    if (is.na(b.fix[i])) b[i] <- cov(log(as.numeric(index[i,avail])),log(as.numeric(abundance[avail])))/var(log(as.numeric(abundance[avail]))) else b[i] <- b.fix[i]
                }
            }else
            {
                if (is.null(b.fix)) b[i] <- 1 else b[i] <- b.fix[i]
            }
            if (is.null(q.fix))
            {
                #q[i] <- exp(mean(log(as.numeric(index[i,avail]))-b[i]*log(as.numeric(abundance[avail]))))
				 q[i] <- sum(as.numeric(index[i,avail])*as.numeric(abundance[avail])^b[i])/sum(as.numeric(abundance[avail])^(2*b[i])) #changed
            }else
            {
                q[i] <- q.fix[i]
            }
            #obj <- c(obj,index.w[i]*sum((log(as.numeric(index[i,avail]))-log(q[i])-b[i]*log(as.numeric(abundance[avail])))^2))
			obj <- c(obj,index.w[i]*sum((as.numeric(index[i,avail])-q[i]*as.numeric(abundance[avail ])^b[i])^2))
        }
    }
    if (est.method=="ml")
    {
      if (use.index[1] != "all") sigma.constraint <- sigma.constraint[use.index]

        if(!(length(sigma.constraint)==nindex))
        {
            stop("length of sigma constraint does not match the number of indices!!!!")#sigma.constraintの長さがindexの本数と異なる場合にはエラーを出して停止。
        }
        Abund <- nn <- sigma <- b <- NULL
        for (i in 1:nindex)
        {
            abundance <- abund.extractor(abund=abund[i], naa, faa, dat, min.age=min.age[i], max.age=max.age[i], link=link[i], base=base[i], af=af[i], catch.prop=catch.prop, sel.def=sel.def, p.m=p.m, omega=omega, scale=scale)
            if (!is.null(remove.abund)) {
              abundance <- abundance - as.numeric(remove.abund[i,])
              abundance[abundance<1e-6] <- 1e-6
            }
            Abund <- rbind(Abund, abundance)
            avail <- which(!is.na(as.numeric(index[i,])))
            nn[i] <- length(avail)
            if (b.est)
            {
                if (is.null(b.fix))
                {
                    b[i] <- cov(log(as.numeric(index[i,avail])),log(as.numeric(abundance[avail])))/var(log(as.numeric(abundance[avail])))
                }else
                {
                    if (is.na(b.fix[i]))
                    {
                        b[i] <- cov(log(as.numeric(index[i,avail])),log(as.numeric(abundance[avail])))/var(log(as.numeric(abundance[avail])))
                    }else
                    {
                        b[i] <- b.fix[i]
                    }
                }
            }else
            {
                if (is.null(b.fix))
                {
                    b[i] <- 1
                }else
                {
                    b[i] <- b.fix[i]
                }
            }
            if (is.null(q.fix))
            {
                q[i] <- exp(mean(log(as.numeric(index[i,avail]))-b[i]*log(as.numeric(abundance[avail]))))
            }else
            {
                q[i] <- q.fix[i]
            }
            #sigma[i] <- sqrt(sum((log(as.numeric(index[i,avail]))-log(q[i])-b[i]*log(as.numeric(abundance[avail])))^2)/nn[i])
            #obj <- c(obj,index.w[i]*(-as.numeric(na.omit(dnorm(log(as.numeric(index[i,avail])),log(q[i])+b[i]*log(as.numeric(abundance[avail])),sigma[i],log=TRUE)))))
        }
        unique.sigma.constraint <- unique(sigma.constraint)
        for(i in 1:length(unique.sigma.constraint))
        {
            index.num <- which(sigma.constraint==unique.sigma.constraint[i])
            sq.error <- 0
            for(j in index.num)
            {
                abundance <- abund.extractor(abund=abund[j], naa, faa, dat, min.age=min.age[j], max.age=max.age[j], link=link[j], base=base[j], af=af[j], catch.prop=catch.prop, sel.def=sel.def, p.m=p.m, omega=omega, scale=scale)
                if (!is.null(remove.abund)) {
                  abundance <- abundance - as.numeric(remove.abund[i,])
                  abundance[abundance<1e-6] <- 1e-6
                }
                avail <- which(!is.na(as.numeric(index[j,])))
                sq.error <- sq.error + sum((log(as.numeric(index[j,avail]))-log(q[j])-b[j]*log(as.numeric(abundance[avail])))^2)
            }
            sigma[index.num] <- sqrt(sq.error/sum(nn[index.num]))
        }
        for (i in 1:nindex)
        {
            abundance <- abund.extractor(abund=abund[i], naa, faa, dat, min.age=min.age[i], max.age=max.age[i], link=link[i], base=base[i], af=af[i], catch.prop=catch.prop, sel.def=sel.def, p.m=p.m, omega=omega, scale=scale)
            if (!is.null(remove.abund)) {
              abundance <- abundance - as.numeric(remove.abund[i,])
              abundance[abundance<1e-6] <- 1e-6
            }
            avail <- which(!is.na(as.numeric(index[i,])))
            obj <- c(obj,index.w[i]*(-as.numeric(na.omit(dnorm(log(as.numeric(index[i,avail])),log(q[i])+b[i]*log(as.numeric(abundance[avail])),sigma[i],log=TRUE)))))
        }
    }

      obj0 <- obj
      obj <- sum(obj)
      convergence <- 1
      saa <- sel.func(faa, def=sel.def)

      if (penalty=="p" && isFALSE(p_by_age)) {

        if (is.null(eta) || eta==-1) {
          obj <- (1-lambda)*obj + lambda*sum(p^beta)
        } else {
          eta.age <- eta.age + 1
          obj <- (1-lambda)*obj + lambda*(eta*sum(p[eta.age]^beta) + (1-eta)*sum(p[-eta.age]^beta))
        }
        }

	  if (penalty=="p" && isTRUE(p_by_age)) {

        if (is.null(eta)) {
          if (is.null(penalty_age)) {stop("please specify penalty_age")}#etaがNULLでpenalty="p"で選択率更新法を採用していて,penaltyを年齢別に与えたいのにpenalty_ageを未指定の場合にはエラーを出して停止。
           penalty_age <- penalty_age + 1
		  obj <- (1-lambda)*obj + lambda*sum(faa[penalty_age,ny]^beta)
        } else {
          if (is.null(no_eta_age)) {stop("please specify no_eta_age")}#etaがNULLでpenalty="p"で選択率更新法を採用していて,penaltyを年齢別に与えたいのにpenalty_ageを未指定の場合にはエラーを出して停止。
           eta.age <- eta.age + 1
		  no_eta_age <- no_eta_age +1
		  obj <- (1-lambda)*obj + lambda*(eta*sum(faa[eta.age,ny]^beta) + (1-eta)*sum(faa[no_eta_age,ny]^beta))
        }
        }

      if (penalty=="f") obj <- (1-lambda)*obj + lambda*sum((abs(faa[1:(na[ny]-1),ny]-apply(faa[1:(na[ny]-1), years %in% tf.year],1,get(stat.tf))))^beta)

      if (penalty=="s") obj <- (1-lambda)*obj + lambda*sum((abs(saa[1:(na[ny]-1),ny]-apply(saa[1:(na[ny]-1), years %in% tf.year],1,get(stat.tf))))^beta)

      if (!is.null(sel.rank)) obj <- obj+1000000*sum((rank(saa[,ny])-sel.rank)^2)

      rownames(Abund) <- 1:nindex
      }
    }
  }
  else {obj <- (p - alpha*faa[na[ny]-1, ny])^2; obj0 <- NA}

  #

    if (isTRUE(out)) {
        # next year

        if (n.add==1 & is.na(naa[1,ny+n.add])){
          new.naa <- forward.calc(faa,naa,M,na,ny+n.add,plus.group=plus.group)
          if (!is.null(f.new) & !is.null(saa.new)) faa[,ny+n.add] <- f.new*saa.new else faa[,ny+n.add] <- 0
          naa[,ny+n.add] <- new.naa
          baa <- naa*waa
          ssb <- baa*maa*exp(-ssb.coef*(faa+M))

          if (is.null(rec.new) & !isTRUE(last.catch.zero)) {
            new.naa[1] <- NA # median((naa[1,]/colSums(ssb))[years %in% rps.year])*sum(ssb[,ny+n.add],na.rm=TRUE)
          }else if(!is.null(rec.new)) {new.naa[1] <- rec.new
          }else if(is.null(rec.new) & isTRUE(last.catch.zero)){new.naa[1] <-NA}

          naa[1,ny+n.add] <- new.naa[1]
          baa[1,ny+n.add] <- naa[1,ny+n.add]*waa[1,ny+n.add]

          if (isTRUE(Pope)) caa[,ny+n.add] <- naa[,ny+n.add]*(1-exp(-faa[,ny+n.add]))*exp(-M[,ny+n.add]/2) else caa[,ny+n.add] <- naa[,ny+n.add]*(1-exp(-faa[,ny+n.add]-M[,ny+n.add]))*faa[,ny+n.add]/(faa[,ny+n.add]+M[,ny+n.add])

          ssb[1,ny+n.add] <- baa[1,ny+n.add]*maa[1,ny+n.add]*exp(-ssb.coef*(faa[1,ny+n.add]+M[1,ny+n.add]))

        if (ssb.lag==1) ssb <- cbind(NA, ssb[,-ncol(ssb)])
        }
        else {
          baa <- naa*waa
          ssb <- baa*maa*exp(-ssb.coef*(faa+M))
          if (ssb.lag==1) ssb <- cbind(NA, ssb[,-ncol(ssb)])
        }


        obj <- list(minimum=obj, minimum.c=obj0, caa=caa, naa=naa, faa=faa, baa=baa, ssb=ssb)
        if (isTRUE(eq.tf.mean)) obj$p <- max(faa[,ny],na.rm=TRUE)

        if (isTRUE(tune)) {
          if (est.method=="ls" ){
            if (use.index[1]=="all") Nindex <- sum(!is.na(index[index.w > 0,])) else Nindex <- sum(!is.na(index[index.w[use.index > 0] > 0,]))
            Sigma2 <- obj$minimum/Nindex
            neg.logLik <- Nindex/2*log(2*pi*Sigma2)+Nindex/2
            obj$q <- q
            obj$b <- b
            obj$sigma <- sqrt(Sigma2)
            obj$convergence <- convergence
            obj$Abund <- Abund
            obj$logLik <- -neg.logLik
          }
		   if (est.method=="ls_nolog" ){
            if (use.index[1]=="all") Nindex <- sum(!is.na(index[index.w > 0,])) else Nindex <- sum(!is.na(index[index.w[use.index > 0] > 0,]))
            Sigma2 <- obj$minimum/Nindex
            #neg.logLik <- Nindex/2*log(2*pi*Sigma2)+Nindex/2
			neg.logLik <-log(obj$minimum)
            obj$q <- q
            obj$b <- b
            obj$sigma <- sqrt(Sigma2)
            obj$convergence <- convergence
            obj$Abund <- Abund
            obj$logLik <- -neg.logLik
          }
          if (est.method=="ml"|est.constraint| !is.null(fixed.index.var)){
            if (est.constraint){
              names(q) <- q.const
              names(b) <- b.const
              names(sigma) <- sigma.constraint
            }
            obj$convergence <- convergence
            obj$q <- q
            obj$b <- b
            obj$sigma <- sigma
            obj$Abund <- Abund
            obj$logLik <- -obj$minimum
          }

        }
    }

    return(obj)   # 目的関数を返す
  }


########################################################################################

  # execution of optimization
#  if (isTRUE(ADMB)){
#    require(R2admb)
#
#    index2 <- as.matrix(t(apply(index,1,function(x) ifelse(is.na(x),0,x))))
#
#    Type <- ifelse(abund=="SSB", 1, ifelse(abund=="B",4,ifelse(abund=="N",3,2)))
#
#    if(is.null(dat$maa.tune)) MAA <- as.matrix(dat$maa) else MAA <- as.matrix(dat$maa.tune)
#    if (is.na(af[1])) af <- rep(0,nindex)
#
#    data2 <- list(A=nrow(dat$caa),Y=ncol(dat$caa),K=length(use.index),Est=ifelse(est.method=="ls",0,1),b_est=as.numeric(b.est),alpha=alpha,lambda=lambda,beta=beta,Type=Type,w=index.w,af=af,CATCH=as.matrix(dat$caa),WEI=as.matrix(dat$waa/scale),MAT=MAA,M=as.matrix(dat$M),CPUE=index2,MISS=ifelse(index2==0,1,0))
#
#    init <- log(p.init)
#
#    write_dat("vpa",data2)
#    write_pin("vpa",init)

#    system("vpa -nohess")

#    summary.p.est <- read_pars("vpa")
#    summary.p.est$estimate <- exp(summary.p.est$coeflist$log_F)
#    summary.p.est$minimum <- -summary.p.est$loglik
#    summary.p.est$gradient <- summary.p.est$maxgrad
#    summary.p.est$code <- 0
#    log.p.hat <- log(summary.p.est$estimate)
#  } else {
  if (isTRUE(TMB)){

  if(isTRUE(TMB) & isTRUE(sel.update)){print("TMB is not supported for sel.update method");stop()}
  if(isTRUE(TMB) & alpha!=1){print("TMB is not supported for alpha!=1");stop()}
  if(isTRUE(TMB) & isFALSE(Pope)){print("TMB is not supported for baranov equation option. only for Pope");stop()}
  if(isTRUE(TMB) & penalty=="f"){print("TMB is not supported for penalty=f. only for penalty=p");stop()}
  if(isTRUE(TMB) & penalty=="s"){print("TMB is not supported for penalty=s. only for penalty=p");stop()}

    index2 <- as.matrix(t(apply(index,1,function(x) ifelse(is.na(x),0,x))))

    Ab_type <- ifelse(abund=="SSB", 1, ifelse(abund=="N", 2, ifelse(abund=="B", 3, 4)))
    sel_def <- ifelse(sel.def=="max",0,ifelse(sel.def=="mean",1,2))
    Ab_type_age <- ifelse(is.na(min.age),0,min.age)
    Ab_type_max_age <- ifelse(is.na(max.age),0,max.age)+1

    if(is.null(dat$maa.tune)) MAA <- as.matrix(dat$maa) else MAA <- as.matrix(dat$maa.tune)
    if (is.na(af[1])) af <- rep(0,nindex)

    if (isTRUE(b.est)) b_fix <- rep(0,nindex) else b_fix <- rep(1,nindex)
    if (!is.null(b.fix)) b_fix <- ifelse(is.na(b.fix),0,b.fix)
    # if (use.index[1] != "all") b_fix <- b_fix[use.index]

    # if (!(length(sigma.constraint)==nindex)) {
    #   stop("length of sigma constraint does not match the number of indices!!!!")#sigma.constraintの長さがindexの本数と異なる場合にはエラーを出して停止。
    # }
    unique.sigma.constraint <- unique(sigma.constraint)
    sigma_constraint <- sigma.constraint
    if (use.index[1] != "all") {
      unique.sigma.constraint <- unique.sigma.constraint[use.index]
      sigma_constraint <- sigma_constraint[use.index]
      }
    nsigma <- length(unique.sigma.constraint)
    for (i in 1:nsigma) {
      pos <- which(sigma_constraint == unique.sigma.constraint[i])
      sigma_constraint[pos] <- i-1
    }
    if (is.null(eta)) eta <- -1.0
    eta_age <- rep(1,length(p.init))
    eta_age[eta.age+1] <- 0
    data2 <- list(Est=ifelse(est.method=="ls",0,1),b_fix=as.numeric(b_fix),
                  alpha=alpha,lambda=lambda,beta=beta,Ab_type=Ab_type,sel_def=sel_def,
                  Ab_type_age=Ab_type_age,Ab_type_max_age=Ab_type_max_age,w=index.w,af=af,CATCH=t(as.matrix(dat$caa)),WEI=t(as.matrix(dat$waa)),MAT=t(MAA),M=t(as.matrix(dat$M)),CPUE=t(index2),MISS=t(ifelse(index2==0,1,0)),Last_Catch_Zero=ifelse(isTRUE(last.catch.zero),1,0),sigma_constraint=sigma_constraint,eta=eta,eta_age=eta_age)

    parameters <- list(
      log_F=log(p.init)
    )

    obj2 <- try(TMB::MakeADFun(data2, parameters, DLL=tmb.file))
    if (class(obj2) == "try-error") {
      stop("Please run use_rvpa_tmb() first!")
    }
    opt <- nlm(obj2$fn, obj2$par, gradient=obj2$gr, hessian=hessian)
    if (sdreport) rep <- TMB::sdreport(obj2)

    summary.p.est <- opt
    # summary.p.est <- list()
    # summary.p.est$estimate <- exp(opt$estimate)
    # summary.p.est$minimum <- -opt$minimum
    # summary.p.est$gradient <- opt$gradient
    # summary.p.est$code <- opt$code
    # log.p.hat <- opt$estimate
    log.p.hat <- summary.p.est$estimate
  } else {
    if (isTRUE(no.est)){
      if (isTRUE(eq.tf.mean)) {
        summary.p.est <- p.est(log(p.init), out=TRUE)
        summary.p.est <- list(estimate=summary.p.est$p, minimum=p.est(log(summary.p.est$p)), gradient=NA, code=NA)
        log.p.hat <- log(summary.p.est$estimate)
      }else{
        summary.p.est <- list(estimate=log(p.init), minimum=p.est(log(p.init)), gradient=NA, code=NA)
        log.p.hat <- summary.p.est$estimate
      }
    }else{
      if (optimizer=="nlm") summary.p.est <- nlm(p.est, log(p.init), hessian=hessian)
      if (optimizer=="nlminb") {
        summary.p.est <- nlminb(log(p.init), p.est, hessian=hessian, lower=Lower, upper=Upper)
        summary.p.est$estimate <- summary.p.est$par
        summary.p.est$minimum <- summary.p.est$objective
        summary.p.est$gradient <- NA
        summary.p.est$code <- summary.p.est$convergence
      }
      log.p.hat <- summary.p.est$estimate
    }
  }

  gradient <- summary.p.est$gradient
  code <- summary.p.est$code
  message.nlminb <- summary.p.est$message

  np <- length(summary.p.est$estimate)

  out <- p.est(log.p.hat, out=TRUE)

  term.f <- exp(log.p.hat)

  #

  if(isTRUE(hessian)) hessian <- summary.p.est$hessian

  naa <- as.data.frame(out$naa)
  faa <- as.data.frame(out$faa)
  baa <- as.data.frame(out$baa)
  ssb <- as.data.frame(out$ssb)
  saa <- as.data.frame(sel.func(faa, def=sel.def))

  if(isTRUE(tune)){
    logLik <- out$logLik
    sigma <- out$sigma
    q <- out$q
    b <- out$b
    convergence <- out$convergence
    message <- message.nlminb
    pred.index <- q*out$Abund^b
  }
 else logLik <- sigma <- q <- b <- convergence <- message <- pred.index <- NULL

Ft <- mean(faa[,ny],na.rm=TRUE)
  Fc.at.age <- apply(faa[,years %in% fc.year,drop=FALSE],1,mean)  # drop=FALSEで，行列のベクトル化を防ぐ
  Fc.mean <- mean(Fc.at.age,na.rm=TRUE)
  Fc.max <- max(Fc.at.age,na.rm=TRUE)

  res <- list(input=arglist, term.f=term.f, np=np, minimum=out$minimum, minimum.c=out$minimum.c, logLik=logLik, gradient=gradient, code=code, q=q, b=b, sigma=sigma, convergence=convergence, message=message, hessian=hessian, Ft=Ft, Fc.at.age=Fc.at.age, Fc.mean=Fc.mean, Fc.max=Fc.max, last.year=last.year, Pope=Pope, ssb.coef=ssb.coef, pred.index=pred.index, wcaa=caa*waa.catch,naa=naa, faa=faa, baa=baa, ssb=ssb, saa=saa)

  invisible(res2 <- list(input=arglist, use.index=use.index, abund=abund, min.age=min.age, max.age=max.age, link=link, base=base, af=af, index.w=index.w, q=q, naa=naa, faa=faa, baa=baa, ssb=ssb, pred.index=pred.index, sigma=sigma, b=b)) #use.indexを考慮し，実際にVPAのチューニングで与えた値

  # print(list(use.index=use.index, abund=abund, min.age=min.age, max.age=max.age, link=link, base=base, af=af, index.w=index.w))

  if (isTRUE(plot) & isTRUE(tune)){
    if(is.null(plot.year)) plot.year <- colnames(naa) %>% as.numeric()
    graph <- try(plot_residual_vpa2(res2, index_name = NULL, plot_year = plot.year)) # plot.yearに対応する引数を追加してください
    if(class(graph)=="try-error"){
      for (i in 1:nindex){
        Y <- years %in% plot.year
        Pred <- (index[i,Y]/q[i])^(1/b[i])
        plot(years[Y], Pred, ylim=range(Pred, out$Abund[i,Y], na.rm=TRUE),col=3,pch=16,xlab="Year",ylab=paste("index", i), main=abund[i])
        lines(years[Y],out$Abund[i,Y],col=2,lwd=2)
      } # for(i) 従来のplot
    } else {
      gridExtra::grid.arrange(graph$year_resid, graph$fitting_Index, graph$abund_Index) # 3つのggplotを並べる
      res <- c(res, list(plot = graph))
    } # 加筆（浜辺）
  }

  if (isTRUE(TMB) & isTRUE(sdreport)) {
    res$obj <- obj2
    res$rep <- rep
    }

  class(res) <- "vpa"
  return(invisible(res))



}

#'
#' VPAの推定値についてprofile likelihood (one parameter)を実施する
#'
#' @param res vpa関数の出力値
#' @encoding UTF-8
#'
#' @export
#'

 profile_likelihood.vpa <- function(res,Alpha=0.95,min.p=1.0E-6,max.p=1,L=20,method="ci"){

   res.c <- res
   res.c$input$no.est <- TRUE
   res.c$input$plot <- FALSE

   like <- function(p,method="ci") {
     res.c$input$p.init <- p

     res1 <- do.call(vpa,res.c$input)

     if (method=="ci") obj <- (-2*(res1$logLik - res$logLik)-qchisq(Alpha,1))^2
#     if (method=="dist") obj <- res1$logLik
     if (method=="dist"){   # 市野川変更
          obj <- list(logLik=res1$logLik,LLs=res1$minimum.c)
     }
     return(obj)
  }

  if (method=="ci"){
    res.lo <- nlminb(start=res$term.f*0.5, like, lower=0.001, upper=0.999*res$term.f, method="ci")
    res.up <- nlminb(start=res$term.f*1.5, like, lower=1.001*res$term.f, upper=Inf, method="ci")
    out <- list(lower=res.lo,upper=res.up,ci=c(res.lo$par, res.up$par))
  }
  if (method=="dist"){
    p0 <- seq(min.p,max.p,len=L)
    tmp <- lapply(p0, like, method="dist")
    out <- list(logLik=sapply(tmp,function(x) x$logLik),
    			LLs = sapply(tmp,function(x) x$LLs))
    out$TLL <- -out$logLik - min(-out$logLik)
    out$RLLs <- sweep(out$LLs,1,apply(out$LLs,1,min),FUN="-")
    out$p0 <- p0
#    out <- p0 # 市野川変更
  }

  return(out)
}

dp.est <- function(p,res,Ref,target="F",beta=1){
    res.c <- res
    res.c$input$no.est <- TRUE
    res.c$input$plot <- FALSE

    res.c$input$p.init <- p

    ny <- length(res.c$faa[1,])
    na <- length(res.c$faa[,ny])

     res1 <- do.call(vpa,res.c$input)

     if (target=="B") out <- -res1$logLik+beta*(sum(res1$baa[,ny])-Ref)^2
     if (target=="F") out <- -res1$logLik+beta*(res1$faa[na,ny]-Ref)^2

     return(out)
}

pl.ci.dp <- function(res,target="F",beta=10^5,Alpha=0.8,lo.p=0.1,up.p=2.0,lo.Ref=0.5,up.Ref=3,method="ci"){
    res.c <- res
    res.c$input$no.est <- TRUE
    res.c$input$plot <- FALSE

     p.est <- function(Ref) optimize(dp.est,c(lo.p,up.p),res=res,Ref=Ref,target=target,beta=beta)$minimum

    ny <- length(res$faa[1,])
    na <- length(res$faa[,ny])

    if (target=="F") Ref0 <- res$faa[na,ny]
    if (target=="B") Ref0 <- sum(res$baa[,ny])

      like <- function(Ref,method="ci") {

        p <- p.est(Ref)

        res.c$input$p.init <- p

        res1 <- do.call(vpa,res.c$input)

  if (method=="ci") obj <- -2*(res1$logLik - res$logLik)-qchisq(Alpha,1)
     if (method=="dist") obj <- res1$logLik
     return(obj)
  }

  if (method=="ci"){
    res.lo <- uniroot(like, lower=Ref0*lo.Ref, upper=Ref0, method="ci")
    res.up <- uniroot(like, lower=Ref0, upper=Ref0*up.Ref, method="ci")
    out <- list(lower=res.lo,upper=res.up,ci=c(res.lo$root, res.up$root))
  }
  if (method=="dist"){
    p0 <- seq(Ref0*lo.Ref,Ref0*up.Ref,len=L)
    out <- sapply(p0, like, method="dist")
  }

  return(out)
}

profile.likelihood.vpa.B <- function(res,Alpha=0.95,min.p=1.0E-6,max.p=1,L=20,method="ci"){

   res.c <- res
   res.c$input$no.est <- TRUE

   like <- function(p,method="ci") {

     Bm <- exp(p)

     p0 <- res.c$term.f

     f1 <- function(p0){
       res.c$input$p.init <- p0
       res1 <- do.call(vpa,res.c$input)
       (sum(res1$baa[,37])-Bm)^2
     }

     res1 <- nlm(f1,p0)

     res.c$input$p.init <- res1$estimate
     res1 <- do.call(vpa,res.c$input)

     if (method=="ci") obj <- (-2*(res1$logLik - res$logLik)-qchisq(Alpha,1))^2
     if (method=="dist") obj <- res1$logLik
     return(obj)
  }

  if (method=="ci"){
    res.lo <- nlminb(start=log(sum(res$baa[,37])*0.5), like, lower=-Inf, upper=log(sum(res$baa[,37])), method="ci")
    res.up <- nlminb(start=log(sum(res$baa[,37])*1.5), like, lower=log(sum(res$baa[,37])), upper=Inf, method="ci")
    out <- list(lower=res.lo,upper=res.up,ci=c(res.lo$par, res.up$par))
  }
  if (method=="dist"){
    p0 <- seq(min.p,max.p,len=L)
    out <- sapply(p0, like, method="dist")
  }

  return(out)
}

#' bootstrapを実施する
#'
#' @param res vpaの出力値
#' @param method "p": パラメトリックブートストラップ。"n": ノンパラメトリックブートストラップ。"r": 残差をスムージングした後ブートストラップt法
#'
#' @encoding UTF-8
#'
#' @seealso
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#' @export

boo.vpa <- function(res,B=5,method="p",mean.correction=FALSE){
  ## method == "p": parametric bootstrap
  ## method == "n": non-parametric bootstrap
  ## method == "r": smoothed residual bootstrap-t

  assertthat::assert_that(method == "p" | method == "n" | method == "r")

  if(is.numeric(res$input$use.index)){
    index <- res$input$dat$index[res$input$use.index,]
    p.index <- res$pred.index
    assertthat::assert_that(nrow(index) == nrow(p.index)) # obsとpredのデータの種類数が違ったら止める
    resid <- log(as.matrix(index))-log(as.matrix(p.index))
  } else { # use.index = "all" (default)
    index <- res$input$dat$index
    p.index <- res$pred.index
    resid <- log(as.matrix(index))-log(as.matrix(p.index))
  }

  R <- nrow(resid)                               # 残差の行数(=用いたデータの数)
  n <- apply(resid,1,function(x) sum(!is.na(x))) # 残差の数
  np <- res$np                                   # パラメータ数
  rs2 <- rowSums(resid^2, na.rm=TRUE)/(n-np)

  res.c <- res

  if(np == 1){
    res.c$input$p.init <- res$term.f[1] # original code
  } else {
    res.c$input$p.init <- c(res$term.f,
                            rev(res$term.f)[1]*res$input$alpha
                            )    # 初期値とターミナルF合わせてる
  }

  b.index <- res$input$dat$index # ブートストラップCPUEの箱を先に作っておく
  Res1 <- list()

  for (b in 1:B){
    print(b)

    for (i in 1:R){
    #use.indexオプションを使っているとき、b.indexの行数とずれるので、b.indexの行番号をjで設定
      j <- ifelse(is.numeric(res$input$use.index), res$input$use.index[i],i)
      if (method=="p") b.index[j,!is.na(index[i,])] <- exp(log(p.index[i,!is.na(index[i,])]) + rnorm(sum(!is.na(index[i,])),0,sd=sqrt(rs2[i])))
      if (method=="n") b.index[j,!is.na(index[i,])] <- exp(log(p.index[i,!is.na(index[i,])]) + sample(resid[i,!is.na(index[i,])],length(index[i,!is.na(index[i,])]),replace=TRUE))
      if (isTRUE(mean.correction)) b.index[j,!is.na(index[i,])] <- b.index[j,!is.na(index[i,])]*exp(-rs2[i]/2)
      if (method=="r") {
        rs.d <- density(resid[i,!is.na(index[i,])])
        rs.db <- sample(rs.d$x,length(index[i,!is.na(index[i,])]),prob=rs.d$y,replace=TRUE)
        sd.j <- sd(rs.db)
        s.rs.b <- rs.db/sd.j
        b.index[j,!is.na(index[i,])] <- exp(log(p.index[i,!is.na(index[i,])]) + s.rs.b*sqrt(rs2[i]))
      }
      if (isTRUE(mean.correction)) b.index[j,!is.na(index[i,])] <- b.index[j,!is.na(index[i,])]*exp(-rs2[i]/2)
    }

    res.c$input$dat$index <- b.index

    res1 <- try(safe_call(vpa,res.c$input, force=TRUE))  # do.callからsafe_callに変更(浜辺'20/06/30)
    if(class(res1)=="try-error"){
      Res1[[b]] <- "try-error"
    }
    else{
      Res1[[b]] <- list(index=b.index,naa=res1$naa,baa=res1$baa,ssb=res1$ssb,faa=res1$faa,saa=res1$saa,
                      Fc.at.age=res1$Fc.at.age,q=res1$q,b=res1$b,sigma=res1$sigma) # 2013.8.20追記(市野川)
    }
  }

  return(Res1)
}


logit <- function(x) log(x/(1-x))

##

cv.est <- function(res,n=5){

   nr <- ifelse(res$input$use.index=="all", 1:nrow(res$input$dat$index), res$input$use.index)
   nc <- ncol(res$input$dat$index)

   obj <- NULL

   for (i in 0:(n-1)){
     res.c <- res

     res.c$input$dat$index[,nc-i] <- NA
     res.c$input$plot <- FALSE
#     res.c$input$p.init <- res$term.f

     res1 <- do.call(vpa,res.c$input)

     if (abs(res1$gradient) < 10^(-3)){
       obj <- c(obj,mean(dnorm(log(res$input$dat$index[nr,nc-i]),log(res1$pred[,nc-i]),res1$sigma,log=TRUE),na.rm=TRUE))
     }
   }

   return(mean(obj,na.rm=TRUE))
}

#' レトロスペクティブ解析の実施
#'
#' @param res VPAの出力
#' @param n 除く年数
#' @param b.fix b推定してる場合にbを固定するかどうか
#' @param remove.maxAgeF Mohn's rhoを計算する際に最高齢のFを除くか（alphaを仮定して計算していることが多いから）
#' @param ssb.forecast Mohn's rhoを計算する際にSSBは1年後を計算するか(last.catch.zero=TRUEのときのみ有効)
#' @param grid.add.ini \code{add.p.ini}をgridで変えて初期値を事前に探索する
#' @param grid.init \code{p.init}でgridを変えて初期値を事前に探索する
#' @param remove.short.index 年数が指定された数字以下になった指標値を除いて計算する(初期設定:-1)
#' @encoding UTF-8
#' @export
#'

retro.est <- function(res,n=5,stat="mean",init.est=FALSE, b.fix=TRUE,
                      remove.maxAgeF=FALSE,ssb.forecast=FALSE,sel.mat=NULL,
                      grid.add.ini = NULL,grid.init = NULL,
                      remove.short.index=-1){
   res.c <- res
   res.c$input$plot <- FALSE
   Res <- list()
   obj.n <- obj.b <- obj.s <- obj.r <- obj.f <- NULL
   A <- nrow(res$faa)

   if (isTRUE(b.fix)){
     res.c$input$b.fix <- res$b
     res.c$input$b.est <- FALSE
   }

   #if (res$input$last.catch.zero) res.c$input$last.catch.zero <- FALSE
   if (ssb.forecast && !(res.c$input$last.catch.zero)) {
     warning("'ssb.forecast' is usable only when 'last.catch.zero=TRUE' and so ignored")
   }


   for (i in 1:n){
     nc <- ncol(res.c$input$dat$caa)

     res.c$input$dat$caa <- res.c$input$dat$caa[,-nc]
     res.c$input$dat$maa <- res.c$input$dat$maa[,-nc]
     res.c$input$dat$maa.tune <- res.c$input$dat$maa.tune[,-nc]
     res.c$input$dat$waa <- res.c$input$dat$waa[,-nc]
     res.c$input$dat$waa.catch <- res.c$input$dat$waa.catch[,-nc]
     res.c$input$dat$M <- res.c$input$dat$M[,-nc]
     res.c$input$dat$catch.prop <- res.c$input$dat$catch.prop[,-nc]

     if(!is.null(res$input$dat$index)){ # チューニングなしVPAにも対応
       # 毎年等しく取り除くのではなく、データごとに1年分取り除くように修正（浜辺07/07）
       label_tmp <- which(is.na(res.c$input$dat$index[,nc]))
       res.c$input$dat$index <- res.c$input$dat$index[,-nc,drop=FALSE]
       res.c$input$dat$index[label_tmp,length(res.c$input$dat$index[1,])] <- NA
     }

     res.c$input$tf.year <- res.c$input$tf.year-1
     res.c$input$fc.year <- res.c$input$fc.year-1
     if (!is.null(res.c$input$tf.mat)) res.c$input$tf.mat <- res.c$input$tf.mat[,-1]
     if (!is.null(res.c$input$remove.abund)) res.c$input$remove.abund <- res.c$input$remove.abund[,-nc]

     if (!is.null(sel.mat)) {
       if (any(dim(sel.mat) != c(nrow(res$saa),n))) {
         stop("Dimension of 'sel.mat' is not appropriate")
       } else {
         res.c$input$sel.f <- sel.mat[,i] #2stepの方のときチューニングなしの時の選択率を行列で使う
       }
     }
     if (isTRUE(init.est)) res.c$input$p.init <- res.c$term.f

     # last.catch.zero = TRUE用に修正
     if (res.c$input$last.catch.zero) {res.c$input$dat$caa[,nc-1] <- 0; Y <- nc-2} else Y <- nc-1

     if (remove.short.index>0) {
         index_n = apply(res.c$input$dat$index,1,function(x) length(x)-sum(is.na(x)))
         use.index = 1:nrow(res.c$input$dat$index)
         if (res.c$input$use.index[1]=="all") {
           use.index = use.index[index_n > remove.short.index]
         } else {
           use.index = intersect(res.c$input$use.index,use.index[index_n > remove.short.index])
         }
         res.c$input$use.index <- use.index
     }

     # res1 <- do.call(vpa,res.c$input) # do.callからsafe_callに変更(浜辺'20/06/30)

     if (!is.null(grid.add.ini)) {
       input_tmp = res.c$input
       res_list = grid.add.ini %>% map(function(x) {
         input_tmp$no.est = TRUE
         input_tmp$add.p.ini <- x
         RES = do.call(vpa,input_tmp)
       })
       nan_v = sapply(1:10, function(i) sum(is.nan(res_list[[i]]$sigma)))
       loglik_v = sapply(1:10, function(i) res_list[[i]]$logLik)
       pos = which(loglik_v == max(loglik_v[nan_v==0]))
       res.c$input$add.p.ini <- grid.add.ini[pos]
     }

     if (!is.null(grid.init)) {
       input_tmp = res.c$input
       res_list = grid.init %>% map(function(x) {
         input_tmp$no.est = TRUE
         input_tmp$p.init <- x
         RES = do.call(vpa,input_tmp)
       })
       nan_v = sapply(1:10, function(i) sum(is.nan(res_list[[i]]$sigma)))
       loglik_v = sapply(1:10, function(i) res_list[[i]]$logLik)
       pos = which(loglik_v == max(loglik_v[nan_v==0]))
       res.c$input$p.init <- grid.init[pos]
     }

     res1 <- safe_call(vpa,res.c$input, force=TRUE) # do.callからsafe_callに変更(浜辺'20/06/30)

     Res[[i]] <- res1

     if ((max(abs(res1$gradient)) < 10^(-3) & !isTRUE(res1$input$ADMB)) | (max(abs(res1$gradient)) > 0 & max(abs(res1$gradient)) < 10^(-3) & isTRUE(res1$input$ADMB)) | (is.na(max(abs(res1$gradient))) & res1$input$optimizer=="nlminb")){
       obj.n <- c(obj.n, (sum(res1$naa[,Y])-sum(res$naa[,Y]))/sum(res$naa[,Y]))
       obj.b <- c(obj.b, (sum(res1$baa[,Y])-sum(res$baa[,Y]))/sum(res$baa[,Y]))
       if (ssb.forecast && res.c$input$last.catch.zero) {
         obj.s <- c(obj.s, (sum(res1$ssb[,Y+1])-sum(res$ssb[,Y+1]))/sum(res$ssb[,Y+1]))
       } else {
         obj.s <- c(obj.s, (sum(res1$ssb[,Y])-sum(res$ssb[,Y]))/sum(res$ssb[,Y]))
       }
       obj.r <- c(obj.r, (res1$naa[1,Y]-res$naa[1,Y])/res$naa[1,Y])
       if (remove.maxAgeF) {
         obj.f <- c(obj.f, (sum(res1$faa[-A,Y])-sum(res$faa[-A,Y]))/sum(res$faa[-A,Y]))
         } else {
           obj.f <- c(obj.f, (sum(res1$faa[,Y])-sum(res$faa[,Y]))/sum(res$faa[,Y]))
         }
     } else {
       obj.n <- c(obj.n, NA)
       obj.b <- c(obj.b, NA)
       obj.s <- c(obj.s, NA)
       obj.r <- c(obj.r, NA)
       obj.f <- c(obj.f, NA)
     }
   }

   mohn <- c(get(stat)(obj.n,na.rm=TRUE),get(stat)(obj.b,na.rm=TRUE),get(stat)(obj.s,na.rm=TRUE),get(stat)(obj.r,na.rm=TRUE),get(stat)(obj.f,na.rm=TRUE))

   names(mohn) <- c("N","B","SSB","R","F")

   return(list(Res=Res,retro.n=obj.n, retro.b=obj.b, retro.s=obj.s, retro.r=obj.r, retro.f=obj.f, mohn=mohn))
}

#最新年の漁獲量が0の場合 (last.zero.catch=0)、最新年のFが0となり、加入量も推定できないため、もう1年前の推定値でMohn's rhoを計算するための関数
## ------------------------------------------------------  ##
# retro.est中にlast.zero.catch=0の場合のif文加えた
# retro.estが問題なければ以下関数は捨てても大丈夫（浜辺07/07）
## ------------------------------------------------------  ##
retro.est2 <- function(res,n=5,stat="mean",init.est=FALSE, b.fix=TRUE){
   res.c <- res
   res.c$input$plot <- FALSE
   Res <- list()
   obj.n <- obj.b <- obj.s <- obj.r <- obj.f <- NULL

   if (isTRUE(b.fix)){
     res.c$input$b.fix <- res$b
     res.c$input$b.est <- FALSE
   }

   for (i in 1:n){
     nc <- ncol(res.c$input$dat$caa)

     res.c$input$dat$caa <- res.c$input$dat$caa[,-nc]
     res.c$input$dat$maa <- res.c$input$dat$maa[,-nc]
     res.c$input$dat$maa.tune <- res.c$input$dat$maa.tune[,-nc]
     res.c$input$dat$waa <- res.c$input$dat$waa[,-nc]
     res.c$input$dat$M <- res.c$input$dat$M[,-nc]
     res.c$input$dat$index <- res.c$input$dat$index[,-nc,drop=FALSE]
     res.c$input$dat$catch.prop <- res.c$input$dat$catch.prop[,-nc]

     res.c$input$tf.year <- res.c$input$tf.year-1
     res.c$input$fc.year <- res.c$input$fc.year-1
     if (!is.null(res.c$input$tf.mat)) res.c$input$tf.mat <- res.c$input$tf.mat[,-1]

     if (isTRUE(init.est)) res.c$input$p.init <- res.c$term.f

     if (res.c$input$last.catch.zero) {res.c$input$dat$caa[,nc-1] <- 0; Y <- nc-2} else Y <- nc-1

     res1 <- do.call(vpa,res.c$input)

     Res[[i]] <- res1

     if ((max(abs(res1$gradient)) < 10^(-3) & !isTRUE(res1$input$ADMB)) | (max(abs(res1$gradient)) > 0 & max(abs(res1$gradient)) < 10^(-3) & isTRUE(res1$input$ADMB)) | (is.na(max(abs(res1$gradient))) & res1$input$optimizer=="nlminb")){
        obj.n <- c(obj.n, (sum(res1$naa[,Y])-sum(res$naa[,Y]))/sum(res$naa[,Y]))
        obj.b <- c(obj.b, (sum(res1$baa[,Y])-sum(res$baa[,Y]))/sum(res$baa[,Y]))
        obj.s <- c(obj.s, (sum(res1$ssb[,Y])-sum(res$ssb[,Y]))/sum(res$ssb[,Y]))
        obj.r <- c(obj.r, (res1$naa[1,Y]-res$naa[1,Y])/res$naa[1,Y])
        obj.f <- c(obj.f, (sum(res1$faa[,Y])-sum(res$faa[,Y]))/sum(res$faa[,Y]))
     } else {
       obj.n <- c(obj.n, NA)
       obj.b <- c(obj.b, NA)
       obj.s <- c(obj.s, NA)
       obj.r <- c(obj.r, NA)
       obj.f <- c(obj.f, NA)
     }
   }

   mohn <- c(get(stat)(obj.n,na.rm=TRUE),get(stat)(obj.b,na.rm=TRUE),get(stat)(obj.s,na.rm=TRUE),get(stat)(obj.r,na.rm=TRUE),get(stat)(obj.f,na.rm=TRUE))

   names(mohn) <- c("N","B","SSB","R","F")

   return(list(Res=Res,retro.n=obj.n, retro.b=obj.b, retro.s=obj.s, retro.r=obj.r, retro.f=obj.f, mohn=mohn))
}


#' VPAの結果を受けて，バブルプロット(横軸が年，縦軸が年齢，バブルの大きさと色で値を示す）を生成する
#'
#' @param res VPAの出力
#' @param target プロットしたい対象　例）"faa"や"naa"
#' @param years プロットしたい年数　デフォルトはNULLで全年．例）１９９７：１９９９とすると１９９７～１９９９年のみプロット
#' @param fix_ratio x軸とy軸のスケールの比．ここが２だと１：１になる
#' @param max_size バブルの最大の大きさ
#' @param legend_position legendの位置（see `theme_SH`）
#' @encoding UTF-8
#' @export
#'

bubble_plot2 <- function(
res,
target="faa",
years=NULL,
fix_ratio=2, #x軸とy軸のスケールの比.ここが2だと１：１
max_size=10, #バブルの最大の大きさ
legend_position = "bottom" # legendの位置
){

res00 <- res
dat <- res00[names(res00)==paste0(target)]

dat <- dat[[1]] %>% as.data.frame() %>% tibble::rownames_to_column("age") %>% pivot_longer(!age, names_to = "year", values_to = "value") %>% mutate_at(vars(year), as.factor)

dat <- dat %>% mutate_at(vars(age), as.factor)
dat$age <- factor(dat$age, levels=rev(levels(dat$age)))

#range <- range[1:(length(range)-1)]

if(!is.null(years)) dat <- subset(dat, year %in% years)

ggplot(dat,aes(x=year, y=age)) +
  geom_point(aes(size=value, fill=value), alpha = 0.75, shape = 21) +
  coord_fixed(ratio=fix_ratio)+
  scale_size_continuous(limits = c(0.0001, 1.0001)*max(dat$value,na.rm=TRUE), range = c(1,max_size))+
  scale_fill_continuous(trans="reverse")+
  labs(x="Year", y="Age",size=paste0(target),fill=paste0(target),fill=guide_colourbar(reverse=TRUE))  +
  guides(color=c("none"))+
  ggtitle(toupper(target))+
  theme_SH(legend.position = legend_position) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
}


#' リッジVPAのλ探索の自動化
#'
#' @param input `vpa`関数の引数をlist形式で与える
#' @param target_retro ペナルティλの選択に用いるmohn's rhoのパラメータ
#' @param n_retro レトロスペクティブ解析で遡る年数。デフォルトは`5`
#' @param b_fix レトロスペクティブ解析内でbを固定するか。デフォルトは`TRUE`
#' @details etaのオプションに対応できているかは怪しい。
#' @encoding UTF-8
#'
#' @export
#'

autocalc_ridgevpa <- function(input,
                              target_retro,
                              n_retro = 5,
                              b_fix   = TRUE,
                              bin     = 0.1 # lambdaとetaのペナルティ探索の幅
){
  assert_that(target_retro %in% c("F", "B", "N","SSB","R"))

  if(is.null(input$eta)){
    search_lambda_vpa <- function(x){
      print(x)
      input$lambda <- x
      tmp          <- do.call(vpa,input)
      tmp_rho      <- retro.est(tmp, n = n_retro, b.fix = b_fix)$mohn
      print(tmp_rho)
      tmp_rho[names(tmp_rho)==target_retro] %>% as.numeric()
    }# search_lambda_vpa

    lambda_set1 <- seq(0,1,0.1)
    res1 <- purrr::map(as.list(lambda_set1), search_lambda_vpa) %>% as.numeric()
    tmp_min <- which(abs(res1) == min(abs(res1)))
    lambda_set2 <- seq(lambda_set1[tmp_min-1],lambda_set1[tmp_min+1],0.01)
    res2 <- purrr::map(as.list(lambda_set2), search_lambda_vpa) %>% as.numeric()
    tmp_min <- which(abs(res2) == min(abs(res2)))

    lambda_mat1 <- data.frame(lambda = lambda_set1, mohn = res1) %>%
      mutate(delta_mohn = abs(mohn)-min(abs(mohn)))
    lambda_mat2 <- data.frame(lambda = lambda_set2, mohn = res2) %>%
      mutate(delta_mohn = abs(mohn)-min(abs(mohn)))

    g2 <- ggplot(lambda_mat2,
                 aes(x=1, y=lambda, fill=delta_mohn)) +
      geom_tile()+
      geom_hline(aes(yintercept = lambda_mat2$lambda[tmp_min]))+
      xlab("")+
      coord_cartesian(xlim = c(0.75,1.25))
  } #if(is.null(eta))

  if(!is.null(input$eta)){
    search_lambda_vpa <- function(x,y){
	 print(x)
	 print(y)
      input$lambda <- x
      input$eta    <- y
      tmp          <- do.call(vpa,input)
      tmp_rho      <- retro.est(tmp, n = n_retro, b.fix = b_fix)$mohn
	 print(tmp_rho)
      tmp_rho[names(tmp_rho)==target_retro] %>% as.numeric()
    }# search_lambda_vpa

    penalty_set <- expand.grid(lambda = seq(0,1,bin), eta = seq(0,1,bin))
    res1 <- purrr::map2(as.list(penalty_set$lambda), as.list(penalty_set$eta), search_lambda_vpa) %>%
      as.numeric()
    penalty_mat <- penalty_set %>%
      mutate(mohn = res1) %>%
      mutate(delta_mohn = abs(mohn)-min(abs(mohn)))

    g2 <- ggplot(penalty_mat %>%
                   mutate(delta_mohn = ifelse(delta_mohn>1,NA,delta_mohn)),
                 aes(x=eta, y=lambda, fill=delta_mohn)) +
      geom_tile()+
      geom_text(aes(label = as.character(round(delta_mohn, 3))))+
      geom_text(data = penalty_mat %>% filter(delta_mohn == 0),
                aes(label = as.character(round(delta_mohn, 3))), color = "red")+
      scale_fill_continuous(trans="reverse", na.value = "grey")

  } #if(is.null(eta))


  return(
    if(is.null(input$eta)){
      list(min_penalty  = lambda_mat2$lambda[tmp_min],
           plot = g2,
           lambda_mat1 = lambda_mat1, lambda_mat2 = lambda_mat2)
    } else {
      list(min_penalty  = penalty_mat %>% filter(delta_mohn == 0),
           plot = g2, penalty_mat = penalty_mat)
    }
  ) # return()

} #autocalc_ridgevpa
