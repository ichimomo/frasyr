source("./tools/generate-testdata/rvpa1.9.4.r")
source("./tools/generate-testdata/future2.1.r")
#source("./tools/generate-testdata/future.r")
source("./tools/generate-testdata/utilities.r")
source("./tools/generate-testdata/stock_recruit.r") #stock_recruit.r 2019/11/26 ver.

caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)

dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
res_vpa_pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                   term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

SRdata_pma <- get.SRdata(res_vpa_pma)

SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
SR.list <- list()

for (i in 1:nrow(SRmodel.list)) {
  SR.list[[i]] <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i],
                         AR = SRmodel.list$AR.type[i], out.AR =SRmodel.list$out.AR[i], hessian = FALSE)
}

SRmodel.list$AICc <- sapply(SR.list, function(x) x$AICc)
SRmodel.list$delta.AIC <- SRmodel.list$AICc - min(SRmodel.list$AICc)
SR.list <- SR.list[order(SRmodel.list$AICc)]  # AICの小さい順に並べたもの
(SRmodel.list <- SRmodel.list[order(SRmodel.list$AICc), ]) # 結果

SRmodel.base <- SR.list[[1]] # AIC最小モデルを今後使っていく

selectedSR <- sprintf("%s.recAR",SRmodel.base$input$SR[1])

# calc.yield curve=FALSEオプションの結果オブジェクトを作成（T/FでestMSYのtraceのオブジェクトサイズが違う）
# future fcurrent ----
res_future_Fcurrent_pma <- future.vpa(res_vpa_pma,
                                      multi=1,
                                      nyear=50, # 将来予測の年数
                                      start.year=2012, # 将来予測の開始年
                                      N=100, # 確率的計算の繰り返し回数
                                      ABC.year=2013, # ABCを計算する年
                                      waa.year=2009:2011, # 生物パラメータの参照年
                                      maa.year=2009:2011,
                                      M.year=2009:2011,
                                      is.plot=TRUE, # 結果をプロットするかどうか
                                      seed=1,
                                      silent=TRUE,
                                      recfunc=eval(parse(text=selectedSR))
                                      , # 再生産関係の関数
                                      # recfuncに対する引数
                                      rec.arg=list(a=SRmodel.base$pars$a,b=SRmodel.base$pars$b,
                                                   rho=SRmodel.base$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                   sd=SRmodel.base$pars$sd,resid=SRmodel.base$resid))

# est MSY----
res_MSY_pma <- est.MSY(res_vpa_pma, # VPAの計算結果
                       res_future_Fcurrent_pma$input, # 将来予測で使用した引数
                       seed=res_future_Fcurrent_pma$input$seed,
                       resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                       N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                       calc.yieldcurve=FALSE,
                       PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                       onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                       B0percent=c(0.2,0.3,0.4),
                       Bempirical=c(round(tail(colSums(res_vpa_pma$ssb),n=1)),
                                    round(max(colSums(res_vpa_pma$ssb))),
                                    24000, # 現行Blimit
                                    SRmodel.base$pars$b) # HSの折れ点
)
save(res_MSY_pma,file = "./inst/extdata/res_MSY_pma_calcyieldcurve_F.rda")

# 確率計算の繰り返し数を10のn乗として各々simN回計算
for(n in 2:4){
#n <-2
    sample.size <- 10^n
# future fcurrent ----
res_future_Fcurrent_pma <- future.vpa(res_vpa_pma,
                                            multi=1,
                                            nyear=58, # 将来予測の年数
                                            start.year=2012, # 将来予測の開始年
                                            N=sample.size, # 確率的計算の繰り返し回数
                                            ABC.year=2013, # ABCを計算する年
                                            waa.year=2009:2011, # 生物パラメータの参照年
                                            maa.year=2009:2011,
                                            M.year=2009:2011,
                                            is.plot=TRUE, # 結果をプロットするかどうか
                                            seed=1,
                                            silent=TRUE,
                                            recfunc=eval(parse(text=selectedSR))
                                            , # 再生産関係の関数
                                            # recfuncに対する引数
                                            rec.arg=list(a=SRmodel.base$pars$a,b=SRmodel.base$pars$b,
                                                         rho=SRmodel.base$pars$rho, # ここではrho=0なので指定しなくてもOK
                                                         sd=SRmodel.base$pars$sd,resid=SRmodel.base$resid))

# iteration num
simN <- 100
# est MSY----
for(i in 1:simN){
res_MSY_pma_test <- est.MSY(res_vpa_pma, # VPAの計算結果
                             res_future_Fcurrent_pma$input, # 将来予測で使用した引数
                             #seed=res_future_Fcurrent_pma$input$seed,
                              seed=i,
                             resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                             N=sample.size, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                             calc.yieldcurve=FALSE,
                             PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                             onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                             B0percent=c(0.2,0.3,0.4),
                             Bempirical=c(round(tail(colSums(res_vpa_pma$ssb),n=1)),
                                          round(max(colSums(res_vpa_pma$ssb))),
                                          24000, # 現行Blimit
                                          SRmodel.base$pars$b) # HSの折れ点
)

#outfile <- sprintf("./tools/generate-testdata/res_MSY_pma_%d.rda",i)
#save(res_MSY_pma,file = outfile)
  assign(sprintf("res_MSY_pma%d",i),res_MSY_pma_test)

}

# 結果meanとsdのオブジェクトサイズ確保
load(system.file("extdata","res_MSY_pma_calcyieldcurve_F.rda",package = "frasyr"))
res_MSY_pma_mean <- res_MSY_pma
res_MSY_pma_mean$summary <- c()
res_MSY_pma_mean$summaryAR <- c()
res_MSY_pma_mean$all.stat <- c()
res_MSY_pma_mean$all.statAR <- c()

res_MSY_pma_sd <- res_MSY_pma
res_MSY_pma_sd$summary <- c()
res_MSY_pma_sd$summaryAR <- c()
res_MSY_pma_sd$all.stat <- c()
res_MSY_pma_sd$all.statAR <- c()

#オブジェクト初期化
for(j in 1:length(names(res_MSY_pma$summary_tb))){
  if(is.character(res_MSY_pma$summary_tb[[j]])){
    res_MSY_pma_mean$summary_tb[[j]] <- res_MSY_pma_sd$summary_tb[[j]] <- res_MSY_pma$summary_tb[[j]]
  }else if(is.logical(res_MSY_pma$summary_tb[[j]])){
    res_MSY_pma_mean$summary_tb[[j]] <- res_MSY_pma_sd$summary_tb[[j]] <- res_MSY_pma$summary_tb[[j]]
  }else{
    res_MSY_pma_mean$summary_tb[[j]] <- res_MSY_pma_sd$summary_tb[[j]] <- res_MSY_pma_mean$summary_tb[[j]]*0
  }
}


for(j in 1:length(res_MSY_pma$F.msy)){
  res_MSY_pma_mean$F.msy[[j]] <- res_MSY_pma_sd$F.msy[[j]] <- 0
}

for(j in 1:length(names(res_MSY_pma$all.stat_tb))){
  if(is.character(res_MSY_pma$all.stat_tb[[j]])){
    res_MSY_pma_mean$all.stat_tb[[j]] <- res_MSY_pma_sd$all.stat_tb[[j]] <- res_MSY_pma$all.stat_tb[[j]]
  }else if(is.logical(res_MSY_pma$all.stat_tb[[j]])){
    res_MSY_pma_mean$all.stat_tb[[j]] <- res_MSY_pma_sd$all.stat_tb[[j]] <- res_MSY_pma$all.stat_tb[[j]]
  }else{
    res_MSY_pma_mean$all.stat_tb[[j]] <- res_MSY_pma_sd$all.stat_tb[[j]] <- res_MSY_pma_mean$all.stat_tb[[j]]*0
  }
}
for(j in 1:length(names(res_MSY_pma$trace))){
  res_MSY_pma_mean$trace[[j]] <- res_MSY_pma_sd$trace[[j]] <- c(0)
}
for(j in 1:length(res_MSY_pma$ssb.ar.mean)){
  res_MSY_pma_mean$ssb.ar.mean[j] <- res_MSY_pma_sd$ssb.ar.mean[j] <- c(0)
}
res_MSY_pma_mean$SPR.msy <- res_MSY_pma_sd$SPR.msy <- 0

#for(i in 1:simN){
#  readlile <- sprintf("./tools/generate-testata/res_MSY_pma%d.rda",i)
#  assign(sprintf("res_vpa_pma%d",i), res_vpa_pma)
#  }

checkitem <-c("summary_tb","all.stat_tb","trace")
for(l in 1:length(checkitem)){
  tmp <- list()
  for(i in 1:simN){
    tmp[[i]] <- as.data.frame(eval(parse(text=paste(sprintf("res_MSY_pma%d$%s",i,checkitem[l])))))
  }

  for(k in 1:length(eval(parse(text=paste(sprintf("res_MSY_pma$%s",checkitem[l])))))){
    each.tmp <- c()
    for(i in 1:simN){
      each.tmp <- cbind(each.tmp,tmp[[i]][k][,])
    }
    #print(sprintf("res_MSY_pma$%s",checkitem[l]))
    if(is.numeric(tmp[[i]][k][,])){
      mean.tmp <- c()
      sd.tmp <- c()
      for(j in 1:length(tmp[[i]][k][,])){
        mean.tmp <- rbind(mean.tmp,mean(each.tmp[j,]))
        sd.tmp <- rbind(sd.tmp, sd(each.tmp[j,]))
      }
      if(l == 1){
        res_MSY_pma_mean$summary_tb[k] <- mean.tmp
        res_MSY_pma_sd$summary_tb[k] <- sd.tmp
      } else if(l == 2){
        res_MSY_pma_mean$all.stat_tb[k] <- mean.tmp
        res_MSY_pma_sd$all.stat_tb[k] <- sd.tmp
      } else if(l == 3){
        res_MSY_pma_mean$trace[k] <- mean.tmp
        res_MSY_pma_sd$trace[k] <- sd.tmp
      }

    }
  }
}

checkitem <-c("F.msy")

for(l in 1:length(checkitem)){
  tmp <- list()
  each.tmp <-c()
  for(i in 1:simN){
    tmp[[i]] <- as.data.frame(eval(parse(text=paste(sprintf("res_MSY_pma%d$%s",i,checkitem[l])))))
    if(is.numeric(tmp[[i]][,])){
    each.tmp <- cbind(each.tmp,tmp[[i]][,])
    }
  }
  mean.tmp <- c()
  sd.tmp <- c()
  for(j in 1:length(tmp[[i]][,])){
    mean.tmp <- rbind(mean.tmp,mean(each.tmp[j,]))
    sd.tmp <- rbind(sd.tmp, sd(each.tmp[j,]))
  }
  res_MSY_pma_mean$F.msy <- mean.tmp
  res_MSY_pma_sd$F.msy <- sd.tmp
}

checkitem <-c("ssb.ar.mean")
for(l in 1:length(checkitem)){
  tmp <- list()
  for(i in 1:simN){
    tmp[[i]] <- as.data.frame(eval(parse(text=paste(sprintf("res_MSY_pma%d$%s",i,checkitem[l])))))
  }
  mean.mat <- c()
  sd.mat <- c()
    for(k in 1:ncol(eval(parse(text=paste(sprintf("res_MSY_pma$%s",checkitem[l])))))){
      each.tmp <-c()
      for(i in 1:simN){
        each.tmp <- cbind(each.tmp,tmp[[i]][k][,])
      }
      mean.tmp <- c()
      sd.tmp <- c()
      for(j in 1:nrow(tmp[[i]][k])){
        mean.tmp <- rbind(mean.tmp,mean(each.tmp[j,]))
        sd.tmp <- rbind(sd.tmp, sd(each.tmp[j,]))
      }
    mean.mat <- cbind(mean.mat,mean.tmp)
    sd.mat <- cbind(sd.mat, sd.tmp)
  }

  res_MSY_pma_mean$ssb.ar.mean <- mean.mat
  res_MSY_pma_sd$ssb.ar.mean <- sd.mat
}

checkitem <-c("SPR.msy")
tmp <- list()
for(i in 1:simN){
  tmp[[i]] <- as.data.frame(eval(parse(text=paste(sprintf("res_MSY_pma%d$%s",i,checkitem)))))
}
each.tmp <-c()
for(i in 1:simN){
  each.tmp <- cbind(each.tmp,tmp[[i]][,])
}

mean.tmp <- mean(each.tmp)
sd.tmp <- sd(each.tmp)

res_MSY_pma_mean$SPR.msy <-mean.tmp
res_MSY_pma_sd$SPR.msy <-sd.tmp

outmeanfile <- sprintf("./inst/extdata/res_MSY_pma_mean_samplesize%d.rda",10^n)
outsdfile <- sprintf("./inst/extdata/res_MSY_pma_sd_samplesize%d.rda",10^n)

save(res_MSY_pma_mean,file=outmeanfile)
save(res_MSY_pma_sd,file=outsdfile)

}

