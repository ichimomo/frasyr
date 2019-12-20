source("./tools/generate-testdata/rvpa1.9.4.r")
source("./tools/generate-testdata/future2.1.r")

caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)

dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
res_vpa_pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                   term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

SRdata_pma <- get.SRdata(res_vpa_pma)

# モデル選択 by AICc----
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

SRmodel.select <- SR.list[[1]] # AIC最小モデルを今後使っていく

# プロットは無視してAICcのみを判断材料としてモデル診断（pmaではHS L2 AR0 (outAR=Tのほうが選ばれているが、Fでも同じはず)）

# ignore plot-----
#plot(SRdata_pma$R ~ SRdata_pma$SSB, cex=2, type = "b",xlab="SSB",ylab="R",
#     main="HS vs. BH vs. RI",ylim=c(0,max(SRdata_pma$R)*1.3),xlim=c(0,max(SRdata_pma$SSB)*1.1))
#points(rev(SRdata_pma$SSB)[1],rev(SRdata_pma$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
#points(resHS$pred$SSB,resHS$pred$R,col=2,type="l",lwd=3)
#points(resBH$pred$SSB,resBH$pred$R,col=3,type="l",lwd=3,lty=2)
#points(resRI$pred$SSB,resRI$pred$R,col=4,type="l",lwd=3,lty=3)
#legend("topleft",
#       legend=c(sprintf("HS %5.2f",resHS$AICc),sprintf("BH %5.2f",resBH$AICc),sprintf("RI %5.2f",resRI$AICc)),
#       lty=1:3,col=2:4,lwd=2,title="AICc",ncol=3)

# 正規性チェック----
srr_check1_pma <- shapiro.test(SRmodel.select$resid)
srr_check2_pma <- ks.test(SRmodel.select$resid,y="pnorm")

# ignore plot ----
#par(mfrow=c(1,2),mar=c(4,4,2,2))
#hist(SRmodel.select$resid,xlab="Residuals",main="Normality test",freq=FALSE)
#X <- seq(min(SRmodel.select$resid)*1.3,max(SRmodel.select$resid)*1.3,length=200)
#points(X,dnorm(X,0,SRmodel.select$pars$sd),col=2,lwd=3,type="l")
#mtext(text=" P value",adj=1,line=-1,lwd=2,font=2)
#mtext(text=sprintf(" SW: %1.3f",check1$p.value),adj=1,line=-2)
#mtext(text=sprintf(" KS: %1.3f",check2$p.value),adj=1,line=-3)
#qqnorm(SRmodel.select$resid2,cex=2)
#qqline(SRmodel.select$resid2,lwd=3)

save(srr_check1_pma, file="./inst/extdata/srr_check1_pma.rda")
save(srr_check2_pma, file="./inst/extdata/srr_check2_pma.rda")

# residual trend & acf ----

# ignore plot----
#par(mfrow=c(1,2),mar=c(4,4,2,2))
#plot(SRdata_pma$year, SRmodel.select$resid2,pch=16,main="",xlab="Year",ylab="Residual")
#abline(0,0,lty=2)
#par(new=T)
#scatter.smooth(SRdata_pma$year, SRmodel.select$resid2, lpars=list(col="red", lwd=2),ann=F,axes=FALSE)

ac_res_pma <- acf(SRmodel.select$resid2,plot=FALSE)
#plot(ac.res_pma,main="",lwd=3)

save(ac_res_pma,file = "./inst/extdata/ac_res_pma.rda")

# boot strap ----
boot_res_pma <- boot.SR(SRmodel.select)

# ignore plot
#par(mfrow=c(2,2),mar=c(4,4,2,2))
#hist(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$a),xlab="",ylab="",main="a")
#abline(v=SRmodel.select$pars$a,col=2,lwd=3)
#abline(v=median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$a)),col=3,lwd=3,lty=2)
#arrows(quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$a),0.1),0,
#       quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$a),0.9),0,
#       col=4,lwd=3,code=3)
#legend("topright",
#       legend=c("Estimate","Median","CI(0.8)"),lty=1:2,col=2:4,lwd=2,ncol=1,cex=1)

#hist(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$b),xlab="",ylab="",main="b")
#abline(v=SRmodel.select$pars$b,col=2,lwd=3)
#abline(v=median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$b)),col=3,lwd=3,lty=2)
#  arrows(quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$b),0.1),0,
#       quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$b),0.9),0,
#       col=4,lwd=3,code=3)

#hist(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$sd),xlab="",ylab="",main="sd")
#abline(v=SRmodel.select$pars$sd,col=2,lwd=3)
#abline(v=median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$sd)),col=3,lwd=3,lty=2)
#arrows(quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$sd),0.1),0,
#       quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$sd),0.9),0,
#       col=4,lwd=3,code=3)

#if (SRmodel.select$input$AR==1) {
#  hist(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$rho),xlab="",ylab="",main="rho")
#  abline(v=SRmodel.select$pars$rho,col=2,lwd=3)
#  abline(v=median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$rho)),col=3,lwd=3,lty=2)
#  arrows(quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$rho),0.1),0,
#         quantile(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$rho),0.9),0,
#         col=4,lwd=3,code=3)
#}
#par(mfrow=c(1,1))

boot_res_pma_median_pars_a <-median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$a))
boot_res_pma_median_pars_b <-median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$b))
boot_res_pma_median_pars_sd <-median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$sd))
boot_res_pma_median_pars_rho <-median(sapply(1:length(boot_res_pma), function(i) boot_res_pma[[i]]$pars$rho))

boot_res_pma_median_pars <- cbind(boot_res_pma_median_pars_a,boot_res_pma_median_pars_b,boot_res_pma_median_pars_sd,boot_res_pma_median_pars_rho)

save(boot_res_pma_median_pars,file = "./inst/extdata/boot_res_pma_median_pars.rda")

# jack knife (関数を定義しているわけではないので無視)----

#jack_res_pma <- lapply(1:length(SRdata_pma$year), function(i){
#  jack_pma <- SRmodel.select
#  jack_pma$input$w[i] <- 0
#  do.call(fit.SR,jack_pma$input)
#})

# ignore plot
#par(mfrow=c(2,2),mar=c(4,4,2,2))
#plot(SRdata_pma$year,sapply(1:length(SRdata_pma$year), function(i) jack_res_pma[[i]]$pars$a),type="b",
#     xlab="Removed year",ylab="",main="a",pch=19)
#abline(SRmodel.select$pars$a,0,lwd=3,col=2)

#plot(SRdata_pma$year,sapply(1:length(SRdata_pma$year), function(i) jack_res_pma[[i]]$pars$b),type="b",
#     xlab="Removed year",ylab="",main="b",pch=19)
#abline(SRmodel.select$pars$b,0,lwd=3,col=2)

#plot(SRdata_pma$year,sapply(1:length(SRdata_pma$year), function(i) jack_res_pma[[i]]$pars$sd),type="b",
#     xlab="Removed year",ylab="",main="sd",pch=19)
#abline(SRmodel.select$pars$sd,0,lwd=3,col=2)

#if (SRmodel.select$input$AR==1){
#  plot(SRdata_pma$year,sapply(1:length(SRdata_pma$year), function(i) jack_res_pma[[i]]$pars$rho),type="b",
#       xlab="Removed year",ylab="",main="rho",pch=19)
#  abline(SRmodel.select$pars$rho,0,lwd=3,col=2)
#}
#par(mfrow=c(1,1))

#jack_res_pma_pars <- c()
#for(i in 1:length(SRdata_pma$year)){
#  jack_res_pma_pars <- rbind(jack_res_pma_pars,jack_res_pma[[i]]$pars)
#}
#save(jack_res_pma_pars,file="./inst/extdata/jack_res_pma_pars.rda")

# profile likelihood (関数を定義しているわけではないので無視)----

#ngrid <- 100
#a.grid <- seq(SRmodel.select$pars$a*0.5,SRmodel.select$pars$a*1.5,length=ngrid)
#b.grid <- seq(min(SRdata_pma$SSB),max(SRdata_pma$SSB),length=ngrid)
#ba.grids <- expand.grid(b.grid,a.grid)
#prof_lik_res_pma <- sapply(1:nrow(ba.grids),function(i) prof.lik(SRmodel.select,a=as.numeric(ba.grids[i,2]),b=as.numeric(ba.grids[i,1])))

#save(prof_lik_res_pma,file = "./inst/extdata/prof_lik_res_pma.rda")

# ignore plot----
#image(b.grid,a.grid,matrix(prof_lik_res_pma,nrow=ngrid),ann=F,col=cm.colors(12),
#      ylim=c(SRmodel.select$pars$a*0.5,SRmodel.select$pars$a*1.5),xlim=c(min(SRdata_pma$SSB),max(SRdata_pma$SSB)))
#par(new=T, xaxs="i",yaxs="i")
#contour(b.grid,a.grid,matrix(prof_lik_res_pma,nrow=ngrid),
#        ylim=c(SRmodel.select$pars$a*0.5,SRmodel.select$pars$a*1.5),xlim=c(min(SRdata_pma$SSB),max(SRdata_pma$SSB)),
#        xlab="b",ylab="a",main="Profile likelihood")
#for(i in 1:length(jack_res_pma)) points(jack_res_pma[[i]]$pars$b,jack_res_pma[[i]]$pars$a,lwd=1,col=1)

#lines(y=as.numeric(quantile(sapply(1:length(boot_res_pma),function(i)boot_res_pma[[i]]$pars$a),c(0.1,0.9))),
#      x=rep(SRmodel.select$pars$b,2),col=4,lwd=2)
#lines(x=as.numeric(quantile(sapply(1:length(boot_res_pma
#                                            ),function(i)boot_res_pma[[i]]$pars$b),c(0.1,0.9))),
#      y=rep(SRmodel.select$pars$a,2),col=4,lwd=2)
#legend("bottomleft",c("Bootstrap CI(0.8)","Jackknife"),lty=1:0,pch=c("","○"),col=c(4,1),lwd=2:1)
