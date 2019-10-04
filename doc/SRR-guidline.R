## ----setup, include=FALSE, warning=FALSE,eval=FALSE,echo=FALSE-----------
#  # ここのチャンクはもう必要ないので，削除してもOKです
#  library(rmdformats)
#  ## Global options
#  options(max.print="75")
#  opts_chunk$set(echo=TRUE,
#                       cache=TRUE,
#                 prompt=FALSE,
#                 tidy=TRUE,
#                 comment=NA,
#                 message=FALSE,
#                 warning=FALSE)

## ---- warning=FALSE------------------------------------------------------
library(frasyr)
data(res_vpa)
SRdata <- get.SRdata(res_vpa) 

resHS <- fit.SR(SRdata,SR="HS",method="L2",AR=0)
resBH <- fit.SR(SRdata,SR="BH",method="L2",AR=0)
resRI <- fit.SR(SRdata,SR="RI",method="L2",AR=0)

plot(SRdata$R ~ SRdata$SSB, cex=2, type = "b",xlab="SSB",ylab="R",
     main="HS vs. BH vs. RI",ylim=c(0,max(SRdata$R)*1.3),xlim=c(0,max(SRdata$SSB)*1.1))
points(rev(SRdata$SSB)[1],rev(SRdata$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
points(resHS$pred$SSB,resHS$pred$R,col=2,type="l",lwd=3)
points(resBH$pred$SSB,resBH$pred$R,col=3,type="l",lwd=3,lty=2)    
points(resRI$pred$SSB,resRI$pred$R,col=4,type="l",lwd=3,lty=3)
legend("topleft",
       legend=c(sprintf("HS %5.2f",resHS$AICc),sprintf("BH %5.2f",resBH$AICc),sprintf("RI %5.2f",resRI$AICc)),
       lty=1:3,col=2:4,lwd=2,title="AICc",ncol=3)

resSR <- resHS #HSを選択

## ----message=FALSE,warning=FALSE-----------------------------------------
resAR1 <- fit.SR(SRdata,SR="HS",method="L2",AR=1)
resL1 <- fit.SR(SRdata,SR="HS",method="L1",AR=0)

plot(SRdata$R ~ SRdata$SSB, cex=2, type = "b",xlab="SSB",ylab="R",
     main="Effects of autocorrelation and L1",ylim=c(0,max(SRdata$R)*1.3),xlim=c(0,max(SRdata$SSB)*1.1))
points(rev(SRdata$SSB)[1],rev(SRdata$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
points(resSR$pred$SSB,resSR$pred$R,col=2,type="l",lwd=3)
points(resAR1$pred$SSB,resAR1$pred$R,col=3,type="l",lwd=3,lty=2)    
points(resL1$pred$SSB,resL1$pred$R,col=4,type="l",lwd=3,lty=3)
legend("topleft",
       legend=c(sprintf("L2&AR0 %5.2f",resSR$AICc),sprintf("L2&AR1 %5.2f",resAR1$AICc),sprintf("L1&AR0 %5.2f",resL1$AICc)),
       lty=1:3,col=2:4,lwd=2,title="AICc",ncol=3)

resSR <- resL1 #L1 normを採用

## ----warning=FALSE-------------------------------------------------------
check1 <- shapiro.test(resSR$resid)
check2 <- ks.test(resSR$resid,y="pnorm",sd=resSR$pars$sd)

par(mfrow=c(1,2),mar=c(4,4,2,2))
hist(resSR$resid,xlab="Residuals",main="Normality test",freq=FALSE)
X <- seq(min(resSR$resid)*1.3,max(resSR$resid)*1.3,length=200)
points(X,dnorm(X,0,resSR$pars$sd),col=2,lwd=3,type="l")
mtext(text=" P value",adj=1,line=-1,lwd=2,font=2)
mtext(text=sprintf(" SW: %1.3f",check1$p.value),adj=1,line=-2)
mtext(text=sprintf(" KS: %1.3f",check2$p.value),adj=1,line=-3)

qqnorm(resSR$resid2,cex=2)
qqline(resSR$resid2,lwd=3)

## ----warning=FALSE-------------------------------------------------------
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(SRdata$year, resSR$resid2,pch=16,main="",xlab="Year",ylab="Residual")
abline(0,0,lty=2)
par(new=T)
scatter.smooth(SRdata$year, resSR$resid2, lpars=list(col="red", lwd=2),ann=F,axes=FALSE)
ac.res <- acf(resSR$resid2,plot=FALSE)
plot(ac.res,main="",lwd=3)

## ----warning=FALSE-------------------------------------------------------
boot.res <- boot.SR(resSR)

par(mfrow=c(2,2),mar=c(4,4,2,2))
hist(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$a),xlab="",ylab="",main="a")
abline(v=resSR$pars$a,col=2,lwd=3)
abline(v=median(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$a)),col=3,lwd=3,lty=2)
arrows(quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$a),0.1),0,
       quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$a),0.9),0,
       col=4,lwd=3,code=3)
legend("topright",
       legend=c("Estimate","Median","CI(0.8)"),lty=1:2,col=2:4,lwd=2,ncol=1,cex=1)

hist(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$b),xlab="",ylab="",main="b")
abline(v=resSR$pars$b,col=2,lwd=3)
abline(v=median(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$b)),col=3,lwd=3,lty=2)
arrows(quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$b),0.1),0,
       quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$b),0.9),0,
       col=4,lwd=3,code=3)

hist(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$sd),xlab="",ylab="",main="sd")
abline(v=resSR$pars$sd,col=2,lwd=3)
abline(v=median(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$sd)),col=3,lwd=3,lty=2)
arrows(quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$sd),0.1),0,
       quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$sd),0.9),0,
       col=4,lwd=3,code=3)

if (resSR$input$AR==1) {
  hist(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$rho),xlab="",ylab="",main="rho")
  abline(v=resSR$pars$rho,col=2,lwd=3)
  abline(v=median(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$rho)),col=3,lwd=3,lty=2)
  arrows(quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$rho),0.1),0,
         quantile(sapply(1:length(boot.res), function(i) boot.res[[i]]$pars$rho),0.9),0,
         col=4,lwd=3,code=3)
}

par(mfrow=c(1,1))
plot(SRdata$R ~ SRdata$SSB, cex=2, type = "b",xlab="SSB",ylab="R",
     main="Residual bootstrap",ylim=c(0,max(SRdata$R)*1.3))
points(rev(SRdata$SSB)[1],rev(SRdata$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
for (i in 1:length(boot.res)) {
  points(boot.res[[i]]$pred$SSB,boot.res[[i]]$pred$R,type="l",lwd=2,col=rgb(0,0,1,alpha=0.1))
}
points(resSR$pred$SSB,resSR$pred$R,col=2,type="l",lwd=3)

## ----warning=FALSE-------------------------------------------------------
jack.res <- lapply(1:length(SRdata$year), function(i){
  jack <- resSR
  jack$input$w[i] <- 0
  do.call(fit.SR,jack$input)
})

par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(SRdata$year,sapply(1:length(SRdata$year), function(i) jack.res[[i]]$pars$a),type="b",
     xlab="Removed year",ylab="",main="a",pch=19)
abline(resSR$pars$a,0,lwd=3,col=2)

plot(SRdata$year,sapply(1:length(SRdata$year), function(i) jack.res[[i]]$pars$b),type="b",
     xlab="Removed year",ylab="",main="b",pch=19)
abline(resSR$pars$b,0,lwd=3,col=2)

plot(SRdata$year,sapply(1:length(SRdata$year), function(i) jack.res[[i]]$pars$sd),type="b",
     xlab="Removed year",ylab="",main="sd",pch=19)
abline(resSR$pars$sd,0,lwd=3,col=2)

if (resSR$input$AR==1){
  plot(SRdata$year,sapply(1:length(SRdata$year), function(i) jack.res[[i]]$pars$rho),type="b",
       xlab="Removed year",ylab="",main="rho",pch=19)
  abline(resSR$pars$rho,0,lwd=3,col=2)
}

par(mfrow=c(1,1))  
plot(SRdata$R ~ SRdata$SSB, cex=2, type = "b",xlab="SSB",ylab="R",
     main="Jackknife estimate",ylim=c(0,max(SRdata$R)*1.3))
points(rev(SRdata$SSB)[1],rev(SRdata$R)[1],col=1,type="p",lwd=3,pch=16,cex=2)
for (i in 1:length(jack.res)) {
  points(jack.res[[i]]$pred$SSB,jack.res[[i]]$pred$R,type="l",lwd=3,col=rgb(0,0,1,alpha=0.1))
}
points(resSR$pred$SSB,resSR$pred$R,col=2,type="l",lwd=3)


## ----warning=FALSE-------------------------------------------------------
ngrid <- 100
a.grid <- seq(resSR$pars$a*0.5,resSR$pars$a*1.5,length=ngrid)
b.grid <- seq(min(SRdata$SSB),max(SRdata$SSB),length=ngrid)
ba.grids <- expand.grid(b.grid,a.grid)
prof.lik.res <- sapply(1:nrow(ba.grids),function(i) prof.lik(resSR,a=as.numeric(ba.grids[i,2]),b=as.numeric(ba.grids[i,1])))

image(b.grid,a.grid,matrix(prof.lik.res,nrow=ngrid),ann=F,col=cm.colors(12),
      ylim=c(resSR$pars$a*0.5,resSR$pars$a*1.5),xlim=c(min(SRdata$SSB),max(SRdata$SSB)))
par(new=T, xaxs="i",yaxs="i")
contour(b.grid,a.grid,matrix(prof.lik.res,nrow=ngrid),
        ylim=c(resSR$pars$a*0.5,resSR$pars$a*1.5),xlim=c(min(SRdata$SSB),max(SRdata$SSB)),
        xlab="b",ylab="a",main="Profile likelihood")
for(i in 1:length(jack.res)) points(jack.res[[i]]$pars$b,jack.res[[i]]$pars$a,lwd=1,col=1)

lines(y=as.numeric(quantile(sapply(1:length(boot.res),function(i)boot.res[[i]]$pars$a),c(0.1,0.9))),
      x=rep(resSR$pars$b,2),col=4,lwd=2)
lines(x=as.numeric(quantile(sapply(1:length(boot.res),function(i)boot.res[[i]]$pars$b),c(0.1,0.9))),
      y=rep(resSR$pars$a,2),col=4,lwd=2)
legend("bottomleft",c("Bootstrap CI(0.8)","Jackknife"),lty=1:0,pch=c("",1),col=c(4,1),lwd=2:1)

