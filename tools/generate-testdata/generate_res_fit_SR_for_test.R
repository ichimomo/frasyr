source("./tools/generate-testdata/rvpa1.9.4.r")
#source("./tools/generate-testdata/future2.1.r")
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

# SRタイプ、L1L2、自己相関タイプごとに異なるオブジェクトへ格納
for (i in 1:nrow(SRmodel.list)) {
  assign(sprintf("SRpma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),SR.list[[i]])

  savefilenameresfres <- sprintf("./inst/extdata/SRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])

  save(list=paste("SRpma_",SRmodel.list$SR.rel[i],"_",SRmodel.list$L.type[i],"_AR", SRmodel.list$AR.type[i],"_outAR",as.numeric(SRmodel.list$out.AR[i]), sep=""),file=savefilenameresfres)

}

SRpma_HS_L1_AR1_outAR0_repoptT <- fit.SR(SRdata_pma, SR = "HS", method = "L1",
                                         AR = 1, out.AR =FALSE, hessian = FALSE, rep.opt = TRUE)

savefilenameresfres <- sprintf("./inst/extdata/SRpma_HS_L1_AR1_outAR0_repoptT.rda")
save(SRpma_HS_L1_AR1_outAR0_repoptT,file=savefilenameresfres)
