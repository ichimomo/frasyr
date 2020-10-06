source("./tools/generate-testdata/rvpa1.9.4.r")
#source("./tools/generate-testdata/future2.1.r")
#source("./tools/generate-testdata/stock_recruit.r") #stock_recruit.r 2019/11/26 ver.
source("./tools/generate-testdata/stock_recruit_pulreq446.R") #from stock_recruit.r based on pull request #446

caa <- read.csv("./inst/extdata/caa_pma.csv",row.names=1)
waa <- read.csv("./inst/extdata/waa_pma.csv",row.names=1)
maa <- read.csv("./inst/extdata/maa_pma.csv",row.names=1)

dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)
res_vpa_pma <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0)

SRdata_pma <- get.SRdata(res_vpa_pma)
year <- as.character(max(res_vpa_pma$input$rec.year))

SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
SR.list <- list()

for (i in 1:nrow(SRmodel.list)) {

    resSR_pma <- fit.SR(SRdata_pma, SR = SRmodel.list$SR.rel[i], method = SRmodel.list$L.type[i], AR = SRmodel.list$AR.type[i], out.AR =SRmodel.list$out.AR[i])

  corSR_pma <- corSR(resSR=resSR_pma)

  assign(sprintf("res_corSR_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),corSR_pma)

  savefilenamerescorSR <- sprintf("./inst/extdata/res_corSR_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])

  save(list=paste("res_corSR_pma_",SRmodel.list$SR.rel[i],"_",SRmodel.list$L.type[i],"_AR", SRmodel.list$AR.type[i],"_outAR",as.numeric(SRmodel.list$out.AR[i]), sep=""),file=savefilenamerescorSR)

  res_steepness_pma <- calc_steepness(SR=SRmodel.list$SR.rel,rec_pars=SR.list[[i]]$pars,M=res_vpa_pma$input$dat$M[,year],waa=res_vpa_pma$input$dat$waa[,year],maa=res_vpa_pma$input$dat$maa[,year])

  assign(sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),res_steepness_pma)

  steepness_file <- sprintf("res_steepness_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])

  save(list=paste("res_steepness_pma_",SRmodel.list$SR.rel[i],"_",SRmodel.list$L.type[i],"_AR", SRmodel.list$AR.type[i],"_outAR",as.numeric(SRmodel.list$out.AR[i]), sep=""),file=steepness_file)

}
