source("./tools/generate-testdata/rvpa1.9.4.r")
source("./tools/generate-testdata/future2.1.r")

load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))
load(system.file("extdata","SRdata_pma.rda",package = "frasyr"))

SRmodel.list <- expand.grid(SR.rel = c("HS","BH","RI"), AR.type = c(0, 1), out.AR=c(TRUE,FALSE), L.type = c("L1", "L2"))
SR.list <- list()

for (i in 1:nrow(SRmodel.list)) {
  infile.name <- sprintf("SRpma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])
  resfitSR <- load(system.file("extdata",infile.name,package = "frasyr"))

  fittedSR <- eval(parse(text=resfitSR))
  fres_pma_recarg_list <-list(a=fittedSR$pars$a,b=fittedSR$pars$b,
                              rho=fittedSR$pars$rho, # ここではrho=0なので指定しなくてもOK
                              sd=fittedSR$pars$sd,resid=fittedSR$resid)

  selectedrecSR <- sprintf("%s.recAR",fittedSR$input$SR[1])

  fres_pma <-future.vpa(res0=res_vpa_pma,
                              multi=1, # res.pma$Fc.at.ageに掛ける乗数
                              nyear=50, # 将来予測の年数
                              start.year=2012, # 将来予測の開始年
                              N=2, # 確率的計算の繰り返し回数
                              ABC.year=2013, # ABCを計算する年
                              waa.year=2009:2011, # 生物パラメータの参照年
                              maa.year=2009:2011,
                              M.year=2009:2011,
                              is.plot=TRUE, # 結果をプロットするかどうか
                              seed=1,
                              silent=TRUE,
                              recfunc=eval(parse(text=selectedrecSR)), # 再生産関係の関数
                              # recfuncに対する引数
                              rec.arg=fres_pma_recarg_list)

  assign(sprintf("fres_pma_%s_%s_AR%d_outAR%d",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i]),fres_pma)

  savefilenameresfres <- sprintf("./inst/extdata/fres_pma_%s_%s_AR%d_outAR%d.rda",SRmodel.list$SR.rel[i],SRmodel.list$L.type[i], SRmodel.list$AR.type[i],SRmodel.list$out.AR[i])

  save(list=paste("fres_pma_",SRmodel.list$SR.rel[i],"_",SRmodel.list$L.type[i],"_AR", SRmodel.list$AR.type[i],"_outAR",as.numeric(SRmodel.list$out.AR[i]), sep=""),file=savefilenameresfres)
}
