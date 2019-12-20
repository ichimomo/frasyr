source("./tools/generate-testdata/future2.1.r")

load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

res_ref_f_pma <- ref.F(res_vpa_pma,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                           waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                           d = 0.001,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                           iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )

save(res_ref_f_pma,file = "./inst/extdata/res_ref_f_pma.rda")
