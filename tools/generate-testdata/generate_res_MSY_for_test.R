source("./tools/generate-testdata/rvpa1.9.4.r")
source("./tools/generate-testdata/future2.1.r")
#source("./tools/generate-testdata/future.r")
source("./tools/generate-testdata/utilities-future-vpa.r")
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

# future fcurrent ----
res_future_Fcurrent_pma <- future.vpa(res_vpa_pma,
                                            multi=1,
                                            nyear=58, # 将来予測の年数
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

save(res_future_Fcurrent_pma,file = "./inst/extdata/res_future_Fcurrent_pma.rda")

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
save(res_MSY_pma,file = "./inst/extdata/res_MSY_pma.rda")


# est MSY (grid) F.msyがみつからな----
res_MSY_pma_grid <- est.MSY(res_vpa_pma, # VPAの計算結果
                       res_future_Fcurrent_pma$input, # 将来予測で使用した引数
                       seed=res_future_Fcurrent_pma$input$seed,
                       resid.year=0, # ARありの場合、最近何年分の残差を平均するかをここで指定する。ARありの設定を反映させたい場合必ずここを１以上とすること（とりあえず１としておいてください）。
                       N=100, # 確率的計算の繰り返し回数=>実際の計算では1000~5000回くらいやってください
                       calc.yieldcurve=TRUE,
                       PGY=c(0.95,0.9,0.6,0.1), # 計算したいPGYレベル。上限と下限の両方が計算される
                       onlylower.pgy=FALSE, # TRUEにするとPGYレベルの上限は計算しない（計算時間の節約になる）
                       optim.method="",
                       B0percent=c(0.2,0.3,0.4),
                       Bempirical=c(round(tail(colSums(res_vpa_pma$ssb),n=1)),
                                    round(max(colSums(res_vpa_pma$ssb))),
                                    24000, # 現行Blimit
                                    SRmodel.base$pars$b) # HSの折れ点
)
save(res_MSY_pma_grid,file = "./inst/extdata/res_MSY_pma_grid.rda")

# base ----
refs_all_pma <- res_MSY_pma$summary

refs_all_pma$RP.definition[refs_all_pma$RP_name=="MSY" & refs_all_pma$AR==FALSE] <- "Btarget0"
refs_all_pma$RP.definition[refs_all_pma$RP_name=="PGY_0.9_lower" & refs_all_pma$AR==FALSE] <- "Blow0"
refs_all_pma$RP.definition[refs_all_pma$RP_name=="PGY_0.6_lower" & refs_all_pma$AR==FALSE] <- "Blimit0"
refs_all_pma$RP.definition[refs_all_pma$RP_name=="PGY_0.1_lower" & refs_all_pma$AR==FALSE] <- "Bban0"

refs_all_pma$RP.definition[refs_all_pma$RP_name=="B0-20%" & refs_all_pma$AR==FALSE] <- "Btarget1"  # たとえばBtargetの代替値としてB020%も候補に残しておきたい場合
refs_all_pma$RP.definition[refs_all_pma$RP_name=="PGY_0.95_lower" & refs_all_pma$AR==FALSE] <- "Btarget2"
refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-19431" & refs_all_pma$AR==FALSE] <- "Bcurrent"
refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-63967" & refs_all_pma$AR==FALSE] <- "Bmax"
refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-24000" & refs_all_pma$AR==FALSE] <- "Blimit1"
refs_all_pma$RP.definition[refs_all_pma$RP_name=="Ben-51882" & refs_all_pma$AR==FALSE] <- "B_HS"

refs_all_pma %>% dplyr::select(RP_name,RP.definition)

refs_base_pma <- refs_all_pma %>%
  dplyr::filter(!is.na(RP.definition)) %>% # RP.definitionがNAでないものを抽出
  arrange(desc(SSB)) %>% # SSBを大きい順に並び替え
  dplyr::select(RP.definition,RP_name,SSB,SSB2SSB0,Catch,Catch.CV,U,Fref2Fcurrent)

save(refs_base_pma,file = "./inst/extdata/refs_base_pma.rda")

