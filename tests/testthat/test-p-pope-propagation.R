library(frasyr)

context("check p.pope propagation: vpa -> ref.F -> MSY -> future_vpa") # ----

# pmaデータ(tools/generate-testdata/generate_res_vpa_for_test.Rと同じレシピ)を使い、
# p.pope=0.3(非デフォルト)のVPA結果を作って、下流の関数まで一貫してp.popeが伝播することを確認する。
# 「p.popeを変えると結果が変わる」ことだけでなく、変わった値が手計算(該当関数の数式を直接書き下したもの)
# と一致することも確認する。

caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"),row.names=1)
waa <- read.csv(system.file("extdata","waa_pma.csv",package="frasyr"),row.names=1)
maa <- read.csv(system.file("extdata","maa_pma.csv",package="frasyr"),row.names=1)
dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5)

res_vpa_p03 <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                   term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0,p.pope=0.3)
res_vpa_p05 <- vpa(dat,fc.year=2009:2011,rec=585,rec.year=2011,tf.year = 2008:2010,
                   term.F="max",stat.tf="mean",Pope=TRUE,tune=FALSE,p.init=1.0,p.pope=0.5)

test_that("vpa(): p.popeがnaaの計算に反映され、res$input$p.popeに記録される",{
  expect_equal(res_vpa_p03$input$p.pope, 0.3)
  expect_equal(res_vpa_p05$input$p.pope, 0.5)

  # 最高齢×最終年のnaaはcaa*exp(M*p.pope)/(1-exp(-faa))という直接式で計算される(rvpa.r内のPopeの終端条件式)
  na <- nrow(res_vpa_p03$naa)
  ny <- ncol(res_vpa_p03$naa)
  caa_term <- res_vpa_p03$input$dat$caa[na,ny]
  M_term   <- res_vpa_p03$input$dat$M[na,ny]
  faa_term <- res_vpa_p03$faa[na,ny]
  manual_naa <- caa_term * exp(M_term*0.3) / (1-exp(-faa_term))
  expect_equal(res_vpa_p03$naa[na,ny], manual_naa)

  # p.pope=0.5(デフォルト)とは異なる結果になる
  expect_false(isTRUE(all.equal(res_vpa_p03$naa, res_vpa_p05$naa)))
})

test_that("ref.F(): res_vpaのinput$p.popeが自動継承され、calc.rel.abund()の手計算と一致する",{
  Fcur  <- res_vpa_p03$Fc.at.age %>% as.numeric()
  waa_y <- res_vpa_p03$input$dat$waa$"2011"
  maa_y <- res_vpa_p03$input$dat$maa$"2011"
  M_y   <- res_vpa_p03$input$dat$M$"2011"

  # res_vpaを渡した場合、p.popeを明示しなくてもres_vpa$input$p.pope(=0.3)が自動的に使われる
  ref_inherit <- ref.F(res_vpa_p03,Fcurrent=NULL,Pope=TRUE,plot=FALSE,pSPR=NULL)
  ref_explicit_p03 <- ref.F(res=NULL,Fcurrent=Fcur,waa=waa_y,maa=maa_y,M=M_y,waa.catch=waa_y,
                            rps.vector=NULL,Pope=TRUE,p.pope=0.3,min.age=0,max.age=Inf,
                            plot=FALSE,pSPR=NULL)
  ref_explicit_p05 <- ref.F(res=NULL,Fcurrent=Fcur,waa=waa_y,maa=maa_y,M=M_y,waa.catch=waa_y,
                            rps.vector=NULL,Pope=TRUE,p.pope=0.5,min.age=0,max.age=Inf,
                            plot=FALSE,pSPR=NULL)

  # 自動継承(res_vpa_p03)と明示指定(p.pope=0.3)のYPR曲線が完全に一致する
  expect_equal(ref_inherit$ypr.spr$ypr, ref_explicit_p03$ypr.spr$ypr)

  # ref.Fのypr.spr(YPR曲線)は、同じF値でcalc.rel.abund()を直接呼んだ手計算と一致する
  sel <- Fcur/max(Fcur)
  Fval <- ref_explicit_p03$ypr.spr$F.range[5]
  manual <- calc.rel.abund(sel=sel,Fr=Fval,na=length(Fcur),M=M_y,waa=waa_y,waa.catch=waa_y,
                           maa=maa_y,min.age=0,max.age=Inf,Pope=TRUE,p.pope=0.3,ssb.coef=0)
  expect_equal(sum(manual$ypr,na.rm=TRUE), ref_explicit_p03$ypr.spr$ypr[5])

  # p.popeが違えばYPRの値も変わる(Popeの式でのみ影響。MはpmaデータではFmaxの位置自体は動かない)
  expect_false(isTRUE(all.equal(ref_explicit_p03$ypr.spr$ypr, ref_explicit_p05$ypr.spr$ypr)))

  # res=NULLで生パラメータを与える場合、p.popeを省略するとPopeと同様エラーになる
  expect_error(ref.F(res=NULL,Fcurrent=Fcur,waa=waa_y,maa=maa_y,M=M_y,waa.catch=waa_y,
                     rps.vector=NULL,Pope=TRUE,min.age=0,max.age=Inf,plot=FALSE,pSPR=NULL))
})

test_that("MSY計算(est_MSYRP_proxy): p.popeがmake_future_data経由で自動継承され、内部のref.F計算と一致する",{
  SRdata <- get.SRdata(res_vpa_p03)
  res_sr <- fit.SR(SRdata, SR="HS", method="L2", AR=0, hessian=FALSE)

  common_args <- list(res_vpa=res_vpa_p03, res_SR=res_sr,
                      nsim=10, nyear=5,
                      future_initial_year_name=2011, start_F_year_name=2012,
                      start_biopar_year_name=2012, start_random_rec_year_name=2012,
                      waa_year=2009:2011, waa_catch_year=2009:2011,
                      maa_year=2009:2011, M_year=2009:2011, faa_year=2009:2011,
                      start_ABC_year_name=2012, silent=TRUE)

  # p.popeを省略 => res_vpa_p03$input$p.pope(=0.3)を自動継承
  data_future_inherit <- do.call(make_future_data, common_args)
  # p.popeを明示的に0.5に上書き
  data_future_p05      <- do.call(make_future_data, c(common_args, list(p.pope=0.5)))

  expect_equal(data_future_inherit$data$p_pope, 0.3)
  expect_equal(data_future_p05$data$p_pope, 0.5)

  msy_args <- list(Fmsy_proxy_candidate=c("Fmax","F0.1","F%spr"), msy_SPR_candidate=30,
                   Blimit_candidate="Bmin", Bban_candidate="0", select_Btarget="F%spr30")
  res_msy_inherit <- do.call(est_MSYRP_proxy, c(list(data_future=data_future_inherit), msy_args))
  res_msy_p05      <- do.call(est_MSYRP_proxy, c(list(data_future=data_future_p05), msy_args))

  # est_MSYRP_proxy内部でref.Fに渡されたp.popeが正しいことを、同じ生物パラメータで
  # 独立にref.F()を手計算して確認する
  lastyear <- dim(data_future_inherit$data$waa_mat[,,1])[[2]]
  tmp <- data_future_inherit$data$waa_mat[,lastyear,1]!=0
  waa_l       <- data_future_inherit$data$waa_mat[tmp,lastyear,1]
  waa_catch_l <- data_future_inherit$data$waa_catch_mat[tmp,lastyear,1]
  maa_l       <- data_future_inherit$data$maa_mat[tmp,lastyear,1]
  M_l         <- data_future_inherit$data$M_mat[tmp,lastyear,1]
  futureF     <- data_future_inherit$data$faa_mat[tmp,lastyear,1]

  manual_refF_p03 <- ref.F(res=NULL,Fcurrent=futureF,waa=waa_l,maa=maa_l,M=M_l,waa.catch=waa_catch_l,
                           rps.vector=NULL,Pope=TRUE,p.pope=0.3,min.age=0,max.age=Inf,
                           plot=FALSE,pSPR=NULL)
  expect_equal(res_msy_inherit$res_refF$Fmax, manual_refF_p03$Fmax)

  # p.pope=0.3(自動継承)とp.pope=0.5(明示)でMSY管理基準値(Fmaxでの漁獲量)が異なる
  fmax_row_i <- which(res_msy_inherit$summary$RP_name=="Fmax")
  fmax_row_5 <- which(res_msy_p05$summary$RP_name=="Fmax")
  expect_false(isTRUE(all.equal(res_msy_inherit$summary$Catch[fmax_row_i],
                                res_msy_p05$summary$Catch[fmax_row_5])))
})

test_that("future_vpa(): p.popeがmake_future_data経由で将来予測の漁獲量計算に反映され、catch_equation()の手計算と一致する",{
  SRdata <- get.SRdata(res_vpa_p03)
  res_sr <- fit.SR(SRdata, SR="HS", method="L2", AR=0, hessian=FALSE)

  common_args <- list(res_vpa=res_vpa_p03, res_SR=res_sr,
                      nsim=10, nyear=5,
                      future_initial_year_name=2011, start_F_year_name=2012,
                      start_biopar_year_name=2012, start_random_rec_year_name=2012,
                      waa_year=2009:2011, waa_catch_year=2009:2011,
                      maa_year=2009:2011, M_year=2009:2011, faa_year=2009:2011,
                      start_ABC_year_name=2012, silent=TRUE)

  data_future_inherit <- do.call(make_future_data, common_args)
  data_future_p05      <- do.call(make_future_data, c(common_args, list(p.pope=0.5)))

  res_future_inherit <- future_vpa(tmb_data=data_future_inherit$data, optim_method="none", multi_init=1)
  res_future_p05      <- future_vpa(tmb_data=data_future_p05$data,      optim_method="none", multi_init=1)

  # future_vpaが実際に計算した漁獲量が、catch_equation()による手計算(p.pope=0.3)と一致する
  manual_catch <- sum(catch_equation(res_future_inherit$naa[,"2013",1],
                                     res_future_inherit$faa[,"2013",1],
                                     res_future_inherit$waa_catch_mat[,"2013",1],
                                     res_future_inherit$M[,"2013",1],
                                     Pope=1,p.pope=0.3))
  expect_equal(sum(res_future_inherit$wcaa[,"2013",1]), manual_catch)

  # p.pope=0.3(自動継承)とp.pope=0.5(明示)で将来の漁獲量が異なる
  expect_false(isTRUE(all.equal(sum(res_future_inherit$wcaa[,"2013",1]),
                                sum(res_future_p05$wcaa[,"2013",1]))))
})
