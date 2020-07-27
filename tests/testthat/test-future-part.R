library(frasyr)

# ref.F test ----
context("check ref.F") 

test_that("ref.F (level 2)",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  res_ref_f_pma_check <- ref.F(res_vpa_pma,Fcurrent=NULL,waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                               waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                               d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,pSPR = seq(10,90,by=10),
                               iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )

  res.ref.f_2times <- ref.F(res_vpa_pma,Fcurrent=res_vpa_pma$Fc.at.age*2,
                            waa=NULL,maa=NULL,M=NULL,waa.catch=NULL,M.year=NULL,
                            waa.year=NULL,maa.year=NULL,rps.year = NULL,max.age = Inf,min.age = 0,
                            d = 0.001,Fem.init = 0.5,Fmax.init = 1.5,F0.1.init = 0.7,
                            pSPR = c(seq(10,90,by=10),res_ref_f_pma_check$currentSPR$perSPR*100),
                            iterlim=1000,plot=TRUE,Pope=FALSE,F.range = seq(from=0,to=2,length=101) )


  times2_check <- as.numeric(table(unlist(res_ref_f_pma_check$summary/res.ref.f_2times$summary[,1:16])))

  expect_equal(times2_check[1],2)
  expect_equal(times2_check[2],31)
  expect_equal(times2_check[3],15)
  expect_equal(round(res.ref.f_2times$summary[3,17],2),0.5)

  #上記引数での計算結果を読み込み
  load(system.file("extdata","res_ref_f_pma.rda",package = "frasyr"))

  #照合内容
  testcontents <-c("sel","max.age","min.age","rps.q","spr.q","Fcurrent","Fmed","Flow","Fhigh","Fmax","F0.1","Fmean","rps.data","FpSPR","summary","ypr.spr","waa","waa.catch","maa","spr0")

  for(i in 1:length(testcontents)){
      tmp1 <- eval(parse(text=paste("res_ref_f_pma$",testcontents[i])))
      tmp2 <- eval(parse(text=paste("res_ref_f_pma_check$",testcontents[i])))
      expect_equivalent(tmp1,tmp2,tolerance=1e-4,label=testcontents[i])
      # %SPRを計算するところで、初期値が変わると1e-4以下の誤差で値がずれるので1e-4をtoleranceに入れる
      # toleranceのつづりが間違っていても誰も教えてくれない（無言でtoleranceを無視する）ため注意
  }

  # 同じ機能を持つcalc_Fratioとの整合性をチェック=> ref.Fとcalc_Fratioは同じ機能を提供
  Fratio_test <- purrr::map_dfr((1:4) * 10, function(x)
    calc_Fratio(res_vpa_pma$Fc.at.age,
                rev(res_vpa_pma$input$dat$waa)[,1],
                rev(res_vpa_pma$input$dat$maa)[,1],
                rev(res_vpa_pma$input$dat$M  )[,1],
                SPRtarget=x,
                waa.catch=NULL,
                Pope=res_vpa_pma$input$Pope,
                return_SPR=TRUE))
  for_test_tmp <- 1/res_ref_f_pma_check$summary[str_c("FpSPR.",1:4 * 10,".SPR")][3,] %>%
    unlist() %>% as.numeric()
  expect_equal(for_test_tmp,Fratio_test$Fratio,tol=0.0001)
  expect_equal(1:4 * 10,Fratio_test$SPR_est,tol=0.0001)  

  MSY_perSPR1 <- calc_future_perSPR(res_future_0.8HCR,res_vpa_example,res_MSY$F.msy)  
  MSY_perSPR2 <- calc_future_perSPR(res_future_0.8HCR,res_vpa_example,res_MSY$F.msy,target.year=2040:2045)
  MSY_perSPR3 <- calc_future_perSPR(res_future_0.8HCR,res_vpa_example,res_MSY$F.msy,target.col=30) 
  expect_equal(MSY_perSPR1,0.2307774,tol=1e-4)
  expect_equal(MSY_perSPR1,MSY_perSPR2,tol=1e-4)      
  expect_equal(MSY_perSPR1,MSY_perSPR3,tol=1e-4)

  # just for run
  MSY_Fratio <- calc_perspr(res_future_0.8HCR,res_vpa_example,res_MSY$F.msy,SPRtarget=0.3)

})

# check SR data ----
context("check get.SRdata")

test_that("get.SRdata (level 2)",{
  load(system.file("extdata","res_vpa_pma.rda",package = "frasyr"))

  SRdata0 <- get.SRdata(R.dat=exp(rnorm(10)),SSB.dat=exp(rnorm(10)))
  SRdata0usingPeriodFrom1990To2000 <- get.SRdata(res_vpa_pma,years=1990:2000)

  #上記引数での計算結果を読み込み
  SRdata0_pma_check <- read.csv(system.file("extdata","future_SRdata0_pma_check.csv",package="frasyr"),row.names=1)
  SRdata0usingPeriodFrom1990To2000_pma_check <- read.csv(system.file("extdata","future_SRdata0usingPeriodFrom1990To2000_pma_check.csv",package="frasyr"),row.names=1)

  #結果の数値を照合

  expect_equal(SRdata0$year, SRdata0_pma_check$year)
  #expect_equal(SRdata0$SSB, SRdata0_pma_check$SSB)
  #expect_equal(SRdata0$R, SRdata0_pma_check$R)

  expect_equal(SRdata0usingPeriodFrom1990To2000$year, SRdata0usingPeriodFrom1990To2000_pma_check$year)
  expect_equal(SRdata0usingPeriodFrom1990To2000$SSB, SRdata0usingPeriodFrom1990To2000_pma_check$SSB)
  expect_equal(SRdata0usingPeriodFrom1990To2000$R, SRdata0usingPeriodFrom1990To2000_pma_check$R)

})

# check future_vpa with sample data ----
context("check future_vpa with sample data") 

test_that("future_vpa function (with sample vpa data) (level 2)",{
  
  data(res_vpa)
  data(res_sr_HSL2)
  
  data_future_test <- make_future_data(res_vpa, # VPAの結果
                                       nsim = 1000, # シミュレーション回数
                                       nyear = 50, # 将来予測の年数
                                       future_initial_year_name = 2017, # 年齢別資源尾数を参照して将来予測をスタートする年
                                       start_F_year_name = 2018, # この関数で指定したFに置き換える最初の年
                                       start_biopar_year_name=2018, # この関数で指定した生物パラメータに置き換える最初の年
                                       start_random_rec_year_name = 2018, # この関数で指定した再生産関係からの加入の予測値に置き換える最初の年
                                       # biopar setting
                                       waa_year=2015:2017, waa=NULL, # 将来の年齢別体重の設定。過去の年を指定し、その平均値を使うか、直接ベクトルで指定するか。以下も同じ。
                                       waa_catch_year=2015:2017, waa_catch=NULL,
                                       maa_year=2015:2017, maa=NULL,
                                       M_year=2015:2017, M=NULL,
                                       # faa setting
                                       faa_year=2015:2017, # currentF, futureFが指定されない場合だけ有効になる。将来のFを指定の年の平均値とする
                                       currentF=NULL,futureF=NULL, # 将来のABC.year以前のFとABC.year以降のFのベクトル 
                                       # HCR setting (not work when using TMB)
                                       start_ABC_year_name=2019, # HCRを適用する最初の年
                                       HCR_beta=1, # HCRのbeta
                                       HCR_Blimit=-1, # HCRのBlimit
                                       HCR_Bban=-1, # HCRのBban
                                       HCR_year_lag=0, # HCRで何年遅れにするか
                                       # SR setting
                                       res_SR=res_sr_HSL2, # 将来予測に使いたい再生産関係の推定結果が入っているfit.SRの返り値
                                       seed_number=1, # シード番号
                                       resid_type="lognormal", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                                       resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                                       bias_correction=TRUE, # バイアス補正をするかどうか
                                       recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                                       # Other
                                       Pope=res_vpa$input$Pope,
                                       fix_recruit=list(year=c(2020,2021),rec=c(1000,2000)),
                                       fix_wcatch=list(year=c(2020,2021),wcatch=c(1000,2000))		 
  ) 

  # backward-resamplingの場合
  data_future_backward <- make_future_data(res_vpa, # VPAの結果
                                           nsim = 1000, # シミュレーション回数
                                           nyear = 50, # 将来予測の年数
                                           future_initial_year_name = 2017, # 年齢別資源尾数を参照して将来予測をスタートする年
                                           start_F_year_name = 2018, # この関数で指定したFに置き換える最初の年
                                           start_biopar_year_name=2018, # この関数で指定した生物パラメータに置き換える最初の年
                                           start_random_rec_year_name = 2018, # この関数で指定した再生産関係からの加入の予測値に置き換える最初の年
                                           # biopar setting
                                           waa_year=2015:2017, waa=NULL, # 将来の年齢別体重の設定。過去の年を指定し、その平均値を使うか、直接ベクトルで指定するか。以下も同じ。
                                           waa_catch_year=2015:2017, waa_catch=NULL,
                                           maa_year=2015:2017, maa=NULL,
                                           M_year=2015:2017, M=NULL,
                                           # faa setting
                                           faa_year=2015:2017, # currentF, futureFが指定されない場合だけ有効になる。将来のFを指定の年の平均値とする
                                           currentF=NULL,futureF=NULL, # 将来のABC.year以前のFとABC.year以降のFのベクトル 
                                           # HCR setting (not work when using TMB)
                                           start_ABC_year_name=2019, # HCRを適用する最初の年
                                           HCR_beta=1, # HCRのbeta
                                           HCR_Blimit=-1, # HCRのBlimit
                                           HCR_Bban=-1, # HCRのBban
                                           HCR_year_lag=0, # HCRで何年遅れにするか
                                           # SR setting
                                           res_SR=res_sr_HSL2, # 将来予測に使いたい再生産関係の推定結果が入っているfit.SRの返り値
                                           seed_number=1, # シード番号
                                           resid_type="backward", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                                           resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                                           backward_duration=5,
                                           bias_correction=TRUE, # バイアス補正をするかどうか
                                           recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                                           # Other
                                           Pope=res_vpa$input$Pope,
                                           fix_recruit=list(year=c(2020,2021),rec=c(1000,2000)),
                                           fix_wcatch=list(year=c(2020,2021),wcatch=c(1000,2000))		 
  )     
  
  # 単なる将来予測の場合(simple)
  res_future_test <- future_vpa(tmb_data=data_future_test$data,
                                optim_method="none", 
                                multi_init = 1) 
  # option fix_recruit、fix_wcatchのチェック
  catch <- apply(res_future_test$wcaa,c(2,3),sum)
  expect_equal(mean(res_future_test$naa[1,"2020",]), 1000)
  expect_equal(mean(catch["2020",]), 1000, tol=0.001)
  expect_equal(mean(catch["2021",]), 2000, tol=0.001)
  
  # 単なる将来予測の場合(backward)
  res_future_backward <- future_vpa(tmb_data=data_future_backward$data, 
                                    optim_method="none", 
                                    multi_init = 1) 
  # option fix_recruit、fix_wcatchのチェック
  expect_equal(mean(res_future_backward$naa[1,"2020",]), 1000)
  catch <- apply(res_future_backward$wcaa,c(2,3),sum)
  expect_equal(mean(catch["2020",]), 1000, tol=0.001)
  expect_equal(mean(catch["2021",]), 2000, tol=0.001)
  
  # MSY計算の場合(MSY estimation)
  res_future_test_R <- future_vpa(tmb_data=data_future_test$data, 
                                  optim_method="R", 
                                  multi_init  = 1,
                                  multi_lower = 0.001, multi_upper = 5,
                                  objective="MSY")
  # [1] 0.5269326
  expect_equal(round(res_future_test_R$multi,3),0.527)
  
  # MSY計算の場合(MSY estimation)
  res_future_test_backward <- future_vpa(tmb_data=data_future_backward$data, 
                                  optim_method="R", 
                                  multi_init  = 1,
                                  multi_lower = 0.001, multi_upper = 5,
                                  objective="MSY")
  expect_equal(round(res_future_test_backward$multi,3),0.525)
  
  
  if(sum(installed.packages()[,1]=="TMB")){
      # res_future_test_tmb <- future_vpa(tmb_data=data_future_test$data,
      #                                  optim_method="tmb", 
      #                                  multi_init  = 1,
      #                                  multi_lower = 0.001, multi_upper = 5,
      #                                  objective="MSY")
      # expect_equal(round(res_future_test_tmb$multi,3),0.527)
  }
  
  
})

# check future_vpa with dummy data ----
context("check future_vpa_function2 with dummy data") 

test_that("future_vpa function (with dummy vpa data) (level 2-3?)",{

  # test-data.handler.Rで作成したVPAオブジェクトを読み込んでそれを使う  
  load("res_vpa_files.rda")
  
  # estimate SR function ----
  # VPA結果がほとんど同じになるres_vpa_base0_nontune, res_vpa_base1_nontune, res_vpa_rec0$nontueの
  # パラメータ推定値は同じになる。
  vpa_list <- tibble::lst(res_vpa_base0_nontune,
                          res_vpa_base1_nontune,
                          res_vpa_pgc0_nontune,
                          res_vpa_rec0_nontune)
  res_sr_list <- list()
  for(i in 1:length(vpa_list)){
      x <- vpa_list[[i]]
      res_sr <- res_sr_list[[i]] <- get.SRdata(x) %>% fit.SR(AR=0, SR="HS")
      # SB、Rが同じにならない（SD>0）ケースは単純テストから除く
      if(res_sr$pars$sd < 0.001){
        const_ssr <- mean(colSums(x$ssb))
        const_R   <- mean(unlist(x$naa[1,]))
        expect_equal(res_sr$pars$sd, 0, tol=0.000001)
        expect_equal(res_sr$pars$b, const_ssr, tol=0.000001)
        expect_equal(const_R/const_ssr, res_sr$pars$a)
  }}

  names(res_sr_list) <- names(vpa_list)
  
  # res_sr_listの２つのパラメータ推定値はほとんど同じだが、res_vpa_rec0のほうは加入年齢が１歳からなので、再生産関係に使うデータ（SRdata）は一点少ない
  expect_equal(length(res_sr_list$res_vpa_base0_nontune$input$SRdata$SSB),
               length(res_sr_list$res_vpa_rec0_nontune$input$SRdata$SSB)+1)           
  
  # future projection with dummy data ----
  
  max_vpa_year <- max(as.numeric(colnames(res_vpa_base0_nontune$naa)))
  bio_year <- rev(as.numeric(colnames(res_vpa_base0_nontune$naa)))[1:3]
  
  target_eq_naa <- 12
  f <- function (x) 4+4*x+4*x^2+4*x^3 - target_eq_naa
  x <- uniroot(f, c(0, 1))$root
  Fvalue <- -log(x)
  data_future_test <- 
    make_future_data(res_vpa_base0_nontune, 
                     nsim = 10,
                     nyear = 20, 
                     future_initial_year_name = max_vpa_year, # 年齢別資源尾数を参照して将来予測をスタートする年
                     start_F_year_name = max_vpa_year+1, # この関数で指定したFに置き換える最初の年
                     start_biopar_year_name=max_vpa_year+1, # この関数で指定した生物パラメータに置き換える最初の年
                     start_random_rec_year_name = max_vpa_year+1, # この関数で指定した再生産関係からの加入の予測値に置き換える最初の年
                     # biopar setting
                     waa_year=bio_year, waa=NULL, # 将来の年齢別体重の設定。過去の年を指定し、その平均値を使うか、直接ベクトルで指定するか。以下も同じ。
                     waa_catch_year=bio_year, waa_catch=NULL,
                     maa_year=bio_year, maa=NULL,
                     M_year=bio_year, M=c(0,0,0,Inf),
                     # faa setting
                     faa_year=2015:2017, # currentF, futureFが指定されない場合だけ有効になる。将来のFを指定の年の平均値とする
                     currentF=rep(Fvalue,4),futureF=rep(Fvalue,4), # 将来のABC.year以前のFとABC.year以降のFのベクトル 
                     # HCR setting (not work when using TMB)
                     start_ABC_year_name=max_vpa_year+2, # HCRを適用する最初の年
                     HCR_beta=1, # HCRのbeta
                     HCR_Blimit=-1, # HCRのBlimit
                     HCR_Bban=-1, # HCRのBban
                     HCR_year_lag=0, # HCRで何年遅れにするか
                     # SR setting
                     res_SR=res_sr_list$res_vpa_base0_nontune, # 将来予測に使いたい再生産関係の推定結果が入っているfit.SRの返り値
                     seed_number=1, 
                     resid_type="lognormal", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                     resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                     bias_correction=TRUE, # バイアス補正をするかどうか
                     recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                     # Other
                     Pope=res_vpa_base0_nontune$input$Pope,
                     fix_recruit=NULL,
                     fix_wcatch=NULL 
    ) 
  
  # simple
  res_future_F0.1 <- future_vpa(tmb_data=data_future_test$data,
                                optim_method="none", 
                                multi_init = 1) 
  # 平衡状態ではtarget_eq_naaと一致する（そのようなFを使っているので）
  expect_equal(mean(colSums(res_future_F0.1$naa[,"2035",])),
               target_eq_naa, tol=0.0001)
  
  
  # simple, MSY
  res_future_MSY <- future_vpa(tmb_data=data_future_test$data,
                               optim_method="R", objective ="MSY",
                               multi_init = 2, multi_lower=0.01) 
  # expect_equal(round(res_future_MSY$multi,4),2.8273)
  
  # F=0
  res_future_F0 <- data_future_test$input %>%
    list_modify(currentF=rep(0,4),
                futureF =rep(0,4)) %>%
    safe_call(make_future_data,.) %>%
    future_vpa(tmb_data=.$data, optim_method="none", multi_init = 1)
  # 平衡状態ではすべて４匹づつになる
  expect_equal(mean(res_future_F0$naa[,as.character(2025:2030),]),4, tol=0.0001)
  
  # specific weight at age, F=0.1
  res_future_F0.1_wcatch <- data_future_test$input %>%
    list_modify(res_vpa = res_vpa_base1_nontune) %>%
    safe_call(make_future_data,.) %>%
    future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)
  
  # waa.catchを別に与えた場合の総漁獲量は倍になっているはず
  expect_equal(sum(res_future_F0.1$wcaa[,,1])*2,
               sum(res_future_F0.1_wcatch$wcaa[,,1]))

  # change plus group (途中でプラスグループが変わるVPA結果でもエラーなく計算できることだけ確認（レベル１）
  res_future_pgc <- data_future_test$input %>%
    list_modify(res_vpa = res_vpa_pgc0_nontune) %>%
    safe_call(make_future_data,.) %>%
    future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)  
  
  # backward resampling 
  res_future_backward <- data_future_test$input %>%
    list_modify(resid_type="backward", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                backward_duration=5) %>%
    safe_call(make_future_data,.) %>%
      future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)

  # 残差がゼロのVPA結果なので、バックワードでも対数正規でも結果は同じ
  expect_equal(res_future_backward$naa[,as.character(2025:2030),],
               res_future_F0.1$naa[,as.character(2025:2030),])
})

# check future_vpa with dummy data ----
context("check future_vpa_function for regime shift") 

test_that("future_vpa function (with dummy vpa data) for regime shift (level 2-3?)",{

  # test-data.handler.Rで作成したVPAオブジェクトを読み込んでそれを使う  
  load("res_vpa_files.rda")
  
  # estimate SR function ----
  # VPA結果がほとんど同じになるres_vpa_base0_nontune, res_vpa_base1_nontune, res_vpa_rec0$nontueの
  # パラメータ推定値は同じになる。
  vpa_list <- tibble::lst(res_vpa_base0_nontune,
                          res_vpa_base1_nontune,
                          res_vpa_pgc0_nontune,
                          res_vpa_rec0_nontune)
  res_sr_list <- list()
  res_sr_list[[1]] <- fit.SRregime(get.SRdata(vpa_list[[1]]),
                                   SR="HS",method="HS",regime.key=c(0,1),
                                   regime.par=c("a","b"),regime.year=2005)

#  names(res_sr_list) <- names(vpa_list)
  
  # 2つのレジームでパラメータ推定値は同じになる
  tmp <- (res_sr_list[[1]]$regime_pars[1,-1]-res_sr_list[[1]]$regime_pars[2,-1]) %>%
      unlist() %>% as.numeric()
  expect_equal(round(tmp,4),c(0,0,0))
  
  # future projection with dummy data ----
  
  max_vpa_year <- max(as.numeric(colnames(res_vpa_base0_nontune$naa)))
  bio_year <- rev(as.numeric(colnames(res_vpa_base0_nontune$naa)))[1:3]
  
  target_eq_naa <- 12
  f <- function (x) 4+4*x+4*x^2+4*x^3 - target_eq_naa
  x <- uniroot(f, c(0, 1))$root
  Fvalue <- -log(x)
  data_future_test <- 
    make_future_data(res_vpa_base0_nontune, 
                     nsim = 10,
                     nyear = 20, 
                     future_initial_year_name = max_vpa_year, # 年齢別資源尾数を参照して将来予測をスタートする年
                     start_F_year_name = max_vpa_year+1, # この関数で指定したFに置き換える最初の年
                     start_biopar_year_name=max_vpa_year+1, # この関数で指定した生物パラメータに置き換える最初の年
                     start_random_rec_year_name = max_vpa_year+1, # この関数で指定した再生産関係からの加入の予測値に置き換える最初の年
                     # biopar setting
                     waa_year=bio_year, waa=NULL, # 将来の年齢別体重の設定。過去の年を指定し、その平均値を使うか、直接ベクトルで指定するか。以下も同じ。
                     waa_catch_year=bio_year, waa_catch=NULL,
                     maa_year=bio_year, maa=NULL,
                     M_year=bio_year, M=c(0,0,0,Inf),
                     # faa setting
                     faa_year=2015:2017, # currentF, futureFが指定されない場合だけ有効になる。将来のFを指定の年の平均値とする
                     currentF=rep(Fvalue,4),futureF=rep(Fvalue,4), # 将来のABC.year以前のFとABC.year以降のFのベクトル 
                     # HCR setting (not work when using TMB)
                     start_ABC_year_name=max_vpa_year+2, # HCRを適用する最初の年
                     HCR_beta=1, # HCRのbeta
                     HCR_Blimit=-1, # HCRのBlimit
                     HCR_Bban=-1, # HCRのBban
                     HCR_year_lag=0, # HCRで何年遅れにするか
                     # SR setting
                     res_SR=res_sr_list[[1]], # 将来予測に使いたい再生産関係の推定結果が入っているfit.SRの返り値
                     seed_number=1, 
                     resid_type="lognormal", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                     resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                     bias_correction=TRUE, # バイアス補正をするかどうか
                     recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                     # Other
                     Pope=res_vpa_base0_nontune$input$Pope,
                     fix_recruit=NULL,
                     fix_wcatch=NULL,regime_shift_option=list(future_regime=1)
    ) 
  
  # simple
  res_future_F0.1 <- future_vpa(tmb_data=data_future_test$data,
                                optim_method="none", 
                                multi_init = 1) 
  # 平衡状態ではtarget_eq_naaと一致する（そのようなFを使っているので）
  expect_equal(mean(colSums(res_future_F0.1$naa[,"2035",])),
               target_eq_naa, tol=0.0001)
  
  
  # simple, MSY
  res_future_MSY <- future_vpa(tmb_data=data_future_test$data,
                               optim_method="R", objective ="MSY",
                               multi_init = 2, multi_lower=0.01) 
  # expect_equal(round(res_future_MSY$multi,4),2.8273)
  
  # F=0
  res_future_F0 <- data_future_test$input %>%
    list_modify(currentF=rep(0,4),
                futureF =rep(0,4)) %>%
    safe_call(make_future_data,.) %>%
    future_vpa(tmb_data=.$data, optim_method="none", multi_init = 1)
  # 平衡状態ではすべて４匹づつになる
  expect_equal(mean(res_future_F0$naa[,as.character(2025:2030),]),4, tol=0.0001)
  
  # specific weight at age, F=0.1
  res_future_F0.1_wcatch <- data_future_test$input %>%
    list_modify(res_vpa = res_vpa_base1_nontune) %>%
    safe_call(make_future_data,.) %>%
    future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)
  
  # waa.catchを別に与えた場合の総漁獲量は倍になっているはず
  expect_equal(sum(res_future_F0.1$wcaa[,,1])*2,
               sum(res_future_F0.1_wcatch$wcaa[,,1]))

  # change plus group (途中でプラスグループが変わるVPA結果でもエラーなく計算できることだけ確認（レベル１）
  res_future_pgc <- data_future_test$input %>%
    list_modify(res_vpa = res_vpa_pgc0_nontune) %>%
    safe_call(make_future_data,.) %>%
    future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)  
  
  # backward resampling 
  res_future_backward <- data_future_test$input %>%
    list_modify(resid_type="backward", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                backward_duration=5) %>%
    safe_call(make_future_data,.) %>%
      future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)

  # 残差がゼロのVPA結果なので、バックワードでも対数正規でも結果は同じ
  expect_equal(res_future_backward$naa[,as.character(2025:2030),],
               res_future_F0.1$naa[,as.character(2025:2030),])

  #--- １年アップデートしたデータを作る
  vpa_list[[2]] <- vpa_list[[1]]
  vpa_list[[2]]$naa$"2018" <- vpa_list[[2]]$naa$"2017"
  vpa_list[[2]]$faa$"2018" <- vpa_list[[2]]$faa$"2017"
  vpa_list[[2]]$ssb$"2018" <- vpa_list[[2]]$ssb$"2017"
  vpa_list[[2]]$baa$"2018" <- vpa_list[[2]]$baa$"2017"    
  vpa_list[[2]]$input$dat$waa$"2018" <- vpa_list[[2]]$input$dat$waa$"2017"
  vpa_list[[2]]$input$dat$maa$"2018" <- vpa_list[[2]]$input$dat$waa$"2017"
  vpa_list[[2]]$input$dat$caa$"2018" <- vpa_list[[2]]$input$dat$caa$"2017"
  vpa_list[[2]]$input$dat$M$"2018"   <- vpa_list[[2]]$input$dat$M$"2017"

  max_vpa_year <- max(as.numeric(colnames(vpa_list[[2]]$naa)))
  bio_year <- rev(as.numeric(colnames(vpa_list[[2]]$naa)))[1:3]
  
  data_future_regime2 <- 
    make_future_data(vpa_list[[2]],
                     nsim = 10,
                     nyear = 20, 
                     future_initial_year_name = max_vpa_year, # 年齢別資源尾数を参照して将来予測をスタートする年
                     start_F_year_name = max_vpa_year+1, # この関数で指定したFに置き換える最初の年
                     start_biopar_year_name=max_vpa_year+1, # この関数で指定した生物パラメータに置き換える最初の年
                     start_random_rec_year_name = max_vpa_year+1, # この関数で指定した再生産関係からの加入の予測値に置き換える最初の年
                     # biopar setting
                     waa_year=bio_year, waa=NULL, # 将来の年齢別体重の設定。過去の年を指定し、その平均値を使うか、直接ベクトルで指定するか。以下も同じ。
                     waa_catch_year=bio_year, waa_catch=NULL,
                     maa_year=bio_year, maa=NULL,
                     M_year=bio_year, M=c(0,0,0,Inf),
                     # faa setting
                     faa_year=2015:2017, # currentF, futureFが指定されない場合だけ有効になる。将来のFを指定の年の平均値とする
                     currentF=rep(Fvalue,4),futureF=rep(Fvalue,4), # 将来のABC.year以前のFとABC.year以降のFのベクトル 
                     # HCR setting (not work when using TMB)
                     start_ABC_year_name=max_vpa_year+2, # HCRを適用する最初の年
                     HCR_beta=1, # HCRのbeta
                     HCR_Blimit=-1, # HCRのBlimit
                     HCR_Bban=-1, # HCRのBban
                     HCR_year_lag=0, # HCRで何年遅れにするか
                     # SR setting
                     res_SR=res_sr_list[[1]], # 将来予測に使いたい再生産関係の推定結果が入っているfit.SRの返り値
                     seed_number=1, 
                     resid_type="lognormal", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                     resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                     bias_correction=TRUE, # バイアス補正をするかどうか
                     recruit_intercept=0, # 移入や放流などで一定の加入がある場合に足す加入尾数
                     # Other
                     Pope=res_vpa_base0_nontune$input$Pope,
                     fix_recruit=NULL,
                     fix_wcatch=NULL,regime_shift_option=list(future_regime=1)
                     )

  # simple
  res_future_F0.1 <- future_vpa(tmb_data=data_future_regime2$data,
                                optim_method="none", 
                                multi_init = 1) 
  # 平衡状態ではtarget_eq_naaと一致する（そのようなFを使っているので）
  expect_equal(mean(colSums(res_future_F0.1$naa[,"2035",])),
               target_eq_naa, tol=0.0001)  

  
})

