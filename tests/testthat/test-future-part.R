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


# utility function definition ----
to_vpa_data <- function(x, label_name){
  vpa_data <- x %>% dplyr::filter(label==label_name) %>% as.data.frame()
  if(label_name!="abund"){  
    rownames(vpa_data) <- str_c("X",vpa_data$year)
    vpa_data <- vpa_data %>% select(-year, -label) 
    vpa_data$value <- vpa_data$fishery <- NULL    
    vpa_data <- as.data.frame(t(vpa_data))
    rownames(vpa_data) <- str_sub(rownames(vpa_data),5,5)
  }
  else{
    vpa_data <- vpa_data %>% 
      select(-starts_with("age_")) %>%
      pivot_wider(names_from=year, values_from=value) 
    abund_name <- vpa_data$fishery
    vpa_data <- as.data.frame(vpa_data)
    vpa_data$label <- vpa_data$fishery <- NULL
    rownames(vpa_data) <- abund_name
    colnames(vpa_data) <- str_c("X",colnames(vpa_data))
  }
  vpa_data  
}

# check future_vpa with dummy data ----
context("check future_vpa_function2 with dummy data") 

test_that("future_vpa function (with dummy vpa data) (level 2-3?)",{
  # read various data ----
  # data with caa=maa=waa=1, M=0
  data_base <- read_csv(system.file("extdata","all_dummy_data_base.csv",package="frasyr")) 
  # data with caa=maa=waa=1, M=0 but plus group have changed 
  data_pgc <- read_csv(system.file("extdata","all_dummy_data_plus_group_change.csv",package="frasyr")) 
  # data with caa=maa=waa=1, M=0 but first recruit age is 1
  data_rec <- read_csv(system.file("extdata","all_dummy_data_rec.csv",package="frasyr")) 

  # create various vpa data ----
  vpadat_base0 <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
                               waa=to_vpa_data(data_base, label_name="waa"),
                               maa=to_vpa_data(data_base, label_name="maa"),
                               M  = 0,
                               index = to_vpa_data(data_base, label_name="abund"),
                               maa.tune = NULL,
                               waa.catch = NULL,
                               catch.prop = NULL)
  # waa.catch is given 
  vpadat_base1 <- data.handler(caa=to_vpa_data(data_base, label_name="caa"),
                               waa=to_vpa_data(data_base, label_name="waa"),
                               maa=to_vpa_data(data_base, label_name="maa"),
                               M  = 0,
                               index = to_vpa_data(data_base, label_name="abund"),
                               maa.tune = NULL,
                               waa.catch = to_vpa_data(data_base, label_name="waa")*2,
                               catch.prop = NULL)
  
  vpadat_pgc0 <- data.handler(caa=to_vpa_data(data_pgc, label_name="caa"),
                              waa=to_vpa_data(data_pgc, label_name="waa"),
                              maa=to_vpa_data(data_pgc, label_name="maa"),
                              M  = 0,
                              index = to_vpa_data(data_pgc, label_name="abund"),
                              maa.tune = NULL,
                              waa.catch = NULL,
                              catch.prop = NULL)
  
  vpadat_rec0 <- data.handler(caa=to_vpa_data(data_rec, label_name="caa"),
                              waa=to_vpa_data(data_rec, label_name="waa"),
                              maa=to_vpa_data(data_rec, label_name="maa"),
                              M  = 0,
                              index = to_vpa_data(data_rec, label_name="abund"),
                              maa.tune = NULL,
                              waa.catch = NULL,
                              catch.prop = NULL)
  
  # vpa (no tuning) ----
  # この結果が4,3,2,2になるのはなんとなくそんな感じする             
  res_vpa_base0_nontune <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_base0_nontune$naa)), 
               c(4,3,2,2))
  
  res_vpa_base1_nontune <- vpa(vpadat_base1, tf.year=2015:2016, last.catch.zero = FALSE, 
                               Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_base1_nontune$naa)), 
               c(4,3,2,2))
  
  res_vpa_pgc0_nontune <- vpa(vpadat_pgc0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(unlist(res_vpa_pgc0_nontune$naa["2017"])), 
               c(3,2,2,NA), tol=0.0001)
  
  res_vpa_rec0_nontune <- vpa(vpadat_rec0, tf.year=2015:2016, last.catch.zero = FALSE, 
                              Pope = TRUE, p.init = 0.5) 
  expect_equal(as.numeric(rowMeans(res_vpa_rec0_nontune$naa)), 
               c(4,3,2,2))
  
  # catch計算用のwaaを２倍にしているbase1データでは漁獲量が倍になる
  expect_equal(res_vpa_base0_nontune$wcaa*2,
               res_vpa_base1_nontune$wcaa)
  
  # vpa (tuning) ----
  # この記述は正しいか不明。計算しただけ。
  res_vpa_base0_tune <- vpa(vpadat_base0, tf.year=2015:2016, last.catch.zero = FALSE, 
                            Pope = TRUE, p.init = 0.5, tune=TRUE, sel.update=TRUE)
  
  # estimate SR function ----
  library(magrittr)
  
  # VPA結果がほとんど同じになるres_vpa_base0_nontune, res_vpa_base1_nontune, res_vpa_rec0$nontueの
  # パラメータ推定値は同じになる。
  res_sr_list <- purrr::map(tibble::lst(res_vpa_base0_nontune,
                                        res_vpa_base1_nontune,
                                        res_vpa_rec0_nontune), 
                            function(x){
                              const_ssr <- mean(colSums(x$ssb))
                              const_R   <- mean(unlist(x$naa[1,]))
                              
                              res_sr <- get.SRdata(x) %>% #%T>% plot_SRdata() %>%
                                fit.SR(AR=0, SR="HS") %>% # %T>% SRplot_gg() 
                              
                              expect_equal(res_sr$pars$sd, 0, tol=0.000001)
                              expect_equal(res_sr$pars$b, const_ssr, tol=0.000001)
                              expect_equal(const_R/const_ssr, res_sr$pars$a)
                              return(res_sr)
                            })
  
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
  
  # backward resampling 
  res_future_backward <- data_future_test$input %>%
    list_modify(resid_type="backward", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                resample_year_range=0, # リサンプリングの場合、残差をリサンプリングする年の範囲
                backward_duration=5) %>%
    safe_call(make_future_data,.) %>%
    future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)
  expect_equal(res_future_backward$naa[,as.character(2025:2030),],
               res_future_F0.1$naa[,as.character(2025:2030),])
})

