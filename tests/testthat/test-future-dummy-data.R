context("check future_vpa_function2 with dummy data") # ダミーデータでの将来予測 ----

# test-data.handler.Rで作成したVPAオブジェクトを読み込んでそれを使う
load("res_vpa_files.rda")
# estimate SR function ----
# VPA結果がほとんど同じになるres_vpa_base0_nontune, res_vpa_base1_nontune, res_vpa_rec0$nontueの
# パラメータ推定値は同じになる。
vpa_list <- tibble::lst(res_vpa_base0_nontune,
                        res_vpa_base1_nontune,
                        res_vpa_pgc0_nontune,
                        res_vpa_rec0_nontune)

test_that("future_vpa function (with dummy vpa data) (level 2-3?)",{

  res_sr_list <- list()
  for(i in 1:length(vpa_list)){
      biopar <- derive_biopar(vpa_list[[i]],derive_year=2017)
      biopar$M[nrow(biopar)] <- 100
      x <- vpa_list[[i]]
      res_sr <- res_sr_list[[i]] <- get.SRdata(x) %>% fit.SR(AR=0, SR="HS",bio_par=biopar)
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
  save(data_future_test, file="data_future_test.rda")

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
  expect_equal(res_future_MSY$multi,1.345,tol=0.001)

  # F=0
  res_future_F0 <- redo_future(data_future_test,
                               list(currentF=rep(0,4),futureF =rep(0,4)),
                               optim_method="none", multi_init = 1)  
  # 平衡状態ではすべて４匹づつになる
  expect_equal(mean(res_future_F0$naa[,as.character(2025:2030),]),4, tol=0.0001)

  # specific weight at age, F=0.1
  res_future_F0.1_wcatch <- redo_future(data_future_test,
                                        list(res_vpa=res_vpa_base1_nontune),
                                        optim_method="none", multi_init = 1)
  # waa.catchを別に与えた場合の総漁獲量は倍になっているはず
  expect_equal(sum(res_future_F0.1$wcaa[,,1])*2,
               sum(res_future_F0.1_wcatch$wcaa[,,1]))

  #plot_futures(res_vpa_base1_nontune, list(res_future_F0.1_wcatch))

  data_future_F0.1_wcatch <- redo_future(data_future_test,
                                         list(res_vpa=res_vpa_base1_nontune),
                                         only_data=TRUE)
  aa <- est_MSYRP(data_future_F0.1_wcatch)
  expect_equal(aa$summary$cB,aa$summary$B*2)
  expect_equal(aa$summary$U,aa$summary$Catch/aa$summary$cB)  

  # change plus group (途中でプラスグループが変わるVPA結果でもエラーなく計算できることだけ確認（レベル１）
  res_future_pgc <- redo_future(data_future_test,
                                list(res_vpa=res_vpa_pgc0_nontune),
                                optim_method="none", multi_init = 1)  

  # backward resampling
  res_future_backward <- redo_future(data_future_test,
                                     list(resid_type="backward", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                                          resample_year_range=as.numeric(colnames(data_future_test$input$res_vpa$naa)), # リサンプリングの場合、残差をリサンプリングする年の範囲
                                          backward_duration=5),
                                optim_method="none", multi_init = 1)     

  # 残差がゼロのVPA結果なので、バックワードでも対数正規でも結果は同じ
  expect_equal(res_future_backward$naa[,as.character(2025:2030),],
               res_future_F0.1$naa[,as.character(2025:2030),])


  # sd>0, ARありの場合のテスト
  res_sr_sd02ar00 <- res_sr_list$res_vpa_base0_nontune
  res_sr_sd02ar00$pars$sd <- 0.1

  data_future_sd02ar00 <- redo_future(data_future_test,
                                      list(nsim=5000, res_SR=res_sr_sd02ar00),
                                      only_data=TRUE)
  res_future_sd02ar00 <- future_vpa(data_future_sd02ar00$data,
                                    optim_method="none", multi_init=0)

  expect_equal(sd(log(res_future_sd02ar00$naa[1,43,])), 0.1, tol=0.01)
  expect_equal(mean(  res_future_sd02ar00$naa[1,43,]) ,   4, tol=0.001)

  # sd=0, AR=0.5
  res_sr_sd00ar10 <- res_sr_list$res_vpa_base0_nontune
  res_sr_sd00ar10$pars$sd <- 0
  res_sr_sd00ar10$pars$rho <- 0.5
  res_vpa_base0_nontune2 <- res_vpa_base0_nontune
  res_vpa_base0_nontune2$naa$"2017"[1] <- 6

  # 通常の将来予測
  res_future_sd00ar10 <- redo_future(data_future_test,
             list(nsim=3, res_SR=res_sr_sd00ar10, res_vpa=res_vpa_base0_nontune2),
             optim_method="none", multi_init=0)
  # 2018年の値
  expect_equal(res_future_sd00ar10$naa[1,"2018",1],
               exp(log(res_future_sd00ar10$naa[1,"2017",1]/4) * 0.5) * 4, tol=0.001)
  # 最終年の値
  expect_equal(mean(res_future_sd00ar10$naa[1,43,]), 4, tol=0.01)

 
  # sd>0
  res_sr_sd01ar10 <- res_sr_list$res_vpa_base0_nontune
  res_sr_sd01ar10$pars$sd <- 0.1
  res_sr_sd01ar10$pars$rho <- 0.5  

  res_future_sd01ar10 <- redo_future(data_future_test,
             list(nsim=1000, res_SR=res_sr_sd01ar10, res_vpa=res_vpa_base0_nontune2),
             optim_method="none", multi_init=0)
  # 2018年の平均
  expect_equal(mean(res_future_sd01ar10$naa[1,"2018",]),
               exp(log(res_future_sd01ar10$naa[1,"2017",1]/4) * 0.5) * 4, tol=0.01)
  # 最終年のSDと平均
  expect_equal(mean(res_future_sd01ar10$naa[1,43,]), 4, tol=0.01)
  expect_equal(sd(log(res_future_sd01ar10$naa[1,43,])), sqrt(1/(1-(0.5^2)) * 0.1 ^2), tol=0.01)

  # fix_recruitオプション+ARあり(ichimomo/frasyr_tool#256を再現)
  res_future <- redo_future(data_future_test,
            list(nsim=1000, res_SR=res_sr_sd01ar10, res_vpa=res_vpa_base0_nontune2,
                 fix_recruit=list(year=2018,rec=6)),
             optim_method="none", multi_init=0)
  # 2019年の平均  
  expect_equal(mean(res_future$naa[1,"2019",]),
               exp(log(res_future$naa[1,"2018",1]/4) * 0.5) * 4, tol=0.01)
  # 最終年のSDと平均
  expect_equal(mean(res_future$naa[1,43,]), 4, tol=0.01)
  expect_equal(sd(log(res_future$naa[1,43,])), sqrt(1/(1-(0.5^2)) * 0.1 ^2), tol=0.01)

  # fix_recruitオプション（加入複数）
  res_future <- redo_future(data_future_test,
            list(nsim=10, res_SR=res_sr_sd01ar10, res_vpa=res_vpa_base0_nontune2,
                 fix_recruit=list(year=c(2018,2019),rec=list(rep(6,10),rep(5,10)))),
            optim_method="none", multi_init=0)
  expect_equal(as.numeric(rowMeans(res_future$naa[1,c("2018","2019"),])),c(6,5))

  res_future <- redo_future(data_future_test,
            list(nsim=10, res_SR=res_sr_sd01ar10, res_vpa=res_vpa_base0_nontune2,
                 fix_recruit=list(year=c(2018),rec=list(rep(5,10)))),
            optim_method="none", multi_init=0)
  expect_equal(as.numeric(mean(res_future$naa[1,c("2018"),])),c(5))  
  
  ## vpa function without plus group (with dummy vpa data)
  vpa_no_plus_input <- vpa_list[[1]]$input %>% list_modify(plus.group=FALSE)
  vpa_no_plus_input$dat$caa[4,] <- vpa_no_plus_input$dat$caa[4,]/2
  vpa_no_plus <- do.call(vpa,vpa_no_plus_input)
  vpa_no_plus$naa %>% apply(1,mean) %>% as.numeric() %>%
    expect_equal(c(4,3,2,1))
  
  expect_equal(detect_plus_group(vpa_no_plus),FALSE)
  expect_equal(detect_plus_group(vpa_list[[1]]),TRUE)  

  # VPAでプラスグループなしオプションを使う場合は将来予測にそのまま引き継がれるはずだが、do.callを使う場合にはそうでない
  # 明示的に指定する
  res_future_no_plus <- redo_future(data_future_test,
                                    list(res_vpa=vpa_no_plus,M=c(0,0,0,0),plus_group=FALSE),
                                    optim_method="none", multi_init=0)
  expect_equal(mean(res_future_no_plus$naa[,as.character(2025:2030),]),4, tol=0.0001)
  rev(rowMeans(res_future_no_plus$SR_mat[,,"ssb"]))[1] %>% as.numeric() %>%
      expect_equal(res_sr_list[[1]]$steepness$SB0)

  BRP <- ref.F(vpa_no_plus,M=c(0,0,0,0.0001),Fcurrent=rep(Fvalue,4))
  res_future_no_plus_F90 <- redo_future(data_future_test,
                                    list(res_vpa=vpa_no_plus,M=c(0,0,0,0),plus_group=FALSE),
                                    optim_method="none",
                                    multi_init=BRP$summary$FpSPR.90.SPR[3])
  
  x <- rev(rowMeans(res_future_no_plus_F90$SR_mat[,,"ssb"]))[1]/rev(rowMeans(res_future_no_plus$SR_mat[,,"ssb"]))[1]  %>%  unlist() %>% as.numeric() %>% as.numeric()
  expect_equal(round(as.numeric(x),3),0.9)
  

  
})

context("check future_vpa_function for regime shift") # ダミーデータ・レジームシフト将来予測 ----

test_that("future_vpa function (with dummy vpa data) for regime shift (level 2-3?)",{

  res_sr_list <- list()
  res_sr_list[[1]] <- fit.SRregime(get.SRdata(vpa_list[[1]]),
                                   SR="HS",method="L2",regime.key=c(0,1),
                                   regime.par=c("a","b"),regime.year=2005)

  # sdが異なるケースもテストしないといけない
  res_sr_list[[2]] <- fit.SRregime(get.SRdata(vpa_list[[1]]),
                                   SR="HS",method="L2",regime.key=c(0,1),
                                   regime.par=c("a","b","sd"),regime.year=2005)
#  res_sr_list[[2]]$pars$sd[2] <- 0.3 # 本当は両方ゼロだがテストのために0.3を入れる
  res_sr_list[[2]]$regime_pars$sd[2] <- 0.3 # 本当は両方ゼロだがテストのために0.3を入れる

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
                     resample_year_range=as.numeric(colnames(res_vpa_base0_nontune$naa)), # リサンプリングの場合、残差をリサンプリングする年の範囲
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

  # simple (different SD)
  data_future_test2 <- data_future_test
  data_future_test2 <- safe_call(make_future_data,
                                 list_modify(data_future_test2$input,res_SR=res_sr_list[[2]]))
  res_future_F2 <- future_vpa(tmb_data=data_future_test2$data,
                              optim_method="none",
                              multi_init = 1)
  # 修正前のプログラムだとここで交互にrand_residが0,乱数,0,乱数となるようになってしまっている
  # boxplot(res_future_F2$SR_mat[as.character(2023:2030),,"rand_resid"],type="b")

  data_future_test3 <- data_future_test2
  data_future_test3 <- safe_call(make_future_data,
                                 list_modify(data_future_test2$input,
                                             regime_shift_option=list(future_regime=0)))
  res_future_F2 <- future_vpa(tmb_data=data_future_test3$data,
                              optim_method="none",
                              multi_init = 1)


  # simple, MSY
  res_future_MSY <- future_vpa(tmb_data=data_future_test$data,
                               optim_method="R", objective ="MSY",
                               multi_init = 2, multi_lower=0.01)
#  expect_equal(mean(res_future_MSY$naa[,as.character(2025:2030),]),4, tol=0.0001)

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
                resample_year_range=as.numeric(colnames(data_future_test$input$res_vpa$naa)), # リサンプリングの場合、残差をリサンプリングする年の範囲
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
                     resample_year_range=as.numeric(colnames(vpa_list[[2]]$naa)), # リサンプリングの場合、残差をリサンプリングする年の範囲
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

test_that("future_vpa function (with dummy vpa data) for regime shift & shepherd ",{

  res_sr_list <- list()
  data(res_vpa_example)
  res_vpa <- res_vpa_example
  res_sr_list[[1]] <- fit.SRregime(get.SRdata(res_vpa),
                                   SR="Shepherd",method="L2",regime.key=c(0,1),
                                   regime.par=c("a","b","sd"),regime.year=2005,gamma=0.5)
  res_sr_list[[1]]$regime_pars$sd[2] <- 0
  
  res_sr_list[[2]] <- fit.SRregime(get.SRdata(res_vpa),
                                   SR="Cushing",method="L2",regime.key=c(0,1),
                                   regime.par=c("a","b","sd"),regime.year=2005)
  res_sr_list[[2]]$regime_pars$sd[2] <- 0

 
  res_sr_list[[3]] <- fit.SRregime(get.SRdata(res_vpa),
                                   SR="BH",method="L2",regime.key=c(0,1),
                                   regime.par=c("a","b","sd"),regime.year=2005)

  res_sr_list[[4]] <- fit.SRregime(get.SRdata(res_vpa),
                                   SR="Shepherd",method="L2",regime.key=c(0,1),
                                   regime.par=c("a","b","sd"),regime.year=2005,gamma=1,
                                   p0=as.numeric(unlist(t(res_sr_list[[3]]$regime_pars[c("a","b")]))))  

  expect_equal(mean(unlist(res_sr_list[[3]]$regime_pars[c("a","b","sd")]/res_sr_list[[4]]$regime_pars[c("a","b","sd")])),
               1,tol=0.0001)
  res_sr_list[[3]]$regime_pars$sd[2] <- 0  
  res_sr_list[[4]]$regime_pars$sd[2] <- 0  
  
  # future projection with dummy data ----
  max_vpa_year <- max(as.numeric(colnames(res_vpa$naa)))
  bio_year <- rev(as.numeric(colnames(res_vpa$naa)))[1:3]
  Fvalue <- 0.1
  data_future_test <-
    make_future_data(res_vpa,
                     nsim = 10,
                     nyear = 20,
                     future_initial_year_name = max_vpa_year, 
                     start_F_year_name = max_vpa_year+1, 
                     start_biopar_year_name=max_vpa_year+1, 
                     start_random_rec_year_name = max_vpa_year+1, 
                     # biopar setting
                     waa_year=bio_year, waa=NULL, 
                     waa_catch_year=bio_year, waa_catch=NULL,
                     maa_year=bio_year, maa=NULL,
                     M_year=bio_year, M=NULL,
                     # faa setting
                     faa_year=2015:2017, 
                     currentF=rep(Fvalue,4),futureF=rep(Fvalue,4), 
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
                     Pope=res_vpa$input$Pope,
                     fix_recruit=NULL,
                     fix_wcatch=NULL,regime_shift_option=list(future_regime=1)
                     )


  # simple
  res_future_F0.1 <- future_vpa(tmb_data=data_future_test$data,
                                optim_method="none",
                                multi_init = 1)
  pars <- res_sr_list[[1]]$regime_pars
  # SSB -> SRF_SH -> recruitment
  pred_rec <- res_future_F0.1$summary %>% dplyr::filter(year>max_vpa_year) %>% select(SSB) %>%
      SRF_SH(a=pars$a[2], b=pars$b[2], gamma=res_sr_list[[1]]$input$gamma) %>% unlist() %>%
      as.numeric()
  # calculated recruitment
  cacl_rec <- res_future_F0.1$summary %>% dplyr::filter(year>max_vpa_year) %>% select(recruit) %>% unlist() %>% as.numeric()
  expect_equal(cacl_rec, pred_rec)  

  # Cushing
  data_future_test2 <- data_future_test
  data_future_test2 <- safe_call(make_future_data,
                                 list_modify(data_future_test2$input,res_SR=res_sr_list[[2]]))
  res_future_F2 <- future_vpa(tmb_data=data_future_test2$data,
                              optim_method="none",
                              multi_init = 1)
  pars <- res_sr_list[[2]]$regime_pars
  # SSB -> SRF_SH -> recruitment
  pred_rec <- res_future_F2$summary %>% dplyr::filter(year>max_vpa_year) %>% select(SSB) %>%
      SRF_CU(a=pars$a[2], b=pars$b[2]) %>% unlist() %>%
      as.numeric()
  # calculated recruitment
  calc_rec <- res_future_F2$summary %>% dplyr::filter(year>max_vpa_year) %>% select(recruit) %>% unlist() %>% as.numeric()
  expect_equal(calc_rec, pred_rec)  
  
  # simple, MSY
  res_future_MSY <- future_vpa(tmb_data=data_future_test$data,
                               optim_method="R", objective ="MSY",
                               multi_init = 2, multi_lower=0.01)

  # backward resampling
  res_future_backward <- data_future_test$input %>%
    list_modify(resid_type="backward", # 加入の誤差分布（"lognormal": 対数正規分布、"resample": 残差リサンプリング）
                resample_year_range=1988:2017, # リサンプリングの場合、残差をリサンプリングする年の範囲
                backward_duration=5) %>%
    safe_call(make_future_data,.) %>%
      future_vpa(tmb_data=.$data, optim_method="none", multi_init=1)

  # 結果の数値チェックはまだない

  MSY_BH <- redo_future(data_future_test, list(res_SR=res_sr_list[[3]]), only_data=TRUE) %>%
      est_MSYRP()
  MSY_SH <- redo_future(data_future_test, list(res_SR=res_sr_list[[4]]), only_data=TRUE) %>%
      est_MSYRP()  
  expect_equal(MSY_BH$summary$SSB/MSY_SH$summary$SSB, c(1,1,1,1),tol=0.001)


})
