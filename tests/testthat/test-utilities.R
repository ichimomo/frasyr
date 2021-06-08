context("Utilities")

result_vpa  <- load_data("../../inst/extdata/res_vpa_pma.rda")
result_msy  <- load_data("../../inst/extdata/res_MSY_pma_pre.rda")
data(res_vpa)
data(res_sr_HSL2)

test_that("make_kobe_ratio", {
  kobe_ratio <- make_kobe_ratio(result_vpa, result_msy)

  expect_is(kobe_ratio, "data.frame")
  expect_equal(colnames(kobe_ratio), c("year", "Fratio", "Bratio"))
  expect_equal(kobe_ratio$year, as.character(1982:2011))
  expect_is(kobe_ratio$Fratio, "numeric")
  expect_is(kobe_ratio$Bratio, "numeric")
})

test_that("pull single table from table list", {
  # pull_var_from_kobeII_table() is not tested yet
  expect_equal(1, 1)
})

test_that("convert table from kobeIItable to the required format", {
  before <- data.frame(beta = c(1, 0.5, 0),
                       y2019 = seq(12345, 11234, length.out = 3),
                       y2020 = seq(22345, 12345, length.out = 3),
                       y2021 = seq(32345, 23456, length.out = 3))

  after_raw <- format_beta_table(before, divide_by = 1)
  expect_equal(after_raw[, 1], c(1, 0.5, 0))
  expect_equal(after_raw[, "y2019"], c(12345, 11790, 11234))
  expect_equal(after_raw[, "y2020"], c(22345, 17345, 12345))
  expect_equal(after_raw[, "y2021"], c(32345, 27900, 23456))

  after_ton <- format_beta_table(before, divide_by = 1000, round = FALSE)
  expect_equal(after_ton[, 1], c(1, 0.5, 0))
  expect_equal(after_ton[, "y2019"], c(12.345, 11.790, 11.234), tolerance = 5e-4)
  expect_equal(after_ton[, "y2020"], c(22.345, 17.345, 12.345), tolerance = 5e-4)
  expect_equal(after_ton[, "y2021"], c(32.345, 27.900, 23.456), tolerance = 5e-4)

  after_ton_rounded <- format_beta_table(before, divide_by = 1000, round = TRUE)
  expect_equal(after_ton_rounded[, 1], c(1, 0.5, 0))
  expect_equal(after_ton_rounded[, "y2019"], c(12, 12, 11))
  expect_equal(after_ton_rounded[, "y2020"], c(22, 17, 12))
  expect_equal(after_ton_rounded[, "y2021"], c(32, 28, 23))
})

test_that("test for HCR function", {
    res_HCR <- HCR_default(matrix(1:10,2,5),matrix(5,2,5),matrix(1,2,5),matrix(0.8,2,5))
    expect_equal(res_HCR, matrix(c(0,0.2,0.4,0.6,rep(0.8,6)),2,5))
})

test_that("calc_future_perSPR accepts list with different length vectors", {
  future_data <- generate_dummy_future_data(result_vpa)

  perspr <- calc_future_perSPR(fout = list(waa       = future_data$data$waa_mat,
                                           maa       = future_data$data$maa_mat,
                                           M         = future_data$data$M_mat,
                                           waa.catch = future_data$data$waa_catch_mat),
                               res_vpa=result_vpa,
                               Fvector = apply_year_colum(result_vpa$faa, 2007:2011),
                               target.year = list(waa       = 2014:2018,
                                                  waa.catch = 2014:2018,
                                                  maa       = 2016:2018,
                                                  M         = 2014:2018))
  expect_is(perspr, "numeric")
})

test_that("test caa.est.mat", {

  expect_catch <- 0.5

  # set usual F => OK
  res <- caa.est.mat(c(1,1,1,1),c(1,1,1,1),c(1,1,1,1),c(0,0,0,0),catch.obs=expect_catch,Pope=TRUE)
  expect_equal(round(sum(res$caa),3),round(expect_catch,3))

  res <- caa.est.mat(c(1,1,1,1),c(1,1,1,1),c(1,1,1,1),c(0,0,0,0),catch.obs=expect_catch,Pope=FALSE)
  expect_equal(round(sum(res$caa),3),round(expect_catch,3))  

  # set very small F => OK
  res <- caa.est.mat(c(1,1,1,1),c(0.00001,0.00001,0.00001,0.00001),c(1,1,1,1),c(0,0,0,0),
                     catch.obs=expect_catch,Pope=TRUE)
  expect_equal(round(sum(res$caa),3),round(expect_catch,3))

  res <- caa.est.mat(c(1,1,1,1),c(0.00001,0.00001,0.00001,0.00001),c(1,1,1,1),c(0,0,0,0),
                     catch.obs=expect_catch,Pope=FALSE)
  expect_equal(round(sum(res$caa),3),round(expect_catch,3))  

  # naa is very small => warning
  expect_warning(caa.est.mat(c(0.01,0.01,0.01,0.01),c(1,1,1,1),c(1,1,1,1),c(0,0,0,0),
                             catch.obs=expect_catch,Pope=TRUE))
  expect_warning(caa.est.mat(c(0.01,0.01,0.01,0.01),c(1,1,1,1),c(1,1,1,1),c(0,0,0,0),
                             catch.obs=expect_catch,Pope=FALSE))  
})

test_that("calc.rel.abund",{

  age.test <- c(1:5)
  Fc.test <- rep(1,length(age.test))
  waa.test <- c(1:5)
  maa.test <- c(0,0.5,1,1,1)
  M.test <- rep(0.5,length(age.test))

  calc_rel_abund_popeT_check <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = TRUE)

  calc_rel_abund_popeF_check <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = FALSE)

  #上記計算内容をエクセルで計算したもの(../../tools/generate-testdata/check.calc.rel.abund.xlsxを数値のみのcsvに変換)を読み込み

  calc_rel_abund <- system.file("extdata", "check_calc_rel_abund.csv", package = "frasyr") %>%
    read.csv(header = T)
  #データ整形
  calc_rel_abund_popeT <- list(calc_rel_abund$rel.abundant,calc_rel_abund$ypr1.popeT,calc_rel_abund$spr)
  names(calc_rel_abund_popeT) <- c("rel.abund","ypr","spr")
  calc_rel_abund_popeF <- list(calc_rel_abund$rel.abundant,calc_rel_abund$ypr1.popeF,calc_rel_abund$spr)
  names(calc_rel_abund_popeF) <- c("rel.abund","ypr","spr")
  # 結果照合
  expect_equal(calc_rel_abund_popeT,calc_rel_abund_popeT_check)
  expect_equal(calc_rel_abund_popeF,calc_rel_abund_popeF_check)

  # 2歳分しか年齢がない場合
  age.test <- c(1:2)
  Fc.test <- rep(1,length(age.test))
  waa.test <- c(1:2)
  maa.test <- c(0,1)
  M.test <- rep(0.5,length(age.test))

  calc_rel_abund_age2 <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = TRUE)

  calc_rel_abund_age2_noplus <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = TRUE, max.age=2)

  # 3歳分あるとき
  age.test <- c(1:3)
  Fc.test <- rep(1,length(age.test))
  waa.test <- c(1:2,2)
  maa.test <- c(0,1,1)
  M.test <- rep(0.5,length(age.test))

  calc_rel_abund_age3 <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = TRUE)
  
  expect_equal(sapply(calc_rel_abund_age3,sum),
               sapply(calc_rel_abund_age2,sum))

  # 2歳分, F=0, M=0
  age.test <- c(1:3)
  Fc.test <- rep(0,length(age.test))
  waa.test <- rep(1,length(age.test))
  maa.test <- rep(1,length(age.test))
  M.test <- rep(0.001,length(age.test))
  calc_rel_abund_age2 <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = TRUE)
  calc_rel_abund_age2_noplus <- calc.rel.abund(sel = Fc.test,Fr=1,na = length(age.test),M = M.test,waa = waa.test, maa = maa.test, Pope = TRUE, max.age=length(age.test),min.age=1)
  expect_equal(sum(calc_rel_abund_age2$rel.abund), 1000.5, tol=0.001)
  expect_equal(calc_rel_abund_age2_noplus$rel.abund,c(1,1,1), tol=0.01)

})

test_that("catch_equation",{

  expect_equal(catch_equation(1,1,1,1), 1*(1-exp(-1))*exp(-0.5)*1)
  expect_equal(catch_equation(1,1,1,1,Pope = F), 1*(1-exp(-1-1))*1/(1+1)*1 )

})

test_that("solv.Feq",{

  age.test <- c(1:5)
  faa.test <- rep(1,length(age.test))
  waa.test <- c(1:5)
  naa.test <- rep(3,length(age.test))
  maa.test <- c(0,0.5,1,1,1)
  M.test <- rep(0.5,length(age.test))
  # Baranov eqをもちいてfaa,M,naaをつかってcaaを求める
  caa.test <- faa.test/(faa.test+M.test) *(1-exp(-faa.test- M.test)) *naa.test

  # tolerance=1e-4　でチェック
  expect_equal(faa.test, solv.Feq(cvec = caa.test,nvec = naa.test,mvec = M.test),tolerance=1e-4)

})

test_that("Generation.Time",{
    Generation.Time(vpares=res_vpa) %>% round(3) %>%
        expect_equal(2.919)
    Generation.Time(vpares=res_vpa, Plus=0) %>% round(3) %>%
        expect_equal(1.91)
    Generation.Time(maa=c(1,1,1),M=c(0,0,0),age=0:2,Plus=0)%>%
        expect_equal(mean(0:2))
    Generation.Time(maa=c(1,1,1),M=c(0,0,0),age=0:2)%>%
        expect_equal(mean(0:(2+19)))    

})

test_that("get.SPR", {

  target.SPR <- 30
  Fmax <- 10
  max.age <- Inf

  byear <- colnames(result_vpa$faa)[1]
  res_ref_F <- ref.F(result_vpa,waa.year=byear,maa.year=byear,M.year=byear,rps.year=2000:2011,pSPR=round(target.SPR), F.range=c(seq(from=0,to=ceiling(max(result_vpa$Fc.at.age,na.rm=T)*Fmax),length=301),max(result_vpa$Fc.at.age,na.rm=T)),plot=FALSE,max.age=max.age)

  SPR_pma_check <- get.SPR(result_vpa)
  naa <- result_vpa$naa
  years <- dimnames(naa)[[2]]

  # 出力ができているかのみ。計算結果の照合はしていない。
  expect_equal(nrow(SPR_pma_check$ysdata),length(years))
  expect_equal(colnames(SPR_pma_check$ysdata),c("perSPR","YPR","SPR","SPR0","F/Ftarget"))

})

test_that("convert_faa_perSPR", {

  res_convert_faa_perSPR <- convert_faa_perSPR(result_vpa,sel_year = 2000:2011, faa_year = 2009:2011)
  # 動いているかだけで、関数の戻り値の型だけチェック
  expect_is(res_convert_faa_perSPR,"numeric")
})

test_that("make_summary_table", {

  matrix_test <- matrix(1:9, nrow=3, ncol=3, byrow = T)

  res_mat_sum_table <- make_summary_table(matrix_test)

  means <- c()
  percent10s <-c()
  percent50s <-c()
  percent80s <-c()
  for(i in 1:nrow(matrix_test)){
    means <- c(means,mean(matrix_test[i,]))
    percent10s <- c(percent10s,quantile(matrix_test[i,],0.1))
    percent50s <- c(percent50s,quantile(matrix_test[i,],0.5))
    percent80s <- c(percent80s,quantile(matrix_test[i,],0.8))
      }
  expect_equal(res_mat_sum_table[,1], means)
  expect_equal(res_mat_sum_table[,2], as.numeric(percent10s))
  expect_equal(res_mat_sum_table[,3], as.numeric(percent50s))
  expect_equal(res_mat_sum_table[,4], as.numeric(percent80s))

})

test_that("outvpa", {

})

test_that("readvpa", {

})

test_that("to_vpa_data", {

})

test_that("get.stat", {

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
  res_future_test <- future_vpa(tmb_data=data_future_test$data,
                                                      optim_method="none",
                                                      multi_init = 1)

  # 計算結果が入っているかだけ
  expect_false(is_null(get.stat(res_future_test,use_new_output = TRUE)))
})

test_that("get.stat4",{
  # 使われていない。
  })

test_that("get.trace",{
    data("res_MSY_HSL2")
    get.trace(res_MSY_HSL2$trace) %>%
        select(ssb.mean) %>% slice(1) %>% as.numeric() %>%
        expect_equal(463754.3, tol=0.1)
})

test_that("make kobeII table", {
  test_that("beta.simulation() works", {
    kobe_data <- beta.simulation(generate_dummy_future_new_object(nsim=2)$input,
                                 beta_vector = seq(0, 1, 0.5),
                                 year.lag = 0,
                                 type = "new")

    expect_is(kobe_data, "data.frame")
    expect_setequal(colnames(kobe_data),
                    c("year", "sim", "value", "stat", "HCR_name", "beta"))

    test_that("make_kobeII_table() works", {
      kobe_table <- make_kobeII_table(kobe_data,
                        load_data("../../inst/extdata/res_vpa_pma.rda"))

      expect_is(kobe_table, "list")
#     ここの出力はフレキシブルに変わるのでテスト対象からとりあえずはずす      
#      expect_setequal(names(kobe_table),
#                      c("catch.mean", "ssb.mean", "ssb.lower10percent",
#                        "ssb.upper90percent", "prob.over.ssbtarget",
#                        "prob.over.ssblimit", "prob.over.ssbban",
#                        "prob.over.ssbmin", "prob.over.ssbmax", "catch.aav",
#                        "kobe.stat", "catch.risk", "bban.risk", "blimit.risk"))
    })
  })
})

test_that("load_folder() loads 'rda's in the given directory", {
  expect_is(load_folder("../../inst/extdata"), "list")
  test_that("each object exists", {
    expect_true(exists("res_MSY"))
    expect_true(exists("res_future_0.8HCR"))

    # これらのオブジェクトはロードされているかのように見えるが、実際はされていない
      # 原因: これらの名前が関数内にハードコードされているため
    expect_failure(expect_true(exists("res_SR")))
    expect_failure(expect_true(exists("kobeII.table")))
    expect_failure(expect_true(exists("model_selection")))
  })
})

test_that("apply_year_colum",{

  waa.year <- dimnames(result_vpa$naa)[[2]]
  trans_waa_dataframe <- as.data.frame(result_vpa$input$dat$waa)

  appl_year_col_check <- apply_year_colum(result_vpa$input$dat$waa,waa.year)

  test_matrix_rows <- c()
  for(i in 1:nrow(trans_waa_dataframe)){
    test_matrix_rows <- c(test_matrix_rows, mean(as.matrix(trans_waa_dataframe[i,])) )
  }
    expect_equal(as.numeric(appl_year_col_check),test_matrix_rows)

})

test_that("convert_df",{
  # 一旦スキップ
})

test_that("convert_2d_future",{
  # 一旦スキップ
})

test_that("derive_biopar",{
    
    a1 <- derive_biopar(res_vpa, derive_year=2000)
    a2 <- derive_biopar(res_vpa, derive_year=2000:2003)
    expect_equal(a1[,1:3],a2[,1:3])
    (a1$faa + a2$faa) %>% round(2) %>% as.numeric() %>%
        expect_equal(c(1.12,3.47,3.82,3.82))

    rm(list=ls())
    a1 <- derive_biopar(res_future_0.8HCR, derive_year=2030)
    a2 <- derive_biopar(res_future_0.8HCR, derive_year=2031:2032)
    expect_equal(a1[,c(1,3,4)],a2[,c(1,3,4)])
    (a1$faa + a2$faa) %>% round(2) %>% as.numeric() %>%
        expect_equal(c(0.22, 0.53, 0.60, 0.60))
    expect_equal(a1$waa, a1$waa.catch)

})


test_that("calc_forward",{
  naa2 <- calc_forward(naa=res_vpa$naa,faa=res_vpa$faa,M=res_vpa$input$dat$M,t=1,plus_age=4,plus_group=TRUE)[,2]
  naa2 %>% expect_equal(res_vpa$naa[,2])
})

