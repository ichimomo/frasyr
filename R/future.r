# basic dynamics ----

#'
#' future_vpaにインプットとして入れる将来予測の空のarrayを生成する関数
#'
#' @param res_vpa vpaの結果 (vpa関数の返り値)
#' @param nsim シミュレーションの繰り返し回数
#' @param nyear 将来予測の実施年数
#' @param plus_age プラスグループとして計算する行（年齢ではないことに注意）。デフォルト値（NULL）ならfuture_initial_year_name年にNA以外の数値が入っている一番大きい行をプラスグループの行とする。plus_ageという名前だが、plus_group=FALSEの場合には、たんに最大年齢になる。年齢の行数よりもプラスグループが小さいような場合（対馬マイワシなど）に対応。
#' @param plus_group プラスグループを考慮するかどうか。与えない場合、res_vpa$input$plus.groupの設定を引き継ぐ。ただし、inputを使ってdo.callする場合などにはうまく調整できないことがあるので、明示的に与えたほうがよいかも。
#' @param future_initial_year_name 将来予測の「初期値となる」年齢別資源尾数を参照する年。この年の年齢別資源尾数を使って翌年の個体群動態が将来予測で決定される
#' @param start_F_year_name 将来予測でF全体にmultiplierを乗じる場合、multiplierを乗じる最初の年
#' @param start_biopar_year_name 生物パラメータを将来の生物パラメータとして設定された値に置き換える年の最初の年
#' @param start_random_rec_year_name 将来の加入を再生産関係から予測する最初の年
#' @param waa_year 将来の年齢別体重を過去の平均値とする場合、過去のパラメータを平均する期間, maa_year, M_yearも同様
#' @param waa 将来の年齢別体重を直接与える場合, maa, M_yearも同様
#' @param waa_fun FALSE: 使わない、TRUE: log(weight)~log(number)の回帰式から将来のweightを予測する,
#' @param waa_fun_name カスタマイズされたwaa_funを使う場合、その関数のオブジェクトの名前
#' @param waa_catch_fun FALSE: 使わない、TRUE: log(weight)~log(number)の回帰式から将来のweightを予測する
#' @param waa_catch_fun_name カスタマイズされたwaa_funを使う場合、その関数のオブジェクトの名前
#' @param maa_fun maturity ~ number の回帰式から将来のmaturityを予測する(暫定的、太平洋マダラでのみ利用)
#' @param start_waafun_year_name 上記の設定がスタートする最初の年。それ以外の年は上で設定されたパラメータが使われる
#' @param faa_year 将来のFを過去の平均値とする場合、平均をとる年を指定する。下のcurrentF, futureFが指定されている場合にはこの設定は無視される
#' @param currentF start_ABC_yar_name以前に使うFのベクトル（いわゆるcurrent F）
#' @param futureF  start_ABC_yar_name以降に使うFのベクトル（いわゆるFmsy）
#' @param start_ABC_year_name HCRを有効にする年
#' @param HCR_beta HCRのbeta
#' @param HCR_Blimit HCRのBlimit
#' @param HCR_Bban HCRのBban
#' @param HCR_year_lag HCRするときにいつのタイミングのssbを参照するか.0の場合、ABC計算年のSSBを参照する。正の値1を入れると1年前のssbを参照する
#' @param HCR_beta_year betaを年によって変える場合。tibble(year=2020:2024, beta=c(1.3,1.2,1.1,1,0.9))　のようにtibble形式で与える。HCR_betaで設定されたbetaは上書きされる。
#' @param HCR_Blimit_year Blimitを年によって変える場合。tibble(year=2020:2024, Blimit=c(1.3,1.2,1.1,1,0.9))　のようにtibble形式で与える。HCR_Blimitで設定されたBlimitは上書きされる。
#' @param HCR_Bban_year Bbanを年によって変える場合。tibble(year=2020:2024, Bban=c(1.3,1.2,1.1,1,0.9))　のようにtibble形式で与える。HCR_Bbanで設定されたBbanは上書きされる。
#' @param HCR_TAC_reserve_rate TACの取り残し率
#' @param HCR_TAC_carry_rate TACの何％まで持ち越せるか
#' @param HCR_TAC_reserve_amount その年の総漁獲可能量に対する獲り残し量
#' @param HCR_TAC_carry_amount 当初TACのうち何トンまで持ち越せるか
#' @param HCR_TAC_upper_CV 漁獲量が前年の漁獲量のHCR_TAC_upper_CV倍と比較し、それよりも変化が大きい場合には前年の漁獲量xHCR_TAC_upper_CVを上限とする。単一の値か、tibble形式 tibble(year=2020:2024, TAC_upper_CV=rep(0.1,5)) で与える
#' @param HCR_TAC_lower_CV 漁獲量が前年の漁獲量のHCR_TAC_lower_CV倍と比較し、それよりも変化が大きい場合には前年の漁獲量xHCR_TAC_lower_CVを下限とする。単一の値か、tibble形式 tibble(year=2020:2024, TAC_lower_CV=rep(0.1,5)) で与える
#' @param Pope 漁獲方程式にPopeの近似式を使うかどうか。与えない場合には、VPAのオプションが引き継がれる
#' @param HCR_function_name デフォルトは"HCR_default" ここを変更(関数名を文字列で与える。関数は別に定義しておく)すると自作のHCRが適用される。その場合、データのほうで定義されているbeta, Blimit, Bbanなど、同じ名前のものは有効
#' @param fix_recruit 将来予測において再生産関係を無視して加入量を一定値で与える場合、その加入の値。list(year=2020, rec=1000)のように与える。
#' @param fix_wcatch 将来予測において漁獲量をあらかじめ決める場合
#' @param res_SR 再生産関係の推定関数 (fit.SR　of fit.SRregime) の返り値
#' @param seed_number 乱数のシードの数
#' @param resid_type 加入量の残差の発生方法。"lognormal":対数正規分布, "resample": リサンプリング, "backward": backward resampling
#' @param bias_correction 将来予測でバイアス補正をするかどうか
#' @param resample_year_range "resampling", "backward"で有効。年の範囲を入れると、対象とした年の範囲で計算される残差を用いる。
#' @param backward_duration "backward"の場合、何年で１ブロックとするか。"backward"で有効。デフォルトは5 。
#' @param recruit_intercept 将来の加入の切片。将来の加入は R=f(ssb) + intercept となる。
#' @param setting_release より詳細な放流シナリオを設定する場合. 将来の放流数は list(number=tibble(value=)) or list(number=tibble(year=))) or list(number=tibble(year=xxx, value=xxx))) として与える。yearのみ単独で与えるのは過去年を与え、その期間の平均放流数を使う。yearとvalueの両方を与える場合には、将来年が想定されており、その年の放流尾数をvalueとして用いる。例えば tibble(year=c(2021:2100), value=c(100,200,rep(300,length(21:2100)-2)) と与えると2021,2022年は100, 200匹の総放流尾数、2023年以降は300匹となる。年は、実際に将来予測を実施するよりも十分長い年数を与えること（長すぎても問題はないので）。将来の添加効率も list(rate=tibble(value=)) or list(rate=tibble(year=))) として与える。この場合のyearは過去年。valueを与えると、与えたvalueのベクトルからランダムサンプリングされる。結果的にsetting_release = list(number=tibble(year=1990:2000), rate=tibble(value=c(0.4, 0.5, 0.3))) みたいな、リストの中に２つのtibbleが入ったちょっとややこしい構造になります。この将来予測へのデータにはres_vpaとres_SR$input$SRdataの両方に放流に関するデータが格納されています。過去の放流データをどちらのソースからとってくるかについても、data_source="SR" or "VPA"で選択してください。setting_releaseで指定されない限り、デフォルトはVPA結果からとってきます。
#' @param model_average_option model averagingをする場合のオプション. SR_matのlistとweightをlist形式で入れる(list(SR_list=list(res_SR1,res_SR2),weight=c(0.5,0.5)))
#' @param regime_shift_option res_SRにfit.SRregimeの返り値を入れた場合に指定する。将来予測で再生産関係のどのフェーズがおこるかを指定する。list(future_regime=将来のregimeの仮定。keyで指定された番号を入れる)
#' @param special_setting list形式で与えるmake_future_dataの返り値のdataと同じ名前の要素について、最後にデータをここで示されたarrayのシミュレーション1回めの値で上書きする。arrayのデータに対してのみ有効。
#'
#' @return 以下の要素からなるリスト
#' \describe{
#' \item{\code{input}}{使用した引数のリスト。\code{do.call(make_future_data, input)}で計算を再現できる}
#' \item{\code{data}}{future_vpa関数に渡すデータのセット}
#' \item{\code{data$naa_mat, data$caa_mat, data$waa_mat, data$waa_catch_mat, data$maa_mat, data$M_mat, data$faa_mat}}{年齢×年数（VPA期間年＋将来予測年）×シミュレーション回数の3次元データ。順に、年齢別年別シミュレーション別の資源尾数、漁獲尾数、体重、漁獲量計算用の体重、成熟率、自然死亡係数、漁獲死亡係数。このうち、資源尾数はVPA期間年までのみデータが入っていて、\code{future_vpa}を実行することによってここに推定値が入る。\code{faa_mat}についても、future_vpa実行時にHCRなどの適用の設定によって適宜書き換えられる。}
#' \item{\code{data$SR_mat}}{
#'
#' }
#' }
#'
#' @export
#' @encoding UTF-8

make_future_data <- function(res_vpa,
                             nsim = 1000, # number of simulation
                             nyear = 50, # number of future year
                             plus_age  = NULL, # if null, equal to row number as plus group
                             plus_group = NULL,
                             future_initial_year_name = 2017,
                             start_F_year_name = 2018,
                             start_biopar_year_name=2018,
                             start_random_rec_year_name = 2018,
                             # biopar setting
                             waa_year, waa=NULL,
                             waa_catch_year, waa_catch=NULL,
                             waa_fun = FALSE,
                             start_waafun_year_name = start_biopar_year_name,
                             waa_fun_name = NA, 
                             waa_catch_fun = FALSE,                             
                             start_waacatchfun_year_name = start_biopar_year_name,
                             waa_catch_fun_name = NA,                              
                             maa_year, maa=NULL,
                             maa_fun = FALSE,
                             start_maafun_year_name = start_biopar_year_name,
                             M_year, M=NULL,
                             # faa setting
                             faa_year=NULL,
                             currentF=NULL,futureF=NULL,
                             # HCR setting (not work when using TMB)
                             start_ABC_year_name=2019,
                             HCR_beta=1,
                             HCR_Blimit=-1,
                             HCR_Bban=-1,
                             HCR_year_lag=0,
                             HCR_beta_year=NULL,
                             HCR_Blimit_year=NULL,
                             HCR_Bban_year=NULL,
                             HCR_TAC_reserve_rate=NA,
                             HCR_TAC_reserve_amount=NA,  #
                             HCR_TAC_carry_rate=NA,
                             HCR_TAC_carry_amount=NA,    #
                             HCR_TAC_upper_CV=NA,
                             HCR_TAC_lower_CV=NA,
                             HCR_function_name="HCR_default",
                             # Other
                             Pope=res_vpa$input$Pope,
                             fix_recruit=NULL, # list(year=2020, rec=1000)
                             fix_wcatch=NULL, # list(year=2020, wcatch=2000)
                             max_F=exp(10),
                             max_exploitation_rate=0.99,
                             # SR setting
                             res_SR=NULL,
                             seed_number=1,
                             scale_ssb=1,
                             scale_R=1,                             
                             resid_type="lognormal", # or resample or backward
                             bias_correction=TRUE,
                             resample_year_range=NA, # only when "resample" or backward
                             backward_duration=5, # only when backward
                             recruit_intercept=0, # number of additional recruitment (immigration or enhancement)
                             setting_release=NULL, # より詳細な放流シナリオを設定する場合
                             model_average_option=NULL,
                             regime_shift_option =NULL,
                             silent=FALSE,
                             # special
                             special_setting=NULL
)
{

  argname <- ls()
  input <- lapply(argname,function(x) eval(parse(text=x)))
  names(input) <- argname

  if(!is.na(HCR_TAC_reserve_rate  )) assertthat::assert_that(min(HCR_TAC_reserve_rate  ) >= 0)
  if(!is.na(HCR_TAC_carry_rate    )) assertthat::assert_that(min(HCR_TAC_carry_rate    ) >= 0)
  if(!is.na(HCR_TAC_reserve_amount)) assertthat::assert_that(min(HCR_TAC_reserve_amount) >= 0)
  if(!is.na(HCR_TAC_carry_amount  )) assertthat::assert_that(min(HCR_TAC_carry_amount  ) >= 0)

  assertthat::assert_that(is.logical(waa_fun),
                          is.logical(waa_catch_fun),
                          is.logical(maa_fun),
                          is.logical(bias_correction))

  if(!is.na(HCR_TAC_reserve_rate) && !is.na(HCR_TAC_reserve_amount)) stop("HCR_TAC_reserve_rateとHCR_TAC_reserve_amountが同時に指定されています（同時には指定できません）")
  if(!is.na(HCR_TAC_carry_rate) && !is.na(HCR_TAC_carry_amount))     stop("HCR_TAC_carry_rateとHCR_TAC_carry_amountが同時に指定されています（同時には指定できません）")

  # define age and year
  nage <- nrow(res_vpa$naa)
  age_name    <- as.numeric(rownames(res_vpa$naa))
  recruit_age <- min(as.numeric(rownames(res_vpa$naa)))

  vpa_nyear           <- ncol(res_vpa$naa)
  future_initial_year <- which(colnames(res_vpa$naa)==future_initial_year_name)
  if(length(future_initial_year)  ==0) stop("future_initial_year_name is invalid.")
  total_nyear         <- future_initial_year + nyear
  allyear_name        <- min(as.numeric(colnames(res_vpa$naa)))+c(0:(total_nyear-1))
  allyear_label       <- c(rep("VPA",future_initial_year),rep("future",nyear))
  start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
  if(length(start_random_rec_year)==0) stop("start_random_rec_year_name is invalid.")

  tmpdata <- tibble(allyear_name, allyear_label) %>%
    group_by(allyear_label) %>%
    summarize(start=min(allyear_name),end=max(allyear_name)) %>%
    arrange(start)
  if(is.null(plus_age)) plus_age <- max(which(!is.na(res_vpa$naa[,future_initial_year])))
  if(is.null(plus_group)) plus_group <- res_vpa$input$plus.group

  if(silent==FALSE){
      print(tmpdata)
      cat("plus.group =",plus_group,"\n")
      cat("Pope =",Pope,"\n")
  }

  # define empty array
  waa_mat <- waa_catch_mat <- M_mat <- maa_mat <- naa_mat <- faa_mat <- caa_mat <- waa_catch_mat <-
    array(0, dim=c(nage, total_nyear, nsim),
          dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim))
  class(waa_mat) <- class(M_mat) <- class(maa_mat) <- class(naa_mat) <- class(faa_mat) <- class(caa_mat) <- class(waa_catch_mat) <- "myarray"
  SR_mat <- array(0, dim=c(total_nyear, nsim, 15),
                  dimnames=list(year=allyear_name, nsim=1:nsim,
                                par=c("a","b","rho", #1-3
                                      "SR_type", # 4
                                      "rand_resid", # 5
                                      "deviance", #6
                                      "recruit","ssb",
                                      "intercept","sd",#9-10
                                      "bias_factor", #11
                                      "gamma", #12 gamma for shepherd
                                      "biomass","cbiomass", # add biomass statistics
                                      "blank5")))

  #HCR_mat <- array(0, dim=c(total_nyear, nsim, 7),
  HCR_mat <- array(0, dim=c(total_nyear, nsim, 11),
                   dimnames=list(year=allyear_name, nsim=1:nsim,
                                 par=c("beta","Blimit","Bban","year_lag", #1-4
                                       "expect_wcatch",# 5 漁獲量。ここにあらかじめ値を入れているとこの漁獲量どおりに漁獲する
                                       # 以下、取り残し用の設定
                                       "TAC_reserve_rate", # 6 全漁獲可能量の何割まで獲り残すか      #
                                       "TAC_carry_rate", # 7 TACの何割まで翌年に持ち越しを許容するか  #
                                       "TAC_reserve_amount", # 8 全漁獲可能量のうち何トンまで獲り残すか  #
                                       "TAC_carry_amount", # 9 何トンまで持ち越しを許容するか             #
                                       "TAC_upper_CV", # 10 漁獲量の変動の上側の上限
                                       "TAC_lower_CV" # 11 漁獲量の変動の下側の下限
                                       )))
  class(SR_mat)  <- "myarray"
  class(HCR_mat) <- "myarray"

  HCR_mat[,,"Blimit"] <- HCR_mat[,,"Bban"] <- -1
#  HCR_mat[,,"beta"] <- HCR_mat[,,"beta_gamma"] <- 1
  HCR_mat[,,"beta"] <- 1
  HCR_mat[,,"TAC_reserve_rate"] <- HCR_mat[,,"TAC_carry_rate"] <- NA
  HCR_mat[,,"TAC_reserve_amount"] <- HCR_mat[,,"TAC_carry_amount"] <- NA
  HCR_mat[,,"TAC_upper_CV"] <- NA
  HCR_mat[,,"TAC_lower_CV"] <- NA

  # fill vpa data
  waa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa)
  maa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$maa)
  M_mat  [,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$M)
  faa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$faa)
  naa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$naa)
  caa_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$caa)

  waa_mat <- make_array(waa_mat, waa, waa_year, start_biopar_year_name)
  maa_mat <- make_array(maa_mat, maa, maa_year, start_biopar_year_name)
  M_mat   <- make_array(M_mat  , M  , M_year  , start_biopar_year_name)

  if(is.null(res_vpa$input$dat$waa.catch)){
    waa_catch_mat <- waa_mat
  }
  else{
    waa_catch_mat[,1:vpa_nyear,] <- as.matrix(res_vpa$input$dat$waa.catch)
    waa_catch_mat <- make_array(waa_catch_mat, waa_catch, waa_catch_year, start_biopar_year_name)
  }

  # set SR parameter
  SR_mat <- set_SR_mat(res_vpa=res_vpa,
                       res_SR=res_SR,
                       SR_mat=SR_mat, seed_number=seed_number,
                       start_random_rec_year_name=start_random_rec_year_name, resid_type=resid_type,
                       resample_year_range=resample_year_range,
                       bias_correction=bias_correction,
                       recruit_intercept=recruit_intercept,
                       recruit_age=recruit_age,
                       setting_release=setting_release,
                       backward_duration=backward_duration,
                       scale_ssb=scale_ssb,
                       scale_R  =scale_R,
                       model_average_option=model_average_option,
                       regime_shift_option=regime_shift_option,
                       fix_recruit=fix_recruit)

  naa_mat[1,,] <- SR_mat[,,"recruit"]

  # set F & HCR parameter
  start_F_year <- which(allyear_name==start_F_year_name)
  start_ABC_year <- which(allyear_name==start_ABC_year_name)
  faa_mat <- make_array(faa_mat, currentF, faa_year, start_F_year_name)
  faa_mat <- make_array(faa_mat, futureF,  faa_year, start_ABC_year_name)

  HCR_mat[start_ABC_year:total_nyear,,"beta"    ] <- HCR_beta
  HCR_mat[start_ABC_year:total_nyear,,"Blimit"  ] <- HCR_Blimit
  HCR_mat[start_ABC_year:total_nyear,,"Bban"    ] <- HCR_Bban
  HCR_mat[start_ABC_year:total_nyear,,"year_lag"] <- HCR_year_lag
  HCR_mat[start_ABC_year:total_nyear,,"TAC_reserve_rate"] <- HCR_TAC_reserve_rate
  HCR_mat[start_ABC_year:total_nyear,,"TAC_carry_rate"  ]   <- HCR_TAC_carry_rate
  HCR_mat[start_ABC_year:total_nyear,,"TAC_reserve_amount"] <- HCR_TAC_reserve_amount #
  HCR_mat[start_ABC_year:total_nyear,,"TAC_carry_amount"]   <- HCR_TAC_carry_amount     #

  assign_HCR_ <- function(HCR_mat, HCR_year, target){
    if(!is.null(HCR_year)){
      assert_that(all(c("year",target)%in%names(HCR_year)))
      tmp <- which(dimnames(HCR_mat)[[1]] %in% as.character(HCR_year$year) )
      tmp2 <- which(as.character(HCR_year$year)%in% dimnames(HCR_mat)[[1]] )
      if(length(tmp)>0){
          tmp3 <- HCR_year[target][1] %>% unlist() %>% as.numeric()
          HCR_mat[tmp,,target] <- tmp3[tmp2]
        }
    }
    return(HCR_mat)
  }

  HCR_mat <- assign_HCR_(HCR_mat, HCR_beta_year,   target="beta")
  HCR_mat <- assign_HCR_(HCR_mat, HCR_Blimit_year, target="Blimit")
  HCR_mat <- assign_HCR_(HCR_mat, HCR_Bban_year,   target="Bban")

  # set upper and lower limit of TAC
  if(!is.list(HCR_TAC_upper_CV)){
    HCR_mat[start_ABC_year:total_nyear,,"TAC_upper_CV"] <- HCR_TAC_upper_CV
  }else{
    HCR_mat <- assign_HCR_(HCR_mat, HCR_TAC_upper_CV,   target="TAC_upper_CV")
  }
  if(!is.list(HCR_TAC_lower_CV)){
    HCR_mat[start_ABC_year:total_nyear,,"TAC_lower_CV"] <- HCR_TAC_lower_CV
  }else{
    HCR_mat <- assign_HCR_(HCR_mat, HCR_TAC_lower_CV,   target="TAC_lower_CV")
  }

  # when fix wcatch
  # VPA期間中のwcatchをfixするかどうか？？
  #    if(isTRUE(Pope)){
  #        tmp <- naa_mat*(1-exp(-faa_mat))*exp(-M_mat/2) * waa_catch_mat
  #    }
  #    else{
  #        tmp <- naa_mat*(1-exp(-faa_mat-M))*faa_mat/(faa_mat+M) * waa_catch_mat
  #    }
  #    HCR_mat[,,"wcatch"] <- apply(tmp,c(2,3),sum)
  HCR_mat[as.character(fix_wcatch$year), ,"expect_wcatch"] <- fix_wcatch$wcatch

  waa_mat[is.na(waa_mat)] <- 0
  waa_catch_mat[is.na(waa_catch_mat)] <- 0
  maa_mat[is.na(maa_mat)] <- 0
  M_mat  [is.na(M_mat)  ] <- 0
  faa_mat[is.na(faa_mat)] <- 0
  naa_mat[is.na(naa_mat)] <- 0
  caa_mat[is.na(caa_mat)] <- 0


  # set data and parameter for TMB
  tmb_data <- list(naa_mat=naa_mat,
                   caa_mat=caa_mat,
                   SR_mat=SR_mat,
                   waa_mat = waa_mat,
                   waa_catch_mat = waa_catch_mat,
                   maa_mat = maa_mat,
                   M_mat = M_mat,
                   faa_mat = faa_mat,
                   Pope = as.numeric(Pope),
                   total_nyear = total_nyear,
                   future_initial_year = future_initial_year,
                   start_ABC_year=start_ABC_year,
                   start_random_rec_year=start_random_rec_year,
                   nsim = nsim,
                   nage = nage,
                   plus_age = plus_age,
                   plus_group = plus_group,
                   recruit_age = recruit_age,
                   max_exploitation_rate=max_exploitation_rate,
                   max_F=max_F,
                   scale_ssb=scale_ssb,
                   scale_R  =scale_R,                   
                   HCR_mat = HCR_mat,
                   obj_stat = 0, # 0: mean, 1:geomean
                   objective = 0, # 0: MSY, 1: PGY, 2: percentB0 or Bempirical
                   obj_value = -1,
                   HCR_function_name=HCR_function_name
  )

  if(isTRUE(waa_fun)){
    waa_rand_mat <- array(0,dim=c(nage,total_nyear,nsim),
                          dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim))
    waa_par_mat <- array(0,dim=c(nage,nsim,3),
                         dimnames=list(age=age_name, nsim=1:nsim, pars=c("sd", "b0", "b1")))
    class(waa_rand_mat) <- class(waa_par_mat) <- "myarray"
    waa_fun_year <- which(allyear_name %in% start_waafun_year_name:max(allyear_name))    
    if(is.na(waa_fun_name)){
      for(a in 1:nage){
        for(i in 1:nsim){
          log.w <- as.numeric(log(waa_mat[a,,i]))
          log.n <- as.numeric(log(naa_mat[a,,i]))
          observed <- log.n>-Inf
          log.w <- log.w[observed]
          log.n <- log.n[observed]
          tmp <- lm(log.w~log.n)
          waa_par_mat[a,i,c("b0","b1")] <- as.numeric(tmp$coef[1:2])
          waa_par_mat[a,i,c("sd")] <- sqrt(mean(tmp$residual^2))
          waa_rand_mat[a,observed,i] <- tmp$residual
          waa_rand_mat[a,waa_fun_year,i] <- rnorm(length(waa_fun_year),-0.5*waa_par_mat[a,i,c("sd")]^2,waa_par_mat[a,i,c("sd")])
        }}
    }
    else{
      dimnames(waa_par_mat)[[3]] <- c("waa_fun_name", "x1","x2")
      waa_par_mat[,,1] <- waa_fun_name
      waa_rand_mat[,waa_fun_year,] <- 1
    }
    tmb_data$waa_rand_mat <- waa_rand_mat
    tmb_data$waa_par_mat <- waa_par_mat
    tmb_data$waa_mat[,waa_fun_year,] <- 0    
  }


  if(isTRUE(waa_catch_fun)){
    waa_catch_rand_mat <- array(0,dim=c(nage,total_nyear,nsim),
                                dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim))
    waa_catch_par_mat <- array(0,dim=c(nage,nsim,3),
                               dimnames=list(age=age_name, nsim=1:nsim, pars=c("sd", "b0", "b1")))
    class(waa_catch_rand_mat) <- class(waa_catch_par_mat) <- "myarray"
    waa_catch_fun_year <- which(allyear_name %in% start_waacatchfun_year_name:max(allyear_name))
    if(is.na(waa_catch_fun_name)){
      for(a in 1:nage){
        for(i in 1:nsim){
          log.w <- as.numeric(log(waa_catch_mat[a,,i]))
          log.n <- as.numeric(log(naa_mat[a,,i]))
          observed <- log.n>-Inf
          log.w <- log.w[observed]
          log.n <- log.n[observed]
          tmp <- lm(log.w~log.n)
          waa_catch_par_mat[a,i,c("b0","b1")] <- as.numeric(tmp$coef[1:2])
          waa_catch_par_mat[a,i,c("sd")] <- sqrt(mean(tmp$residual^2))
          waa_catch_rand_mat[a,observed,i] <- tmp$residual
          
          # 特別なwaa_funを使っていないけどwaa_catch_funとcaa_funの両方を使っている場合
          # 両者の乱数をSDで調整してから一致させる
          if(waa_fun==TRUE && is.na(waa_fun_name)){
            waa_catch_rand_mat[a,waa_catch_fun_year,i]  <- waa_rand_mat[a,waa_catch_fun_year,i] *
              waa_catch_par_mat[a,i,c("sd")]/waa_par_mat[a,i,c("sd")]
            warning("暫定的にwaa_funとwaa_catch_funで同じ年・年齢の場合には同じ正規化された残差を使うことにします。このオプションには今後より妥当な設定についての検討が必要です（漁獲量計算用の体重と資源量計算用の体重が同じ場合には大丈夫です）\n")            
          }
          else{
            waa_catch_rand_mat[a,waa_catch_fun_year,i] <- rnorm(length(waa_catch_fun_year),-0.5*waa_catch_par_mat[a,i,c("sd")]^2,waa_catch_par_mat[a,i,c("sd")])
          }
        }}
    }
    else{ # when using specific waa function
      dimnames(waa_catch_par_mat)[[3]] <- c("waa_catch_fun_name", "x1","x2")
      waa_catch_par_mat[,,1] <- waa_catch_fun_name
      waa_catch_rand_mat[,waa_catch_fun_year,] <- 1
    }
    tmb_data$waa_catch_rand_mat <- waa_catch_rand_mat
    tmb_data$waa_catch_par_mat <- waa_catch_par_mat
    tmb_data$waa_catch_mat[,waa_catch_fun_year,] <- 0    
  }

  if(isTRUE(maa_fun)){
    maa_rand_mat <- array(0,dim=c(nage,total_nyear,nsim),
                          dimnames=list(age=age_name, year=allyear_name, nsim=1:nsim))
    maa_par_mat <- array(0,dim=c(nage,nsim,5),
                         dimnames=list(age=age_name, nsim=1:nsim, pars=c("sd", "b0", "b1","min","max")))
    class(maa_rand_mat) <- class(maa_par_mat) <- "myarray"
    maa_fun_year <- which(allyear_name %in% start_maafun_year_name:max(allyear_name))

    for(a in 1:nage){
      for(i in 1:nsim){
        data_tmp <- data.frame(naa=naa_mat[a,,i], maa=maa_mat[a,,i])
        observed <- naa_mat[a,,i]>0
        tmp <- lm(maa~naa, data=data_tmp[observed,])
        maa_par_mat[a,i,c("b0","b1")] <- as.numeric(tmp$coef[1:2])
        maa_par_mat[a,i,c("sd")] <- sqrt(mean(tmp$residual^2))
        maa_par_mat[a,i,c("min")] <- min(maa_mat[a,,i])
        maa_par_mat[a,i,c("max")] <- max(maa_mat[a,,i])
        maa_rand_mat[a,observed,i] <- tmp$residual
        maa_rand_mat[a,maa_fun_year,i] <- rnorm(length(maa_fun_year),-0.5*maa_par_mat[a,i,c("sd")]^2,maa_par_mat[a,i,c("sd")])
      }}
    tmb_data$maa_rand_mat <- maa_rand_mat
    tmb_data$maa_par_mat <- maa_par_mat
    tmb_data$maa_mat[,maa_fun_year,] <- 0
  }

  if(!is.null(special_setting)){
    set_name <- names(special_setting)
    for(i in seq_len(length(set_name))){
      tmb_data[[which(set_name[[i]]==tmb_data)[[1]]]][] <- special_setting[[i]][]
    }}

  return(tibble::lst(data=tmb_data,input=input))
}

#' 将来予測の実施関数
#'
#' @param tmb_data make_future_dataの返り値。将来の生物パラメータや再生産関係のシナリオを年齢×年×シミュレーション回数で指定した様々なarrayが含まれる。
#' @param optim_method "none"の場合、通常の将来予測を実施. "R": 以下のobj関係の設定とあわせてMSYなどを探索する. "tmb"もあるが、限定した設定でしか使えない
#' @param SPR_target 目標とする\%SPR。NULL以外の値の場合、過去〜将来のそれぞれの年・シミュレーションが、目標とするF\%SPRに対して何倍にあたるか(F/Ftarget)を計算して、HCR_realizedの"Fratio"に入れる。HCRが生きている年については"beta_gamma"と一致するはず。
#' @param max_F 漁獲量一定方策を実施する際のF at ageの最大値の上限（将来的にはmake_future_data関数に入れたい)
#' @param max_exploitation_rate 漁獲量一定方策を実施する際のMを考慮した上での漁獲率の上限（将来的にはmake_future_data関数に入れたい)
#' @param do_MSE 簡易MSEを実施するか
#' @param MSE_input_data 簡易MSEを実施する場合、ABC計算するための将来予測を実施するための設定ファイル
#' @param MSE_nsim 簡易MSEを実施する場合、ABC計算するための将来予測の繰り返し回数。ここを1にすると、決定論的な将来予測の漁獲量が用いられる。
#' @param MSE_sd 簡易MSEをする場合の加入変動の大きさ。ここをゼロにすると決定論的な将来予測の値を将来の漁獲量として用いる。その場合MSE_nsimは自動的に２に設定される。単純なモデルの場合、ここがゼロでも多分問題ない。モデル平均を使っている場合にはちゃんとした簡易MSEをすること。
#' @param objective MSY:MSYの推定、PGY:PGYの値をobj_valueに入れる、percentB0:B0パーセント、何％にするかはobj_valueで指定, SSB:obj_valueで指定した特定の親魚資源量に一致するようにする
#' @param obj_stat 目的関数を計算するときに利用する計算方法（"mean"だと平均、"median"だと中央値、"geomean"だと幾何平均）
#' @param obj_value 目的とする値
#'
#' @export
#' @encoding UTF-8

future_vpa <- function(tmb_data,
                       optim_method="none", # or "R" or "none" or "tmb"
                       multi_init = 1,
                       multi_lower = 0.0001,
                       multi_upper = 10,
                       objective ="MSY", # or PGY, percentB0, Bempirical
                       obj_value = 0,
                       obj_stat  ="mean",
                       do_MSE=NULL,
                       MSE_input_data=NULL,
                       MSE_nsim=NULL,
                       MSE_sd=NULL,
                       MSE_catch_exact_TAC=FALSE,
                       compile=FALSE,
                       output_format="new",
                       attach_input=TRUE,
                       SPRtarget=NULL,
                       calc_SPR_year_name=NULL){

  argname <- ls()
  input <- lapply(argname,function(x) eval(parse(text=x)))
  names(input) <- argname

  # conversion of estimation option
  tmb_data$objective <-
    dplyr::case_when(objective=="MSY"   ~ 0,
                     objective=="PGY"   ~ 1,
                     objective=="SSB"   ~ 2)
  tmb_data$obj_stat <-
    dplyr::case_when(obj_stat=="mean"    ~ 0,
                     obj_stat=="geomean" ~ 1,
                     obj_stat=="median"  ~ 2)
  tmb_data$obj_value <- obj_value

  if(optim_method=="tmb"){
    if(!is.null(tmb_data$waa_par_mat) || !is.null(tmb_data$maa_par_mat) || !is.null(tmb_data$waa_catch_par_mat)){
      cat("Warning: waa_fun and maa_fun option cannot be used in TMB. The optimization is conducted by R..\n")
      optim_method <- "R"
  }}

  if(optim_method=="tmb"){

    #        if(!is.null(rec_new) | !is.null(wcatch_fix)){
    #            stop("rec_new or wcatch_fix option cannot be used in \"tmb\" option\n")
    #        }
    tmb_data$HCR_function_name <- NULL

    # comple & load cpp file
    use_rvpa_tmb(TmbFile = "est_MSY_tmb",
                 CppDir = system.file("executable",package="frasyr"),
                 RunDir = getwd(), overwrite=compile)

    objAD <- TMB::MakeADFun(tmb_data, list(x=log(multi_init)), DLL="est_MSY_tmb")
    msy_optim <- nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                        lower=list(x=log(multi_lower)),
                        upper=list(x=log(multi_upper)))#,contol=nlminb_control)

    multi <- as.numeric(exp(msy_optim$par))
    msy <- exp(-msy_optim$objective)
    ad_report <- objAD$report()

    res_future <- naming_adreport(tmb_data, ad_report)
    res_future$multi <- multi
  }

  #--- R
  if(optim_method=="R" | optim_method=="none"){

    tmb_data$do_MSE <- do_MSE
    tmb_data$MSE_input_data <- MSE_input_data
    tmb_data$MSE_nsim <- MSE_nsim
    tmb_data$MSE_sd <- MSE_sd
    tmb_data$MSE_catch_exact_TAC <- MSE_catch_exact_TAC

    R_obj_fun <- function(x, tmb_data, what_return="obj"){
      tmb_data$x <- x
      tmb_data$what_return <- what_return
      obj <- safe_call(future_vpa_R, tmb_data)
      return(obj)
    }

    if(optim_method=="R"){
      # nlminbを使うとrangeによって壁に当たることがある。optimizeを使ったほうがよさそうなので、optimizeに変更
      # 一方、tmbはnlminbでないとうまくいかないらしい
      #            msy_optim <- nlminb(start=0, objective=R_obj_fun, tmb_data=tmb_data,
      #                                lower=list(x=log(multi_lower)), upper=list(x=log(multi_upper)))
      #            tmb_data$x <- msy_optim$par
      msy_optim <- optimize(R_obj_fun,lower=log(multi_lower),
                            upper=log(multi_upper), tmb_data=tmb_data)
      tmb_data$x <- msy_optim$minimum

    }
    else{
      tmb_data$x <- log(multi_init)
    }
    tmb_data$what_return <- "stat"
    res_future <- safe_call(future_vpa_R, tmb_data)
  }

  res_future$input <- input

  if(output_format!="new"){
    res_future <- format_to_old_future(res_future)
  }
  else{
    class(res_future) <- "future_new"
  }
  if(!isTRUE(attach_input)) res_future$input <- NULL

  # modify results
  # remove HCR control parameter before start_ABC_year
  res_future$HCR_mat[1:(tmb_data$start_ABC_year-1),,] <- NA
  # add Fratio if needed
  if(!is.null(SPRtarget)){
    if(is.null(calc_SPR_year_name)){
       calc_SPR_year <- 1:dim(res_future$faa)[[2]]
    }
    else{
       calc_SPR_year <- which(dimnames(res_future$faa)[[2]] %in% calc_SPR_year_name)
    }

    for(i in calc_SPR_year){
      for(j in 1:dim(res_future$faa)[[3]]){
        if(j>2 &&
           all(res_future$faa[,i,j]==res_future$faa[,i,j-1]) &&
           all(res_future$waa[,i,j]==res_future$waa[,i,j-1]) &&
           all(res_future$waa_catch_mat[,i,j:(j-1)]==res_future$waa_catch_mat[,i,j:(j)]) &&
           all(res_future$maa[,i,j]==res_future$maa[,i,j-1]) &&
           all(res_future$M  [,i,j]==res_future$M  [,i,j-1])){
          res_future$HCR_realized[i,j,"Fratio"] <- res_future$HCR_realized[i,j-1,"Fratio"]
        }
        else{
          tmp <- res_future$naa[,i,j]>0
          res_future$HCR_realized[i,j,"Fratio"] <-
            calc_Fratio(faa=res_future$faa[tmp,i,j],
                        waa=res_future$waa[tmp,i,j],
                        maa=res_future$maa[tmp,i,j],
                        M  =res_future$input$tmb_data$M_mat[tmp,i,j],
                        waa.catch=res_future$waa_catch_mat[tmp,i,j],
                        SPRtarget=SPRtarget,
                        plus_group=tmb_data$plus_group)
        }
      }}}

  res_future$summary <- derive_future_summary(res_future)
  return(res_future)

  # 足りないもの
  # 推定結果の簡易的グラフ(MSY_est_plot)
}

#' @export

future_vpa_R <- function(naa_mat,
                         caa_mat,
                         SR_mat,
                         waa_mat,
                         waa_catch_mat,
                         M_mat,
                         maa_mat,
                         faa_mat,
                         Pope,
                         total_nyear,
                         future_initial_year,
                         start_ABC_year, # the year for estimating multiplier F & conduct HCR
                         start_random_rec_year,
                         nsim,
                         nage,
                         plus_age,
                         plus_group,
                         recruit_age,
                         obj_stat,
                         objective,
                         obj_value,
                         x,
                         scale_ssb=1,
                         scale_R=1,
                         what_return="obj",
                         HCR_mat,
                         HCR_function_name,
                         max_F=exp(10),
                         max_exploitation_rate=0.99,
                         do_MSE=NULL,
                         MSE_input_data=NULL,
                         MSE_nsim = NULL,
                         MSE_sd = NULL,
                         MSE_catch_exact_TAC=FALSE,
                         waa_par_mat  = NULL, # option for waa_fun
                         waa_rand_mat = NULL,
                         waa_catch_par_mat  = NULL, # option for waa_catch_fun
                         waa_catch_rand_mat = NULL,                         
                         maa_par_mat  = NULL, # option for maa_fun
                         maa_rand_mat = NULL
){

  options(deparse.max.lines=10)

  # setting for density dependent biological parameter
  is_waa_fun <- !is.null(waa_par_mat)
  if(is_waa_fun) when_waa_fun <- apply(waa_mat[,,1],2,sum)==0
  is_waa_catch_fun <- !is.null(waa_catch_par_mat)
  if(is_waa_catch_fun) when_waa_catch_fun <- apply(waa_catch_mat[,,1],2,sum)==0
  is_maa_fun <- !is.null(maa_par_mat)
  if(is_maa_fun) when_maa_fun <- apply(maa_mat[,,1],2,sum)==0

  # setting for specific function for waa_fun
  if(is_waa_fun && dimnames(waa_par_mat)[[3]][1]=="waa_fun_name"){
    update_waa_mat <- get(waa_par_mat[1,1,"waa_fun_name"])
  }
  if(is_waa_catch_fun && dimnames(waa_catch_par_mat)[[3]][1]=="waa_catch_fun_name"){
    update_waa_catch_mat <- get(waa_catch_par_mat[1,1,"waa_catch_fun_name"])
  }    

  HCR_function <- get(HCR_function_name)
  allyear_name <- as.numeric(dimnames(SR_mat)[[1]])

  argname <- ls()
  tmb_data <- lapply(argname,function(x) eval(parse(text=x)))
  names(tmb_data) <- argname

  HCR_realized_name <- c("wcatch", "beta_gamma", "Fratio","reserved_catch","original_ABC","original_ABC_plus")
  HCR_realized <- array(0,dim=c(dim(HCR_mat)[[1]],dim(HCR_mat)[[2]],length(HCR_realized_name)),
                        dimnames=list(dimnames(HCR_mat)[[1]],
                                      dimnames(HCR_mat)[[2]],
                                      HCR_realized_name))
  class(HCR_realized) <- "myarray"

  if(isTRUE(do_MSE)){
    MSE_seed <- MSE_input_data$input$seed_number + 1
    if(!is.null(MSE_sd) && MSE_sd==0){
      MSE_nsim <- 2
      if(MSE_input_data$input$resid_type=="resampling"|MSE_input_data$input$resid_type=="backward")
        warning("MSE_sd=0 option cannot be used in resampling option")
    }
    if(!is.null(MSE_nsim)) MSE_input_data$input$nsim <- MSE_nsim
    if( is.null(MSE_nsim)) MSE_nsim <- MSE_input_data$input$nsim
    SR_MSE <- SR_mat
    SR_MSE[,,"recruit"] <- SR_MSE[,,"ssb"] <- 0
    dimnames(SR_MSE)$par[12] <- "real_true_catch"
    dimnames(SR_MSE)$par[13] <- "pseudo_true_catch"

    # max_F, max_exploitation_rateはそのままMSEに引き継ぐとしたけどやめる
    #
    # max_F_MSE <- max_F; max_exploitation_rate_MSE <- max_exploitation_rate
  }

  F_mat <- N_mat <-  naa_mat
  F_mat[] <- N_mat[] <-  0
  class(F_mat) <- class(N_mat) <- "myarray"

  N_mat <- naa_mat
  spawner_mat <- apply(N_mat * waa_mat * maa_mat, c(2,3) , sum)

  F_mat[,1:(start_ABC_year-1),] <- faa_mat[,1:(start_ABC_year-1),]
  F_mat[,start_ABC_year:total_nyear,] <- faa_mat[,start_ABC_year:total_nyear,] * exp(x)

  ## この時点でもwaa_funを入れる必要がある

  for(t in future_initial_year:total_nyear){
    # 再生産関係から加入を推定するため、０歳以外のSSB計算用の体重と成熟率を更新
    if(is_waa_fun)
      waa_mat[,t,]       <- update_waa_mat(t=t,waa=waa_mat,rand=waa_rand_mat,naa=N_mat,
                                     pars_b0=waa_par_mat[,,"b0"],pars_b1=waa_par_mat[,,"b1"])
#    if(is_waa_catch_fun)
#      waa_catch_mat[,t,] <- update_waa_catch_mat(t=t,waa=waa_catch_mat,rand=waa_catch_rand_mat,naa=N_mat,
#                                                 pars_b0=waa_catch_par_mat[,,"b0"],pars_b1=waa_catch_par_mat[,,"b1"])    
    if(is_maa_fun) maa_mat[,t,] <- update_maa_mat(maa=maa_mat[,t,],rand=maa_rand_mat[,t,],naa=N_mat[,t,],
                                                  pars_b0=maa_par_mat[,,"b0"],pars_b1=maa_par_mat[,,"b1"],
                                                  min_value=maa_par_mat[,,"min"],max_value=maa_par_mat[,,"max"])
    spawner_mat[t,] <- colSums(N_mat[,t,,drop=F] * waa_mat[,t,,drop=F] * maa_mat[,t,,drop=F])
    spawner_mat[t,spawner_mat[t,]<0.0001] <-  0.0001

    if(t>=start_random_rec_year){
      spawn_t <- t-recruit_age
      # 加入を再生産関係からの予測値とする場合
      if(all(N_mat[1,t,]==0)){
        N_mat[1,t,] <- purrr::pmap_dbl(tibble(x=SR_mat[t,,"SR_type"],
                                              ssb=spawner_mat[spawn_t,]*scale_ssb,
                                              a=SR_mat[t,,"a"],b=SR_mat[t,,"b"],gamma=SR_mat[t,,"gamma"]),
                                       function(x,ssb,a,b,gamma){
                                         fun <- list(SRF_HS,SRF_BH,SRF_RI,SRF_SH,SRF_CU,SRF_BHS)[[x]];
                                         fun(ssb,a,b,gamma)
                                       })*scale_R
        N_mat[1,t,] <- N_mat[1,t,]*exp(SR_mat[t,,"deviance"]) + SR_mat[t,,"intercept"]
      }else{
        # fix_recruitですでに加入尾数が入っていて、自己相関ありの場合
        # make_future_dataの段階では対応するSSBがいくつかわからないので、SSBが計算された段階で
        # 残差を計算しなおす必要がある -> new_deviance
        if(!all(N_mat[1,t,]==0) && all(SR_mat[t-1,,"rho"]>0)){
          rec_predict <- purrr::pmap_dbl(tibble(x=SR_mat[t,,"SR_type"],
                                              ssb=spawner_mat[spawn_t,]*scale_ssb,
                                              a=SR_mat[t,,"a"],b=SR_mat[t,,"b"],gamma=SR_mat[t,,"gamma"]),
                                         function(x,ssb,a,b,gamma){
                                           fun <- list(SRF_HS,SRF_BH,SRF_RI,SRF_SH,SRF_CU)[[x]];
                                           fun(ssb,a,b,gamma)
                                         })*scale_R
          new_deviance <- log(N_mat[1,t,]) - log(rec_predict)
          SR_mat[t,,"deviance"] <- new_deviance

          # t年の残差をもとに、将来予測年のdevianceを置き換える
          if(t<total_nyear){ # replace recruit deviance
            for(t_tmp in t:(total_nyear-1)){
              SR_mat[t_tmp+1, ,"deviance"] <- SR_mat[t_tmp, ,"deviance"]*SR_mat[t_tmp,,"rho"] + SR_mat[t_tmp, ,"rand_resid"]
            }
            # そんで、バイアス補正もする
            SR_mat[(t+1):total_nyear,,"deviance"] <- SR_mat[(t+1):total_nyear,,"deviance"] - SR_mat[(t+1):total_nyear,,"bias_factor"]
          }
        }
      } # close else in [all(N_mat[1,t,]==0)]

      if(is.na(N_mat[1,t,1])) stop("Error: Recruitment cannot be estimated correctly...")
      # 加入が入力されたあとで０歳のSSB計算用の体重と漁獲量計算用の体重を更新する
      if(is_waa_fun)
        waa_mat[1,t,]       <- update_waa_mat(t=t,waa=waa_mat,rand=waa_rand_mat,naa=N_mat,
                                             pars_b0=waa_par_mat[,,"b0"],pars_b1=waa_par_mat[,,"b1"])[1,]
      if(is_waa_catch_fun)
        waa_catch_mat[,t,] <- update_waa_catch_mat(t=t,waa=waa_catch_mat,rand=waa_catch_rand_mat,naa=N_mat,
                                                   pars_b0=waa_catch_par_mat[,,"b0"],pars_b1=waa_catch_par_mat[,,"b1"])
    }

    if(t>=start_ABC_year){
      # harvest control rule
      ssb_tmp <- spawner_mat[cbind(t-HCR_mat[t,,"year_lag"],1:nsim)]
      HCR_realized[t,,"beta_gamma"] <- HCR_function(ssb_tmp, HCR_mat[t,,"Blimit"],
                                                    HCR_mat[t,,"Bban"], HCR_mat[t,,"beta"],allyear_name[t])
      F_mat[,t,] <- sweep(F_mat[,t,],2,HCR_realized[t,,"beta_gamma"],FUN="*")
    }

   if(isTRUE(do_MSE) && t>=start_ABC_year){
     MSE_input_data$input$silent <- TRUE
     # ここでmax_Fの設定を上書きするようにしていたけど、それを廃止
     # MSE_input_dataそのままの設定を使うようにする
#     MSE_input_data$input$max_F <- max_F_MSE
     #     MSE_input_data$input$max_exploitation_rate <- max_exploitation_rate_MSE
     MSE_dummy_data <- safe_call(make_future_data,MSE_input_data$input)$data
     MSE_dummy_data <- MSE_dummy_data %>%
        purrr::list_modify(future_initial_year   = t-2,
                           start_random_rec_year = t-1,
                           start_ABC_year        = t,
                           total_nyear           = t,
                           x                     = 0, # = 1 in normal scale
                           what_return           = "stat",
                           do_MSE                = FALSE)
      for(i in 1:nsim){
        #                cat(t,"-",i,"\n")
        MSE_dummy_data$naa_mat[] <-  N_mat[,,i] # true dynamics
        MSE_dummy_data$naa_mat[,(t-1):t,] <-  0 # estiamted as future
        MSE_dummy_data$faa_mat[,1:(t-1),] <-  F_mat[,1:(t-1),i] # we know true F even in future
        MSE_dummy_data$faa_mat[,t,] <-  MSE_input_data$data$faa[,t,i] # alpha in ABC year depends on future SSB
        MSE_dummy_data$waa_mat[] <-  waa_mat[,,i] # in case
        MSE_dummy_data$waa_catch_mat[] <-  waa_catch_mat[,,i] # in case
        MSE_dummy_data$maa_mat[] <-  maa_mat[,,i] # in case
        MSE_dummy_data$M_mat[]   <-  M_mat[,,i] # in case
        # MSEのシミュレーション内では繰越設定はオフにする（自然な条件下でABCがいくつになるか知りたいので。繰越がある場合はこのあとで漁獲量が調整されるので大丈夫）
        MSE_dummy_data$HCR_mat[,,"TAC_reserve_rate"] <- NA
        MSE_dummy_data$HCR_mat[,,"TAC_carry_rate"] <- NA
        MSE_dummy_data$HCR_mat[,,"TAC_reserve_amount"] <- NA #
        MSE_dummy_data$HCR_mat[,,"TAC_carry_amount"] <- NA   #
        # 同様にTACの変動の上限設定もオフにする
        MSE_dummy_data$HCR_mat[,,"TAC_upper_CV"] <- NA
        MSE_dummy_data$HCR_mat[,,"TAC_lower_CV"] <- NA

        # TACどおりに漁獲すると将来予測でも仮定して将来予測する!!
        if(MSE_catch_exact_TAC==TRUE) MSE_dummy_data$HCR_mat[t-1,,"expect_wcatch"] <- HCR_mat[t-1,i,"expect_wcatch"]

        # 漁獲量に上限設定があってそれが厳しい場合に上限を予測値から決定しないといけない
        if(sum(HCR_mat[t,i,"expect_wcatch"]>0)>0) MSE_dummy_data$HCR_mat[t,,"expect_wcatch"] <- HCR_mat[t,i,"expect_wcatch"]

        for(k in 1:MSE_nsim){
          MSE_dummy_data$SR_mat[,k,]  <- SR_mat[,i,]
          MSE_dummy_data$SR_mat[,k,"ssb"]  <- spawner_mat[,i] # true ssb
          MSE_dummy_data$SR_mat[,k,"recruit"]  <- N_mat[1,,i] # true recruit
        }
        # re-calculate past deviance and produce random residual in future
        if(!is.null(MSE_sd) && MSE_sd==0){
          MSE_input_data$input$res_SR$pars$sd[] <- 0
          if(!is.null(MSE_input_data$input$res_SR$regime_pars)) MSE_input_data$input$res_SR$regime_pars$sd[] <- 0
        }
        MSE_dummy_data$SR_mat <-
          set_SR_mat(res_vpa   = NULL, # past deviande is calculated by true ssb
                     res_SR    = MSE_input_data$input$res_SR,
                     SR_mat    = MSE_dummy_data$SR_mat,
                     seed_number=MSE_seed,
                     start_random_rec_year_name = dimnames(naa_mat)[[2]][t-1],
                     recruit_age = recruit_age,
                     scale_ssb=scale_ssb,
                     scale_R  =scale_R,                     
                     resid_type                 = MSE_input_data$input$resid_type,
                     resample_year_range        = dimnames(naa_mat)[[2]][1]:dimnames(naa_mat)[[2]][t-2],
                     backward_duration          = MSE_input_data$input$backward_duration,
                     bias_correction            = MSE_input_data$input$bias_correction,
                     recruit_intercept          = MSE_input_data$input$recruit_intercept,
                     setting_release            = MSE_input_data$input$setting_release,
                     model_average_option       = MSE_input_data$input$model_average_option,
                     regime_shift_option        = MSE_input_data$input$regime_shift_option,
                     fix_recruit                = MSE_input_data$input$fix_recruit)
        # fix_recruitの年だけ加入を固定する（めんどくさいけど、今後要改善）
        if(!is.null(MSE_input_data$input$fix_recruit)){
          fix_year <- which(allyear_name%in%MSE_input_data$input$fix_recruit$year)
          MSE_dummy_data$naa_mat[1,fix_year,] <- MSE_dummy_data$SR_mat[fix_year,,"recruit"]
        }
        res_tmp <- safe_call(future_vpa_R,MSE_dummy_data) # do future projection
        HCR_mat[t,i,"expect_wcatch"] <- mean(apply(res_tmp$wcaa[,t,],2,sum)) # determine ABC in year t here
        SR_MSE[t,i,"recruit"] <- mean(res_tmp$naa[1,t,])
        SR_MSE[t,i,"ssb"]     <- mean(res_tmp$SR_mat[t,,"ssb"])
        # MSEでなく本来のFで漁獲していたらどうなっていたかの計算
        if(Pope==1){
          SR_MSE[t,i,"real_true_catch"] <- sum(N_mat[,t,i]*(1-exp(-F_mat[,t,i]))*exp(-M_mat[,t,i]/2) * waa_catch_mat[,t,i])
        }
        else{
          SR_MSE[t,i,"real_true_catch"] <- sum(N_mat[,t,i]*(1-exp(-F_mat[,t,i]-M_mat[,t,i]))*F_mat[,t,i]/(F_mat[,t,i]+M_mat[,t,i]) * waa_catch_mat[,t,i])
        }

        MSE_seed <- MSE_seed+1
      }
    }

    # TAC carry over setting
    if((sum(!is.na(HCR_mat[t,,"TAC_reserve_rate"]  ))>0) ||
       (sum(!is.na(HCR_mat[t,,"TAC_reserve_amount"]))>0)  ){
      if(sum(HCR_mat[t,,"expect_wcatch"])==0){
        # non-MSE
        HCR_realized[t,,"original_ABC"] <-
          catch_equation(N_mat[,t,],F_mat[,t,],waa_catch_mat[,t,],M_mat[,t,],Pope=Pope) %>% colSums()
      }
      else{
        # MSE
        HCR_realized[t,,"original_ABC"] <- HCR_mat[t,,"expect_wcatch"]
      }
      HCR_realized[t,,"original_ABC_plus"] <- HCR_realized[t,,"original_ABC"] + HCR_realized[t,,"reserved_catch"]
      if(all(!is.na(HCR_mat[t,,"TAC_reserve_rate"]))){
        HCR_mat[t,,"expect_wcatch"] <- HCR_realized[t,,"original_ABC_plus"] * (1-HCR_mat[t,,"TAC_reserve_rate"])
      } #
      if(all(!is.na(HCR_mat[t,,"TAC_reserve_amount"]))){
        tmpcatch <- HCR_realized[t,,"original_ABC_plus"] - HCR_mat[t,,"TAC_reserve_amount"]
        HCR_mat[t,,"expect_wcatch"] <- ifelse(tmpcatch<0, 0.01, tmpcatch)
      } #expect_wcatchをゼロにすると不具合がありそうなので、微小値（0.01）を与える

      if(t<total_nyear){
        if(all(!is.na(HCR_mat[t,,"TAC_carry_rate"]))){
          max_carry_amount <- HCR_mat[t,,"TAC_carry_rate"]*HCR_realized[t,,"original_ABC"]
        } #
        if(all(!is.na(HCR_mat[t,,"TAC_carry_amount"]))){
          max_carry_amount <- HCR_mat[t,,"TAC_carry_amount"]
        }                             #
        ABC_reserve_amount <- HCR_realized[t,,"original_ABC"] - HCR_mat[t,,"expect_wcatch"]
        ABC_reserve_amount[ABC_reserve_amount<0] <- 0
        HCR_realized[t+1,,"reserved_catch"] <- cbind(max_carry_amount, ABC_reserve_amount) %>%
          apply(1,min)
      }
    }

    # 漁獲量の変動の上限・下限設定をする場合
    if(t>=start_ABC_year && sum(!is.na(HCR_mat[t,,"TAC_upper_CV"]))){
      # expect_wcatchが全部空だったらexpect catchを計算して入れる
      if(all(HCR_mat[t,,"expect_wcatch"]==0)){
        HCR_mat[t,,"expect_wcatch"] <- catch_equation(N_mat[,t,],F_mat[,t,],waa_catch_mat[,t,],M_mat[,t,],Pope=Pope) %>% colSums()
      }
      # CVよりも小さい・大きかったら上限値にexpect_wcatchを置き換える
      HCR_mat[t,, "expect_wcatch"] <-
        set_upper_limit_catch(HCR_realized[t-1,,"wcatch"], HCR_mat[t,,"expect_wcatch"], HCR_mat[t,,"TAC_upper_CV"])
    }

    if(t>=start_ABC_year && sum(!is.na(HCR_mat[t,,"TAC_lower_CV"]))){
      # expect_wcatchが全部空だったらexpect catchを計算して入れる
      if(all(HCR_mat[t,,"expect_wcatch"]==0)){
        HCR_mat[t,,"expect_wcatch"] <- catch_equation(N_mat[,t,],F_mat[,t,],waa_catch_mat[,t,],M_mat[,t,],Pope=Pope) %>% colSums()
      }
      HCR_mat[t,, "expect_wcatch"] <-
        set_lower_limit_catch(HCR_realized[t-1,,"wcatch"], HCR_mat[t,,"expect_wcatch"], HCR_mat[t,,"TAC_lower_CV"])
    }

    if(sum(HCR_mat[t,,"expect_wcatch"])>0){
      F_max_tmp <- apply(faa_mat[,t,],2,max) # betaを乗じる前のもとのfaaを用いる
      # F_mat[,t,]がものすごく小さい値になっている場合、そのままで入れるとFへの乗数が壁に当たる場合があるので最大１で正規化する
      saa.tmp <- sweep(faa_mat[,t,],2,F_max_tmp,FUN="/")
      fix_catch_multiplier <- purrr::map_dbl(which(F_max_tmp>0),
                                             function(x) caa.est.mat(N_mat[,t,x],saa.tmp[,x],#F_mat[,t,x],#saa.tmp[,x],
                                                                     waa_catch_mat[,t,x],M_mat[,t,x],
                                                                     HCR_mat[t,x,"expect_wcatch"],
                                                                     max_exploitation_rate=max_exploitation_rate,
                                                                     max_F=max_F,
                                                                     Pope=Pope)$x)
      F_mat[,t,which(F_max_tmp>0)] <- sweep(saa.tmp[,which(F_max_tmp>0)],2, fix_catch_multiplier, FUN="*")
      HCR_realized[t,which(F_max_tmp>0),"beta_gamma"] <- HCR_realized[t,which(F_max_tmp>0),"beta_gamma"] *
        fix_catch_multiplier / F_max_tmp[which(F_max_tmp>0)]
    }

    if(t<total_nyear){
      # forward calculation
      for(iage in 1:(plus_age-1)) {
        N_mat[iage+1,t+1,] <- N_mat[iage,t,]*exp(-M_mat[iage,t,]-F_mat[iage,t,])
      }
      if(plus_group == TRUE) N_mat[plus_age,t+1,] <- N_mat[plus_age,t+1,] + N_mat[plus_age,t,]*exp(-M_mat[plus_age,t,]-F_mat[plus_age,t,])
      # waaとmaaの更新(ここ、不要では？)
      ## if(is_waa_fun)
      ##   waa_mat[,t,]       <- update_waa_mat(t=t,waa=waa_mat,rand=waa_rand_mat,naa=N_mat,
      ##                                        pars_b0=waa_par_mat[,,"b0"],pars_b1=waa_par_mat[,,"b1"])
      ## if(is_waa_catch_fun)
      ##   waa_catch_mat[,t,] <- update_waa_catch_mat(t=t,waa=waa_catch_mat,rand=waa_catch_rand_mat,naa=N_mat,
      ##                                              pars_b0=waa_catch_par_mat[,,"b0"],pars_b1=waa_catch_par_mat[,,"b1"])    
      ## if(is_maa_fun) maa_mat[,t,] <- update_maa_mat(maa=maa_mat[,t,],rand=maa_rand_mat[,t,],naa=N_mat[,t,],
      ##                                               pars_b0=maa_par_mat[,,"b0"],pars_b1=maa_par_mat[,,"b1"],
      ##                                               min_value=maa_par_mat[,,"min"],max_value=maa_par_mat[,,"max"])
    }
    HCR_realized[t,,"wcatch"] <- catch_equation(N_mat[,t,],F_mat[,t,],waa_catch_mat[,t,],M_mat[,t,],Pope=Pope) %>% colSums()
  }

  if(Pope==1){
    wcaa_mat <- N_mat*(1-exp(-F_mat))*exp(-M_mat/2) * waa_catch_mat
  }
  else{
    wcaa_mat <- N_mat*(1-exp(-F_mat-M_mat))*F_mat/(F_mat+M_mat) * waa_catch_mat
  }
  HCR_realized[,,"wcatch"] <- apply(wcaa_mat,c(2,3),sum)

  if(isTRUE(do_MSE)){
    F_pseudo_mat <- MSE_input_data$data$faa
    beta_gamma <- HCR_function(spawner_mat,
                              MSE_input_data$data$HCR_mat[,,"Blimit"],
                              MSE_input_data$data$HCR_mat[,,"Bban"],
                              MSE_input_data$data$HCR_mat[,,"beta"],
                              allyear_name[t])
    F_pseudo_mat[] <- sweep(F_pseudo_mat,c(2,3),beta_gamma,FUN="*")

    if(Pope==1){
      wcaa_tmp <- N_mat*(1-exp(-F_pseudo_mat))*exp(-M_mat/2) * waa_catch_mat
    }
    else{
      wcaa_tmp <- N_mat*(1-exp(-F_pseudo_mat-M_mat))*F_pseudo_mat/
        (F_pseudo_mat+M_mat) * waa_catch_mat
    }
    SR_MSE[,,"pseudo_true_catch"] <- apply(wcaa_tmp, c(2,3), sum)
  }

  if(objective<2){
    last_catch <- colSums(wcaa_mat[,total_nyear,])
    if(obj_stat==0) obj <- mean(last_catch)
    if(obj_stat==1) obj <- geomean(last_catch)
  }
  else{
    if(obj_stat==0) obj <- mean(spawner_mat[total_nyear,])
    if(obj_stat==1) obj <- geomean(spawner_mat[total_nyear,])
    if(obj_stat==2) obj <- median(spawner_mat[total_nyear,])
  }


  {if(objective==0) obj <- -log(obj) # MSY case
    else{
      obj <- (log(obj/obj_value))^2 # PGY, B0% case
    }
  }

  if(what_return=="obj")  return(obj)
  if(what_return=="stat"){
    tmb_data$SR_mat[,,"ssb"]  <- spawner_mat
    tmb_data$SR_mat[,,"recruit"]   <- N_mat[1,,]
    tmb_data$SR_mat[,,"biomass"]   <- apply(N_mat*waa_mat,c(2,3),sum)
    tmb_data$SR_mat[,,"cbiomass"]  <- apply(N_mat*waa_catch_mat,c(2,3),sum)  
    res <- list(naa=N_mat, wcaa=wcaa_mat, faa=F_mat, SR_mat=tmb_data$SR_mat,maa=maa_mat,
                HCR_mat=HCR_mat,HCR_realized=HCR_realized,multi=exp(x),waa=waa_mat, waa_catch_mat=waa_catch_mat)
    if(isTRUE(do_MSE)) res$SR_MSE <- SR_MSE
    return(res)
  }

}

# utility functions for stock recruitment relationship ----

#'
#' 将来予測用の再生産関係の設定を行う関数
#'
#' 再生産関係をres_SRで与えると、res_vpaを見ながら残差を再計算したのち、start_random_rec_year_name以降の加入のdeviationを計算しSR_mat[,,"deviance"]に入れる。
#'
#' @param res_vpa VPAの推定結果
#' @param res_SR 再生産関係の推定結果
#' @param SR_mat 将来予測用の再生産関係パラメータが格納する３次元行列
#' @param seed_number シード番号
#' @param start_random_rec_year_name ランダム加入を仮定する最初の年
#' @param resid_type 残差の発生パターン；対数正規分布は"lognormal"、単純リサンプリングは"resampling"、backward-resamplingは"backward"
#' @param resample_year_range "resampling", "backward"で有効。年の範囲を入れると、対象とした年の範囲で計算される残差を用いる。
#' @param backward_duration "backward"の場合、何年で１ブロックとするか。"backward"で有効。デフォルトは5 。
#' @param model_average_option model averagingをする場合のオプション. SR_matのlistとweightをlist形式で入れる(list(SR_list=list(res_SR1,res_SR2),weight=c(0.5,0.5))). 上で設定されたres_SRは使われない.
#' @param regime_shift_option レジームシフトを仮定する場合のオプション. この場合, res_SRにはfit.SRregimeの結果オブジェクトを入れる. オプションの設定は list(future_regime=将来のregimeの仮定。keyで指定された番号を入れる)
#' @param setting_release 詳細はmake_future_dataへ
#'
#' @export
#' @encoding UTF-8

set_SR_mat <- function(res_vpa=NULL,
                       start_random_rec_year_name,
                       SR_mat,
                       res_SR,
                       seed_number,
                       resid_type="lognormal",
                       bias_correction=TRUE,
                       resample_year_range=NA,
                       backward_duration=5,
                       recruit_intercept=0,
                       setting_release=NULL,
                       recruit_age=0,
                       scale_ssb=1,
                       scale_R=1,                       
                       model_average_option=NULL,
                       regime_shift_option=NULL,
                       fix_recruit=NULL
){

  nsim <- dim(SR_mat)[[2]]
  allyear_name <- dimnames(SR_mat)[[1]]
  start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
  random_rec_year_period <- (start_random_rec_year):length(allyear_name)

  # check arguments in fix_recruit
  assertthat::assert_that(sum(!fix_recruit$year %in% allyear_name)==0)
  if(!resid_type%in%c("lognormal","resample","backward")) stop("resid_type is invalid.")

  if(!is.null(model_average_option)){
    weight <- arrange_weight(model_average_option$weight,nsim)
    SR_mat <- average_SR_mat(res_vpa     = res_vpa,
                             res_SR_list = model_average_option$SR_list,
                             range_list  = weight,
                             SR_mat      = SR_mat,
                             seed_number = seed_number+1,
                             start_random_rec_year_name=start_random_rec_year_name,
                             resid_type  = resid_type,
                             recruit_intercept=recruit_intercept,
                             resample_year_range=resample_year_range,
                             recruit_age = recruit_age,
                             regime_shift_option = regime_shift_option,
                             bias_correction = bias_correction
    )
  }

  if(is.null(model_average_option)){
    # define SR function
    if(res_SR$input$SR=="HS"){
      SR_mat[,,"SR_type"] <- 1
      SRF <- SRF_HS
    }
    if(res_SR$input$SR=="BH"){
      SR_mat[,,"SR_type"] <- 2
      SRF <- SRF_BH
    }
    if(res_SR$input$SR=="RI"){
      SR_mat[,,"SR_type"] <- 3
      SRF <- SRF_RI
    }
    if(res_SR$input$SR=="Shepherd" | res_SR$input$SR=="SH"){
      SR_mat[,,"SR_type"] <- 4
      SRF <- SRF_SH
    }
    if(res_SR$input$SR=="Cushing" | res_SR$input$SR=="CU"){
      SR_mat[,,"SR_type"] <- 5
      SRF <- SRF_CU
    }
    if(res_SR$input$SR=="BHS"){
      SR_mat[,,"SR_type"] <- 6
      SRF <- SRF_BHS
    }            

    # define SR parameter
    if(is.null(regime_shift_option)){ # when no-regime shift
      SR_mat[,,"a"] <- res_SR$pars$a
      SR_mat[,,"b"] <- res_SR$pars$b
      SR_mat[,,"sd"] <- res_SR$pars$sd
    }
    else{ # when regime_shift
      regime_data <- res_SR$regime_resid %>%
        left_join(res_SR$regime_pars, by="regime") %>%
        bind_cols(res_SR$input$SRdata)
      SR_mat[as.character(regime_data$year),,"a"] <- regime_data$a
      SR_mat[as.character(regime_data$year),,"b"] <- regime_data$b
      SR_mat[as.character(regime_data$year),,"sd"] <- regime_data$sd
      future_regime_par <- res_SR$regime_pars %>% dplyr::filter(regime==regime_shift_option$future_regime)
      SR_mat[random_rec_year_period,,"a"] <- future_regime_par$a
      SR_mat[random_rec_year_period,,"b"] <- future_regime_par$b
      SR_mat[random_rec_year_period,,"sd"] <- future_regime_par$sd

      # どのレジームに属するか明示的に指定されないケースが出てくる。
      # そういう場合は将来予測のレジームのパラメータを使うようにする→要改善
      missing_year <- which(!(1:dim(SR_mat)[[1]] %in% c(which(allyear_name %in% as.character(regime_data$year)),random_rec_year_period)))

      if(length(missing_year)>0){
        SR_mat[missing_year,,"a"] <- future_regime_par$a
        SR_mat[missing_year,,"b"] <- future_regime_par$b
        SR_mat[missing_year,,"sd"] <- future_regime_par$sd
      }

      res_SR$pars$rho <- 0
    }

    # set gamma parameter (暫定. regimeありでもなしでも同じgammaがinput$gammaで与えられている場合の特殊ケース)
    SR_mat[,,"gamma"] <- ifelse(!is.null(res_SR$input$gamma), res_SR$input$gamma, NA)
#    if(res_SR$input$SR!="Shepherd"&&res_SR$input$SR!="SH") SR_mat[,,"gamma"] <- NA
      
    SR_mat[,,"rho"] <- res_SR$pars$rho
    SR_mat[random_rec_year_period,,"intercept"] <- recruit_intercept # future intercept

    # ここ大事なところ。res_vpaが別に与えられている（VPA結果が更新されている）場合には、SR_matのssbとRについても更新された値に置き換える。過去の加入の残差はここで計算するため、将来の自己相関を計算したりするときに影響する
    if(!is.null(res_vpa)){
      SR_mat[1:(start_random_rec_year-1),,"ssb"] <- as.numeric(colSums(res_vpa$ssb,na.rm=T))[1:(start_random_rec_year-1)]
      SR_mat[1:(start_random_rec_year-1),,"recruit"] <- as.numeric(res_vpa$naa[1,1:(start_random_rec_year-1)])
    }

    ## 放流関連の設定. res_SRとres_vpaのどちらの情報を使うかまず判断する
    data_SR_release <- res_SR$input$SRdata  
    if(!is.null(res_vpa)){
      if(!is.null(setting_release) && is.null(setting_release$data_source) ||
           (!is.null(setting_release) && setting_release$data_source=="VPA")){
        data_SR_release <- get.SRdata(res_vpa)
      }}
    
    # 過去データの入力
    if("release" %in% str_sub(names(data_SR_release),1,7))
      if("release_alive" %in% names(data_SR_release)) SR_mat[as.character(data_SR_release$year),,"intercept"] <- data_SR_release$release_alive
      if("release" %in% names(data_SR_release)) SR_mat[as.character(data_SR_release$year),,"intercept"] <- data_SR_release$release        

    ## 未来の放流関連の設定
    if(!is.null(setting_release)){
      assert_that(!is.null(data_SR_release$release_ratealive), TRUE)
      assert_that(!is.null(data_SR_release$release_alive),     TRUE)
      assert_that(ncol(setting_release$number)<3 , TRUE)
      assert_that(ncol(setting_release$rate  )==1, TRUE)      

      tmp <- SR_mat[random_rec_year_period,,"intercept"]
      # for value_average
      if(ncol(setting_release$number)==1){ # yearとvalueの両方が入っているか・いないか
        if("year" %in% colnames(setting_release$number)){
          tmp[] <- data_SR_release %>%
            dplyr::filter(year %in% setting_release$number$year) %>%
            select(release_all) %>% unlist() %>% mean(na.rm=TRUE)
        }
        if("value" %in% colnames(setting_release$number)){
          tmp[] <- mean(setting_release$number$value, na.rm=TRUE)
        }
      }else{
        # exclude data without random-recruitment assumption
        setting_release$number <- setting_release$number %>%
            dplyr::filter(year>=start_random_rec_year_name)
        tmp[dimnames(tmp)[[1]] %in% setting_release$number$year,] <- setting_release$number$value[setting_release$number$year %in% dimnames(tmp)[[1]]]
      }

      if(!is.null(setting_release$rate$year)){
        rate_value <- data_SR_release %>%
            dplyr::filter(year %in% setting_release$number$year) %>%
            select(release_ratealive) %>% unlist()
      }
      else{
        rate_value <- setting_release$rate$value
      }
      tmp[] <- tmp[] * sample(rate_value, dim(tmp)[[1]] * dim(tmp)[[2]], replace=TRUE) 

      SR_mat[random_rec_year_period,,"intercept"] <- tmp
    }
        
    recruit_range <- (recruit_age+1):(start_random_rec_year-1)
    ssb_range     <- 1:(start_random_rec_year-1-recruit_age)

    # re-culcurate past recruitment deviation
    # intercept=release fish
    SR_mat[recruit_range,,"deviance"] <- SR_mat[recruit_range,,"rand_resid"] <-
        log(SR_mat[recruit_range,,"recruit"]-SR_mat[recruit_range,,"intercept"]) -
        log(scale_R*SRF(SR_mat[ssb_range,,"ssb"]*scale_ssb,SR_mat[recruit_range,,"a"],SR_mat[recruit_range,,"b"],SR_mat[recruit_range,,"gamma"]))

      # define future recruitment deviation
      set.seed(seed_number)

      if(resid_type=="lognormal"){
        if(isTRUE(bias_correction)){
          #            sd_with_AR <- sqrt(res_SR$pars$sd^2/(1-res_SR$pars$rho^2))
          #            bias_factor <- 0.5* sd_with_AR^2
          sd_with_AR <- sqrt(SR_mat[,,"sd"]^2/(1-SR_mat[,,"rho"]^2))
          SR_mat[,,"bias_factor"] <- 0.5 * sd_with_AR^2
          SR_mat[-random_rec_year_period,,"bias_factor"] <- 0
        }
        else{
          #            bias_factor <- 0
          SR_mat[,,"bias_factor"] <- 0
        }

        tmp_SR <- t(SR_mat[random_rec_year_period,,"rand_resid"])
        tmp_SR[] <- rnorm(nsim*length(random_rec_year_period), mean=0,
                          sd=t(SR_mat[random_rec_year_period,,"sd"]))
        SR_mat[random_rec_year_period,,"rand_resid"] <- t(tmp_SR)

        for(t in random_rec_year_period){
          SR_mat[t, ,"deviance"] <- SR_mat[t-1, ,"deviance"]*SR_mat[t,,"rho"] + SR_mat[t, ,"rand_resid"]
        }
        SR_mat[random_rec_year_period,,"deviance"] <- SR_mat[random_rec_year_period,,"deviance"] - SR_mat[random_rec_year_period,,"bias_factor"]
      }

      if(resid_type=="resample" | resid_type=="backward"){
        # 推定された残差をそのまま使う
        #if(resample_year_range==0){
        #  resample_year_range <- sort(res_SR$input$SRdata$year[res_SR$input$w==1])
        #}

        sampled_residual <- SR_mat[as.character(resample_year_range),,"rand_resid"]
        if(isTRUE(bias_correction)){
          #            bias_factor <- log(colMeans(exp(sampled_residual)))
          SR_mat[random_rec_year_period,,"bias_factor"] <- rep(log(colMeans(exp(sampled_residual))),
                                                               each=length(random_rec_year_period))
        }
        else{
          #            bias_factor <- rep(0,ncol(sampled_residual))
          SR_mat[random_rec_year_period,,"bias_factor"] <- 0
        }
        for(i in 1:ncol(sampled_residual)){
          if(resid_type=="resample"){
            SR_mat[random_rec_year_period,i,"rand_resid"] <- sample(sampled_residual[,i], length(random_rec_year_period), replace=TRUE)
          }
          if(resid_type=="backward"){
            SR_mat[random_rec_year_period,i,"rand_resid"] <- sample_backward(sampled_residual[,i], length(random_rec_year_period), backward_duration)
          }
          SR_mat[random_rec_year_period,i,"deviance"] <- SR_mat[random_rec_year_period,i,"rand_resid"]-SR_mat[random_rec_year_period,i,"bias_factor"]
        }
      }
  }

  # when fix recruitment
  if(!is.null(fix_recruit)){
    if(!is.list(fix_recruit$rec)){
      # scalar
      SR_mat[as.character(fix_recruit$year),,"recruit"] <- fix_recruit$rec
    }
    else{
      # vector
      for(i in seq_len(length(fix_recruit$year))){
        if(length(fix_recruit$rec[[i]])!=dim(SR_mat)[[2]]) stop("invalid length of recruit")
        SR_mat[as.character(fix_recruit$year[i]),,"recruit"] <- as.numeric(unlist(fix_recruit$rec[i]))
        
      }
    }
  }

  return(SR_mat)
}

#' @export
SRF_HS <- function(x,a,b,gamma){
    ifelse(x>b,b*a,x*a)
}

#' @export
SRF_BH <- function(x,a,b,gamma){
    a*x/(1+b*x)
}

#' @export
SRF_RI <- function(x,a,b,gamma){
    a*x*exp(-b*x)
}

#' @export
SRF_SH <- function(x,a,b,gamma){
    a*x/(1+(b*x)^gamma)
}

#' @export
SRF_CU <- function(x,a,b,gamma){
    a*x^b
}

#' @export
SRF_BHS <- function(x,a,b,gamma){
    res <- ifelse(x < b, a*b*(x/b)^{1-((x)/b)^gamma}, a*b)
    return(res)    
}

#'
#' 将来予測用の三次元行列（年齢×年×シミュレーション）を与えられたら, pars.yearで指定された期間のパラメータを平均するか、parで指定されたパラメータを、year_replace_future以降の年で置き換える
#'
#' @param d3_mat 将来予測用の３次元行列
#' @param pars 置き換えるべき生物パラメータ
#' @param pars.year この期間の生物パラメータを平均して、将来のパラメータとする
#' @param year_replace_future 生物パラメータを置き換える最初の年
#' @encoding UTF-8
#' @export
#'

make_array <- function(d3_mat, pars, pars.year, year_replace_future){
  if(length(dim(pars))==3){
    return(pars)
  }
  else{
    years <- dimnames(d3_mat)[[2]]
    if(is.null(pars)){
      pars.future <- rowMeans(d3_mat[,years%in%pars.year,1])
    }
    else{
      if(length(pars)==dim(d3_mat)[[1]]) pars.future <- pars
      else stop("length of parameter is different from what is expected.")
    }
    d3_mat[,which(year_replace_future==years):length(years),] <- pars.future
    if(sum(year_replace_future==years)==0) stop("year_replace_future is invalid.")

    return(d3_mat)
  }
}

#'
#' モデルのAICから計算される重みに応じて、複数のモデルからシミュレーション結果を抜き出す
#'
#' @param weight 複数モデルのAkaike weightのベクトル
#' @param nsim シミュレーション回数
#' @export
arrange_weight <- function(weight, nsim){
  weight <- weight / sum(weight)
  weight <- round(cumsum(weight) * nsim)
  weight2 <- c(1,weight[-length(weight)]+1)
  purrr::map(seq_len(length(weight)),function(x) weight2[x]:weight[x])
}

#'
#' モデル平均的な再生産関係を与える
#'
#' @param res_vpa VPAの推定結果
#' @param res_SR_list 再生産関係の推定結果のリスト
#' @param range_list
#' @param SR_mat 将来予測用の再生産関係パラメータを格納する３次元行列
#' @param seed_number シード番号
#' @param start_random_rec_year_name ランダム加入を仮定する最初の年
#' @param resid_type 残差の発生パターン；対数正規分布は"lognormal"、単純リサンプリングは"resampling"
#' @param resample_year_range 0の場合、推定に使ったデータから計算される残差を用いる。年の範囲を入れると、対象とした年の範囲で計算される残差を用いる
#'
#' @export
#' @encoding UTF-8

average_SR_mat <- function(res_vpa,
                           res_SR_list,
                           range_list=list(1:500,501:1000),
                           SR_mat,
                           seed_number,
                           start_random_rec_year_name,
                           recruit_age,
                           resid_type="lognormal",
                           resample_year_range=0,
                           regime_shift_option=NULL,
                           recruit_intercept=0,
                           bias_correction=TRUE){

  allyear_name <- dimnames(SR_mat)[[1]]
  start_random_rec_year  <- which(allyear_name==start_random_rec_year_name)
  random_rec_year_period <- (start_random_rec_year):length(allyear_name)

  for(i in seq_len(length(res_SR_list))){
    SR_mat_tmp <- set_SR_mat(res_vpa=res_vpa,
                             start_random_rec_year_name,
                             SR_mat=SR_mat,
                             res_SR=res_SR_list[[i]],
                             seed_number=seed_number+i,
                             resid_type=resid_type,
                             resample_year_range=resample_year_range,
                             recruit_age=recruit_age,
                             recruit_intercept=recruit_intercept,
                             bias_correction=bias_correction,
                             regime_shift_option=regime_shift_option)
    SR_mat[,as.character(range_list[[i]]),] <-
      SR_mat_tmp[,range_list[[i]],]
  }

  return(SR_mat)
}

#' 過去にさかのぼってブロックサンプリングをおこなう
#'
#' @param residual リサンプリングする残差
#' @param n 将来にわたって何年分のリサンプリング残差を作るか
#' @param duration 1ブロックの年の長さ
#'
#' @exmaples
#'
#' set.seed(1)
#' res <- sample_backward(rep(1:5,each=5), 30, 5)
#' apply(matrix(res,5,6),2,min)
#'
#' @export
#' @encoding UTF-8
#'

sample_backward <- function(residual, n, duration){
  residual_rev <- rev(residual)
  nblock <- floor(length(residual_rev)/duration)
  block <- (1:(nblock))*(duration)
  block2 <- block-(duration)+1
  block[length(block)] <- length(residual_rev)
  block.list <- purrr::map(seq_len(length(block)),function(x) residual_rev[block2[x]:block[x]])
  # calculate sampling probability in the case of tail block (different length)
  block.probability <- sapply(block.list,length)
  block.probability <- block.probability/sum(block.probability)
  resid_future <- ceiling(1:n/duration)
  resid_future <- ifelse(resid_future>nblock,nblock,resid_future)

  for(i in 1:n){
    # block is resampled according to the data length of each block
    if(i%%duration==1){
      block_select <- sample(x=1:resid_future[i],size=1,prob=block.probability[1:resid_future[i]])
    }
    # iごとに指定されたブロックをリサンプリング＝重複あり
    resid_future[i] <- sample(x=block.list[[block_select]],size=1)
  }
  return(resid_future)
}


#' @export
print.myarray <- function(x) cat("array :", dim(x),"\n")


# utility functions for other parts -----

#'
#' @export
#' @encoding UTF-8

naming_adreport <- function(tmb_data, ad_report){
  #    ssb <- ad_report$spawner_mat
  wcaa <- ad_report$catch_mat
  naa <- ad_report$N_mat
  faa <- ad_report$F_mat
  dimnames(wcaa) <- dimnames(naa) <- dimnames(faa) <-
    dimnames(tmb_data$naa)
  class(wcaa) <- class(naa) <- class(faa) <- "myarray"

  tmb_data$SR_mat[,,"ssb"] <- ad_report$spawner_mat
  tmb_data$SR_mat[,,"recruit"] <- ad_report$N_mat[1,,]

  return(list(wcaa=wcaa, naa=naa, faa=faa,
              SR_mat        = tmb_data$SR_mat,
              HCR_mat       = tmb_data$HCR_mat,
              waa           = tmb_data$waa_mat,
              waa_catch_mat = tmb_data$waa_catch_mat))
}


#'
#' @export
#' @encoding UTF-8

trace_future <- function(tmb_data,
                         trace.multi=c(seq(from=0,to=0.9,by=0.1),1,
                                       seq(from=1.1,to=2,by=0.1),3:5,7,20,100),
                         ncore=0){
  #    R_obj_fun <- function(x, tmb_data, what_return="obj"){
  #        tmb_data$x <- x
  #        tmb_data$what_return <- what_return
  #        obj <- do.call(future_vpa_R, tmb_data)
  #        return(obj)
  #    }

  if(ncore==0){
    trace <- purrr::map_dfr(trace.multi, function(x)
      #       R_obj_fun(x=x, tmb_data = tmb_data, what_return="stat") %>%
      future_vpa(tmb_data,
                 optim_method="none",
                 multi_init=x,
                 multi_lower=x,
                 multi_upper=x) %>%
        get.stat(use_new_output=TRUE) %>% as_tibble())
  }
  else{
    library(foreach)
    #        library(doParallel)
    cl <- parallel::makeCluster(ncore, type="FORK")
    doParallel::registerDoParallel(cl)
    trace <- foreach::foreach(x=trace.multi,.combine="rbind")%dopar%{
      future_vpa(tmb_data,
                 optim_method="none",
                 multi_init=x,
                 multi_lower=x,
                 multi_upper=x) %>%
        get.stat(use_new_output=TRUE) %>% as_tibble()
    }
    parallel::stopCluster(cl)
  }

  return(trace)
}

#'
#' @export
#' @encoding UTF-8

get_summary_stat <- function(all.stat){

  sumvalue <- all.stat %>% as_tibble %>%
    mutate(SSB2SSB0=all.stat$ssb.mean/all.stat$ssb.mean[2]) %>%
    select(RP_name,ssb.mean,SSB2SSB0,biom.mean,cbiom.mean, U.mean,catch.mean,catch.CV,Fref2Fcurrent)
  colnames(sumvalue) <- c("RP_name","SSB","SSB2SSB0","B","cB","U","Catch","Catch.CV","Fref/Fcur")

  sumvalue <- bind_cols(sumvalue,all.stat[,substr(colnames(all.stat),1,1)=="F"])

  sumvalue$RP.definition <- NA
  sumvalue$RP.definition[sumvalue$RP_name=="MSY"]           <- "Btarget0"
  sumvalue$RP.definition[sumvalue$RP_name=="PGY_0.6_lower"] <- "Blimit0"
  sumvalue$RP.definition[sumvalue$RP_name=="PGY_0.1_lower"] <- "Bban0"
  sumvalue <- sumvalue %>% select(1,ncol(sumvalue),2:(ncol(sumvalue)-1))

  Fvector <- select(sumvalue,num_range("F",0:40))

  #    sumvalue$perSPR <- NA
  #    for(i in 1:nrow(Fvector)){
  #        sumvalue$perSPR[i] <- calc_perspr(input.list[[1]],Fvector[i,])
  #    }

  tibble::lst(sumvalue,Fvector)

}

#'
#' @export
#' @encoding UTF-8

format_to_old_future <- function(fout){
  fout_old <- fout[c("naa","faa","multi","input","waa")]
  #    fout_old$waa       <- fout$input$tmb_data$waa_mat
  fout_old$waa.catch <- fout$waa_catch_mat
  if(!is.null(fout$maa)){
      fout_old$maa       <- fout$maa #input$tmb_data$maa_mat
  }else{
      fout_old$maa       <- fout$input$tmb_data$maa_mat
  }
  fout_old$M         <- fout$input$tmb_data$M_mat
  fout_old$vssb      <- apply(fout$naa * fout_old$waa * fout_old$maa, c(2,3), sum, na.rm=T)
  fout_old$vbiom_catch <- apply(fout$naa * fout_old$waa.catch, c(2,3),sum, na.rm=T)  
  fout_old$vbiom     <- apply(fout$naa * fout_old$waa, c(2,3),sum, na.rm=T)
  fout_old$vwcaa     <- apply(fout$wcaa,c(2,3),sum, na.rm=T)
  fout_old$currentF  <- fout$faa[,fout$input$tmb_data$start_ABC_year-1,1]
  fout_old$futureF   <- fout$faa[,fout$input$tmb_data$start_ABC_year,1]
  fout_old$finalmeanF<- fout$faa[,dim(fout$faa)[[2]],] %>% apply(1,mean) # newly define
  fout_old$caa       <- fout$wcaa/fout_old$waa
  fout_old$multi     <- fout$multi
  fout_old$recruit   <- fout$SR_mat[,,"recruit"]
  if(!is.null(fout$HCR_realized)){
      fout_old$beta_gamma     <- fout$HCR_realized[,,"beta_gamma"]
      fout_old$alpha     <- fout$HCR_realized[,,"beta_gamma"]
      fout_old$Fratio     <- fout$HCR_realized[,,"Fratio"]
  }
  else{
    fout_old$beta_gamma     <- fout$HCR_mat[,,"beta_gamma"]
      fout_old$alpha     <- fout$HCR_mat[,,"beta_gamma"]
      fout_old$Fratio     <- fout$HCR_mat[,,"Fratio"]
  }
  fout_old$wcaa <- fout$wcaa
  return(fout_old)
}


#'
#' do.callのsafe版
#'
#' do.callで与えたリストの中にfuncで定義されていないものが混じっていた場合に、実際にdo.callを呼び出す前にerorrを吐いて関数をストップさせる。非常に大きいオブジェクトを与えていながらdo.callで上記の場面でエラーが出ると、デバッグモードで長時間Rがフリーズするのを避けるため。force=TRUEにすると、func内で定義されていない引数はリストから除外してdo.callを実行する.
#'
#' @export
#' @encoding UTF-8

safe_call <- function(func,args,force=FALSE,...){
  argname <- names(formals(func))
  check_argument <- names(args) %in% argname

  # make_future_dataでの引数追加への対応
  is_make_future_data <- sum("start_random_rec_year_name"==names(args))
  if(is_make_future_data){ # あとから追加された引数のリスト
    non_defined_arg <- names(args)[check_argument==FALSE]
    if("maa_fun" %in% non_defined_arg) args$maa_fun <- FALSE
    if("waa_catch_fun" %in% non_defined_arg) args$waa_catch_fun <- FALSE    
    if("start_waacatchfun_year_name" %in% non_defined_arg) args$start_waacatchfun_year_name <- args$start_biopar_year
    if("waa_fun_name" %in% non_defined_arg) args$waa_fun_name <- NA
    if("waa_catch_fun_name" %in% non_defined_arg) args$waa_catch_fun_name <- NA
    check_argument <- names(args) %in% argname
    force <- TRUE
  }

  if(sum(check_argument==FALSE)>0){
    if(force==FALSE){
      stop(paste(names(args)[check_argument==FALSE]), " is not used in func\n")
    }
    else{
      args <- args[check_argument==TRUE]
    }
  }
  return(do.call(func,args,...))
}

#'
#' @export
#' @encoding UTF-8

HCR_default <- function(ssb, Blimit, Bban, beta, year_name){
  beta_gamma <- beta
  tmp <- ssb < Blimit
  beta_gamma[tmp] <- beta[tmp]*(ssb[tmp]-Bban[tmp])/(Blimit[tmp]-Bban[tmp])
  beta_gamma[beta_gamma < 0] <- 0
  return(beta_gamma)
}


#'
#'
#'

if(0){
  update_waa_mat(waa=matrix(c(rep(0,5),10),2,3),
                 naa=matrix(1,2,3),
                 rand=matrix(0,2,3),
                 pars=array(c(0,0,log(1),log(1),2,2),dim=c(2,3),                                                            dimnames=list(age=1:2, pars=c("sd", "b0", "b1"))))
  #   pars
  #age sd b0 b1
  #  1  0  0  2
  #  2  0  0  2
}

update_waa_mat <- update_waa_catch_mat <-function(t,waa,rand,naa,pars_b0,pars_b1){
  waa_tmp <- exp(pars_b0+pars_b1*log(naa[,t,])+rand[,t,])
  replace_tmp <- waa[,t,]==0 & naa[,t,]>0  
  waa[,t,][replace_tmp] <- waa_tmp[replace_tmp] # ここでwaa=0のところにだけ数値を入れるので、もともと数値が入っていたら置き換わらない
  waa[,t,]
}

update_maa_mat <- function(maa,rand,naa,pars_b0,pars_b1,min_value,max_value){
  maa_tmp <- pars_b0+pars_b1*naa+rand
  maa_tmp[maa_tmp <= min_value] <- min_value[maa_tmp <= min_value]
  maa_tmp[maa_tmp >= max_value] <- max_value[maa_tmp >= max_value]
  is_maa_zero <- apply(maa,2,sum)==0
  maa[,is_maa_zero] <- maa_tmp[,is_maa_zero]
  maa[naa==0] <- 0
  maa
}

#' @export
#' @encoding UTF-8

get_wcatch <- function(res){
    if(class(res)=="future_new")  return(apply(res$wcaa,c(2,3),sum))
    if(class(res)=="vpa" || "tune" %in% names(res$input)){
        if(!is.null(res$input$dat$waa.catch)) return(colSums(res$input$dat$caa * res$input$dat$waa.catch,na.rm=T))
        if(is.null(res$input$dat$waa.catch)) return(colSums(res$input$dat$caa * res$input$dat$waa,na.rm=T))
    }
}

#' @export
get_ssb <- function(res){
    if(class(res)=="future_new")  return(res$SR_mat[,,"ssb"])
    if(class(res)=="vpa" || "tune" %in% names(res$input))  return(colSums(res$ssb, na.rm=TRUE))
}

#' @export
get_U <- function(res){
    if(class(res)=="future_new")  return(res$HCR_realized[,,"wcatch"]/apply(res$naa * res$waa_catch,c(2,3),sum))
    if(class(res)=="vpa" || "tune" %in% names(res$input)){
        wcatch <- get_wcatch(res)
        biomass <- colSums(res$naa * res$input$dat$waa,na.rm=T)
        return(wcatch/biomass)
    }
}



#' @export
#' @encoding UTF-8
#'

calc_forward <- function(naa,M,faa,t, plus_age, plus_group = TRUE){
  if(length(dim(naa))==3){
    for(iage in 1:(plus_age-1)) {
      naa[iage+1,t+1,] <- naa[iage,t,]*exp(-M[iage,t,]-faa[iage,t,])
    }
    if(plus_group == TRUE) naa[plus_age,t+1,] <- naa[plus_age,t+1,] + naa[plus_age,t,]*exp(-M[plus_age,t,]-faa[plus_age,t,])
  }
  else{
    for(iage in 1:(plus_age-1)) {
      naa[iage+1,t+1] <- naa[iage,t]*exp(-M[iage,t]-faa[iage,t])
    }
    if(plus_group == TRUE) naa[plus_age,t+1] <- naa[plus_age,t+1] + naa[plus_age,t]*exp(-M[plus_age,t]-faa[plus_age,t])
  }
  return(naa)
}

set_upper_limit_catch <- function(catch_previous_year, catch_current_year, upper_limit){
  upper_catch <- catch_previous_year * upper_limit
  is_over_upper_catch  <- catch_current_year > upper_catch
  catch_current_year[is_over_upper_catch] <- upper_catch[is_over_upper_catch]
  catch_current_year
}

set_lower_limit_catch <- function(catch_previous_year, catch_current_year, lower_limit){
  lower_catch <- catch_previous_year * lower_limit
  is_under_lower_catch  <- catch_current_year < lower_catch
  catch_current_year[is_under_lower_catch] <- lower_catch[is_under_lower_catch]
  catch_current_year
}


#'
#' future_vpaを使ってMSY管理基準値などを計算するwrapper関数
#'
#' @param data_future make_future_dataの返り値
#' @param candidate_PGY PGYの計算候補
#' @param candidate_B0 b0の計算候補
#' @param candidate_Babs Babsの計算候補
#' @param trace_multi （このベクトル）×（管理基準値として計算されるF）の平行状態の資源状態などを計算する。HCR_catchのプロットするときに、この値をもっと細かく設定することが必要になってくるかも
#' @param multi_upper_PGY PGYを計算するときの探索範囲の上限。デフォルトは10だが、うまくいかない場合にはこの数字を少し変えてみるとよいかも。
#'
#' @export
#' @encoding UTF-8
#'

est_MSYRP <- function(data_future, ncore=0, optim_method="R", compile_tmb=FALSE, candidate_PGY=c(0.1,0.6),
                      only_lowerPGY="lower", candidate_B0=-1, candidate_Babs=-1, candidate_Fbase=-1,
                      calc_yieldcurve=TRUE,
                      trace_multi=c(0.9,0.925,0.95,0.975,1.025,1.05,1.075),
                      select_Btarget=0, select_Blimit=0, select_Bban=0,
                      multi_upper_PGY=10){

  res_vpa_MSY <- data_future$input$res_vpa
  res_SR_MSY <-  data_future$input$res_SR
  # F=0からssbがゼロになるまでFを順次大きくしたtraceを実行する
  trace.multi <- unique(sort(c(0.001,seq(from=0,to=2,by=0.1),10,100)))
  trace_pre <- frasyr::trace_future(data_future$data, trace.multi=trace.multi, ncore=ncore)
  B0stat <- trace_pre %>% dplyr::filter(fmulti==0) %>% mutate(RP_name="B0")
  trace.multi2 <- unique(range(trace.multi[trace_pre$ssb.mean>0.001]))

  f_range <- range(trace.multi[which.max(trace_pre$catch.mean)+c(-1,1)])
  f_range[f_range==0] <- 0.0001

  # 以降、初期値はそれを参考に決める
  res_future_MSY <- future_vpa(tmb_data = data_future$data,
                               optim_method=optim_method,
                               multi_init=mean(f_range),
                               multi_lower=f_range[1],
                               multi_upper=f_range[2],
                               compile=compile_tmb)

    MSYstat <- res_future_MSY %>% get.stat(use_new_output=TRUE) %>%
      mutate(RP_name="MSY")

      # 他管理基準値を推定するためのオブジェクトを作っておく
    obj_mat <- NULL
    if(candidate_PGY[1]>0){

        obj_mat <- bind_rows(obj_mat,
                             tibble(RP_name    = str_c("PGY",candidate_PGY,"lower",sep="_"),
                                    obj_value  = candidate_PGY * MSYstat$catch.mean,
                                    optim_method=optim_method,
                                    multi_init = res_future_MSY$multi*1.2,
                                    multi_lower= res_future_MSY$multi,
                                    multi_upper= multi_upper_PGY,
                                    objective="PGY"
                                    ))

        if(only_lowerPGY=="both"){
            obj_mat2 <- tibble(RP_name    = str_c("PGY",candidate_PGY,"upper",sep="_"),
                               obj_value  = candidate_PGY * MSYstat$catch.mean,
                               optim_method=optim_method,
                               multi_init = res_future_MSY$multi*0.5,
                               multi_upper= res_future_MSY$multi,
                               multi_lower= 0.001,
                               objective="PGY")
            obj_mat <- bind_rows(obj_mat, obj_mat2)
        }
    }

    if(candidate_B0[1]>0){
        fssb.range <- trace.multi[trace_pre$ssb.mean>0.1]
        obj_mat <- bind_rows(obj_mat,
                             tibble(RP_name    = str_c("B0-",candidate_B0*100,"%"),
                                    obj_value  = candidate_B0 * B0stat$ssb.mean,
                                    optim_method=optim_method,
                                    multi_init = mean(fssb.range),
                                    multi_upper= max (fssb.range),
                                    multi_lower= 0.001,
                                    objective="SSB"
                                    ))
    }

    if(candidate_Babs[1]>0){
        fssb.range <- trace.multi[trace_pre$ssb.mean>0.1]
        obj_mat <- bind_rows(obj_mat,
                             tibble(RP_name    = str_c("Ben-",candidate_Babs,""),
                                    obj_value  = candidate_Babs,
                                    optim_method=optim_method,
                                    multi_init = mean(fssb.range),
                                    multi_upper= max (fssb.range),
                                    multi_lower= 0.001,
                                    objective="SSB"
                                    ))
    }

    if(candidate_Fbase[1]>0){
        obj_mat <- bind_rows(obj_mat,
                             tibble(RP_name    = str_c("Fbase-",round(candidate_Fbase,3),""),
                                    obj_value  = 0,
                                    multi_init = candidate_Fbase,
                                    multi_upper= candidate_Fbase,
                                    multi_lower= candidate_Fbase,
                                    objective="MSY",
                                    optim_method="none"
                                    ))
    }

    # obj_matをまとめてmapで回す
    other_RP_stat <- NULL
    if(!is.null(obj_mat)){

        other_RP_stat <-
            purrr::map_dfr(1:nrow(obj_mat),
                           function(x){
                               res <- future_vpa(tmb_data     = data_future$data,
                                                 optim_method = obj_mat$optim_method[x],
                                                 multi_init   = obj_mat$multi_init[x],
                                                 multi_lower  = obj_mat$multi_lower[x],
                                                 multi_upper  = obj_mat$multi_upper[x],
                                                 compile      = FALSE,
                                                 objective    = obj_mat$objective[x],
                                                 obj_value    = obj_mat$obj_value[x],
                                                 obj_stat     = "mean") %>%
                                   get.stat(use_new_output=TRUE)})

        other_RP_stat <- bind_cols(other_RP_stat, select(obj_mat, RP_name))
        print(bind_cols(obj_mat[,1:2], select(other_RP_stat,catch.mean, ssb.mean)))
    }

    all.stat <- bind_rows(MSYstat, B0stat, other_RP_stat)
    sum.stat <- get_summary_stat(all.stat)

    if(calc_yieldcurve==TRUE){
        # update trace
        trace.multi2 <- c(sum.stat$sumvalue$"Fref/Fcur",trace.multi2)
        trace.multi2 <- trace.multi2[trace.multi2>0] %>%
            purrr::map(function(x) x * trace_multi) %>%
            unlist() %>% sort() %>% unique()
        diff.trace <- diff(log(trace.multi2))
        trace.multi2[which(mean(diff.trace)<diff.trace)]
        trace.multi2 <- c(trace.multi2,
                          trace.multi2[which(mean(diff.trace)<diff.trace)] +
                          diff.trace[which(mean(diff.trace)<diff.trace)]/2) %>%
            sort()
        trace_pre2 <- trace_future(data_future$data,
                                   trace.multi=trace.multi2, ncore=ncore)
        trace_pre <- bind_rows(trace_pre,trace_pre2)
        trace_pre <- trace_pre[!duplicated(trace_pre$ssb.mean),]
    }

    res_MSY <- lst(all.stat=all.stat, summary=sum.stat$sumvalue,
                   Fvector=sum.stat$Fvector,input=res_future_MSY$input,
                   input_data=data_future$input,
                   trace=trace_pre, res_vpa=res_vpa_MSY, res_SR=res_SR_MSY)

    res_MSY$summary$perSPR <-
        purrr::map_dbl(1:dim(res_MSY$Fvector)[1],
                   function(x)
                       calc_perspr(fout=format_to_old_future(res_future_MSY),
                                   res_vpa=res_vpa_MSY,Fvector=res_MSY$Fvector[x,]))

    # define RP.definition for Btarget
    if(select_Btarget!=0){
        if(select_Btarget<0){
            print(select(res_MSY$summary,-Catch.CV))
            select_Btarget <- readline("Enter row number to be Btarget: ")
            select_Btarget <- as.integer(select_Btarget)
        }
        res_MSY$summary$RP.definition[1] <- NA
        res_MSY$summary$RP.definition[select_Btarget] <- "Btarget0"
    }
    # define RP.definition for Blimit
    if(select_Blimit!=0){
        if(select_Blimit<0){
            print(select(res_MSY$summary,-Catch.CV))
            select_Blimit <- readline("Enter row number to be Blimit: ")
            select_Blimit <- as.integer(select_Blimit)
        }
        res_MSY$summary$RP.definition[which(res_MSY$summary$RP.definition=="Blimit0")] <- NA
        res_MSY$summary$RP.definition[select_Blimit] <- "Blimit0"
        
    }
    # define RP.definition for Bban
    if(select_Bban!=0){
        if(select_Bban<0){
            print(select(res_MSY$summary,-Catch.CV,-RP.definition))
            select_Bban <- readline("Enter row number to be Bban: ")
            select_Bban <- as.integer(select_Bban)
        }
        res_MSY$summary$RP.definition[which(res_MSY$summary$RP.definition=="Bban0")] <- NA
        res_MSY$summary$RP.definition[select_Bban] <- "Bban0"
    }

  res_MSY$res_future_MSY <- res_future_MSY
  res_MSY$data_future_MSY <- data_future
  return(res_MSY)

}

#' @export
#' 

est_MSYRP_proxy <- function(data_future,
                            Fmsy_proxy_candidate=c("Fmax","F0.1","F%spr"),
                            msy_SPR_candidate=c(30,40),
                            Blimit_candidate=c("Bmin","10%B0","Babs"),
                            Babs_value = 1000,
                            Bban_candidate=c("0","0.1Blimit","0.2Blimit"), # not calculate all statistics
                            select_Btarget="F%spr30",
                            select_Blimit="Bmin",
                            select_Bban="0",
                            F.range = seq(from=0,to=2,length=101)
                            ){

  # waa_fun, maa_fun=TRUEの場合にはFbaseの管理基準値は計算できない
  assertthat::assert_that(data_future$input$waa_fun==FALSE)
  assertthat::assert_that(data_future$input$maa_fun==FALSE)
  assertthat::assert_that(data_future$input$waa_catch_fun==FALSE)  

  # candidateとして受け付けるもの以外ははねる
  assertthat::assert_that(all(Fmsy_proxy_candidate %in% c("Fmax","F0.1","F%spr")))
  assertthat::assert_that(all(Blimit_candidate     %in% c("Bmin","10%B0","Babs")))
  assertthat::assert_that(all(Bban_candidate       %in% c("0","0.1Blimit","0.2Blimit")))
  
  Babs_vector <- numeric()
  if("Bmin" %in% Blimit_candidate) Babs_vector <- c(Babs_vector, "Bmin" = min(colSums(data_future$input$res_vpa$ssb, na.rm=TRUE)))
  if("Babs" %in% Blimit_candidate){
    names(Babs_value) <- str_c("Babs",1:length(Babs_value))
    Babs_vector <- c(Babs_vector, Babs_value)
  }
  if("10%B0" %in% Blimit_candidate){
    candidate_B0 <- c("10%B0"=0.1)
  }
  else{
    candidate_B0 <- -1
  }

  lastyear <- dim(data_future$data$waa_mat[,,1])[[2]]
  waa <- data_future$data$waa_mat[,lastyear,1]
  tmp <- waa!=0
  waa <- waa[tmp]
  maa <- data_future$data$maa_mat[tmp,lastyear,1]
  M   <- data_future$data$M_mat  [tmp,lastyear,1]
  waa.catch <- data_future$data$waa_catch_mat[tmp,lastyear,1]
  futureF <- data_future$data$faa_mat[tmp,lastyear,1]  

  # Fmsy proxyを計算
  age_name <- as.numeric(rownames(data_future$input$res_vpa$naa))
  age_name <- age_name[tmp]

  res_refF <- ref.F(res      = NULL,
                    Fcurrent = futureF,
                    waa      = waa,
                    maa      = maa,
                    M        = M,
                    waa.catch  = waa.catch,
                    rps.vector = NULL, 
                    Pope     = data_future$input$Pope,
                    min.age  = min(age_name),
                    max.age  = ifelse(data_future$input$res_vpa$input$plus.group==FALSE, max(age_name), Inf),
                    pSPR     = msy_SPR_candidate,plot=FALSE,
                    F.range = F.range)

  f_proxy_vector <- numeric()

  if("Fmax" %in% Fmsy_proxy_candidate) f_proxy_vector <- c(f_proxy_vector, Fmax=res_refF$summary$Fmax[3])
  if("F0.1" %in% Fmsy_proxy_candidate) f_proxy_vector <- c(f_proxy_vector, F0.1=res_refF$summary$F0.1[3])  
  if("F%spr" %in% Fmsy_proxy_candidate){
    if(length(msy_SPR_candidate)>1) label <- str_c("FpSPR.",msy_SPR_candidate,".SPR")
    if(length(msy_SPR_candidate)==1) label <- str_c("X",msy_SPR_candidate,".SPR")
    tmp <- res_refF$summary[3,label]
    names(tmp) <- str_c("F%spr", msy_SPR_candidate)
    f_proxy_vector <- c(f_proxy_vector, tmp)  
  }    

  f_proxy_vector <- f_proxy_vector %>% unlist()
  
  #一気にest_MSYRPで計算
  res_MSY <- est_MSYRP(data_future=data_future, ncore=0, optim_method="R",
                       candidate_PGY=-1,
                       candidate_Babs=Babs_vector, candidate_B0=candidate_B0,
                       candidate_Fbase=f_proxy_vector,
                       calc_yieldcurve=FALSE)

  allRP <- c(candidate_B0, Babs_vector,f_proxy_vector)
  allRP <- allRP[allRP>0]

  # dynamic Bmsyを計算
  # res_MSY$dynamic_Bmsy <- calc_dmsy()
  res_MSY$summary$RP_name[-1:-2] <- names(allRP)

  # selct reference points
  assertthat::assert_that(select_Btarget %in% res_MSY$summary$RP_name)
  assertthat::assert_that(select_Blimit  %in% res_MSY$summary$RP_name)

  res_MSY$summary$RP.definition[] <- NA
  res_MSY$summary$RP.definition[res_MSY$summary$RP_name==select_Btarget] <- "Btarget0"
  res_MSY$summary$RP.definition[res_MSY$summary$RP_name==select_Blimit]  <- "Blimit0"

  Blimit <- derive_RP_value(res_MSY$summary, "Blimit0")$SSB
  # Bbanは数値のみ入れる
  for(i in 1:length(Bban_candidate)){
    bban_summary <- tibble(RP_name = Bban_candidate[i],
                           SSB     = case_when(
                             Bban_candidate[i]=="0.1Blimit" ~ 0.1*Blimit,
                             Bban_candidate[i]=="0.2Blimit" ~ 0.2*Blimit,
                             Bban_candidate[i]=="0"          ~ 0)) %>%
      mutate(SSB2SSB0=SSB/res_MSY$summary$SSB[2])
    res_MSY$summary <-bind_rows(res_MSY$summary,bban_summary)
  }
  assertthat::assert_that(select_Bban    %in% res_MSY$summary$RP_name)  
  res_MSY$summary$RP.definition[res_MSY$summary$RP_name==select_Bban]    <- "Bban0"

  res_MSY$res_refF <- res_refF

  return(res_MSY)
}
