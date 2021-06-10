
#' VPAの感度分析のための関数
#'
#' @param res VPAの結果のオブジェクト
#' @param what_replace 感度分析の対象。"M"（自然死亡係数）, "maa"（年齢別成熟率）, "waa"（年齢別体重）, "waa.catch"（漁獲物の年齢別体重）, "alpha"（+グループA歳とA-1歳の漁獲係数の比）, "tuning"（最終年の漁獲係数の扱い）, "lambda"（リッジVPAのペナルティ）,"est.method"（推定方法）, "b"（hyperstability/depletion）から選べる。
#' @param value 感度分析の対象の値。ベースケースの結果に一定割合をかけたい場合は\code{numeric}型（"M", "waa", "waa.catch", "alpha"で使用可）、ベースケースの結果とは異なる表を入れたい場合は\code{list}型に\code{matrix}型を（"maa", "waa", "waa.catch"で使用可）入れてください。"tuning"を選んだ場合、資源量指数を用いていないVPA結果なら、tf.yearの期間を入れること。
#' @param what_plot 作図したい項目を選べる。引数を与えない場合（\code{NULL}）、全て（SSB, biomass, U, catch, Recruitment, fish_number, fishing_mortality, weight, maturity, catch_number）をプロットする。
#' @param ncol 作図の列数。標準で5列なので、\code{what_plot}の数が5以下の場合は適宜変えた方がよい。
#'
#' @return 返ってくる値:
#'     \code{result} 感度分析の結果が\code{list}型式で得られる。
#'     \code{graph} 図が得られる。
#'
#' @details "M"（自然死亡係数）:
#' \code{value}にはベースケースの結果にかけたい割合を\code{numeric}型で記述のこと。
#'
#' @details "maa"（年齢別成熟率）:
#' \code{value}にはベースケースと同じ行列数の\code{matrix}型で記述のこと。
#'
#' @details "waa"（年齢別体重）, "waa.catch"（漁獲物の年齢別体重）:
#' \code{value}にはベースケースの結果にかけたい割合を\code{numeric}型、またはベースケースと同じ行列数の\code{matrix}型で記述のこと。
#'
#' @details "alpha"（+グループA歳とA-1歳の漁獲係数の比）:
#' \code{value}には仮定したい比率を\code{numeric}型で記述のこと。
#'
#' @details "tuning"（最終年の漁獲係数の扱い）:
#' resの計算方法に応じて\code{value}に入れる値が変わる。
#' チューニングしていないVPAの場合、感度分析も同様にチューニングしないVPAを異なるターミナルFの仮定で行う。\code{value}は仮定したい\code{tf.year}を\code{list}型で記述のこと。
#' 全F推定の場合、感度分析では選択率更新法を行う。\code{value}は仮定したい\code{tf.year}を\code{list}型で記述のこと。
#' 選択率更新法の場合、感度分析では全F推定を行う。\code{value}は入力不要である。
#' リッジVPAの場合、感度分析では異なるlambdaでの解析を行う。\code{value}は仮定したい\code{tf.year}を\code{list}型で記述のこと。
#'
#' @details "lambda"（リッジVPAのペナルティ）: \code{value}は与えたいlambdaを\code{numeric}型で記述のこと。
#'
#'
#' @details "est.method"（推定方法）: \code{value}は入力不要である。
#'
#' @details "b"(hyperstability/depletion):
#' bを推定した場合、bを推定しない(b=1)の感度分析を行う。
#' bを1以外の値で固定した、あるいはbを考慮しなかった場合、\code{value="b.est"}とするとbの推定を、仮定したいbを\code{list}で与えるとbを固定して考慮出来る。
#'
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' vpa計算について:  \code{\link{vpa}}
#' 作図について: \code{\link{plot_vpa}}
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#'
#' @encoding UTF-8
#'
#' @export

# helpは↑のように、"#`"のあとに一定の書式に従って記述するのでこちらも書いてみてください。少しだけ試しに書いてみました。@export, @encoding UTF-8とかはそのままにしておいてください

# author: Kohei Hamabe

do_sensitivity_vpa <- function(res,
                               what_replace,
                               value,
                               what_plot = NULL,
                               ncol = 5,
                               plot_year = NULL,
                               scale_value = NULL
){

  res_vpa.s <- list() ; lab.tmp <- numeric()
  res$input$plot <- FALSE # 背後で沢山plotするのを防ぐ

  if(what_replace == "M"){
    if(!class(value) == "numeric")stop(paste('Set values which shape is only "numeric"'))
    for(i in 1:length(value)){
      input0 <- res$input
      input0$dat$M <- input0$dat$M *value[i]
      res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp[i] <- paste("Sensitivity M= x", value[i], sep = "")
    } # for

    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "waa") {

    if(class(value) == "numeric") {
      ## 任意の比率(ex. 1.5倍)で入れた場合
      for(i in 1:length(value)){
        input0 <- res$input
        input0$dat$waa <- input0$dat$waa *value[i]
        res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
        lab.tmp[i] <- paste("Sensitivity waa= x", value[i], sep = "")
      } # for

    } else if (class(value) == "list"){
      ## 任意のwaaをリスト形式で入れた場合
      for(i in 1:length(value)){
        if(! dim(value[[i]])[2] == dim(res$input$dat$waa)[2]){
          stop(paste("Dimentional of value data [[", i, "]] were different...", sep = ""))
        } else {
          next
        }
      } # for(stop())

      colnames.tmp <- colnames(res$input$dat$waa) ; rownames.tmp <- rownames(res$input$dat$waa)

      for(i in 1:length(value)){
        input0 <- res$input
        input0$dat$waa <- value[[i]]
        colnames(input0$dat$waa) <- colnames.tmp ; rownames(input0$dat$waa) <- rownames.tmp
        res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
        lab.tmp[i] <- paste("Sensitivity waa case:", i, sep = "")
      } # for

    } else {
      stop(paste('You have to input values which shape is "numeric" or "list"', sep = ""))
    } # if(データの型)

    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "waa.catch") {
    if(is.null(res$input$dat$waa.catch)) stop(paste("waa.catch is NULL !!"))
    if(class(value) == "numeric") {
      for(i in 1:length(value)){
        input0 <- res$input
        input0$dat$waa.catch <- input0$dat$waa.catch *value[i]
        res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
        lab.tmp[i] <- paste("Sensitivity waa.catch= x", value[i], sep = "")
      } # for

    } else if (class(value) == "list"){
      ## 任意のwaaをリスト形式で入れた場合
      for(i in 1:length(value)){
        if(! dim(value[[i]])[2] == dim(res$input$dat$waa)[2]){
          stop(paste("Dimentional of value data [[", i, "]] were different...", sep = ""))
        } else {
          next
        }
      } # for(stop())

      colnames.tmp <- colnames(res$input$dat$waa.catch) ; rownames.tmp <- rownames(res$input$dat$waa.catch)

      for(i in 1:length(value)){
        input0 <- res$input
        input0$dat$waa.catch <- value[[i]]
        colnames(input0$dat$waa.catch) <- colnames.tmp ; rownames(input0$dat$waa.catch) <- rownames.tmp
        res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
        lab.tmp[i] <- paste("Sensitivity waa.catch case:", i, sep = "")
      } # for

    } else {
      stop(paste('You have to input values which shape is "numeric" or "list"', sep = ""))
    } # if(データの型)

    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "maa") {

    if(! class(value) == "list"){
      stop(paste('For "maa", the shape of values were only "list"', sep = ""))
    }

    for(i in 1:length(value)){
      if (! dim(value[[i]])[2] == dim(res$input$dat$maa)[2]){
        stop(paste("Dimentional of value data [[", i, "]] were different...", sep = ""))
      } else {
        next
      }
    } # for(stop())

    colnames.tmp <- colnames(res$input$dat$maa) ; rownames.tmp <- rownames(res$input$dat$maa)

    for(i in 1:length(value)){
      input0 <- res$input
      input0$dat$maa <- value[[i]]
      colnames(input0$dat$maa) <- colnames.tmp
      rownames(input0$dat$maa) <- rownames.tmp
      res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp[i] <- paste("Sensitivity maa case:", i, sep = "")
    } # for

    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "alpha") {
    for(i in 1:length(value)){
      input0 <- res$input
      input0$alpha <- value[i]
      res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp[i] <- paste("Sensitivity alpha=", value[i], sep = "")
    } # for

    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "tuning") {
    if (res$input$tune == FALSE){
      # non-tuned VPAを行う
      if(!class(value) == "list")stop(paste('For "tuning", the shape of values were only "list" in non-tuned VPA', sep = ""))
      for(i in 1:length(value)){
        input0 <- res$input
        input0$tune <- FALSE
        input0$term.F <- "max"
        input0$sel.update <- FALSE
        input0$tf.year <- value[[i]]
        res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
        lab.tmp[i] <- paste("tf.year = ", value[[i]][1], ":", value[[i]][length(value[[i]])], sep = "")
      }
    } else if(res$input$sel.update == TRUE){  # resは選択率更新
      # 全F推定を試みる
      input0 <- res$input
      input0$sel.update <- FALSE
      input0$term.F <- "all"
      res_vpa.s[[1]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp <- "ALL_F.est"
      #   } else if(!res$input$lambda == 0){  # resはRidge VPA
      #       # lambdaを変えて感度分析
      #      if(!class(value) == "numeric")stop(paste('For "tuning", the shape of values were only "numeric" in ridge VPA', sep = ""))
      #      for(i in 1:length(value)){
      #        input0 <- res$input
      #        input0$lambda <- value[i]
      #        res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
      #        lab.tmp[i] <- paste("lambda = ", value[i], sep = "")
      #      }
    } else  {
      # 選択率更新法を行う
      if(!class(value) == "list")stop(paste('Set tf.years in the values, and the shape was only "list"'))
      for(i in 1:length(value)){
        input0 <- res$input
        input0$sel.update <- TRUE
        input0$term.F <- "max"
        input0$tf.year <- value[[i]]
        res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
        lab.tmp[i] <- paste("sel.update, tf.year = ", value[[i]][1], ":", value[[i]][length(value[[i]])], sep = "")
      }
    }
    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "lambda") {
    if(!class(value) == "numeric")stop(paste('For "lambda", the shape of values were only "numeric" in ridge VPA', sep = ""))
    for(i in 1:length(value)){
      input0 <- res$input
      input0$lambda <- value[i]
      res_vpa.s[[i]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp[i] <- paste("lambda = ", value[i], sep = "")
    }

    #-------------------------------------------------------------------------------------#

  } else if(what_replace == "b") {
    if(res$input$b.est == TRUE){
      input0 <- res$input
      input0$b.est <- FALSE
      res_vpa.s[[1]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp <- "b = 1"
    } else if(class(value) == "numeric") {
      if(!length(value) == length(res$b)) stop(paste0("Length of b was different !!"))
      input0 <- res$input
      input0$b.fix <- value
      res_vpa.s[[1]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp <- paste("b.fix", sep = "")
    } else if(class(value) == "list") {
      for(j in 1:length(value)){
        if(!length(value[[j]]) == length(res$b)) stop(paste0("Length of b was different !!"))
        input0 <- res$input
        input0$b.fix <- value[[j]]
        res_vpa.s[[j]] <- safe_call(vpa, input0, force=TRUE)
        lab.tmp[j] <- paste("b.fix: pattern", j, sep = "")
      }
    } else {
      if(!value == "b.est")stop(paste('Input the "b.est" or some values(numeric) in "value="!!', sep = ""))
      input0 <- res$input
      input0$b.fix <- NULL
      input0$b.est <- TRUE
      res_vpa.s[[1]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp <- "b.est: TRUE"
    }

    #-------------------------------------------------------------------------------------#

  } else if(what_replace == "est.method") {
    if(res$input$est.method == "ls"){
      input0 <- res$input
      input0$est.method <- "ml"
      #input0$abund <- value    # $abundで資源量指数の対応するものを明記。デフォだと長さがずれてしまう#
      #input0$sigma.const<- input0$sigma.constraint <- 1:length(input0$abund)
      # valueをlistにして、分散の傾斜をオプションでいじれてもいいかも
      res_vpa.s[[1]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp <- "est.method: ml"
    } else {
      input0 <- res$input
      input0$est.method <- "ls"
      res_vpa.s[[1]] <- safe_call(vpa, input0, force=TRUE)
      lab.tmp <- "est.method: ls"
    }

    #-------------------------------------------------------------------------------------#

  } else {
    stop(paste('You have to select the what_replace from "M", "waa", "waa.catch", "maa", "alpha", "tuning", "b" or "est.method"', sep = ""))
  }

  names(res_vpa.s) <- lab.tmp
  # 結果のラベルに名前を付ける
  if(!(is.null(what_plot)))what.plot <- factor(what_plot, levels = as.character(what_plot))
  g1 <- plot_vpa(c(list(Base=res), res_vpa.s), what.plot = what_plot,
                 ncol = ncol, plot_year = plot_year, scale_value = scale_value)

  return(list(result = res_vpa.s, graph = g1))
} # function(do_sensitivity_vpa)





#' VPAのレトロスペクティブ解析結果を自動生成する関数
#'
#' @param res VPAの結果のオブジェクト
#' @param n_retro レトロスペクティブ解析でさかのぼる年数
#' @param b_reest bをレトロスペクティブ解析中で再推定するか
#' @param what_plot 作図したい項目を選べる。\code{NULL}の場合、全て（SSB, biomass, U, catch, Recruitment, fish_number, fishing_mortality, weight, maturity, catch_number）をプロットする。
#' @param ncol 作図の列数。標準で5列なので、\code{what_plot}の数が5以下の場合は適宜変えた方がよい。
#' @param remove_maxAgeF Mohn's rhoを計算する際に最高齢のFを除くか（alphaを仮定して計算していることが多いから）
#' @param ssb_forecast Mohn's rhoを計算する際にSSBは1年後を計算するか(last.catch.zero=TRUEのときのみ有効)
#' @param res_step1 2段階法のレトロ解析をやる場合の1段階目の\code{vpa}オブジェクト
#' @return 返ってくる値:
#'     \code{result} 感度分析の結果が\code{list}型式で得られる。
#'     \code{mohn_rho}
#'     \code{graph} 図が得られる。
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' レトロスペクティブ解析について:  \code{\link{retro.est}}
#' 作図について: \code{\link{plot_vpa}}
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

do_retrospective_vpa <- function(res,
                                 n_retro = 5,
                                 b_reest = FALSE,
                                 what_plot = c("SSB", "biomass", "Recruitment",
                                               "fish_number", "fishing_mortality"),
                                 plot_year = NULL,
                                 ncol = 3,
                                 remove_maxAgeF = FALSE,
                                 ssb_forecast = FALSE,
                                 res_step1 = NULL,
                                 scale_value = NULL
){

  if(b_reest == TRUE && res$input$b.est == FALSE)message(paste('b was not estimated in your vpa model'))
  # vpa内でbの推定をしていないにもかかわらず、b_reestがtrueで入力された場合
  # 推定結果(bを推定している)は得られるが、メッセージを出す

  if (!is.null(res_step1)) { #二段階法の場合のレトロ    i
    retro_step_one <- retro.est(res_step1, n = n_retro)
    yy <- ifelse(res$input$last.catch.zero,2,1)
    sel_mat <- sapply(1:n_retro, function(i) rev(retro_step_one$Res[[i]]$saa)[,yy])
    res_retro <- retro.est(res, n = n_retro, b.fix = !b_reest, remove.maxAgeF=remove_maxAgeF, ssb.forecast=ssb_forecast,sel.mat=sel_mat)
  } else {
    res_retro <- retro.est(res, n = n_retro, b.fix = !b_reest, remove.maxAgeF=remove_maxAgeF, ssb.forecast=ssb_forecast)
  }
  dat_graph <- list()
  for(i in 1:n_retro) dat_graph[[i]] <- res_retro$Res[[i]]

  dat_graph <- c(list(res), dat_graph) # Base case(全データで解析)の追加（浜辺07/08）
  if(res$input$last.catch.zero){ # last.catch.zero=Tの場合、最終年のプロットはしない（Mohn's rhoとずれるから）（浜辺07/08）
    names(dat_graph) <- rev(colnames(res$ssb))[2:(n_retro+2)]
  } else {
    names(dat_graph) <- rev(colnames(res$ssb))[1:(n_retro+1)]  # 図にinputされる結果に名前をつける
  }

  # 図にMohn's rhoの重ね書き用rho data from 市野川さん
  rho_data <- tibble(index = names(res_retro$mohn), value = res_retro$mohn) %>%
    left_join(tibble(index = c("N", "B", "SSB", "R", "F"),
                     stat = c("fish_number", "biomass" ,"SSB" ,"Recruitment" ,"fishing_mortality"))) %>%
    #mutate(y=0, x=as.numeric(min(colnames(res_retro[[1]][[1]]$naa))))
    mutate(y = 0,
           x = if(is.null(plot_year))as.numeric(min(colnames(res_retro[[1]][[1]]$naa))) else plot_year[1])
  if(!length(what_plot) == 5) rho_data <- rho_data[match(what_plot, rho_data$stat),]
  # ここもしかすると長さが5でもwhat_plotのデフォルトと一致しないとエラー出るかも
  # そういった変数についてはレトロして見る需要は少ないのだろうけど

  g1 <- plot_vpa(dat_graph,
                 what.plot = factor(what_plot, levels = as.character(what_plot)),
                 ncol = ncol, plot_year = plot_year, scale_value = scale_value) +
    geom_label(data = rho_data,
               mapping = aes(x = x, y = y, label = str_c("rho=", round(value,2))),
               vjust="inward", hjust="inward")

  # bの結果出力用 from 市野川さん
  b_tmp <- purrr::map_dfc(res_retro$Res, function(x) as_tibble(x$b))
  colnames(b_tmp) <- rev(colnames(res$ssb))[1:n_retro]

  return(list(result = res_retro, mohn_rho = res_retro$mohn, graph = g1, b_res = b_tmp))
} # function(do_retrospective_vpa)





#' vpa結果の収束診断のための関数
#'
#' この関数ではJitter analysisを行います。
#' 乱数で生成した初期値を使ってvpa計算を再度行い、推定結果が同一になるかを確認することで、収束の良さを測る診断法です。
#'
#' @param res VPAの結果のオブジェクト
#' @param n_ite 繰り返し試行数。
#' @param sd_jitter 乱数の分散（logスケール）
#' @param what_plot 作図したい年齢を与える。\code{term.F = "all"}出ない限り、最高齢（+グループ）の初期値の結果のみを図示する。NULLで全てのterm.fを図示、ある特定の年齢のterm.Fを図示したい場合は\code{numeric}で与える必要がある。
#' @param TMB 繰り返し計算をTMBで行う
#'
#' @return 返ってくる値:
#'     \code{initial_value} 乱数で生成された初期値が得られる。
#'     \code{p_name} 与えた初期値（最終年の記載年のF）の名前が確認できる。
#'     \code{value} 各試行の推定値、尤度、収束診断、ヘッセ行列が得られる。
#'     \code{graph} 推定値と尤度について、それぞれ図が得られる。
#'
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

do_estcheck_vpa <- function(res, n_ite = 20, sd_jitter = 1, what_plot = NULL, TMB = FALSE){
  # resの中身の診断
  if(sum(diag(res$hessian))==sum(abs(diag(res$hessian)))){
    message(paste("In your VPA result, Hessian successfully having positive definite!!"))
  } else {
    stop(paste("In your VPA result, Hessian is not positive..."))
  }
  if(res$convergence == 1){
    message(paste("In your VPA result, Successful convergence!!"))
  } else {
    stop(paste("Your VPA result was not converged..."))
  }

  if(length(res$term.f) == 1){
    name_tmp <- "max"
  } else {
    name_tmp <- c(0:(length(res$term.f)-2), "max")
  }

  # 作図用引数をあらかじめチェック <= 解析後にエラーはかわいそう
  if(is.null(what_plot)){
    plot_name <- name_tmp
  } else if(what_plot == "max"){
    plot_name <- "max"
  } else if(what_plot == "numeric"){
    plot_name <- what_plot
  } else {
    stop(paste('what_plot was input age class in numeric, "max", or NULL'))
  }

  #init_list <- purrr::map(res$term.f,
  #                        function(x)exp(log(x) + rnorm(n_ite, 0, sd_jitter))
  #)
  init_list <- list()
  for (i in 1:length(res$term.f)) {
    init_list[[i]] <- seq(log(0.001), log(4), length = n_ite) %>%
      exp() %>% sample(n_ite)
  }
  value_tmp <- Finit <- Fest <- ite_tmp <- ll_tmp <- list()
  Hes_check <- Conv_check <- numeric()

  if(TMB == TRUE){
    #use_rvpa_tmb()
    res$input$TMB <- TRUE
  }

  for (i in 1:n_ite) {
    input0 <- res$input
    init_tmp <- numeric()
    for(j in 1:length(res$term.f)){
      init_tmp[j] <- init_list[[j]][i]
    }  # for(j)
    input0$p.init <- init_tmp
    tmp <- try(safe_call(vpa, input0, force=TRUE))
    if(class(tmp) == "try-error"){
      value_tmp[[i]] <- NA
      ite_tmp[[i]] <- rep(i, length(res$term.f))
      ll_tmp[[i]] <- NA#rep(res$logLik, length(res$logLik))
      Finit[[i]] <- init_tmp
      Fest[[i]] <- NA
    } else {
      value_tmp[[i]] <- list(p_est = tmp$term.f,
                             logLik = tmp$logLik,
                             covergence = tmp$convergence,
                             hessian = tmp$hessian,
                             gradient = tmp$gradient)
      ite_tmp[[i]] <- rep(i, length(res$term.f))
      #      ll_tmp[[i]] <- rep(res$logLik, length(res$term.f))
      ll_tmp[[i]] <- rep(tmp$logLik, length(res$term.f))
      Finit[[i]] <- init_tmp
      Fest[[i]] <- tmp$term.f
    }
    if(sum(diag(res$hessian)) == sum(abs(diag(res$hessian)))){
      Hes_check[i] <- 0
    } else {
      Hes_check[i] <- 1
    }
    Conv_check[i] <- res$convergence
    message(paste('Iteration',i,'has done ...', sep = " "))
  } # for(n_ite)

  names(init_list) <- name_tmp
  # 作図用データフレーム(tidyデータの作成)
  d_tmp <- data.frame(ite = rep(1:n_ite, each=length(res$term.f)),
                      age = rep(name_tmp, n_ite),
                      initial = unlist(Finit),
                      estimated = unlist(Fest),
                      likelihood = unlist(ll_tmp),
                      result_est = rep(res$term.f, n_ite),
                      result_lk = rep(res$logLik, length(unlist(ll_tmp)))
  )
  #est_res <- data.frame(age = name_tmp, estimated = res$term.f)

  yvalue <- max(res$term.f)*2
  g1 <- ggplot(data = d_tmp[d_tmp$age == plot_name,]) +
    geom_segment(aes(x=0, xend = 4, y = result_est, yend = result_est), color = "red", size = 1.3)+
    geom_point(aes(x = initial, y = estimated), size = 5) +
    ylim(c(0, yvalue)) +
    facet_wrap( ~ age) +
    xlab("initial value") +
    theme_SH(base_size = 14)
  g2 <- ggplot(data = d_tmp[d_tmp$age == "max",]) +
    geom_segment(aes(x=0, xend = 4, y = result_lk, yend = result_lk), color = "red", size = 1.3)+
    geom_point(aes(x = initial, y = likelihood), size = 5) +
    ylab("log Likelihood") + xlab("initial value of F of age Max") +
    theme_SH(base_size = 14)

  # Hessianの結果をメッセージで返す
  if(sum(Hes_check) == 0){
    message(paste("Hessian successfully having positive definite for all iterations !!"))
  } else {
    message(paste("All Hessian are not positive ..."))
  }
  # 収束結果をメッセージで返す
  if(sum(Conv_check) == length(Conv_check)){
    message(paste("Successful convergence for all iterations !!"))
  } else {
    lab_tmp <- which(!Conv_check == 1)
    message(paste('Iterations in ', lab_tmp, ' were not convergence ...'))
  }

  maxlike <- max(sapply(value_tmp, function(x) x$logLik))
  cat("Maximum likelihood in jitter analysis is: ",maxlike ,"\n")
  cat("Likelihood with estimated parameters is: ", res$logLik, "\n")

  return(list(initial_value = init_list, #初期値の乱数
              p_name = name_tmp, # 初期値の名前
              value = value_tmp, # 推定値と尤度のリスト
              graph = list(estimated = g1, likelihood = g2)
  ))
} # function(do_estcheck_vpa)





#' 残差プロットの関数
#'
#'
#' @param res VPAの結果のオブジェクト
#' @param index_name 作図時に各指標の名前
#'
#' @return 返ってくる値:
#'     \code{year_resid} 残差の時系列プロット
#'     \code{fitting_CPUE} 観測、予測CPUEのフィッティングの図
#'     \code{abund_CPUE} 資源量と資源量指数の当てはまりの図。
#'     \code{year_sd_resid} 標準化残差の時系列プロット
#'     \code{gg_data} 作図に用いたデータ（tidy形式）
#'
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

plot_residual_vpa <- function(res, index_name = NULL, plot_smooth = TRUE, plot_year = FALSE, plot_scale = FALSE, resid_CI=TRUE, plotAR=FALSE){
  if(is.numeric(res$input$use.index)){
    assertthat::assert_that(length(res$input$dat$index[,1]) >= length(res$input$use.index))
    used_index <- res$input$dat$index[res$input$use.index,]
  } else if(res$input$use.index == "all") {
    used_index <- res$input$dat$index
  } else {
    assertthat::assert_that(is.numeric(res$input$use.index)|res$input$use.index=="all")
  }
  # x軸の範囲
  if(is.numeric(plot_year)) xlim_year <- c(min(plot_year), max(plot_year)) else xlim_year <- c(min(as.numeric(colnames(res_vpa_estb$naa))), max(as.numeric(colnames(res_vpa_estb$naa))))

  d_tmp <- matrix(NA,
                  nrow = length(used_index[1,]),
                  ncol = length(used_index[,1])*8+4)
  d_tmp[,1] <- as.numeric(colnames(used_index))
  d_tmp[,2:(1+length(res$q))] <- as.numeric(t(used_index))
  d_tmp[,(2+length(res$q))] <- as.numeric(apply(res$naa, 2, sum))
  d_tmp[,(3+length(res$q))] <- as.numeric(apply(res$baa, 2, sum))
  d_tmp[,(4+length(res$q))] <- as.numeric(apply(res$ssb, 2, sum))


  q_tmp <- b_tmp <- sig_tmp <- numeric()
  name_tmp1 <- name_tmp2 <- name_tmp3 <- name_tmp4 <- name_tmp5 <- numeric()

  for(i in 1:length(res$q)){
    if(length(res$input$min.age)==1) min_age_tmp <- res$input$min.age[1] else min_age_tmp <- res$input$min.age[i]
    if(length(res$input$min.age)==1) max_age_tmp <- res$input$max.age[1] else max_age_tmp <- res$input$max.age[i]

    resid_tmp <- log(d_tmp[,i+1]) - log(res$pred.index[i,]) # 対数残差
    sd_resid_tmp <- resid_tmp/sd(resid_tmp, na.rm = TRUE) # 対数残差の標準化残差

    #abund.extractor関数で書き換え #catch.prop引数は不要か
    d_tmp[,(i+length(res$q)*1+4)] <- abund.extractor(abund = res$input$abund[i], naa = res$naa, faa = res$faa,
                                                     dat = res$input$dat,
                                                     min.age = res$input$min.age[i], max.age = res$input$max.age[i],
                                                     link = res$input$link, base = res$input$base, af = res$input$af,
                                                     sel.def = res$input$sel.def, p.m=res$input$p.m,
                                                     omega=res$input$omega, scale=1) #res$input$scale)
                                                    #res$ssbはスケーリングしていない結果が出ている(2021/06/09KoHMB)

    d_tmp[,(i+length(res$q)*2+4)] <- res$pred.index[i,] # q*N^B計算結果
    d_tmp[,(i+length(res$q)*3+4)] <- resid_tmp
    d_tmp[,(i+length(res$q)*4+4)] <- sd_resid_tmp
    d_tmp[,(i+length(res$q)*5+4)] <- rep(res$q[i], length(d_tmp[,1]))
    d_tmp[,(i+length(res$q)*6+4)] <- if(res$input$est.method=="ml"){
      rep(res$sigma[i], length(d_tmp[,1])) # ML
    } else {
      rep(res$sigma[1], length(d_tmp[,1])) # LS
    }
    d_tmp[,(i+length(res$q)*7+4)] <- rep(res$b[i], length(d_tmp[,1]))
    if(i >= 10){
      name_tmp1[i] <- paste0("obs_Index",i)
      name_tmp2[i] <- paste0("predabund_Index",i)
      name_tmp3[i] <- paste0("pred_Index",i)
      name_tmp4[i] <- paste0("resid_Index",i)
      name_tmp5[i] <- paste0("sd.resid_Index",i)
      q_tmp[i] <- paste0("q_Index",i)
      b_tmp[i] <- paste0("b_Index",i)
      sig_tmp[i] <- paste0("sigma_Index",i)
    } else {
      name_tmp1[i] <- paste0("obs_Index0",i)
      name_tmp2[i] <- paste0("predabund_Index0",i)
      name_tmp3[i] <- paste0("pred_Index0",i)
      name_tmp4[i] <- paste0("resid_Index0",i)
      name_tmp5[i] <- paste0("sd.resid_Index0",i)
      q_tmp[i] <- paste0("q_Index0",i)
      b_tmp[i] <- paste0("b_Index0",i)
      sig_tmp[i] <- paste0("sigma_Index0",i)
    }
  }
  d_tmp <- as.data.frame(d_tmp)
  names(d_tmp) <- c("year", name_tmp1, "abundance",  "biomass", "ssb", name_tmp2,
                    name_tmp3, name_tmp4, name_tmp5, q_tmp, sig_tmp, b_tmp)

  d_tidy <- tidyr::pivot_longer(d_tmp, col = c(-year, -abundance, -biomass, -ssb),
                                names_to = c(".value", "Index_Label"),
                                names_sep = "_",
                                values_drop_na = TRUE
  )

  if(!is.null(index_name)){
    if(!length(index_name) == length(res$q)) stop(paste0("Length of index_name was different to the number of indices"))
    d_tidy$Index_Label <- rep(index_name, nrow(d_tmp))
  }

  # 自己相関係数推定
  rho.numeric <- signif.numeric <- numeric()
  for(i in 1:length(unique(d_tidy$Index_Label))){
    calc.data <- d_tidy[d_tidy$Index_Label==unique(d_tidy$Index_Label)[i],]
    calc.data <- calc.data[!is.na(calc.data$resid),]
    #ar.res <-  ar(calc.data$resid,aic=FALSE,order.max=1,demean=FALSE,method = "mle")
    acf.res <- acf(calc.data$resid, type = "correlation", plot = FALSE, demean = FALSE)
    rho.numeric[i] <- acf.res$acf[,,1][2]
    signif.numeric[i] <- if((qnorm(0.025)/sqrt(acf.res$n.used) > acf.res$acf[,,1][2])|
                            (qnorm(0.975)/sqrt(acf.res$n.used) < acf.res$acf[,,1][2])) "*" else ""
  }

  # 残差プロットのy軸の上下限
  thred_resid <- range(d_tidy$resid,na.rm = TRUE) %>% abs() %>% max()
  if(plot_scale) thred_y <- c(NA,NA) else thred_y <- c(-thred_resid,thred_resid)
  y.posi_rho_data <- which(!range(d_tidy$resid,na.rm = TRUE)==thred_y)[1]
  thred_sd.resid <- range(d_tidy$sd.resid,na.rm = TRUE) %>% abs() %>% max()
  if(plot_scale) sd.thred_y <- c(NA,NA) else sd.thred_y <- c(-thred_sd.resid,thred_sd.resid)
  sd.y.posi_rho_data <- which(!range(d_tidy$sd.resid,na.rm = TRUE)==sd.thred_y)[1]

  # 残差プロットに追加する観測誤差と自己相関係数のtidy data
  rho_data <- tibble(Index_Label = unique(d_tidy$Index_Label), sigma = res$sigma, ar1 = rho.numeric, signif = signif.numeric) %>%
    mutate(y = thred_y[y.posi_rho_data], x = xlim_year[1], y.sd = sd.thred_y[sd.y.posi_rho_data])

  if(isTRUE(resid_CI)){
    g1 <- ggplot(d_tidy) +
      geom_ribbon(aes(x = year, ymin = -qnorm(0.025)*sigma, ymax = qnorm(0.025)*sigma), alpha=0.05)+
      geom_ribbon(aes(x = year, ymin = -qnorm(0.1)*sigma, ymax = qnorm(0.1)*sigma), alpha=0.1)+
      geom_point(aes(x=year, y=resid, colour = Index_Label), size = 2) +
      geom_smooth(aes(x=year, y=resid, colour = Index_Label), lwd = 0.5, se=FALSE, lty=2) +
      facet_wrap(~Index_Label, scale = if(plot_scale) "free" else "fixed")+
      geom_hline(yintercept = 0, size = 1)+
      xlab("Year") +
      xlim(xlim_year) +
      ylab("log(Residual)") +
      theme_SH(base_size = 14)+
      geom_label(data = rho_data,
                 mapping = aes(x = x, y = if(plot_scale)min(d_tidy$resid) else y,
                               label = str_c("sigma=", round(sigma,2),
                                             ", rho=", round(ar1,2), signif)),
                 vjust="inward", hjust="inward")
    g1_sd <- ggplot(d_tidy) +
      geom_ribbon(aes(x = year, ymin = -qnorm(0.025), ymax = qnorm(0.025)), alpha=0.05)+
      geom_ribbon(aes(x = year, ymin = -qnorm(0.1), ymax = qnorm(0.1)), alpha=0.1)+
      geom_point(aes(x=year, y=sd.resid, colour = Index_Label), size = 2) +
      geom_smooth(aes(x=year, y=sd.resid, colour = Index_Label), lwd = 0.5, se=FALSE, lty=2) +
      facet_wrap(~Index_Label, scale = if(plot_scale) "fixed" else "free")+
      geom_hline(yintercept = 0, size = 1)+
      xlab("Year") +
      xlim(xlim_year) +
      ylab("log(Residual)") +
      theme_SH(base_size = 14)+
      geom_label(data = rho_data,
                 mapping = aes(x = x, y = if(plot_scale)min(d_tidy$sd.resid) else y.sd,
                               label = str_c("sigma=", round(sigma,2),
                                             ", rho=", round(ar1,2), signif)),
                 vjust="inward", hjust="inward")
  } else {
    g1 <- ggplot(d_tidy) +
      geom_point(aes(x=year, y=resid, colour = Index_Label), size = 2) +
      geom_smooth(aes(x=year, y=resid, colour = Index_Label), lwd = 0.5, se=FALSE, lty=2) +
      facet_wrap(~Index_Label, scale = if(plot_scale) "free" else "fixed")+
      geom_hline(yintercept = 0, size = 1)+
      xlab("Year") +
      xlim(xlim_year) +
      ylab("log(Residual)") +
      theme_SH(base_size = 14)+
      geom_label(data = rho_data,
                 mapping = aes(x = x, y = if(plot_scale)min(d_tidy$resid) else y,
                               label = str_c("sigma=", round(sigma,2),
                                             ", rho=", round(ar1,2), signif)),
                 vjust="inward", hjust="inward")
    g1_sd <- ggplot(d_tidy) +
      geom_point(aes(x=year, y=sd.resid, colour = Index_Label), size = 2) +
      geom_smooth(aes(x=year, y=sd.resid, colour = Index_Label), lwd = 0.5, se=FALSE, lty=2) +
      facet_wrap(~Index_Label, scale = if(plot_scale) "fixed" else "free")+
      geom_hline(yintercept = 0, size = 1)+
      xlab("Year") +
      xlim(xlim_year) +
      ylab("log(Residual)") +
      theme_SH(base_size = 14)+
      geom_label(data = rho_data,
                 mapping = aes(x = x, y = if(plot_scale)min(d_tidy$sd.resid) else y.sd,
                               label = str_c("sigma=", round(sigma,2),
                                             ", rho=", round(ar1,2), signif)),
                 vjust="inward", hjust="inward")
  }

  g2 <- ggplot(d_tidy) +
    geom_point(aes(x=year, y=obs, colour = Index_Label), size = 2) +
    geom_line(aes(x=year, y=pred, colour = Index_Label), size = 1) +
    facet_wrap(~Index_Label, scale="free") +
    xlim(xlim_year) + ylim(0, NA) +
    ylab("Abundance index") +
    xlab("Year") +
    theme_SH(base_size = 14)

  # 資源量と指数の（非）線形性のプロット
  Lab_tmp <- unique(d_tidy$Index_Label)
  predIndex_g3 <- predabund_g3 <- list()
  for(i in 1:length(Lab_tmp)){
    tmp_data <- d_tidy[d_tidy$Index_Label == Lab_tmp[i],]
    predIndex_g3[[i]] <- with(tmp_data,
                              seq(#min(tmp_data$pred, na.rm = T),
                                0, max(tmp_data$pred, na.rm = T), length=100))
    predabund_g3[[i]] <- (as.numeric(predIndex_g3[[i]])/res$q[i])^(1/res$b[i])
    tmp <- str_split(res$input$abund[i], "") %>% unlist()
    if(sum(tmp == "N") == 0) predabund_g3[[i]] <- predabund_g3[[i]]*res$input$scale
  }
  # 横軸に資源量（指数に合わせてSSBやNだったり）、縦軸に予測CPUEを
  # 線が描けるように、横軸100刻みほどでデータがある
  ab_Index_tmp <- data.frame(Index_Label = rep(unique(d_tidy$Index_Label), each = 100),
                             X = unlist(predabund_g3),
                             Y = unlist(predIndex_g3))

  g3 <- ggplot(d_tidy) +
    geom_point(aes(x=predabund, y=obs, color = Index_Label), size = 2) +
    geom_line(aes(x = X, y = Y), data = ab_Index_tmp, color = "red") +
    facet_wrap(~Index_Label, scales = "free") +
    xlab("Abundance / Biomass / SSB") +
    ylab("Abundance index") +
    xlim(c(0, NA)) + ylim(c(0, NA)) +
    theme_SH(base_size = 14)

  return(list(year_resid = g1,
              fitting_Index = g2,
              abund_Index = g3,
              year_sd_resid = g1_sd,
              gg_data = d_tidy))
} # function(plot_residual_vpa)





#' ジャックナイフ法を行う関数
#'
#' モデル診断法の一つです。影響力の強いデータや外れ値の検出に用いることができます。
#'
#' @param res VPAの結果のオブジェクト
#' @param method データの抜き方。デフォルトは指数ごとに抜く("index")。他に各指数を各点ごとに抜く場合は"all"。
#' @param what_plot 作図したい項目を選べる。引数を与えない場合（\code{NULL}）、全て（SSB, biomass, U, catch, Recruitment, fish_number, fishing_mortality, weight, maturity, catch_number）をプロットする。
#' @param ncol 作図の列数。標準で5列なので、\code{what_plot}の数が5以下の場合は適宜変えた方がよい。
#'
#' @return 返ってくる値:
#'     \code{JKplot_vpa} ジャックナイフ法を行った結果のプロット。
#'
#' @details method :
#' \code{"index"}を行うことで、どの指数シリーズの影響が大きいのか見ることができます（例：漁業CPUEより調査船調査の方が推定の影響が大きい）。
#' \code{"all"}を行うことで、どの指数シリーズのどの年の影響が大きいのか見ることができます。ただし組み合わせが多くなるので、結果が得られるまでに時間を要します。
#'
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' 作図について: \code{\link{plot_vpa}}
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

do_jackknife_vpa <- function(res,
                             method = "index",
                             what_plot = NULL,
                             ncol = 5,
                             plot_year = NULL,
                             scale_value = NULL
){

  if(res$input$use.index == "all"){
    used_index <- res$input$dat$index
  } else {
    used_index <- res$input$dat$index[res$input$use.index,]
  } # 7月7日加筆（浜辺）vpa関数の引数use.index対策

  year <- as.numeric(colnames(res$input$dat$index))
  res_list <- list()
  abund_tmp <- ssb_tmp <- biom_tmp <- tf_tmp <- list()

  if(method == "index"){
    if(length(used_index[,1]) == 1) stop(paste0('The number of indicies is only 1 !!'))

    if(res$input$use.index == "all"){
      name_tmp <- rep(NA, length = length(row.names(res$input$dat$index)))
      for(i in 1:length(name_tmp)){
        input0 <- res$input
        input0$dat$index <- res$input$dat$index[-i,]
        input0$abund <- input0$abund[-i]
        input0$plot <- FALSE
        input0$sigma.const <- input0$sigma.const[-i]
        input0$sigma.constraint <- input0$sigma.constraint[-i]
        res_tmp <- safe_call(vpa, input0, force=TRUE)  # vpa関数の実行

        res_list[[i]] <- res_tmp
        abund_tmp[[i]] <- apply(res_tmp$naa,2,sum)
        ssb_tmp[[i]] <- apply(res_tmp$ssb,2,sum)
        biom_tmp[[i]] <- apply(res_tmp$baa,2,sum)
        tf_tmp[[i]] <- res_tmp$term.f

        if(i <= 9){
          name_tmp[i] <- paste0('Removed index0',i)
        } else {
          name_tmp[i] <- paste0('Removed index',i)
        }
      } #for(i) データの種類について
    } else {
      # use.indexに指定がある場合用のif文分岐の追加
      ## ------------------------------------------------ ##
      # ここエラー出ないようにコンサバにコーディングしてます
      # 2021年度までにはここ修正加えたい
      ## ------------------------------------------------ ##
      name_tmp <- rep(NA, length = length(row.names(used_index)))
      for(i in 1:length(name_tmp)){
        input0 <- res$input
        input0$use.index <- input0$use.index[-i]
        #input0$dat$index <- res$input$dat$index[-i,]
        input0$abund <- input0$abund[-i]
        input0$plot <- FALSE
        input0$sigma.const <- input0$sigma.const[-i]
        input0$sigma.constraint <- input0$sigma.constraint[-i]
        res_tmp <- safe_call(vpa, input0, force=TRUE)  # vpa関数の実行

        res_list[[i]] <- res_tmp
        abund_tmp[[i]] <- apply(res_tmp$naa,2,sum)
        ssb_tmp[[i]] <- apply(res_tmp$ssb,2,sum)
        biom_tmp[[i]] <- apply(res_tmp$baa,2,sum)
        tf_tmp[[i]] <- res_tmp$term.f

        if(i <= 9){
          name_tmp[i] <- paste0('Removed index0',i)
        } else {
          name_tmp[i] <- paste0('Removed index',i)
        }
      } #for(i) データの種類について
    }

  } else if(method == "all"){ ####-----------------------------------------------------####

    if(res$input$use.index == "all"){
      name_tmp <- rep(NA, length = length(res$input$dat$index[!is.na(res$input$dat$index)]))
      for(i in 1:(dim(res$input$dat$index)[1])){
        index_label <- which(is.na(res$input$dat$index[i,])==FALSE)
        year_tmp <- as.numeric(colnames(res$input$dat$index[i,]))[index_label]

        for(j in 1:length(index_label)){
          input0 <- res$input
          index_tmp <- as.numeric(res$input$dat$index[i,])
          index_tmp[index_label[j]] <- NA
          input0$dat$index[i,] <- index_tmp
          input0$plot <- FALSE
          res_tmp <- safe_call(vpa, input0, force=TRUE)  # vpa関数の実行

          if(i == 1){
            res_list[[j]] <- res_tmp
            abund_tmp[[j]] <- apply(res_tmp$naa,2,sum)
            ssb_tmp[[j]] <- apply(res_tmp$ssb,2,sum)
            biom_tmp[[j]] <- apply(res_tmp$baa,2,sum)
            name_tmp[j] <- paste0('Removed index0',i," ",year_tmp[j])
            #tmp <- length(year)*(j-1)+1
            #tf_mat[tmp:tmp+length(year),] <- res_tmp$term.f
            tf_tmp[[j]] <- res_tmp$term.f
          } else if(i <= 9) {
            next_label <- which(is.na(name_tmp))[1]
            res_list[[next_label]] <- res_tmp
            abund_tmp[[next_label]] <- apply(res_tmp$naa,2,sum)
            ssb_tmp[[next_label]] <- apply(res_tmp$ssb,2,sum)
            biom_tmp[[next_label]] <- apply(res_tmp$baa,2,sum)
            name_tmp[next_label] <- paste0('Removed index0',i,' ',year_tmp[j])
            tf_tmp[[next_label]] <- res_tmp$term.f
          } else {
            next_label <- which(is.na(name_tmp))[1]
            res_list[[next_label]] <- res_tmp
            abund_tmp[[next_label]] <- apply(res_tmp$naa,2,sum)
            ssb_tmp[[next_label]] <- apply(res_tmp$ssb,2,sum)
            biom_tmp[[next_label]] <- apply(res_tmp$baa,2,sum)
            name_tmp[next_label] <- paste0('Removed index',i,' ',year_tmp[j])
            tf_tmp[[next_label]] <- res_tmp$term.f
          }
        } #for(j) 各データの時系列について
        #res_tmp2[[i]] <- res_tmp
      } #for(i) データの種類について

    } else {
      ## ------------------------------------------------ ##
      # ここエラー出ないようにコンサバにコーディングしてます
      # 2021年度までにはここ修正加えたい
      ## ------------------------------------------------ ##
      name_tmp <- rep(NA, length = length(used_index[!is.na(used_index)]))
      for(i in 1:(dim(used_index)[1])){
        index_label <- which(is.na(used_index[i,])==FALSE)
        year_tmp <- as.numeric(colnames(used_index[i,]))[index_label]
        use.index_tmp <- res$input$use.index[i]

        for(j in 1:length(index_label)){
          input0 <- res$input
          index_tmp <- as.numeric(res$input$dat$index[use.index_tmp,])
          index_tmp[index_label[j]] <- NA
          input0$dat$index[use.index_tmp,] <- index_tmp
          input0$plot <- FALSE
          res_tmp <- safe_call(vpa, input0, force=TRUE)  # vpa関数の実行

          if(i == 1){
            res_list[[j]] <- res_tmp
            abund_tmp[[j]] <- apply(res_tmp$naa,2,sum)
            ssb_tmp[[j]] <- apply(res_tmp$ssb,2,sum)
            biom_tmp[[j]] <- apply(res_tmp$baa,2,sum)
            name_tmp[j] <- paste0('Removed index0',i," ",year_tmp[j])
            #tmp <- length(year)*(j-1)+1
            #tf_mat[tmp:tmp+length(year),] <- res_tmp$term.f
            tf_tmp[[j]] <- res_tmp$term.f
          } else if(i <= 9) {
            next_label <- which(is.na(name_tmp))[1]
            res_list[[next_label]] <- res_tmp
            abund_tmp[[next_label]] <- apply(res_tmp$naa,2,sum)
            ssb_tmp[[next_label]] <- apply(res_tmp$ssb,2,sum)
            biom_tmp[[next_label]] <- apply(res_tmp$baa,2,sum)
            name_tmp[next_label] <- paste0('Removed index0',i,' ',year_tmp[j])
            tf_tmp[[next_label]] <- res_tmp$term.f
          } else {
            next_label <- which(is.na(name_tmp))[1]
            res_list[[next_label]] <- res_tmp
            abund_tmp[[next_label]] <- apply(res_tmp$naa,2,sum)
            ssb_tmp[[next_label]] <- apply(res_tmp$ssb,2,sum)
            biom_tmp[[next_label]] <- apply(res_tmp$baa,2,sum)
            name_tmp[next_label] <- paste0('Removed index',i,' ',year_tmp[j])
            tf_tmp[[next_label]] <- res_tmp$term.f
          }
        } #for(j) 各データの時系列について
        #res_tmp2[[i]] <- res_tmp
      } #for(i) データの種類について
    }

  } else {                  ####-----------------------------------------------------####
    stop(paste0('Method must be "index" or "all" ! '))
  }

  # plot_vpaですっきり作図！
  names(res_list) <- name_tmp
  gg <- plot_vpa(c(list(Base=res), res_list), what.plot = what_plot,
                 ncol = ncol, plot_year = plot_year, scale_value = scale_value)

  ## ----------------------------------------------------------------- ##
  ## 問題なさそうなら、この区間は消して問題ない
  tf_name <- numeric()
  age_tmp <- as.numeric(rownames(res$input$dat$caa))
  if(length(res$term.f)==1){
    tf_name <- paste0("tf_age", age_tmp[length(age_tmp)])
  } else {
    for(i in 1:length(res$term.f)){
      if(i == length(res$term.f)){
        tf_name[i] <- paste0("tf_age", age_tmp[length(age_tmp)])
      } else {
        tf_name[i] <- paste0("tf_age", age_tmp[i])
      }
    }
  }# if(tf_name)
  tf_tmp2 <- #purrr::map(, function(x) rep(x, length(year))) %>%
    unlist(tf_tmp) %>%
    matrix(nrow = length(res$term.f)) %>%
    t()
  colnames(tf_tmp2) <- tf_name

  if(method == "index"){
    d_tidy_par <- data.frame(tf_tmp2,
                             JK = name_tmp
                             #map(strsplit(name_tmp," "), ~ first(.)) %>%  unlist()
    ) %>%
      pivot_longer(col = c(-JK),
                   names_to = c(".value", "age"),
                   names_sep = "_",
                   values_drop_na = TRUE)
  } else {
    d_tidy_par <- data.frame(tf_tmp2,
                             Removed_index = substr(name_tmp, 1, 15) ,
                             Removed_year = paste(substr(name_tmp, 1, 7), substr(name_tmp, 17, 20))
    ) %>%
      pivot_longer(col = c(-Removed_year, -Removed_index),
                   names_to = c(".value", "age"),
                   names_sep = "_",
                   values_drop_na = TRUE)
  }

  #inputしたvpa結果のターミナルF推定値の作図用のtidy dataに成形
  result_tf <- matrix(res$term.f, ncol = length(res$term.f)) %>% as.data.frame()
  colnames(result_tf) <- tf_name
  result_tf <- pivot_longer(result_tf,
                            cols = tf_name,
                            names_to = c(".value", "age"),
                            names_sep = "_",
                            values_drop_na = TRUE)
  # 作図
  if(method == "all"){
    g4 <- ggplot(data = d_tidy_par) +
      geom_point(aes(x = age, y = tf, col = JK, shape = JK)) +
      facet_wrap(~ Removed_index) +
      theme_SH(legend.position = "top", base_size = 14) +
      scale_shape_manual(values = scale_value[-1])
  } else {
    g4 <- ggplot(data = d_tidy_par) +
      geom_point(aes(x = age, y = tf, col = JK, shape = JK)) +
      theme_SH(legend.position = "top", base_size = 14) +
      scale_shape_manual(values = scale_value[-1])
  }
  ## ----------------------------------------------------------------- ##

  return(list(JKplot_vpa = gg,
              JKplot_par = g4
              # 一応オブジェクト残しているが、JKplot_vpaで事足りるな...
  ))
} # function(do_jackknife_vpa)





#' ブートストラップ法による信頼区間の算出と作図を自動で行う関数
#'
#' @param res VPAの結果のオブジェクト
#' @param B_ite ブートストラップ計算の数。デフォルトで1000。
#' @param B_method ブートストラップの方法。デフォルトではノンパラメトリックブートストラップ。
#' @param ci_range 信頼区間の幅。デフォルトでは0.95（95％信頼区間）
#'
#' @return 返ってくる値:
#'     \code{plot} 親魚重量、資源尾数、資源重量それぞれについて信頼区間のプロットが得られる。
#'     \code{res_boot} ブートストラップ法の結果が得られる。信頼区間の算出に利用可。
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' ブートストラップ法について:  \code{\link{boo.vpa}}
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

plot_resboot_vpa <- function(res, B_ite = 1000, B_method = "p", ci_range = 0.95){
  if(ci_range >= 1) stop(paste0('"ci_range" must be less than 1'))
  res$input$plot <- FALSE
  res_boo <- boo.vpa(res, B = B_ite, method = B_method)

  year <- res_boo[[1]]$index %>% colnames() %>% as.numeric()
  ssb_mat <- abund_mat <- biomass_mat <- matrix(NA, nrow = B_ite, ncol = length(year))
  cor_mat <- matrix(NA, nrow = B_ite,
                    ncol = length(res_boo[[1]]$Fc.at.age)+#length(res_boo[[1]]$q)+
                      length(res_boo[[1]]$b)#+length(res_boo[[1]]$sigma
                    +2)
  for(i in  1:B_ite){
    tmp <- res_boo[[i]]
    ssb_mat[i,] <- colSums(tmp$ssb)
    abund_mat[i,] <- as.numeric(tmp$naa[1,])
    biomass_mat[i,] <- colSums(tmp$baa)
    cor_mat[i, 1:length(tmp$Fc.at.age)] <- tmp$Fc.at.age
    #    cor_mat[i, (length(tmp$Fc.at.age)+1):
    #              (length(tmp$Fc.at.age)+length(tmp$q))] <- tmp$q
    cor_mat[i, (length(tmp$Fc.at.age)+1):
              (length(tmp$Fc.at.age)+length(tmp$b))] <- tmp$b
    #    cor_mat[i, (length(tmp$Fc.at.age)+length(tmp$q)+length(tmp$b)+1):
    #              (length(tmp$Fc.at.age)+length(tmp$q)+length(tmp$b)+length(tmp$sigma))] <- tmp$sigma
    cor_mat[i, length(tmp$Fc.at.age)+length(tmp$b)+1] <- last(colSums(tmp$ssb))
    cor_mat[i, length(tmp$Fc.at.age)+length(tmp$b)+2] <- last(tmp$naa[1,])
  }

  PB_value <- c((1-ci_range)/2, 0.5, 1-(1-ci_range)/2)
  d_ssb <- t(apply(ssb_mat, 2, quantile, probs = PB_value, na.rm = T))
  d_abund <- t(apply(abund_mat, 2, quantile, probs = PB_value, na.rm = T))
  d_biomass <- t(apply(biomass_mat, 2, quantile, probs = PB_value, na.rm = T))
  colnames(d_ssb) <- c("Lower", "SSB", "Upper")
  colnames(d_abund) <- c("Lower", "Abundance", "Upper")
  colnames(d_biomass) <- c("Lower", "Biomass", "Upper")
  d_ssb <- cbind.data.frame(Year = year, d_ssb)
  d_abund <- cbind.data.frame(Year = year, d_abund)
  d_biomass <- cbind.data.frame(Year = year, d_biomass)

  g1 <- ggplot(d_ssb, aes(x = Year, y = SSB))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    ylim(c(0, NA)) +
    theme_SH()

  g2 <- ggplot(d_abund, aes(x = Year, y = Abundance))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    ylab("Recruitment") +
    ylim(c(0, NA)) +
    theme_SH()

  g3 <- ggplot(d_biomass, aes(x = Year, y = Biomass))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    ylim(c(0, NA)) +
    theme_SH()


  row_damy <- apply(cor_mat, 1, function(x) if(sum(is.nan(x))>=1)0 else 1)
  if(class(which(row_damy == 0)) == "numeric") cor_mat <- cor_mat[-which(row_damy == 0),]
  cor_mat <- as.data.frame(cor_mat)
  names(cor_mat) <- c(paste0("term.F_age",1:length(tmp$Fc.at.age)-1),
                      #                      paste0("q",1:length(tmp$q)),
                      paste0("b",1:length(tmp$b)),
                      #                      paste0("sigma",1:length(tmp$sigma)),
                      paste0("SSB_last"),
                      paste0("Recruitment_last")
  )
  cor_mat2 <- cor(cor_mat)
  g4 <- GGally::ggpairs(cor_mat) + theme_SH()

  return(list(plot_ssb = g1,
              plot_rec = g2,
              plot_biomass = g3,
              res_boot = res_boo,
              cor_mat = cor_mat2,
              plot_cor = g4
  ))

} # function(plot_resboot_vpa)





#' 年齢別漁獲尾数の不確実性評価の関数
#'
#' 年齢別漁獲尾数を乱数で再生成した後、VPA計算を行う関数
#'
#' @param res VPAの結果のオブジェクト
#' @param B_ite ブートストラップ計算の数。デフォルトで1000。
#' @param B_cv 乱数生成の変動係数。デフォルトは0.2。
#' @param ci_range 信頼区間の幅。デフォルトでは0.95（95％信頼区間）
#'
#' @return 返ってくる値:
#'     \code{plot} 親魚重量、資源尾数、資源重量それぞれについて信頼区間のプロットが得られる。
#'     \code{caa_boot_sample} 年齢別漁獲尾数のブートストラップ標本が得られる。
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
#'
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

do_caaboot_vpa <-  function(res, B_ite = 1000, B_cv = 0.2, ci_range = 0.95){
  year <- colnames(res$input$dat$caa) %>% as.numeric()
  age <- rownames(res$input$dat$caa) %>% as.numeric()
  caa_base <- res$input$dat$caa %>% unlist() %>% as.numeric()
  caa_boot <- list()
  for(i in 1:length(caa_base)) caa_boot[[i]] <- exp(log(caa_base[i])+rnorm(B_ite, -0.5*B_cv, B_cv))

  name_tmp <- list() ; tmp <- numeric()
  for(i in 1:length(year)){
    for(j in 1:length(age)){
      tmp[j] <- paste0('age',age[j],'_',year[i])
    }
    name_tmp[[i]] <- tmp
  }
  names(caa_boot) <- unlist(name_tmp)

  input0 <- res$input ; input0$plot <- FALSE
  tmp <- numeric()
  ssb_mat <- abund_mat <- biomass_mat <- matrix(NA, ncol = length(year), nrow = B_ite)
  for(i in 1:B_ite){
    for(j in 1:length(caa_boot)) tmp[j] <- caa_boot[[j]][i]
    caa_tmp <- matrix(tmp, ncol = length(year)) %>% as.data.frame()
    colnames(caa_tmp) <- year
    rownames(caa_tmp) <- age
    input0$dat$caa <- caa_tmp
    res_tmp <- try(safe_call(vpa, input0, force=TRUE))
    if(class(res_tmp) == "try-error"){
      message(paste('Iteration',i,'was errored ...', sep = " "))
      ssb_mat[i,] <- rep(NA, length(year))
      abund_mat[i,] <- rep(NA, length(year))
      biomass_mat[i,] <- rep(NA, length(year))
    } else {
      ssb_mat[i,] <- colSums(res_tmp$ssb)
      abund_mat[i,] <- as.numeric(res_tmp$naa[1,])
      biomass_mat[i,] <- colSums(res_tmp$baa)
      message(paste('Iteration',i,'has done ...', sep = " "))
    }
  }

  PB_value <- c((1-ci_range)/2, 0.5, 1-(1-ci_range)/2)
  d_ssb <- t(apply(ssb_mat, 2, quantile, probs = PB_value, na.rm = T))
  d_abund <- t(apply(abund_mat, 2, quantile, probs = PB_value, na.rm = T))
  d_biomass <- t(apply(biomass_mat, 2, quantile, probs = PB_value, na.rm = T))
  colnames(d_ssb) <- c("Lower", "SSB", "Upper")
  colnames(d_abund) <- c("Lower", "Abundance", "Upper")
  colnames(d_biomass) <- c("Lower", "Biomass", "Upper")
  d_ssb <- cbind.data.frame(Year = year, d_ssb)
  d_abund <- cbind.data.frame(Year = year, d_abund)
  d_biomass <- cbind.data.frame(Year = year, d_biomass)

  g1 <- ggplot(d_ssb, aes(x = Year, y = SSB))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    ylim(c(0, NA)) +
    theme_SH()

  g2 <- ggplot(d_abund, aes(x = Year, y = Abundance))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    ylab("Recruitment") +
    geom_line(size = 1.5)+
    ylim(c(0, NA)) +
    theme_SH()

  g3 <- ggplot(d_biomass, aes(x = Year, y = Biomass))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    ylim(c(0, NA)) +
    theme_SH()

  return(list(plot_ssb = g1,
              plot_rec = g2,
              plot_biomass = g3,
              caa_boot_sample = caa_boot
  ))
} # function(do_caaboot_vpa)
