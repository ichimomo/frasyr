
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
#' リッジVPAの場合、感度分析では異なるlambdaでの解析を行う。\code{value}は仮定したい\code{tf.year}を\code{numeric}型で記述のこと。
#'
#' @details "lambda"（リッジVPAのペナルティ）: \code{value}は与えたいlambdaを\code{numeric}型で記述のこと。
#'
#'
#' @details "est.method"（推定方法）: \code{value}は入力不要である。
#'
#' @details "b"(hyperstability/depletion):
#' bを推定した場合、bを推定しない(b=1)の感度分析を行う。
#' bを1以外の値で固定した、あるいはbを考慮しなかった場合、\code{value="b.est"}とするとbの推定を、仮定したいbを\code{numeric}で与えるとbを固定して考慮出来る。
#'
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @seealso
#' vpa計算について:  \code{\link{vpa}}
#' 作図について: \code{\link{plot_vpa}}
#'
#' @examples # WikiかVignetteにつながるようにする？
#'
#' @encoding UTF-8
#'
#' @export

# helpは↑のように、"#`"のあとに一定の書式に従って記述するのでこちらも書いてみてください。少しだけ試しに書いてみました。@export, @encoding UTF-8とかはそのままにしておいてください

# author: Kohei Hamabe

do_sensitivity_vpa <- function(res, what_replace, value, what_plot = NULL, ncol = 5){
  res_vpa.s <- list() ; lab.tmp <- numeric()
  res$input$plot <- FALSE # 背後で沢山plotするのを防ぐ

  if(what_replace == "M"){
    if(!class(value) == "numeric")stop(paste('Set values which shape is only "numeric"'))
    for(i in 1:length(value)){
      input0 <- res$input
      input0$dat$M <- input0$dat$M *value[i]
      res_vpa.s[[i]] <- do.call(vpa, input0)
      lab.tmp[i] <- paste("Sensitivity M= x", value[i], sep = "")
    } # for

    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "waa") {

    if(class(value) == "numeric") {
      ## 任意の比率(ex. 1.5倍)で入れた場合
      for(i in 1:length(value)){
        input0 <- res$input
        input0$dat$waa <- input0$dat$waa *value[i]
        res_vpa.s[[i]] <- do.call(vpa, input0)
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
        res_vpa.s[[i]] <- do.call(vpa, input0)
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
        res_vpa.s[[i]] <- do.call(vpa, input0)
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
        res_vpa.s[[i]] <- do.call(vpa, input0)
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
      res_vpa.s[[i]] <- do.call(vpa, input0)
      lab.tmp[i] <- paste("Sensitivity maa case:", i, sep = "")
    } # for

    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "alpha") {
    for(i in 1:length(value)){
      input0 <- res$input
      input0$alpha <- value[i]
      res_vpa.s[[i]] <- do.call(vpa, input0)
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
        res_vpa.s[[i]] <- do.call(vpa, input0)
        lab.tmp[i] <- paste("tf.year = ", value[[i]][1], ":", value[[i]][length(value[[i]])], sep = "")
      }
    } else if(res$input$sel.update == TRUE){  # resは選択率更新
      # 全F推定を試みる
      input0 <- res$input
      input0$sel.update <- FALSE
      input0$term.F <- "all"
      res_vpa.s[[1]] <- do.call(vpa, input0)
      lab.tmp <- "ALL_F.est"
      #   } else if(!res$input$lambda == 0){  # resはRidge VPA
      #       # lambdaを変えて感度分析
      #      if(!class(value) == "numeric")stop(paste('For "tuning", the shape of values were only "numeric" in ridge VPA', sep = ""))
      #      for(i in 1:length(value)){
      #        input0 <- res$input
      #        input0$lambda <- value[i]
      #        res_vpa.s[[i]] <- do.call(vpa, input0)
      #        lab.tmp[i] <- paste("lambda = ", value[i], sep = "")
      #      }
    } else  {
      # 選択率更新法を行う
      if(!class(value) == "list")stop(paste('Set tf.years in the values, and the shape was only "list"'))
      for(i in 1:length(value)){
        input0 <- res$input
        input0$sel.update <- TRUE
        input0$tf.year <- value[[i]]
        res_vpa.s[[i]] <- do.call(vpa, input0)
        lab.tmp[i] <- paste("sel.update, tf.year = ", value[[i]][1], ":", value[[i]][length(value[[i]])], sep = "")
      }
    }
    #-------------------------------------------------------------------------------------#

  } else if (what_replace == "lambda") {
    if(!class(value) == "numeric")stop(paste('For "lambda", the shape of values were only "numeric" in ridge VPA', sep = ""))
    for(i in 1:length(value)){
      input0 <- res$input
      input0$lambda <- value[i]
      res_vpa.s[[i]] <- do.call(vpa, input0)
      lab.tmp[i] <- paste("lambda = ", value[i], sep = "")
    }

    #-------------------------------------------------------------------------------------#

  } else if(what_replace == "b") {
    if(res$input$b.est == TRUE){
      input0 <- res$input
      input0$b.est <- FALSE
      res_vpa.s[[1]] <- do.call(vpa, input0)
      lab.tmp <- "b = 1"
    } else if(class(value) == "numeric") {
      for(i in 1:length(value)){
        input0 <- res$input
        input0$b.fix <- value[i]
        res_vpa.s[[i]] <- do.call(vpa, input0)
        lab.tmp[i] <- paste("b.fix = ", value[i], sep = "")
      }
    } else {
      if(!value == "b.est")stop(paste('Input the "b.est" or some values(numeric) in "value="!!', sep = ""))
      input0 <- res$input
      input0$b.est <- TRUE
      res_vpa.s[[1]] <- do.call(vpa, input0)
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
      res_vpa.s[[1]] <- do.call(vpa, input0)
      lab.tmp <- "est.method: ml"
    } else {
      input0 <- res$input
      input0$est.method <- "ls"
      res_vpa.s[[1]] <- do.call(vpa, input0)
      lab.tmp <- "est.method: ls"
    }

    #-------------------------------------------------------------------------------------#

  } else {
    stop(paste('You have to select the what_replace from "M", "waa", "waa.catch", "maa", "alpha", "tuning", "b" or "est.method"', sep = ""))
  }

  names(res_vpa.s) <- lab.tmp
  # 結果のラベルに名前を付ける
  g1 <- plot_vpa(c(list(Base=res), res_vpa.s),
                 what.plot = what_plot, ncol = ncol)
  return(list(result = res_vpa.s, graph = g1))
} # function(do_sensitivity_vpa)





#' VPAのレトロスペクティブ解析結果を自動生成する関数
#'
#' @param res VPAの結果のオブジェクト
#' @param n_retro レトロスペクティブ解析でさかのぼる年数
#' @param b_reest bをレトロスペクティブ解析中で再推定するか
#' @param what_plot 作図したい項目を選べる。引数を与えない場合（\code{NULL}）、全て（SSB, biomass, U, catch, Recruitment, fish_number, fishing_mortality, weight, maturity, catch_number）をプロットする。
#' @param ncol 作図の列数。標準で5列なので、\code{what_plot}の数が5以下の場合は適宜変えた方がよい。
#'
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
#'
#' @examples # GitHubにつながるようにする？
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

do_retrospective_vpa <- function(res, n_retro = 5, b_reest = FALSE, what_plot = NULL, ncol = 5){

  if(b_reest == TRUE && res$input$b.est == FALSE)message(paste('b was not estimated in your vpa model'))
  # vpa内でbの推定をしていないにもかかわらず、b_reestがtrueで入力された場合
  # 推定結果(bを推定している)は得られるが、メッセージを出す

  res_retro <- retro.est(res, n = n_retro, b.fix = !b_reest)
  dat_graph <- list()
  for(i in 1:n_retro) dat_graph[[i]] <- res_retro$Res[[i]]
  names(dat_graph) <- rev(colnames(res$ssb))[1:n_retro]  # 図にinputされる結果に名前をつける
  g1 <- plot_vpa(dat_graph, what.plot = what_plot, ncol = ncol)
  return(list(result = res_retro, mohn_rho = res_retro$mohn, graph = g1))
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
#' @encoding UTF-8
#'
#' @export


do_estcheck_vpa <- function(res, n_ite = 100, sd_jitter = 1, what_plot = NULL, TMB = FALSE){
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

  init_list <- purrr::map(res$term.f,
                          function(x)exp(log(x) + rnorm(n_ite, 0, sd_jitter))
  )
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
    tmp <- try(do.call(vpa, input0))
    if(class(tmp) == "try-error") next
    value_tmp[[i]] <- list(p_est = tmp$term.f,
                           logLik = tmp$logLik,
                           covergence = tmp$convergence,
                           hessian = tmp$hessian,
                           gradient = tmp$gradient)
    ite_tmp[[i]] <- rep(i, length(res$term.f))
    ll_tmp[[i]] <- rep(res$logLik, length(res$logLik))
    Finit[[i]] <- init_tmp
    Fest[[i]] <- tmp$term.f
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
                      likelihood = unlist(ll_tmp)
  )

  g1 <- ggplot(data = d_tmp[d_tmp$age == plot_name,]) +
    geom_point(aes(x = initial, y = estimated), size = 5) +
    facet_wrap( ~ age) +
    theme_SH(base_size = 14)
  g2 <- ggplot(data = d_tmp[d_tmp$age == plot_name,]) +
    geom_point(aes(x = initial, y = likelihood), size = 5) +
    facet_wrap( ~ age) +
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
  return(list(initial_value = init_list, #初期値の乱数
              p_name = name_tmp, # 初期値の名前
              value = value_tmp, #　推定値と尤度のリスト
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
#' @encoding UTF-8
#'
#' @export

plot_residual_vpa <- function(res, index_name = NULL, plot_scale = FALSE){
  d_tmp <- matrix(NA,
                  nrow = length(res$input$dat$index[1,]),
                  ncol = length(res$input$dat$index[,1])*8+4)
  d_tmp[,1] <- as.numeric(colnames(res$input$dat$index))
  d_tmp[,2:(1+length(res$q))] <- as.numeric(t(res$input$dat$index))
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
    if(res$input$abund[i] == "N"){
      d_tmp[,(i+length(res$q)*1+4)] <- (res$pred.index[i,]/res$q[i])^(1/res$b[i])
    } else {
      # 資源重量の場合、リスケーリング分かけて補正
      d_tmp[,(i+length(res$q)*1+4)] <- (res$pred.index[i,]/res$q[i])^(1/res$b[i])*res$input$scale
    }
    d_tmp[,(i+length(res$q)*2+4)] <- res$pred.index[i,]　# q*N^B計算結果
    d_tmp[,(i+length(res$q)*3+4)] <- resid_tmp
    d_tmp[,(i+length(res$q)*4+4)] <- sd_resid_tmp
    d_tmp[,(i+length(res$q)*5+4)] <- rep(res$q[i], length(d_tmp[,1]))
    d_tmp[,(i+length(res$q)*6+4)] <- rep(res$sigma[i], length(d_tmp[,1]))
    d_tmp[,(i+length(res$q)*7+4)] <- rep(res$b[i], length(d_tmp[,1]))
    if(i >= 10){
      name_tmp1[i] <- paste0("obs_CPUE",i)
      name_tmp2[i] <- paste0("predabund_CPUE",i)
      name_tmp3[i] <- paste0("pred_CPUE",i)
      name_tmp4[i] <- paste0("resid_CPUE",i)
      name_tmp5[i] <- paste0("sd.resid_CPUE",i)
      q_tmp[i] <- paste0("q_CPUE",i)
      b_tmp[i] <- paste0("b_CPUE",i)
      sig_tmp[i] <- paste0("sigma_CPUE",i)
    } else {
      name_tmp1[i] <- paste0("obs_CPUE0",i)
      name_tmp2[i] <- paste0("predabund_CPUE0",i)
      name_tmp3[i] <- paste0("pred_CPUE0",i)
      name_tmp4[i] <- paste0("resid_CPUE0",i)
      name_tmp5[i] <- paste0("sd.resid_CPUE0",i)
      q_tmp[i] <- paste0("q_CPUE0",i)
      b_tmp[i] <- paste0("b_CPUE0",i)
      sig_tmp[i] <- paste0("sigma_CPUE0",i)
    }
  }
  d_tmp <- as.data.frame(d_tmp)
  names(d_tmp) <- c("year", name_tmp1, "abundance",  "biomass", "ssb", name_tmp2,
                    name_tmp3, name_tmp4, name_tmp5, q_tmp, sig_tmp, b_tmp)

  d_tidy <- tidyr::pivot_longer(d_tmp, col = c(-year, -abundance, -biomass, -ssb),
                                names_to = c(".value", "CPUE_Label"),
                                names_sep = "_",
                                values_drop_na = TRUE
  )

  if(!is.null(index_name)){
    if(!length(index_name) == length(res$q)) stop(paste0("Length of index_name was different to the number of indices"))
    d_tidy$CPUE_Label <- rep(index_name, length(d_tmp[,1]))
  }

  if(plot_scale) scale_tmp <- "fixed" else scale_tmp <- "free"

  g1 <- ggplot(d_tidy) +
    geom_point(aes(x=year, y=resid, colour = CPUE_Label), size = 2) +
    geom_smooth(aes(x=year, y=resid, colour = CPUE_Label), size = 1.5) +
    facet_wrap(~CPUE_Label, scale=scale_tmp)+
    geom_hline(yintercept = 0, size = 1)+
    theme_SH(base_size = 14)
  g1_sd <- ggplot(d_tidy) +
    geom_point(aes(x=year, y=sd.resid, colour = CPUE_Label), size = 2) +
    geom_smooth(aes(x=year, y=sd.resid, colour = CPUE_Label), size = 1.5) +
    facet_wrap(~CPUE_Label, scale=scale_tmp) +
    geom_hline(yintercept = 0, size = 1)+
    theme_SH(base_size = 14)

  g2 <- ggplot(d_tidy) +
    geom_point(aes(x=year, y=obs, colour = CPUE_Label), size = 2) +
    geom_line(aes(x=year, y=pred, colour = CPUE_Label), size = 1) +
    facet_wrap(~CPUE_Label, scale=scale_tmp) +
    theme_SH(base_size = 14)

  # 資源量と指数の（非）線形性のプロット
  Lab_tmp <- unique(d_tidy$CPUE_Label)
  predCPUE_g3 <- predabund_g3 <- list()
  for(i in 1:length(Lab_tmp)){
    tmp_data <- d_tidy[d_tidy$CPUE_Label == Lab_tmp[i],]
    #Y <- which(is.na(tmp_data$predabund) == is.na(tmp_data$obs))
    predabund_g3[[i]] <- with(tmp_data,
                              seq(min(tmp_data$predabund, na.rm = T),
                                  max(tmp_data$predabund, na.rm = T), length=100))
    predCPUE_g3[[i]] <- with(tmp_data,
                             seq(min(tmp_data$pred, na.rm = T),
                                 max(tmp_data$pred, na.rm = T), length=100))
  }
  # 横軸に資源量（指数に合わせてSSBやNだったり）、縦軸に予測CPUEを
  # 線が描けるように、横軸100刻みほどでデータがある
  ab_CPUE_tmp <- data.frame(CPUE_Label = rep(unique(d_tidy$CPUE_Label), each = 100),
                            X = unlist(predabund_g3),
                            #X = (unlist(CPUE_g3)/unlist(q_g3))^(1/unlist(b_g3))*res$input$scale,
                            Y = unlist(predCPUE_g3))
  #Y = unlist(pred_g3))

  g3 <- ggplot(d_tidy) +
    geom_point(aes(x=predabund, y=obs, color = CPUE_Label), size = 2) +
    geom_line(aes(x = X, y = Y), data = ab_CPUE_tmp, color = "red") +
    facet_wrap(~CPUE_Label, scales = scale_tmp) +
    xlab("Abundance / Biomass / SSB") +
    ylab("") +
    theme_SH(base_size = 14)

  return(list(year_resid = g1,
              fitting_CPUE = g2,
              abund_CPUE = g3,
              year_sd_resid = g1_sd,
              gg_data = d_tidy))
} # function(plot_residual_vpa)





#' ジャックナイフ法を行う関数
#'
#' モデル診断法の一つです。影響力の強いデータや外れ値の検出に用いることができます。
#'
#' @param res VPAの結果のオブジェクト
#'
#' @return 返ってくる値:
#'     \code{JKplot_abund} ジャックナイフ法を行った資源尾数の時系列プロット。
#'     \code{JKplot_ssb} ジャックナイフ法を行った親魚重量の時系列プロット。
#'     \code{JKplot_biomass} ジャックナイフ法を行った資源重量の時系列プロット。
#'     \code{JKplot_par} ジャックナイフ法を行ったターミナルFの時系列プロット。
#'     \code{data1} 資源尾数の時系列プロットに用いたtidyデータ。
#'     \code{data2} ターミナルFのプロットに用いたtidyデータ。
#'
#'
#' @author 濵邉昂平, 市野川桃子
#'
#' @encoding UTF-8
#'
#' @export

do_jackknife_vpa <- function(res){
  year <- as.numeric(colnames(res$input$dat$index))
  res_list <- list()
  abund_tmp <- ssb_tmp <- biom_tmp <- tf_tmp <- list()
  name_tmp <- rep(NA, length = length(res$input$dat$index[!is.na(res$input$dat$index)]))*3
  #tf_mat <- matrix(NA, ncol = length(res$term.f), nrow = length(name_tmp)*length(year))
  for(i in 1:(dim(res$input$dat$index)[1])){
    #for(i in 1:(dim(res$input$dat$index)[1])){
    index_label <- which(is.na(res$input$dat$index[i,])==FALSE)
    year_tmp <- as.numeric(colnames(res$input$dat$index[i,]))[index_label]

    for(j in 1:length(index_label)){
      #for(j in 1:2){
      input0 <- res$input
      index_tmp <- as.numeric(res$input$dat$index[i,])
      index_tmp[index_label[j]] <- NA
      input0$dat$index[i,] <- index_tmp
      input0$plot <- FALSE
      res_tmp <- do.call(vpa, input0)  # vpa関数の実行

      if(i == 1){
        res_list[[j]] <- res_tmp
        abund_tmp[[j]] <- apply(res_tmp$naa,2,sum)
        ssb_tmp[[j]] <- apply(res_tmp$ssb,2,sum)
        biom_tmp[[j]] <- apply(res_tmp$baa,2,sum)
        name_tmp[j] <- paste0('CPUE0',i," ",year_tmp[j])
        #tmp <- length(year)*(j-1)+1
        #tf_mat[tmp:tmp+length(year),] <- res_tmp$term.f
        tf_tmp[[j]] <- res_tmp$term.f
      } else if(i <= 9) {
        next_label <- which(is.na(name_tmp))[1]
        res_list[[next_label]] <- res_tmp
        abund_tmp[[next_label]] <- apply(res_tmp$naa,2,sum)
        ssb_tmp[[next_label]] <- apply(res_tmp$ssb,2,sum)
        biom_tmp[[next_label]] <- apply(res_tmp$baa,2,sum)
        name_tmp[next_label] <- paste0('CPUE0',i,' ',year_tmp[j])
        #tf_mat[next_label:(next_label+length(year)-1),] <- res_tmp$term.f
        tf_tmp[[next_label]] <- res_tmp$term.f
      } else {
        next_label <- which(is.na(name_tmp))[1]
        res_list[[next_label]] <- res_tmp
        abund_tmp[[next_label]] <- apply(res_tmp$naa,2,sum)
        ssb_tmp[[next_label]] <- apply(res_tmp$ssb,2,sum)
        biom_tmp[[next_label]] <- apply(res_tmp$baa,2,sum)
        name_tmp[next_label] <- paste0('CPUE',i,' ',year_tmp[j])
        #tf_mat[next_label:(next_label+length(year)-1),] <- res_tmp$term.f
        tf_tmp[[next_label]] <- res_tmp$term.f
      }
    } #for(j) 各データの時系列について
    #res_tmp2[[i]] <- res_tmp
  } #for(i) データの種類について

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

  d_tidy_abund <- data.frame(year = rep(year[1]:year[length(year)], length(abund_tmp)),
                             abundance = unlist(abund_tmp),
                             ssb = unlist(ssb_tmp),
                             biomass = unlist(biom_tmp)
  ) %>%
    cbind.data.frame(JK_index = map(strsplit(name_tmp," "), ~ first(.)) %>%
                       unlist()  %>%  rep(each = length(year))) %>%
    cbind.data.frame(JK_year = map(strsplit(name_tmp," "), ~ last(.)) %>%
                       unlist()  %>%  rep(each = length(year)))


  #inputしたvpa結果を作図用のtidy dataに成形
  result <- data.frame(year = year,
                       abundance = apply(res$naa, 2, sum),
                       ssb = apply(res$ssb, 2, sum),
                       biomass = apply(res$baa, 2, sum)
  )

  g1 <- ggplot(d_tidy_abund) +
    geom_line(data = result, aes(x = year, y = abundance), size = 2)+
    geom_line(aes(x = year, y = abundance, color = JK_year)) +
    facet_wrap(~ JK_index, scale="free_y")+
    theme_SH()

  g2 <- ggplot(d_tidy_abund) +
    geom_line(data = result, aes(x = year, y = ssb), size = 2)+
    geom_line(aes(x = year, y = ssb, color = JK_year)) +
    facet_wrap(~ JK_index, scale="free_y")+
    theme_SH()

  g3 <- ggplot(d_tidy_abund) +
    geom_line(data = result, aes(x = year, y = biomass), size = 2)+
    geom_line(aes(x = year, y = biomass, color = JK_year)) +
    facet_wrap(~ JK_index, scale="free_y")+
    theme_SH()


  d_tidy_par <- data.frame(tf_tmp2,
                           JK_index = map(strsplit(name_tmp," "), ~ first(.)) %>%  unlist(),
                           JK_year = map(strsplit(name_tmp," "), ~ last(.)) %>%  unlist()
  ) %>%
    pivot_longer(col = c(-JK_year, -JK_index),
                 names_to = c(".value", "age"),
                 names_sep = "_",
                 values_drop_na = TRUE)
  result_tf <- matrix(res$term.f, ncol = length(res$term.f)) %>% as.data.frame()
  colnames(result_tf) <- tf_name
  result_tf <- pivot_longer(result_tf,
                            cols = tf_name,
                            names_to = c(".value", "age"),
                            names_sep = "_",
                            values_drop_na = TRUE)

  g4 <- ggplot(data = d_tidy_par) +
    geom_point(aes(x = JK_year, y = tf, col = age)) +
    #for(i in 1:length(result_tf$age)) geom_hline(yintercept =  result_tf$tf[i]) +
    facet_wrap(~ JK_index) +
    theme_SH()

  return(list(JKplot_abund = g1,
              JKplot_ssb = g2,
              JKplot_biomass = g3,
              JKplot_par = g4,
              data1 = d_tidy_abund,
              data2 = d_tidy_par
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
#'
#' @examples # GitHubにつながるようにする？
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
  for(i in  1:B_ite){
    tmp <- res_boo[[i]]
    ssb_mat[i,] <- colSums(tmp$ssb)
    abund_mat[i,] <- colSums(tmp$naa)
    biomass_mat[i,] <- colSums(tmp$baa)
  }

  PB_value <- c((1-ci_range)/2, 0.5, 1-(1-ci_range)/2)
  d_ssb <- t(apply(ssb_mat, 2, quantile, probs = PB_value))
  d_abund <- t(apply(abund_mat, 2, quantile, probs = PB_value))
  d_biomass <- t(apply(biomass_mat, 2, quantile, probs = PB_value))
  colnames(d_ssb) <- c("Lower", "SSB", "Upper")
  colnames(d_abund) <- c("Lower", "Abundance", "Upper")
  colnames(d_biomass) <- c("Lower", "Biomass", "Upper")
  d_ssb <- cbind.data.frame(Year = year, d_ssb)
  d_abund <- cbind.data.frame(Year = year, d_abund)
  d_biomass <- cbind.data.frame(Year = year, d_biomass)

  g1 <- ggplot(d_ssb, aes(x = Year, y = SSB))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    theme_SH()

  g2 <- ggplot(d_biomass, aes(x = Year, y = Abundance))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    theme_SH()

  g3 <- ggplot(d_biomass, aes(x = Year, y = Biomass))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    theme_SH()

  return(list(plot_ssb = g1,
              plot_abund = g2,
              plot_biomass = g3,
              res_boot = res_boo
  ))

} # function(plot_resboot_vpa)





#' 年齢別漁獲尾数の不確実性評価の関数
#'
#' 年齢別漁獲尾数を乱数で再生成した後、VPA計算を行う関数
#'
#' @param res VPAの結果のオブジェクト
#' @param B_ite ブートストラップ計算の数。デフォルトで1000。
#' @param B_sd 乱数生成の標準誤差。
#' @param ci_range 信頼区間の幅。デフォルトでは0.95（95％信頼区間）
#'
#' @return 返ってくる値:
#'     \code{plot} 親魚重量、資源尾数、資源重量それぞれについて信頼区間のプロットが得られる。
#'     \code{caa_boot_sample} 年齢別漁獲尾数のブートストラップ標本が得られる。
#'
#' @author 濵邉昂平, 市野川桃子
#'
#'
#' @examples GitHubにつながるようにする？
#'
#' @encoding UTF-8
#'
#' @export

# author: Kohei Hamabe

do_caaboot_vpa <-  function(res, B_ite = 1000, B_sd = 1, ci_range = 0.95){
  year <- colnames(res$input$dat$caa) %>% as.numeric()
  age <- rownames(res$input$dat$caa) %>% as.numeric()
  caa_base <- res$input$dat$caa %>% unlist() %>% as.numeric()
  caa_boot <- list()
  for(i in 1:length(caa_base)) caa_boot[[i]] <- exp(log(caa_base[i])+rnorm(B_ite, -0.5*B_sd, B_sd))

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
    res_tmp <- do.call(vpa, input0)
    ssb_mat[i,] <- colSums(res_tmp$ssb)
    abund_mat[i,] <- colSums(res_tmp$naa)
    biomass_mat[i,] <- colSums(res_tmp$baa)
  }

  PB_value <- c((1-ci_range)/2, 0.5, 1-(1-ci_range)/2)
  d_ssb <- t(apply(ssb_mat, 2, quantile, probs = PB_value))
  d_abund <- t(apply(abund_mat, 2, quantile, probs = PB_value))
  d_biomass <- t(apply(biomass_mat, 2, quantile, probs = PB_value))
  colnames(d_ssb) <- c("Lower", "SSB", "Upper")
  colnames(d_abund) <- c("Lower", "Abundance", "Upper")
  colnames(d_biomass) <- c("Lower", "Biomass", "Upper")
  d_ssb <- cbind.data.frame(Year = year, d_ssb)
  d_abund <- cbind.data.frame(Year = year, d_abund)
  d_biomass <- cbind.data.frame(Year = year, d_biomass)

  g1 <- ggplot(d_ssb, aes(x = Year, y = SSB))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    theme_SH()

  g2 <- ggplot(d_biomass, aes(x = Year, y = Abundance))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    theme_SH()

  g3 <- ggplot(d_biomass, aes(x = Year, y = Biomass))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "blue")+
    geom_line(size = 1.5)+
    theme_SH()

  return(list(plot_ssb = g1,
              plot_abund = g2,
              plot_biomass = g3,
              caa_boot_sample = caa_boot
  ))

}


