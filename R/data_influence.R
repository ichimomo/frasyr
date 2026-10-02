#' 推定資源量がどのデータにどれだけ依存しているかを感度解析で調べる
#'
#' データ (漁獲尾数 CAA と資源量指数) の 1 セルずつを対数スケールで ±h 動かして推定し直し、
#' 推定資源量 N の各セル (年齢 x 年) の弾力性 d log N / d log(データ) を中心差分で求める。
#' VPA (frasyr::vpa の結果) と SAM (frasam::sam の結果) の両方に使える。
#'
#' 再推定は `fit$input` をそのまま使う (VPA のリッジの λ も、最適化の出発点 p.init も元のまま)。
#' 出発点を元の解にすると、データを少し動かしても目的関数の変化が許容誤差より小さく、
#' 最適化がその場で止まって弾力性が 0 になることがあるため、出発点は変えない。
#' SAM (sam2) と far_vpa (TMB 版 VPA, rvpa_tmb) の DLL は、呼び出し側で読み込んでおくこと
#' (`use_sam_tmb()` / `use_rvpa_tmb()`)。
#'
#' @param fit   `frasyr::vpa()` または `frasam::sam()` の結果 (`$input` を持つもの)
#' @param h     動かす幅 (対数スケール)。既定 0.01 (±1%)
#' @param sources 動かすデータ。`"caa"` と `"index"` のどちらか、または両方
#' @param cores 並列数 (`parallel::mclapply`)
#'
#' @return list:
#'   `elasticity` = long tibble (source, data_row, data_year, n_age, n_year, elasticity)。
#'   source は "caa" (data_row = 年齢) / "index" (data_row = 指数の番号)。
#'   `base_naa` = 再推定した基準の N、`base_diff` = 元の fit の N との最大相対差 (再現確認)、
#'   `failed` = 再推定に失敗したデータセル、`model`, `h`
#'
#' @examples
#' \dontrun{
#' infl <- data_influence(res$res_hat$vdat[[1]]$vpa, cores = 8)
#' dep  <- data_dependence(infl, caa_sd = 0.3, index_sd = 0.1)
#' plot_data_dependence(dep)
#' plot_data_influence(infl, n_age = 0, n_year = "last")
#' }
#' @export
data_influence <- function(fit, h = 0.01, sources = c("caa", "index"), cores = 1) {
  sources <- match.arg(sources, several.ok = TRUE)
  is_sam <- inherits(fit, "sam")
  input  <- fit$input
  if (is.null(input) || is.null(input$dat))
    stop("fit$input$dat がありません。vpa() / sam() の結果を渡してください。", call. = FALSE)

  refit <- function(dat) {
    inp <- input
    inp$dat <- dat
    m <- if (is_sam) {
      inp$silent <- TRUE
      try(suppressWarnings(safe_do_call(sam, inp)), silent = TRUE)
    } else {
      # force = TRUE: vpa() の引数にない項目 (保存された input に入っている dots など) を外して呼ぶ
      try(suppressWarnings(frasyr::safe_call(frasyr::vpa, inp, force = TRUE)), silent = TRUE)
    }
    if (inherits(m, "try-error") || is.null(m$naa)) return(NULL)
    naa <- as.matrix(m$naa)
    if (any(!is.finite(naa)) || any(naa <= 0)) return(NULL)
    naa
  }

  base <- refit(input$dat)
  if (is.null(base)) stop("元のデータでの再推定に失敗しました。", call. = FALSE)
  base_diff <- max(abs(base / as.matrix(fit$naa) - 1))

  # 動かすデータセル (正で有限の値のみ。指数の欠測は飛ばす)
  cells <- dplyr::bind_rows(lapply(sources, function(src) {
    x <- as.matrix(input$dat[[src]])
    idx <- which(is.finite(x) & x > 0, arr.ind = TRUE)
    tibble::tibble(source = src, i = idx[, 1], j = idx[, 2])
  }))

  one_cell <- function(k) {
    src <- cells$source[k]; i <- cells$i[k]; j <- cells$j[k]
    bump <- function(s) {
      dat <- input$dat
      x <- as.matrix(dat[[src]]); x[i, j] <- x[i, j] * exp(s * h)
      dat[[src]][i, j] <- x[i, j]      # data.frame / matrix どちらでも元の型を保つ
      refit(dat)
    }
    up <- bump(1); dn <- bump(-1)
    if (is.null(up) || is.null(dn)) return(NULL)
    (log(up) - log(dn)) / (2 * h)
  }
  el <- parallel::mclapply(seq_len(nrow(cells)), one_cell, mc.cores = cores)

  ok <- !vapply(el, is.null, logical(1))
  nA <- nrow(base); nY <- ncol(base)
  years  <- colnames(base) %||% as.character(seq_len(nY))
  data_years <- lapply(sources, function(src) colnames(as.matrix(input$dat[[src]])) %||% years)
  names(data_years) <- sources

  long <- dplyr::bind_rows(lapply(which(ok), function(k) {
    tibble::tibble(source = cells$source[k], data_row = cells$i[k],
                   data_year = data_years[[cells$source[k]]][cells$j[k]],
                   n_age = rep(seq_len(nA) - 1, times = nY),
                   n_year = rep(years, each = nA),
                   elasticity = as.vector(el[[k]]))
  }))

  list(elasticity = long, base_naa = base, base_diff = base_diff,
       failed = cells[!ok, ], model = if (is_sam) "sam" else "vpa", h = h)
}


#' 推定資源量の誤差の分散を、CAA 由来と資源量指数由来に分ける
#'
#' [data_influence()] の弾力性 e に、各データの対数スケールの誤差 SD を掛けて二乗和をとる
#' (デルタ法、データの誤差は互いに独立と仮定)。
#' var_source(N[a, y]) = sum_{そのデータ源のセル} (e * sd)^2
#'
#' @param infl     [data_influence()] の返り値
#' @param caa_sd   CAA の誤差 SD。スカラー、年齢別ベクトル、または 年齢 x 年 の行列
#' @param index_sd 指数の誤差 SD。スカラー、指数別ベクトル、または 指数 x 年 の行列
#'
#' @return tibble (n_age, n_year, sd_caa, sd_index, share_index)。
#'   share_index = 指数由来の分散 / 全分散 (1 に近いほど CPUE で決まっている)
#' @export
data_dependence <- function(infl, caa_sd = 1, index_sd = 1) {
  sd_of <- function(sdv, row, col_idx) {
    if (is.matrix(sdv)) sdv[cbind(row, col_idx)]
    else if (length(sdv) == 1) rep(sdv, length(row))
    else sdv[row]
  }
  e <- infl$elasticity |>
    dplyr::mutate(col_idx = match(.data$data_year, unique(.data$data_year)),
                  sd = dplyr::if_else(.data$source == "caa",
                                      sd_of(caa_sd, .data$data_row, .data$col_idx),
                                      sd_of(index_sd, .data$data_row, .data$col_idx)),
                  v = (.data$elasticity * .data$sd)^2)
  e |>
    dplyr::summarise(v = sum(.data$v), .by = c("n_age", "n_year", "source")) |>
    tidyr::pivot_wider(names_from = "source", values_from = "v", values_fill = 0,
                       names_prefix = "v_") |>
    dplyr::mutate(sd_caa   = sqrt(.data$v_caa),
                  sd_index = sqrt(.data$v_index),
                  share_index = .data$v_index / (.data$v_caa + .data$v_index)) |>
    dplyr::select("n_age", "n_year", "sd_caa", "sd_index", "share_index")
}


#' 推定資源量の各セルが資源量指数にどれだけ依存しているかのヒートマップ
#'
#' @param dep   [data_dependence()] の返り値 (モデル比較なら `model` 列を足して行結合したもの)
#' @param value 塗る量。`"share_index"` (指数由来の分散の割合) / `"sd_caa"` / `"sd_index"`
#' @param facet facet に使う列 (例 `"model"`)。`NULL` なら facet しない
#' @return ggplot オブジェクト
#' @export
plot_data_dependence <- function(dep, value = c("share_index", "sd_caa", "sd_index"),
                                 facet = if ("model" %in% names(dep)) "model" else NULL) {
  value <- match.arg(value)
  lab <- c(share_index = "指数由来の\n分散の割合", sd_caa = "CAA 由来の\nSD (log N)",
           sd_index = "指数由来の\nSD (log N)")[[value]]
  p <- ggplot2::ggplot(dep, ggplot2::aes(factor(.data$n_year, levels = unique(.data$n_year)),
                                         .data$n_age, fill = .data[[value]])) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_continuous(breaks = unique(dep$n_age), expand = c(0, 0)) +
    ggplot2::labs(x = "年", y = "年齢", fill = lab,
                  caption = "推定資源量 N の各セルについて、データを ±h 動かしたときの弾力性から計算 (デルタ法)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 6),
                   plot.caption = ggplot2::element_text(size = 7, colour = "grey30"))
  # 小さい = 薄い、大きい = 濃い (白 -> 薄青 -> 青 -> 濃青)
  p <- p + ggplot2::scale_fill_gradientn(
    colours = c("#ffffff", "#c6dbef", "#4292c6", "#08306b"),
    limits  = if (value == "share_index") c(0, 1) else c(0, NA))
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet]]), ncol = 1)
  p
}


#' 推定資源量の 1 セルが、どのデータセルにどれだけ依存しているかのヒートマップ
#'
#' @param infl   [data_influence()] の返り値
#' @param n_age  対象の年齢 (0 始まり)
#' @param n_year 対象の年 (列名)。`"last"` なら終端年
#' @return ggplot オブジェクト (塗り = 弾力性 d log N / d log データ)
#' @export
plot_data_influence <- function(infl, n_age = 0, n_year = "last") {
  e <- infl$elasticity
  if (identical(n_year, "last")) n_year <- utils::tail(unique(e$n_year), 1)
  d <- e[e$n_age == n_age & e$n_year == as.character(n_year), ]
  d$source <- factor(d$source, levels = c("caa", "index"), labels = c("CAA (行 = 年齢)", "指数 (行 = 指数の番号)"))
  lim <- max(abs(d$elasticity), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(factor(.data$data_year, levels = unique(.data$data_year)),
                                  .data$data_row, fill = .data$elasticity)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                                  midpoint = 0, limits = c(-lim, lim)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$source), ncol = 1) +
    ggplot2::labs(x = "データの年", y = NULL, fill = "弾力性",
                  title = sprintf("N(年齢 %s, %s 年) のデータ依存", n_age, n_year)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 6))
}
