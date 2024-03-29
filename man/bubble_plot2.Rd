% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvpa.r
\encoding{UTF-8}
\name{bubble_plot2}
\alias{bubble_plot2}
\title{VPAの結果を受けて，バブルプロット(横軸が年，縦軸が年齢，バブルの大きさと色で値を示す）を生成する}
\usage{
bubble_plot2(
  res,
  target = "faa",
  years = NULL,
  fix_ratio = 2,
  max_size = 10,
  legend_position = "bottom"
)
}
\arguments{
\item{res}{VPAの出力}

\item{target}{プロットしたい対象　例）"faa"や"naa"}

\item{years}{プロットしたい年数　デフォルトはNULLで全年．例）１９９７：１９９９とすると１９９７～１９９９年のみプロット}

\item{fix_ratio}{x軸とy軸のスケールの比．ここが２だと１：１になる}

\item{max_size}{バブルの最大の大きさ}

\item{legend_position}{legendの位置（see `theme_SH`）}
}
\description{
VPAの結果を受けて，バブルプロット(横軸が年，縦軸が年齢，バブルの大きさと色で値を示す）を生成する
}
