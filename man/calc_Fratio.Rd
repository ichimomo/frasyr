% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.r
\encoding{UTF-8}
\name{calc_Fratio}
\alias{calc_Fratio}
\title{calculate F/Ftarget based on F_\%SPR multiplier}
\usage{
calc_Fratio(
  faa,
  waa,
  maa,
  M,
  SPRtarget = 30,
  waa.catch = NULL,
  Pope = TRUE,
  return_SPR = FALSE,
  plus_group = TRUE,
  max.age = Inf
)
}
\arguments{
\item{faa}{F at age}

\item{waa}{weight at age}

\item{maa}{maturity at age}

\item{M}{natural morality at age}

\item{SPRtarget}{target SPR (NULLの場合には最適化しない)}

\item{return_SPR}{return SPR as well as Fratio}

\item{plus_group}{プラスグループを考慮するかどうか}

\item{max.age}{SPR計算を打ち切る最大の年。デフォルトはInf

もともとのF at ageの最大がexp(-7)よりも小さい場合にはFratio=0となる。一方で、F at ageをすごく大きくしても指定されたSPRを実現できないような場合のFratioの上限値を50とする。}
}
\description{
calculate F/Ftarget based on F_\%SPR multiplier
}
