% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.r
\encoding{UTF-8}
\name{caa.est.mat_wrong}
\alias{caa.est.mat_wrong}
\title{年齢別生物パラメータとFと漁獲量を与えると与えた漁獲量と一致するFへの乗数を返す (optimを使わないでやろうとしたけどだめだったやつ）}
\usage{
caa.est.mat_wrong(
  naa,
  saa,
  waa,
  M,
  catch.obs,
  Pope,
  max_exploitation_rate = 0.99,
  max_F = exp(10)
)
}
\arguments{
\item{naa}{numbers at age}

\item{saa}{selectivity at age}

\item{waa}{weight at age}

\item{M}{natural mortality at age}

\item{max_exploitation_rate}{潜在的に漁獲できる漁獲量＜入力した漁獲量の場合、潜在的に漁獲できる漁獲量の何％まで実際に漁獲するか}

\item{max_F}{F at ageの最大値となる値の上限をどこにおくか}
}
\description{
これだと漁獲量をぴったりに返すFのvectorは得られるが、もとの選択率に一致しない
}
