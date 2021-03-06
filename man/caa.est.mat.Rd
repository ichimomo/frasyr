% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.r
\encoding{UTF-8}
\name{caa.est.mat}
\alias{caa.est.mat}
\title{年齢別生物パラメータとFと漁獲量を与えると与えた漁獲量と一致するFへの乗数を返す}
\usage{
caa.est.mat(
  naa,
  saa,
  waa,
  M,
  catch.obs,
  Pope,
  set_max1 = TRUE,
  max_exploitation_rate = 0.99,
  max_F = exp(10)
)
}
\arguments{
\item{max_exploitation_rate}{潜在的に漁獲できる漁獲量＜入力した漁獲量の場合、潜在的に漁獲できる漁獲量の何％まで実際に漁獲するか}

\item{max_F}{F at ageの最大値となる値の上限をどこにおくか}
}
\description{
年齢別生物パラメータとFと漁獲量を与えると与えた漁獲量と一致するFへの乗数を返す
}
