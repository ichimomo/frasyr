% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics_vpa.R
\encoding{UTF-8}
\name{plot_resboot_vpa}
\alias{plot_resboot_vpa}
\title{ブートストラップ法による信頼区間の算出と作図を自動で行う関数}
\usage{
plot_resboot_vpa(res, B_ite = 1000, B_method = "p", ci_range = 0.95)
}
\arguments{
\item{res}{VPAの結果のオブジェクト}

\item{B_ite}{ブートストラップ計算の数。デフォルトで1000。}

\item{B_method}{ブートストラップの方法。デフォルトではノンパラメトリックブートストラップ。}

\item{ci_range}{信頼区間の幅。デフォルトでは0.95（95％信頼区間）}
}
\value{
返ってくる値:
    \code{plot} 親魚重量、資源尾数、資源重量それぞれについて信頼区間のプロットが得られる。
    \code{res_boot} ブートストラップ法の結果が得られる。信頼区間の算出に利用可。
}
\description{
ブートストラップ法による信頼区間の算出と作図を自動で行う関数
}
\seealso{
ブートストラップ法について:  \code{\link{boo.vpa}}
https://ichimomo.github.io/frasyr/articles/Diagnostics-for-VPA.html
}
\author{
濵邉昂平, 市野川桃子
}
