
# correlation among parameters about SR
#' 再生産関係の推定パラメータの相関を出力する関数
#'
#' @inheritParams fit.SR
#' @inheritParams fit.SRregime
#' @param resSR \code{fit.SR}または\code{fit.SRregime}のオブジェクト
#' @return 以下の要素からなるリスト
#' \describe{
#' \item{\code{hessian}}{ヘッセ行列}
#' \item{\code{cov}}{推定されたパラメータの分散共分散行列}
#' \item{\code{cor}}{推定されたパラメータの相関行列}
#' }
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                 out.AR = TRUE)
#' corRes = corSR(resSR)
#' corRes$cor
#' }
#' @encoding UTF-8
#' @export
corSR = function(resSR){
  if (!resSR$input$hessian) {
    resSR$input$hessian <- TRUE
    resSR$input$p0 = resSR$opt$par
    if (class(resSR) == "fit.SR") resSR = do.call(fit.SR, resSR$input)
    if (class(resSR) == "fit.SRregime") resSR = do.call(fit.SRregime, resSR$input)
  }
  hessian = resSR$opt$hessian
  cov = solve(hessian)
  cor = stats::cov2cor(cov)
  return (list(hessian=hessian,cov=cov,cor=cor))
}

# steepness
#' Steepness (h) と関連するパラメータ (SB0,R0,B0)を計算する関数
#'
#' @param SR "HS", "BH", "RI"のいずれか
#' @param rec_pars 再生産関係のパラメータで\code{rec_pars$a},\code{rec_pars$b}で使用する
#' @param M 年齢別自然死亡係数 (ベクトルで与えるか、年齢共通の場合\code{M=0.4}のようにしてもよい)
#' @param waa （親魚量の）年齢別体重
#' @param maa 年齢別親魚量
#' @param plus_group 最高齢がプラスグループかどうか
#' @return 以下の要素からなるデータフレーム
#' \describe{
#' \item{\code{SPR0}}{F=0のときのSPR(この逆数がreplacement lineの傾き)}
#' \item{\code{SB0}}{F=0のときの親魚量}
#' \item{\code{R0}}{F=0のときの加入量}
#' \item{\code{B0}}{F=0のときの資源量}
#' \item{\code{h}}{steepness: BHかRIのときは0.2×SB0のときの加入量がh×R0, HSのときはh=1-b/SB0}}
#' }
#' @examples
#' \dontrun{
#' data(res_vpa)
#' SRdata <- get.SRdata(res_vpa)
#' resSR <- fit.SR(SRdata, SR = c("HS","BH","RI")[1],
#'                 method = c("L1","L2")[2], AR = 1,
#'                out.AR = TRUE)
#' rec_pars = resSR$pars
#' year <- "2017"
#' M = res_vpa$input$dat$M[,year]
#' waa = res_vpa$input$dat$waa[,year]
#' maa = res_vpa$input$dat$maa[,year]
#' Res_h = calc_steepness(SR=SR,rec_pars=rec_pars,M=M,waa=waa,maa=maa,plus_group=TRUE)
#' Res_h
#' }
#' @encoding UTF-8
#' @export
calc_steepness = function(SR="HS",rec_pars,M,waa,maa,plus_group=TRUE) {
  if (length(M)==1) {
    M = rep(M,length(waa))
  }
  if (length(waa) != length(maa) || length(M) != length(maa)) {
    stop("The lengths of 'waa' and 'maa' must be equal")
  }
  NAA0 = 1
  for (i in 1:(length(waa)-1)) {
    NAA0 = c(NAA0,rev(NAA0)[1]*exp(-M[i]))
  }
  if (plus_group) NAA0[length(NAA0)] = rev(NAA0)[1]/(1-exp(-1*rev(M)[1]))
  BAA0 = NAA0*waa
  SSB0 = BAA0*maa
  SPR0 = sum(SSB0) #get.SRRと一致 (testに使える)

  # 再生産関係とy=(1/SPR0)*xの交点を求める
  rec_a = rec_pars$a
  rec_b = rec_pars$b
  if (SR == "HS") {
    R0 = rec_pars$a * rec_pars$b
    SB0 = R0*SPR0
    if (SB0<rec_b) {
      warning("Virgin equilibrium does not exist!")
    }
    h = (SB0-rec_b)/SB0
  }
  if (SR == "BH") {
    SB0 = (rec_a*SPR0-1)/rec_b
    R0 = SB0/SPR0
    h = (rec_a*0.2*SB0/(1+rec_b*0.2*SB0))/R0
  }
  if (SR == "RI") {
    SB0 = (1/rec_b)*log(rec_a*SPR0)
    R0 = SB0/SPR0
    h = (rec_a*0.2*SB0*exp(-rec_b*0.2*SB0))/R0
  }

  B0 = sum(R0*BAA0)
  Res = data.frame(SPR0 = SPR0, SB0 = SB0, R0 = R0, B0 = B0, h = h)
  return(Res)
}
