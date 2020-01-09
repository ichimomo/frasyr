
#' Ridge VPAの計算をTMBで実行するためにcppファイルのコンパイル等をする関数
#' 
#' @param TmbFile Cppファイルの名前
#' @param CppDir Cppファイルが格納されているディレクトリ
#' @param RunDir 実行するディレクトリ
#' @param overwrite Cppファイルをコピーするさい、すでに同じ名前のファイルがあるときにファイルを上書きするかどうか
#' 
#' @encoding UTF-8
#' 
#' @examples
#' \dontrun{
#' use_rvpa_tmb()
#' }
#' 
#' @export

use_rvpa_tmb <- function(TmbFile = "rvpa_tmb",
                         CppDir = system.file("executable",package="frasyr"),
                         RunDir = getwd(),
                         overwrite = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Please install TMB package!")
  }
file.copy( from=paste0(CppDir,"/",TmbFile,".cpp"), to=paste0(RunDir,"/",TmbFile,".cpp"), overwrite=overwrite)
  TMB::compile( paste0(TmbFile,".cpp") )
  dyn.load(TMB::dynlib(TmbFile))
}
