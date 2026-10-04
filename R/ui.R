#' @export
load_data <- function(fname) {
  request_from_developers <- stringr::str_detect(fname, "/inst/extdata")
  if (request_from_developers) {
    fname <- system.file("extdata", extract_filename(fname), package = "frasyr")
  }
  
  res <- switch(return_file_type(fname),
                "csv" = read.vpa(fname),
                "rda" = get(load(fname)),
                "rds" = readRDS(fname),
                stop("Unknown filetype: ", fname, call. = FALSE))

  # 古いバージョンのVPAでresにcaaが入っていない場合に補完
  if(inherits(res,"vpa") || "tune" %in% names(res$input)){
    if(is.null(res$caa)) res$caa <- res$input$dat$caa   # VPAは漁獲量が誤差なしの入力
  }

  if(inherits(res,"sam")){
    # samでもwaa, wcaaはresに保持していないため追加
    if(is.null(res$wcaa)){    
      waa_catch_tmp <- res$input$dat$waa.catch
      if(is.null(waa_catch_tmp)){
        waa_catch_tmp <- res$input$dat$waa
      }
      res$wcaa <- as.data.frame(res$caa * waa_catch_tmp)
      warning("SAMのwcaa (年齢別漁獲重量)を予測caaと観測waaから計算しました")                      
    }
    # 将来予測に対応させるためPope設定を追加
    if(is.null(res$input$Pope)) res$input$Pope <- FALSE # SAMはBaranovの式(frasam R/sam.R:664)    
  }  

  return(res)
  
}

#' Retrieve function argument settings to reuse
#'
#' @param result objects created by make_future_data()
#' @export
retrieve_input <- function(result, new_sd = NULL) {
  assertthat::assert_that(
    assertthat::has_name(result, "input"),
    assertthat::has_name(result$input, "model_average_option")
  )

  retrieved <- result$input

  if (is.null(new_sd)) return(retrieved)

  if (is.numeric(new_sd) == FALSE) stop("'new_sd' should be numeric")

  if(class(retrieved$res_SR)=="fit.SRregime"){
    assertthat::has_name(result$input$res_SR$regime_pars, "sd")
    retrieved$res_SR$regime_pars$sd[] <- new_sd
  }
  
  if(class(retrieved$res_SR)=="fit.SR"){
    assertthat::has_name(result$input$res_SR$pars, "sd")
    retrieved$res_SR$pars$sd <- new_sd      
  }
  
  is_model_averaged <- class(retrieved$res_SR) == "list" && !is.null(retrieved$res_SR$input)
  if(is_model_averaged){
    for(i in 1:length(retrieved$res_SR)){
      retrieved$res_SR[[i]]$pars$sd <- new_sd
    }}

  force(retrieved)
}
