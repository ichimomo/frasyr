#' @export
load_data <- function(fname) {
  request_from_developers <- stringr::str_detect(fname, "/inst/extdata")
  if (request_from_developers) {
    fname <- system.file("extdata", extract_filename(fname), package = "frasyr")
  }
  ftype <- return_file_type(fname)
  if (ftype == "csv") {
    return(read.vpa(fname))
  } else if (ftype == "rda") {
    get(load(fname))
  } else {
    stop("Unknown filetype", call. = TRUE)
  }
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
