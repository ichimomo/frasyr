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
    assertthat::has_name(result$input$res_SR$pars, "sd"),
    assertthat::has_name(result$input, "model_average_option")
  )

  retrieved <- result$input

  if (is.null(new_sd)) return(retrieved)

  if (is.numeric(new_sd) == FALSE) stop("'new_sd' should be numeric")

  retrieved$res_SR$pars$sd <- new_sd

  force(retrieved)
}
