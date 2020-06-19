load_extdata <- function(fname) {
  invisible(get(load(system.file("extdata", fname, package = "frasyr"))))
}
return_file_type <- function(fname) {
  if (stringr::str_detect(fname, "csv$")) {
    "csv"
  } else if (stringr::str_detect(fname, "rda$")) {
    "rda"
  } else {
    stop("Unknown file type", fname)
  }
}

#' Extract value at age from vpadata
#'
#' @param vpadata Object return from \code(vpa)
#' @param x Character one of 'f', 'c', 'w', 'n', 's' or 'wc'
#' @param mean_by If non-null, return mean value:
#' \describe{
#'   \item{"age"}{Return values averaged by age}
#'   \item{"year"}{Return values averaged by year}
#' }
extract_xaa <- function(vpadata, x, year, mean_by = NULL) {
  xaa       <- paste0(x, "aa")
  vars      <- vpadata[[xaa]]
  extracted <- vars[colnames(vars) %in% as.character(year)]

  if (is.null(mean_by)) return(extracted)

  if (mean_by == "year") {
      colMeans(extracted)
  } else if (mean_by == "age") {
      rowMeans(extracted)
  } else {
      stop("'mean_by' sould be either 'age' or 'year'.")
  }
}

extract_year_from <- function(vpadata = NULL, msydata = NULL, var = "faa") {

  use_vpa <- !is.null(vpadata) & is.null(msydata)
  use_msy <- !is.null(msydata) & is.null(vpadata)

  if (use_vpa) {
      year <- colnames(vpadata[[var]])
  } else if (use_msy) {
      stop("Not implemented")
  } else {
      stop("Give me either 'vpadata' or 'msydata'")
  }
  as.numeric(year)
}

select_from_tail <- function(vec, relatives) {

  if (any(relatives >= 0)) stop("'relative' should be negative")

  sort(rev(vec)[-(relatives)])
}

extract_filename <- function(path) {
  stringr::str_extract(path, "(?<=/)(\\w|\\.)+\\.[a-z]+$")
}
