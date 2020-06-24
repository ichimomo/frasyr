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
extract_x <- function(vpadata, x, year, mean_by = NULL) {
  if (x == "perSPR") {
    perspr <- get.SPR(vpadata)$ysdata[x]
    return(perspr[rownames(perspr) == year, ])
  }
  vars      <- vpadata[[x]]
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

add_age_suffix <- function(agevec) {
  ages <- paste0(agevec, "歳")
  ages[length(ages)] <- paste0(ages[length(ages)], "以上")
  ages
}

#' Give age names starts from zero and ends with plus group
#'
#' @param x object to name
#' give_agename(1:3)
#' give_agename(-1:2)
#' give_agename(c(1.234, 2.345))
give_agename <- function(x) {
  UseMethod("give_agename")
}

give_agename.numeric <- function(x) {
  ages <- 0:(length(x) - 1)
  names(x) <- add_age_suffix(ages)
  force(x)
}

wrap_by_paren <- function(..., collapse = ", ") {
  paste0("(", paste(..., sep = "", collapse = collapse), ")")
}

extract_fmsy <- function(result_msy, refpoint_name = "Btarget0", mean = FALSE) {
  # ad-hoc function
  assertthat::assert_that(
    assertthat::has_name(result_msy, "Fvector"),
    assertthat::has_name(result_msy, "summary"),
    assertthat::has_name(result_msy$summary, "RP.definition")
  )
  fmsy_by_age <- result_msy$Fvector %>%
    dplyr::slice(which(result_msy$summary$RP.definition == "Btarget0")) %>%
    as.numeric()

  if (mean) return(mean(fmsy_by_age))

  force(give_agename(fmsy_by_age))
}

convert_unit <- function(tons, to = "千トン", add_unit = FALSE, round = NULL) {
  value <- switch(
    to,
    "百トン" = tons / 100,
    "千トン" = tons / 1000,
    stop("'to' should be either '百トン' or '千トン'"))

  if (!is.null(round)) {
    value <- round(value, digits = round)
  }

  if (add_unit) return(paste(value, to))

  force(value)
}
