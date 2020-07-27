#' Util functions for select SRmodels

#' Pull data frame from list
#'
#' @param list List object such as VPA result
#' @param dfname Name of data frame such as "baa" in VPA result
pull_df_from_list <- function(list, dfname) {
  list[dfname] %>%
    as.data.frame()
}

#' Get fitted SR data with model names
#'
#' @inheritParams fit.SR
#' @inheritParams get.SRdata
make_SR_dframe <- function(SR, method, vpares) {
  fit.SR(SRdata = get.SRdata(vpares = vpares), SR = SR, method = method) %>%
    pull_df_from_list(dfname = "pred") %>%
    dplyr::mutate(name = paste0(SR, "_", method))
}

#' Draw fitted SR line(s) over observed stock-recruitment plots
#'
#' @param fitted Fitted SR data created by fit.SR()
#' @param observed Observed data in format data.frame(year, SSB, R)
#' @param show.year If TRUE, show year label on each points
#' @export
draw_SRline <- function(fitted, observed, show.year) {
  label <- ""
  if (show.year == TRUE) label <- observed$year
  ggplot2::ggplot(data = observed,
                        ggplot2::aes(x = SSB,
                                     y = R,
                                     label = label)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel() +
    ggplot2::geom_line(data = fitted,
                       inherit.aes = FALSE,
                       ggplot2::aes(x = pred.SSB,
                                    y = pred.R,
                                    group = name,
                                    color = name))
}

#' Fit and draw SR line(s) over observed stock-recruitment plots
#'
#' @inheritParams get.SRdata
#' #' @inheritParams draw_SRline
#' @param SR Vector of multiple SR model described in \code{\link{fit.SR}}
#' @param method Vector of multiple method described in \code{\link{fit.SR}}
#' @examples
#' \dontrun{
#' draw_SRline(vpares = res_vpa,
#'             SR = c("HS", "RI", "BH"),
#'             method = c("L1", "L2"))
#' }
#' @export
fit_draw_SRline <- function(vpares, SR, method, show.year = FALSE) {
  rawdata <- get.SRdata(vpares = vpares, return.df = TRUE)
  models  <- expand.grid(model = SR,
                         method = method,
                         stringsAsFactors = FALSE)

  purrr::pmap_df(list(SR = models$model,
                      method = models$method),
                 make_SR_dframe,
                 vpares = vpares) %>%
    draw_SRline(observed = rawdata, show.year = show.year)
}
