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

#' Draw fitted SR line over observed stock-recruitment plots
#'
#' @param fitted Fitted SR data created by fit.SR()
#' @param observed Observed data in format data.frame(year, SSB, R)
#' @param show.year If TRUE, show year label on each points
draw_SRline_ <- function(fitted, observed, show.year) {
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
