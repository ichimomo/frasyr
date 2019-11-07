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
