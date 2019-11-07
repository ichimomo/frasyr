#' Util functions for select SRmodels

#' Pull data frame from list
#'
#' @param list List object such as VPA result
#' @param dfname Name of data frame such as "baa" in VPA result
pull_df_from_list <- function(list, dfname) {
  list[dfname] %>%
    as.data.frame()
}
