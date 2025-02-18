#'
#' create future projection data based on bootstrap results
#' 

make_future_data_boot <- function(vpa_boot, SR_boot, ...){
  future_data_list <- list()
  for(i in 1:length(vpa_boot)){
    future_data_list[[i]] <- make_future_data(res_vpa=vpa_boot[[i]], res_SR=SR_boot[[i]], ...)
  }
  future_data <- unlist_future_data(future_data_lise)
}
