#'
#' create future projection data based on bootstrap results
#'
#' @export

make_future_data_boot <- function(res_vpa, vpa_boot, SR_boot, seed=1, ...){
  future_data_list <- list()
  for(i in 1:length(vpa_boot)){
    vpa_boot[[i]]$input <- res_vpa$input
    future_data_list[[i]] <- make_future_data(res_vpa=vpa_boot[[i]], res_SR=SR_boot[[i]],
                                              seed_number=seed+i, ...)
  }
  future_data <- unlist_future_data(future_data_list)
  
}
