generate_dummy_future_data <- function(result_vpa) {
  years  <- extract_year_from(vpadata = result_vpa)
  lastyr <- tail(years, 1)
  make_future_data(result_vpa,
                   nsim = 1000,
                   nyear = 50,
                   future_initial_year_name = lastyr,
                   start_F_year_name = lastyr + 1,
                   start_biopar_year_name = lastyr + 1,
                   start_random_rec_year_name = lastyr + 1,
                   waa_year = (lastyr - 2):lastyr,
                   waa = NULL,
                   waa_catch_year = (lastyr - 2):lastyr,
                   waa_catch = NULL,
                   maa_year = (lastyr - 2):lastyr,
                   maa = NULL,
                   M_year = (lastyr - 2):lastyr,
                   M = NULL,
                   faa_year = (lastyr - 2):lastyr,
                   currentF = NULL,
                   futureF = NULL,
                   start_ABC_year_name = lastyr + 2,
                   HCR_beta = 1,
                   HCR_Blimit = -1,
                   HCR_Bban = -1,
                   HCR_year_lag = 0,
                   res_SR = res_sr_HSL2,
                   seed_number = 1,
                   resid_type = "lognormal",
                   resample_year_range = 0,
                   bias_correction = TRUE,
                   recruit_intercept = 0,
                   Pope = result_vpa$input$Pope)
}

generate_dummy_future_new_object <- function() {
  data(res_vpa)
  data(res_sr_HSL2)
  dummy_yr    <- 2015:2017
  future_data <- make_future_data(res_vpa,
                   res_SR = res_sr_HSL2,
                   M_year = dummy_yr,
                   maa_year = dummy_yr,
                   waa_catch_year = dummy_yr,
                   waa_year = dummy_yr)
  force(future_vpa(future_data$data))
}


expect_df <- function(tbl) {
  expect_is(tbl, "data.frame")
}
