make_stock_table <- function(result_vpa, result_msy,
                             yr_pre_abc, unit = "百トン") {
  return_ <- function() {
    rbind(recent_five_years_(),
          future_start_year_())
  }

  yr_newest_recent    <- yr_pre_abc[length(yr_pre_abc)]
  yr_oldest_recent    <- yr_pre_abc[1]
  recent_five_years_ <- function() {
    data.frame(Year     = yr_oldest_recent:yr_newest_recent,
               Biomass  = get_x_from_vpa_("biomass"),
               SSB      = get_x_from_vpa_("SSB"),
               Catch    = get_x_from_vpa_("catch"),
               `F/Fmsy` = f_per_fmsy_()) %>%
      dplyr::mutate(HarvestRate = round(Catch / Biomass * 100, 0))
  }
  future_start_year_ <- function() {
    # not implemented
    data.frame(Year = yr_newest_recent + 1,
               Biomass = NA,
               SSB     = NA,
               Catch   = NA,
               `F/Fmsy` = NA,
               HarvestRate = NA)
  }

  get_x_from_vpa_ <- function(x) {
    tbl <- result_vpa %>%
      convert_vpa_tibble() %>%
      dplyr::filter(dplyr::between(year, yr_oldest_recent, yr_newest_recent))
    if (x == "fishing_mortality") {
      tbl %>%
        dplyr::filter(stat == x) %>%
        dplyr::group_by(year) %>%
        dplyr::summarize(output = mean(value)) %>%
        dplyr::pull(output)
    } else {
      tbl %>%
        dplyr::filter(stat == x) %>%
        dplyr::group_by(year) %>%
        dplyr::summarize(output = convert_unit(value, to = unit)) %>%
        dplyr::pull(output) %>%
        round(0)
    }
  }
  f_per_fmsy_ <- function() {
    f    <- get_x_from_vpa_("fishing_mortality")
    fmsy <-  extract_fmsy(result_msy, mean = TRUE)
    force(round(f / fmsy, 2))
  }

  return_()
}
