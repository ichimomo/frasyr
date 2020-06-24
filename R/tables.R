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

make_table <- function(...) {
  UseMethod("make_table")
}

make_table.fit.SR <- function(result_sr) {
  data.frame(kankei     = result_sr$input$SR,
             saitekika  = result_sr$input$method,
             jikosoukan = result_sr$input$AR,
             result_sr$pars) %>%
    magrittr::set_colnames(
      c("再生産関係式", "最適化法", "自己相関", "a", "b", "S.D.", "rho")
    )
}

table1 <- function(...) {
  adhoc_table(number = "one", ...)
}


table2 <- function(...) {
  adhoc_table(number = "two", ...)
}

adhoc_table <- function(result_vpa, yrs_preabc, number, sbtarget = NULL, fmsy = NULL) {
  return_ <- function() {
    switch(number,
           "one" = {
             rbind(sb_latest_(),
                   f_latest_())
           },
           "two" = {
             rbind(pspr_latest_(),
                   pspr_recent_())
           },
           "four" = {
             rbind(sb_latest_over_target_(),
                   f_latest_over_msy_())
           })
  }
  oldest_recentyr <- yrs_preabc[1]
  yr_latest       <- yrs_preabc[length(yrs_preabc)]
  sb_latest_ <- function() {
    make_row(key     = paste0("SB", yr_latest),
             value   = extract_from_vpa_("ssb", mean_by = "year"),
             remarks = paste0(yr_latest, "年漁期の親魚量"))
  }
  f_latest_  <- function(numeric = FALSE) {
    value <- extract_from_vpa_("faa")
    if (numeric) {
      value <- mean(unlist(value))
    } else {
      value <- format_x_at_age(value)
    }
    make_row(key   = paste0("F", yr_latest),
             value = value)
  }
  pspr_latest_ <- function() {
    make_row(key     = paste0("%SPR", wrap_by_paren(yr_latest)),
             value   = perspr_from_latest_(to = yr_latest),
             remarks = paste0(yr_latest, "年漁期の%SPR"))
  }
  pspr_recent_ <- function() {
    years <- paste0(oldest_recentyr, "--", yr_latest)
    make_row(key     = paste0("%SPR", wrap_by_paren(years)),
             value   = perspr_from_latest_(to = oldest_recentyr),
             remarks = paste0("現状", wrap_by_paren(years, "年漁期"),
                               "の漁獲圧に対応する%SPR"))
  }
  sb_latest_over_target_ <- function() {
    make_row(key     = paste0(sb_latest_()$項目, "/ SBtarget(SBmsy)"),
             value   = sb_latest_()$値 / sbtarget,
             remarks = paste0("目標管理基準値（最大持続生産量を実現する親魚量）に対する",
                               yr_latest,
                               "年漁期の親魚量の比"))
  }
  f_latest_over_msy_ <- function() {
    if (length(fmsy) > 1) {
      fmsy <- mean(fmsy)
    }
    make_row(key     = paste0(f_latest_()$項目, "/ Fmsy"),
             value   = f_latest_(numeric = TRUE)$値 / fmsy,
             remarks = paste0("最大持続生産量を実現する漁獲圧に対する",
                              yr_latest,
                              "年漁期の漁獲圧の比"))
  }
  perspr_from_latest_ <- function(to) {
    purrr::map_dbl(yr_latest:to,
                   extract_x,
                   x       = "perSPR",
                   vpadata = result_vpa) %>%
      mean()
  }
  extract_from_vpa_ <- function(what, mean_by = NULL) {
    extract_x(vpadata = result_vpa, what, yr_latest, mean_by = mean_by)
  }
  return_()
}

make_row <- function(key, value, remarks = "") {
  if (length(value) > 1) stop("Only one-length vector is allowed")

  rm_rownames_ <- function(x) {
    row.names(x) <- NULL
    force(x)
  }

  data.frame("項目" = key,
             "値"   = value,
             "備考" = remarks,
             stringsAsFactors = FALSE) %>%
    rm_rownames_()
}

format_x_at_age <- function(df, round = 2) {
  xvec <- unlist(df) %>% round(round)
  ages <- (seq_len(length(xvec)) - 1) %>%
    add_age_suffix()
  force(paste0(wrap_by_paren(ages), " = ", wrap_by_paren(xvec)))
}
