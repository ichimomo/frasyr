## #' Make table for stock assessment result 
## #'
## #' !this function is replaced by make_stock_table2 so that currently not used
## #' 
## #' @param result_vpa Obejct returned from \code{vpa()}
## #' @param result_msy Object MSY result created previous SC meeting...?
## #' @param result_future Object returnd from \code{future_vpa()}
## #' @param yr_future_start Year when future projection starts
## #' @export

## make_stock_table <- function(result_vpa, result_msy, result_future,
##                              yr_future_start, yr_abc, faa_pre_abc, unit = "百トン") {
##   return_ <- function() {
##     rbind(recent_four_years_(),
##           future_start_year_(),
##           abc_year_())
##   }

##   n_years_to_display  <- 6
##   yr_newest_recent    <- yr_abc - 2
##   yr_oldest_recent    <- yr_abc - (n_years_to_display - 1)

##   recent_four_years_ <- function() {

##     f_per_fmsy_ <- function(year) {
##       make_kobe_ratio(result_vpa, result_msy) %>%
##         dplyr::select(year, Fratio) %>%
##         dplyr::filter(year %in% as.character(yr_oldest_recent:yr_newest_recent)) %>%
##         dplyr::pull(Fratio) %>%
##         round(2)
##     }

##     pull_x_from_vpa_result_ <- function(x) {
##       result_vpa %>%
##         convert_vpa_tibble() %>%
##         dplyr::filter(dplyr::between(year, yr_oldest_recent, yr_newest_recent),
##                       stat == x) %>%
##         dplyr::group_by(year) %>%
##         dplyr::summarize(output = convert_unit(value, to = unit)) %>%
##         dplyr::pull(output) %>%
##         round(0)
##     }

##     data.frame(Year     = yr_oldest_recent:yr_newest_recent,
##                Biomass  = pull_x_from_vpa_result_(x = "biomass"),
##                cBiomass  = pull_x_from_vpa_result_(x = "cbiomass"),               
##                SSB      = pull_x_from_vpa_result_(x = "SSB"),
##                Catch    = pull_x_from_vpa_result_(x = "catch"),
##                `F/Fmsy` = f_per_fmsy_()) %>%
##       dplyr::mutate(HarvestRate = round(Catch / Biomass * 100, 0))
##   }

##   future_start_year_ <- function() {

##     fratio <- mean(faa_pre_abc) / extract_fmsy(result_msy, mean = TRUE)

##     data.frame(Year = yr_future_start,
##                Biomass  = calc_x_from_future_result_(x = "biomass", yr = yr_future_start),
##                cBiomass  = calc_x_from_future_result_(x = "cbiomass", yr = yr_future_start),               
##                SSB      = calc_x_from_future_result_(x = "ssb",     yr = yr_future_start),
##                Catch    = calc_x_from_future_result_(x = "catch",   yr = yr_future_start),
##                `F/Fmsy` = round(fratio, 2)) %>%
##       dplyr::mutate(HarvestRate = round(Catch / Biomass * 100, 0))
##   }

##   abc_year_ <- function() {
##     data.frame(Year = yr_abc,
##                Biomass  = calc_x_from_future_result_(x = "biomass", yr = yr_abc),
##                cBiomass  = calc_x_from_future_result_(x = "cbiomass", yr = yr_abc),               
##                SSB      = calc_x_from_future_result_(x = "ssb"    , yr = yr_abc),
##                Catch    = "-",
##                `F/Fmsy` = "-",
##                HarvestRate = "-")
##   }

##   calc_x_from_future_result_ <- function(x, yr) {
##     result_future %>%
##       extract_value.future_new(what = x,
##                                year = yr,
##                                unit = unit) %>%
##       dplyr::pull(average) %>%
##       round(0)
##   }


##   return_()
## }

#' Make table for stock assessment result 2
#'
#' @param result_future0.8 Object returnd from adopted future projection
#' @param yr_future_ABC Year when ABC is calculated
#' @param Fcurrent_per_Fmsy Fcurrent/Fmsy
#' @param Fcurrent_year_range
#' 
#' @export

make_stock_table2 <- function(result_future0.8,
                              yr_future_ABC,
                              Fcurrent_per_Fmsy,
                              Fcurrent_year_range) {

  years_to_display  <- (yr_future_ABC-5):yr_future_ABC

  if(!"cbiomass" %in% names(result_future0.8$summary)) result_future0.8$summary$cbiomass <- apply(result_future0.8$waa_catch_mat * result_future0.8$naa, c(2,3), sum) %>% apply(1,mean)

  stock_table <- result_future0.8$summary %>%
    dplyr::filter(year %in% years_to_display) %>%
    select(year, biomass, cbiomass, SSB, catch, Fratio) %>%
    mutate(harvest_rate=catch/cbiomass) %>%
    rename("F_year/F_msy"=Fratio) %>%
    mutate(year_label=case_when(year==yr_future_ABC ~ "ABC year",
                                year%in%Fcurrent_year_range ~ "Fcurrent year")) %>%
    mutate(year=as.character(year)) 

  stock_table <- bind_rows(stock_table,
                           tibble(year=str_c(range(Fcurrent_year_range), collapse="-"),
                                  "F_year/F_msy"=Fcurrent_per_Fmsy))

  return(stock_table)
}

#' Make table depending on the class of given object
#'
#' @inheritParams make_table.fit.SR
#' @export
make_table <- function(...) {
  UseMethod("make_table")
}

#' @param result_sr Return of \code{fit.SR}
#' @export
make_table.fit.SR <- function(result_sr) {
  data.frame(kankei     = result_sr$input$SR,
             saitekika  = result_sr$input$method,
             jikosoukan = result_sr$input$AR,
             result_sr$pars) %>%
    magrittr::set_colnames(
      c("再生産関係式", "最適化法", "自己相関", "a", "b", "S.D.", "rho")
    )
}

#' @inheritParams make_table.fit.SR Return of \code{fit.SR}
#' @export
make_table.fit.SRregime <- function(result_sr) {
  data.frame(kankei     = result_sr$input$SR,
             saitekika  = result_sr$input$method,
             result_sr$pars) %>%
    magrittr::set_colnames(
      c("再生産関係式", "最適化法", "a", "b", "S.D.")
    )
}

#' Make table of latest SB and F
#'
#' @inheritParams adhoc_table
#' @export
table1 <- function(...) {
  adhoc_table(number = "one", ...)
}

#' Make tables of \%SPRs
#' @inheritParams adhoc_table
#' @export
table2 <- function(...) {
  adhoc_table(number = "two", ...)
}

#' Make tables of latest SB and F relative to MSY values
#' @inheritParams adhoc_table
#' @export
table4 <- function(...) {
  adhoc_table(number = "four", ...)
}

#' Dispacher function for making table
#'
#' @inheritParams make_stock_table
#' @param yrs_preabc Vector of future years preceding ABC year
#' @param number Table number
#' @param sbtarget Value of SB target
#' @param fmsy Value of Fmsy
adhoc_table <- function(result_vpa, yrs_preabc, number, sbtarget = NULL, fmsy = NULL,
                        data_future = NULL, yr_biopar = NULL, result_msy = NULL) {
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
  yr_latest <- tail(extract_year_from(result_vpa), 1)
  sb_latest_ <- function() {
    make_row(key     = paste0("SB", yr_latest),
             value   = colSums(extract_from_vpa_("ssb"), na.rm = TRUE),
             remarks = paste0(yr_latest, "年漁期の親魚量"))
  }
  f_latest_  <- function(numeric = FALSE) {
    value <- extract_from_vpa_("faa")
    if (numeric) {
      value <- mean(unlist(value), na.rm = TRUE)
    } else {
      value <- format_x_at_age(value)
    }
    make_row(key   = paste0("F", yr_latest),
             value = value)
  }
  pspr_latest_ <- function() {
    make_row(key     = paste0("%SPR", wrap_by_paren(yr_latest)),
             value   = mean(extract_x(vpadata = result_vpa,
                                      x       = "perSPR",
                                      year    = yr_latest)),
             remarks = paste0(yr_latest, "年漁期の%SPR"))
  }
  pspr_recent_ <- function() {
    oldest_recentyr <- yrs_preabc[1]
    latest_recentyr <- yrs_preabc[length(yrs_preabc)]
    years <- paste0(oldest_recentyr, "--", latest_recentyr)
    pspr  <- calc_future_perSPR(
      res_vpa     = result_vpa,
      fout        = list(waa       = data_future$data$waa_mat,
                         maa       = data_future$data$maa_mat,
                         M         = data_future$data$M_mat,
                         waa.catch = data_future$data$waa_catch_mat),
      Fvector     = apply_year_colum(result_vpa$faa, latest_recentyr:oldest_recentyr),
      target.year = yr_biopar
    )
    make_row(key     = paste0("%SPR", wrap_by_paren(years)),
             value   = pspr * 100,
             remarks = paste0("現状", wrap_by_paren(years, "年漁期"),
                              "の漁獲圧に対応する%SPR"))
  }
  sb_latest_over_target_ <- function() {
    make_row(key     = paste0(sb_latest_()[["\u5024"]], "/ SBtarget(SBmsy)"),
             value   = sb_latest_()[["\u5024"]] / sbtarget,
             remarks = paste0("目標管理基準値（最大持続生産量を実現する親魚量）に対する",
                               yr_latest,
                               "年漁期の親魚量の比"))
  }
  f_latest_over_msy_ <- function() {
    value <- make_kobe_ratio(result_vpa, result_msy) %>%
      dplyr::filter(year == yr_latest) %>%
      dplyr::pull(Fratio)
    make_row(key     = paste0(f_latest_()[["\u5024"]], "/ Fmsy"),
             value   = value,
             remarks = paste0("最大持続生産量を実現する漁獲圧に対する",
                              yr_latest,
                              "年漁期の漁獲圧の比"))
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

#' Make summary table of MSY estimation
#'
#' @param result_msy Object of msy result
#' @export
make_msytable <- function(result_msy) {
  return_ <- function() {
    rbind(sbrows_(),
          msyrows_())
  }
  sbrows_ <- function() {
    rbind(
      make_sbrow_(name    = "target",
                  remarks = "目標管理基準値。最大持続生産量MSYを実現する親魚量（SBmsy）"),
      make_sbrow_(name    = "limit",
                  remarks = "限界管理基準値。MSYの60%の漁獲量が得られる親魚量（SB0.6msy）"),
      make_sbrow_(name    = "ban",
                  remarks = "禁漁水準。 MSYの10%の漁獲量が得られる親魚量（SB0.1msy）")
    )
  }
  make_sbrow_ <- function(name, remarks) {
    ssb <- derive_RP_value(result_msy$summary, paste0("B", name, "0"))$SSB
    make_row(key     = paste0("SB", name),
             value   = convert_unit(ssb, to ="千トン", round = 0, add_unit = TRUE),
             remarks = remarks)
  }
  msyrows_ <- function() {
    perspr <- derive_RP_value(result_msy$summary, "Btarget0")$perSPR
    msy    <- derive_RP_value(result_msy$summary, "Btarget0")$Catch
    rbind(
      make_row(key   = "Fmsy",
               value = extract_fmsy(result_msy) %>% format_x_at_age()),
      make_row(key     = "%SPR (Fmsy)",
               value   = paste0(round(perspr * 100, 1), "%"),
               remarks = "Fmsy に対応する %SPR"),
      make_row(key     = "MSY",
               value   =  convert_unit(msy, to = "千トン",
                                       round = 0,
                                       add_unit = TRUE),
               remarks = "最大持続生産量 MSY")
    )
  }
  return_()
}

#' Make summary table of ABC estimation
#'
#' @param kobe_table Return of \code{make_kobeII_table}
#' @param result_future Return of \code{future_vpa}
#' @param beta Value of beta to extract from Kobe table
#' @param year ABC year
#' @param faa_pre F at age of \code{yr_preabc}
#' @param faa_after F at age under HCR
#' @param yr_preabc Vector of future years preceding ABC year
#' @export
make_abctable <- function(kobe_table, result_future, beta, year, faa_pre, faa_after, yr_preabc) {
    return_ <- function() {
    df <- data.frame(abc_(),
                     ssb_mean_(),
                     f_over_recentf_(),
                     harvest_rate_())
    colnames(df) <- c(paste0(year, "年漁期のABC(千トン)"),
                      paste0(year, "年漁期の親魚量予測平均値(千トン)"),
                      paste0("現状の漁獲圧に対する比(F",year,"/F", yr_preabc[1], "--", yr_preabc[length(yr_preabc)], ")"),
                      paste0(year, "年漁期の漁獲割合（%）"))
    force(df)
  }

  abc_ <- function() {
    extract_from_kobe_table(kobe_table,
                            beta = beta,
                            what = "catch.mean",
                            year = year,
                            unit = "千トン")
  }
  ssb_mean_  <- function(){
    extract_from_kobe_table(kobe_table,
                            beta = beta,                              
                            what = "ssb.mean",
                            year = year,
                            unit = "千トン")
  }
  f_over_recentf_ <- function() {
    round(mean(faa_pre, na.rm = TRUE) / mean(faa_after), 2)
  }
  harvest_rate_ <- function() {
    biomass_abcyear <- extract_value(from = result_future,
                                     what = "biomass",
                                     year = year)[["average"]]
    round(abc_() / biomass_abcyear * 100, 0)
  }
  return_()
}

#' @export
summary_of_summary <- function(tbl_msy, tbl1, tbl2, tbl4) {
  characterize_valuecol_ <- function(tb) {
    tb[["\u5024"]] <- as.character(tb[["\u5024"]])
    force(tb)
  }
  extract_yr_from_tbl_ <- function() {
    force(stringr::str_extract(tbl1[["\u9805\u76ee"]][1], "[0-9]{4}"))
  }
  name1 <- "管理基準値とMSYに関係する値"
  name2 <- paste0(extract_yr_from_tbl_(), "漁期の親魚量と漁獲圧")
  name3 <- "MSYを実現する水準に対する比率"

  dplyr::bind_rows(
    make_row(key = name1, value = ""),
    characterize_valuecol_(tbl_msy),
    characterize_valuecol_(tbl1),
    make_row(key = name2, value = ""),
    characterize_valuecol_(tbl2),
    make_row(key = name3, value = ""),
    characterize_valuecol_(tbl4))
}


#' Export tables to csv
#'
#' @param to Name of file
#' @param ... Table objects
#'
#' \dontrun{
#' export_tables(to = "hoge.csv",
#'               table1, table2, table3)
#' }
#' @export
export_tables <- function(to, ...) {

  write_tables_to_csv_ <- function() {
    initialize_csv_()
    list2csv_()
  }

  initialize_csv_ <- function() {
    readr::write_excel_csv(data.frame(this_is_dummy = ""),
                           file      = to,
                           append    = FALSE,
                           col_names = FALSE)
  }

  list2csv_ <- function() {

    add_blank_line_ <- function(df) {
      blank <- ""
      suppressWarnings(
        rbind(dplyr::mutate_all(df, as.character()),
              blank,
              stringsAsFactors = FALSE)
      )
    }

    purrr::map(purrr::map(list(...), add_blank_line_),
               readr::write_excel_csv,
               file      = to,
               append    = TRUE,
               col_names = TRUE,
               delim     = ",") %>%
      invisible()
  }

  write_tables_to_csv_()
}
