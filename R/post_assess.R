#' Find a project root by walking up parent directories
#'
#' This helper starts from a directory and walks upward until it finds a
#' directory that contains all requested marker files or directories. It is used
#' by post-assessment scripts to locate a project root without depending on a
#' fixed working directory.
#'
#' @param path Character scalar. Directory to start searching from. Defaults to
#'   the current working directory.
#' @param markers Character vector. File or directory names that must exist
#'   directly under the project root.
#' @param winslash Character scalar passed to [normalizePath()] on Windows.
#' @param mustWork Logical. If `TRUE`, throw an error when no matching parent
#'   directory is found. If `FALSE`, return `NA_character_`.
#'
#' @return A normalized path to the first parent directory containing all
#'   `markers`, or `NA_character_` when `mustWork = FALSE` and no root is found.
#' @export
#'
#' @examples
#' \dontrun{
#' find_project_root(markers = "future_outputs")
#' }
find_project_root <- function(path = getwd(),
                              markers = "future_outputs",
                              winslash = "/",
                              mustWork = TRUE) {
  if (!is.character(path) || length(path) != 1 || is.na(path)) {
    stop("`path` must be a non-missing character scalar.", call. = FALSE)
  }
  if (!is.character(markers) || length(markers) == 0 || any(is.na(markers))) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }

  path <- normalizePath(path, winslash = winslash, mustWork = TRUE)
  repeat {
    if (all(file.exists(file.path(path, markers)))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      if (isTRUE(mustWork)) {
        stop(
          "Could not find a project root containing: ",
          paste(markers, collapse = ", "),
          call. = FALSE
        )
      }
      return(NA_character_)
    }
    path <- parent
  }
}

#' Convert an age-year matrix to a long tibble
#'
#' This helper converts a matrix whose rows are ages and columns are years into
#' a long table with `age`, `year`, and `value` columns. It is useful for
#' post-assessment summaries of age-specific quantities such as numbers at age,
#' fishing mortality, weight, and maturity.
#'
#' @param mat Matrix-like object. Rows are ages and columns are years.
#' @param stat Optional character scalar added as a `stat` column.
#' @param age_names Optional age labels. Defaults to `rownames(mat)`, or
#'   zero-based row indices when row names are missing.
#' @param year_names Optional year labels. Defaults to `colnames(mat)`, or
#'   one-based column indices when column names are missing.
#'
#' @return A tibble with columns `age`, `year`, `value`, and optionally `stat`.
#' @export
as_age_year_tibble <- function(mat,
                               stat = NULL,
                               age_names = NULL,
                               year_names = NULL) {
  mat <- as.matrix(mat)
  if (length(dim(mat)) != 2) {
    stop("`mat` must be a matrix-like object with two dimensions.", call. = FALSE)
  }

  if (is.null(age_names)) {
    age_names <- rownames(mat)
  }
  if (is.null(year_names)) {
    year_names <- colnames(mat)
  }
  if (is.null(age_names) || any(!nzchar(age_names))) {
    age_names <- seq_len(nrow(mat)) - 1
  }
  if (is.null(year_names) || any(!nzchar(year_names))) {
    year_names <- seq_len(ncol(mat))
  }
  if (length(age_names) != nrow(mat)) {
    stop("`age_names` must have the same length as `nrow(mat)`.", call. = FALSE)
  }
  if (length(year_names) != ncol(mat)) {
    stop("`year_names` must have the same length as `ncol(mat)`.", call. = FALSE)
  }

  out <- tibble::tibble(
    age = rep(age_names, times = ncol(mat)),
    year = rep(year_names, each = nrow(mat)),
    value = as.vector(mat)
  ) |>
    dplyr::mutate(
      age = as.numeric(age),
      year = as.numeric(year)
    )

  if (!is.null(stat)) {
    out <- dplyr::mutate(out, stat = stat)
  }
  out
}

#' Predict recruitment from stock-recruitment parameters
#'
#' This helper evaluates stock-recruitment relationships using the parameter
#' coding stored in `SR_mat`. Currently, Hockey-stick (`SR_type = 1`),
#' Beverton-Holt (`SR_type = 2`), and Ricker (`SR_type = 3`) are supported.
#'
#' @param ssb Numeric vector of spawning stock biomass.
#' @param sr_type Numeric vector or scalar. `1` = HS, `2` = BH, `3` = RI.
#' @param a,b Numeric vectors or scalars of stock-recruitment parameters.
#'
#' @return Numeric vector of predicted recruitment.
#' @export
predict_recruitment_from_sr_pars <- function(ssb, sr_type, a, b) {
  common_length <- max(length(ssb), length(sr_type), length(a), length(b))
  ssb <- rep(ssb, length.out = common_length)
  sr_type <- rep(sr_type, length.out = common_length)
  a <- rep(a, length.out = common_length)
  b <- rep(b, length.out = common_length)

  out <- rep(NA_real_, common_length)
  is_hs <- sr_type == 1
  is_bh <- sr_type == 2
  is_ri <- sr_type == 3
  supported <- is_hs | is_bh | is_ri
  if (any(is.na(supported)) || any(!supported)) {
    stop("Only HS (1), BH (2), and RI (3) stock-recruitment relationships are implemented.", call. = FALSE)
  }

  out[is_hs] <- ifelse(ssb[is_hs] > b[is_hs], b[is_hs] * a[is_hs], ssb[is_hs] * a[is_hs])
  out[is_bh] <- a[is_bh] * ssb[is_bh] / (1 + b[is_bh] * ssb[is_bh])
  out[is_ri] <- a[is_ri] * ssb[is_ri] * exp(-b[is_ri] * ssb[is_ri])
  out
}

#' Calculate effective recruitment multipliers from SR_mat
#'
#' This helper compares realized recruitment stored in `SR_mat[, , "recruit"]`
#' with recruitment predicted from `SR_type`, `a`, `b`, and `ssb`. It is useful
#' for approximating a stochastic future projection with a deterministic mean
#' trajectory while preserving the effective recruitment level used in the
#' original projection.
#'
#' @param future_obj A future result object containing `SR_mat`.
#' @param years Numeric or character vector of years to summarize.
#'
#' @return A tibble with one row per year and columns `year`,
#'   `mean_predicted_recruitment`, `mean_recruitment`, and
#'   `effective_recruitment_multiplier`.
#' @export
get_sr_recruitment_multiplier <- function(future_obj, years) {
  year_names <- as.character(years)
  sr_mat <- future_obj$SR_mat
  required_pars <- c("SR_type", "a", "b", "ssb", "recruit")
  if (is.null(sr_mat) || !all(required_pars %in% dimnames(sr_mat)[[3]])) {
    stop("`future_obj$SR_mat` must contain SR_type, a, b, ssb, and recruit.", call. = FALSE)
  }
  if (!all(year_names %in% dimnames(sr_mat)[[1]])) {
    stop("All `years` must be present in `future_obj$SR_mat`.", call. = FALSE)
  }

  purrr::map_dfr(year_names, function(year_name) {
    sr_year <- sr_mat[year_name, , , drop = FALSE][1, , ]
    pred_recruit <- predict_recruitment_from_sr_pars(
      ssb = as.numeric(sr_year[, "ssb"]),
      sr_type = as.numeric(sr_year[, "SR_type"]),
      a = as.numeric(sr_year[, "a"]),
      b = as.numeric(sr_year[, "b"])
    )
    recruit <- as.numeric(sr_year[, "recruit"])
    multiplier <- recruit / pred_recruit
    multiplier <- multiplier[is.finite(multiplier)]
    tibble::tibble(
      year = as.numeric(year_name),
      mean_predicted_recruitment = mean(pred_recruit, na.rm = TRUE),
      mean_recruitment = mean(recruit, na.rm = TRUE),
      effective_recruitment_multiplier = mean(multiplier, na.rm = TRUE)
    )
  })
}

#' Extract age-specific assessment components from future results
#'
#' This helper extracts age- and year-specific quantities from a list of future
#' projection results and returns them in a long table. It is intended for
#' post-assessment comparisons across assessment years.
#'
#' @param future.list A list of future result objects.
#' @param assessment.years Numeric or character vector of assessment years. Its
#'   length must match `future.list`.
#' @param component Character scalar. One of `"selectivity"`, `"F"`,
#'   `"weight"`, `"maturity"`, or `"number"`.
#' @param value_fun Function used to summarise across simulations. Defaults to
#'   [mean()].
#' @param estimate_label,prediction_label Labels used for years before and after
#'   the assessment year.
#'
#' @return A tibble with columns `age`, `year`, `value`, `assess_year`, `type`,
#'   and `Assess`.
#' @export
get_assessment_age_component_data <- function(future.list,
                                              assessment.years,
                                              component = c("selectivity", "F", "weight", "maturity", "number"),
                                              value_fun = mean,
                                              estimate_label = NULL,
                                              prediction_label = "Predicted") {
  component <- match.arg(component)
  if (length(future.list) != length(assessment.years)) {
    stop("`future.list` and `assessment.years` must have the same length.", call. = FALSE)
  }
  if (!is.function(value_fun)) {
    stop("`value_fun` must be a function.", call. = FALSE)
  }

  if (is.null(estimate_label)) {
    estimate_label <- switch(
      component,
      selectivity = "Estimate",
      F = "Estimate",
      number = "Estimated",
      weight = "Observed",
      maturity = "Observed"
    )
  }

  component_tables <- lapply(seq_along(future.list), function(i) {
    future_obj <- future.list[[i]]
    assess_year <- as.numeric(assessment.years[i])
    source_array <- switch(
      component,
      selectivity = future_obj$faa,
      F = future_obj$faa,
      weight = future_obj$waa,
      maturity = future_obj$input$tmb_data$maa_mat,
      number = future_obj$naa
    )
    mat <- switch(
      component,
      selectivity = {
        selectivity_array <- array(NA_real_, dim = dim(future_obj$faa), dimnames = dimnames(future_obj$faa))
        for (j in seq_len(dim(future_obj$faa)[3])) {
          max_f <- apply(future_obj$faa[, , j, drop = FALSE][, , 1], 2, max, na.rm = TRUE)
          selectivity_array[, , j] <- sweep(future_obj$faa[, , j], 2, max_f, FUN = "/")
        }
        apply(selectivity_array, c(1, 2), value_fun, na.rm = TRUE)
      },
      F = apply(future_obj$faa, c(1, 2), value_fun, na.rm = TRUE),
      weight = apply(future_obj$waa, c(1, 2), value_fun, na.rm = TRUE),
      maturity = {
        if (is.null(future_obj$input$tmb_data$maa_mat)) {
          stop("`maa_mat` was not found in a future object.", call. = FALSE)
        }
        apply(future_obj$input$tmb_data$maa_mat, c(1, 2), value_fun, na.rm = TRUE)
      },
      number = apply(future_obj$naa, c(1, 2), value_fun, na.rm = TRUE)
    )

    source_year_names <- dimnames(source_array)[[2]]
    if (is.null(source_year_names) || length(source_year_names) != ncol(mat)) {
      source_year_names <- NULL
    }

    as_age_year_tibble(mat, year_names = source_year_names) |>
      dplyr::mutate(
        assess_year = assess_year
      )
  })

  dplyr::bind_rows(component_tables) |>
    dplyr::mutate(
      type = ifelse(year < assess_year, estimate_label, prediction_label),
      Assess = factor(assess_year)
    )
}

#' Plot age-specific assessment components across assessment years
#'
#' @param future.list A list of future result objects.
#' @param assessment.years Numeric or character vector of assessment years.
#' @param component Character scalar. One of `"selectivity"`, `"F"`,
#'   `"weight"`, `"maturity"`, or `"number"`.
#' @param minyear,maxyear Optional numeric limits for years to draw.
#' @param reference Optional data frame with columns `age` and `value` to overlay
#'   as reference points.
#' @param reference_shape_col Optional column name in `reference` mapped to
#'   point shape.
#' @param reference_shape Shape used when `reference_shape_col = NULL`.
#' @param ncol Number of facet columns.
#' @param ylab Y-axis label. If `NULL`, a component-specific label is used.
#' @param y_scale Either `"identity"` or `"sqrt"`.
#' @param show_points Logical. If `TRUE`, points are added to component lines.
#' @param linetype_values Character vector passed to
#'   [ggplot2::scale_linetype_manual()].
#' @param base_size Base font size passed to [ggplot2::theme_bw()].
#' @param value_fun Function used to summarise across simulations.
#' @param estimate_label,prediction_label Labels used for years before and after
#'   the assessment year.
#'
#' @return A ggplot object. The long table used for plotting is attached as the
#'   `"plot_data"` attribute.
#' @export
plot_assessment_age_component <- function(future.list,
                                          assessment.years,
                                          component = c("selectivity", "F", "weight", "maturity", "number"),
                                          minyear = NULL,
                                          maxyear = NULL,
                                          reference = NULL,
                                          reference_shape_col = NULL,
                                          reference_shape = 4,
                                          ncol = 3,
                                          ylab = NULL,
                                          y_scale = c("identity", "sqrt"),
                                          show_points = component != "selectivity",
                                          linetype_values = c("dashed", "solid"),
                                          base_size = 11,
                                          value_fun = mean,
                                          estimate_label = NULL,
                                          prediction_label = "Predicted") {
  component <- match.arg(component)
  y_scale <- match.arg(y_scale)

  plot_data <- get_assessment_age_component_data(
    future.list = future.list,
    assessment.years = assessment.years,
    component = component,
    value_fun = value_fun,
    estimate_label = estimate_label,
    prediction_label = prediction_label
  )

  if (!is.null(minyear)) {
    plot_data <- dplyr::filter(plot_data, year >= minyear)
  }
  if (!is.null(maxyear)) {
    plot_data <- dplyr::filter(plot_data, year <= maxyear)
  }

  if (is.null(ylab)) {
    ylab <- switch(
      component,
      selectivity = "selectivity",
      F = "F at age",
      weight = "weight at age",
      maturity = "maturity at age",
      number = "number at age"
    )
  }

  g <- ggplot2::ggplot(plot_data, ggplot2::aes(x = age, y = value)) +
    ggplot2::geom_path(ggplot2::aes(colour = Assess, linetype = type), linewidth = 0.5) +
    ggplot2::facet_wrap(ggplot2::vars(year), ncol = ncol) +
    ggplot2::scale_linetype_manual(values = linetype_values) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::ylab(ylab)

  if (isTRUE(show_points)) {
    g <- g +
      ggplot2::geom_point(ggplot2::aes(colour = Assess, shape = type), linewidth = 0.5) +
      ggplot2::scale_shape_manual(values = c(16, NA))
  }

  if (!is.null(reference)) {
    if (!all(c("age", "value") %in% names(reference))) {
      stop("`reference` must contain `age` and `value` columns.", call. = FALSE)
    }
    if (!is.null(reference_shape_col)) {
      reference$.reference_shape <- reference[[reference_shape_col]]
      g <- g + ggplot2::geom_point(
        data = reference,
        mapping = ggplot2::aes(x = age, y = value, shape = .reference_shape),
        inherit.aes = FALSE
      )
    } else {
      g <- g + ggplot2::geom_point(
        data = reference,
        ggplot2::aes(x = age, y = value),
        inherit.aes = FALSE,
        shape = reference_shape
      )
    }
  }

  if (y_scale == "sqrt") {
    g <- g + ggplot2::scale_y_sqrt(limits = c(0, NA))
  } else {
    g <- g + ggplot2::ylim(0, NA)
  }

  attr(g, "plot_data") <- plot_data
  g
}
