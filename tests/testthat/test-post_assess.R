test_that("find_project_root finds a parent directory with the marker", {
  root <- tempfile("frasyr-root-")
  nested <- file.path(root, "a", "b", "c")
  dir.create(file.path(root, "future_outputs"), recursive = TRUE)
  dir.create(nested, recursive = TRUE)

  expect_equal(
    find_project_root(nested, markers = "future_outputs"),
    normalizePath(root, winslash = "/", mustWork = TRUE)
  )
})

test_that("find_project_root can require multiple markers", {
  root <- tempfile("frasyr-root-")
  nested <- file.path(root, "analysis")
  dir.create(file.path(root, "future_outputs"), recursive = TRUE)
  dir.create(file.path(root, "script"), recursive = TRUE)
  dir.create(nested, recursive = TRUE)

  expect_equal(
    find_project_root(nested, markers = c("future_outputs", "script")),
    normalizePath(root, winslash = "/", mustWork = TRUE)
  )
})

test_that("find_project_root returns NA or errors when no root is found", {
  root <- tempfile("frasyr-root-")
  dir.create(root, recursive = TRUE)

  expect_true(is.na(find_project_root(root, markers = "future_outputs", mustWork = FALSE)))
  expect_error(
    find_project_root(root, markers = "future_outputs", mustWork = TRUE),
    "Could not find a project root"
  )
})

test_that("find_project_root validates inputs", {
  expect_error(find_project_root(NA_character_), "`path`")
  expect_error(find_project_root(getwd(), markers = character()), "`markers`")
  expect_error(find_project_root(getwd(), markers = NA_character_), "`markers`")
})

test_that("as_age_year_tibble converts a matrix to long format", {
  mat <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(age = c("0", "1"), year = c("2001", "2002", "2003"))
  )

  out <- as_age_year_tibble(mat, stat = "Number")

  expect_s3_class(out, "tbl_df")
  expect_equal(names(out), c("age", "year", "value", "stat"))
  expect_equal(out$age, c(0, 1, 0, 1, 0, 1))
  expect_equal(out$year, c(2001, 2001, 2002, 2002, 2003, 2003))
  expect_equal(out$value, 1:6)
  expect_true(all(out$stat == "Number"))
})

test_that("as_age_year_tibble validates labels", {
  mat <- matrix(1:4, nrow = 2)

  expect_equal(as_age_year_tibble(mat)$age, c(0, 1, 0, 1))
  expect_equal(as_age_year_tibble(mat)$year, c(1, 1, 2, 2))
  expect_error(as_age_year_tibble(mat, age_names = "0"), "`age_names`")
  expect_error(as_age_year_tibble(mat, year_names = "2001"), "`year_names`")
})

test_that("predict_recruitment_from_sr_pars supports HS, BH, and RI", {
  ssb <- c(50, 200, 100)
  expect_equal(
    predict_recruitment_from_sr_pars(ssb, sr_type = 1, a = 2, b = 100),
    c(100, 200, 200)
  )
  expect_equal(
    predict_recruitment_from_sr_pars(100, sr_type = 2, a = 2, b = 0.01),
    100
  )
  expect_equal(
    predict_recruitment_from_sr_pars(100, sr_type = 3, a = 2, b = 0.01),
    200 * exp(-1)
  )
  expect_error(predict_recruitment_from_sr_pars(100, sr_type = 4, a = 2, b = 1), "Only HS")
})

test_that("get_sr_recruitment_multiplier summarizes SR_mat recruit", {
  sr_mat <- array(
    NA_real_,
    dim = c(2, 2, 5),
    dimnames = list(
      year = c("2019", "2020"),
      iter = c("1", "2"),
      par = c("SR_type", "a", "b", "ssb", "recruit")
    )
  )
  sr_mat[, , "SR_type"] <- 1
  sr_mat[, , "a"] <- 2
  sr_mat[, , "b"] <- 100
  sr_mat["2019", , "ssb"] <- c(50, 200)
  sr_mat["2020", , "ssb"] <- c(100, 120)
  sr_mat["2019", , "recruit"] <- c(150, 300)
  sr_mat["2020", , "recruit"] <- c(200, 400)

  out <- get_sr_recruitment_multiplier(list(SR_mat = sr_mat), years = 2019:2020)

  expect_s3_class(out, "tbl_df")
  expect_equal(out$year, c(2019, 2020))
  expect_equal(out$mean_predicted_recruitment, c(150, 200))
  expect_equal(out$mean_recruitment, c(225, 300))
  expect_equal(out$effective_recruitment_multiplier, c(mean(c(1.5, 1.5)), mean(c(1, 2))))
})

test_that("get_assessment_age_component_data extracts age components", {
  make_future <- function(multiplier = 1) {
    arr <- array(
      multiplier * seq_len(2 * 3 * 2),
      dim = c(2, 3, 2),
      dimnames = list(age = 0:1, year = 2001:2003, iter = 1:2)
    )
    list(
      naa = arr,
      faa = arr / 100,
      waa = arr / 10,
      input = list(tmb_data = list(
        maa_mat = arr / 1000,
        M_mat = arr / 1000
      ))
    )
  }

  out <- get_assessment_age_component_data(
    future.list = list(make_future(), make_future(2)),
    assessment.years = c(2002, 2003),
    component = "number"
  )

  expect_s3_class(out, "tbl_df")
  expect_true(all(c("age", "year", "value", "assess_year", "type", "Assess") %in% names(out)))
  expect_equal(sort(unique(out$assess_year)), c(2002, 2003))
  expect_true(all(out$type %in% c("Estimated", "Predicted")))
})

test_that("plot_assessment_age_component returns a ggplot with attached data", {
  arr <- array(
    seq_len(2 * 3 * 2),
    dim = c(2, 3, 2),
    dimnames = list(age = 0:1, year = 2001:2003, iter = 1:2)
  )
  future_obj <- list(
    naa = arr,
    faa = arr / 100,
    waa = arr / 10,
    input = list(tmb_data = list(maa_mat = arr / 1000))
  )

  g <- plot_assessment_age_component(
    future.list = list(future_obj),
    assessment.years = 2002,
    component = "F",
    minyear = 2001,
    maxyear = 2002
  )

  expect_s3_class(g, "ggplot")
  expect_true(all(attr(g, "plot_data")$year %in% 2001:2002))
})
