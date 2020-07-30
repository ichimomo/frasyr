context("Tools for developers")

test_that("load_extdata() enables simple loading", {
  msy <- load_extdata("res_MSY_pma.rda")

  expect_setequal(
    names(msy),
    c("summary", "summaryAR", "summary_tb", "F.msy", "all.stat", "all.statAR",
      "all.stat_tb", "trace", "input.list", "ssb.ar.mean", "SPR.msy")
  )

  vpa <- load_extdata("res_vpa_pma.rda")

  expect_setequal(
    names(vpa),
    c("input", "term.f", "np", "minimum", "minimum.c", "logLik", "gradient",
      "code", "q", "b", "sigma", "convergence", "message", "hessian", "Ft",
      "Fc.at.age", "Fc.mean", "Fc.max", "last.year", "Pope", "ssb.coef",
      "pred.index", "wcaa", "naa", "faa", "baa", "ssb", "saa")
  )
})


test_that("return_file_type() works", {
  expect_equal(return_file_type("foo.csv"), "csv")
  expect_equal(return_file_type("foo.rda"), "rda")

  expect_error(return_file_type("foo.csvv"), "Unknown file type")
  expect_error(return_file_type("foo.rdaa"), "Unknown file type")
})

context("Extract something from object")

vpadata <- load_data("../../inst/extdata/res_vpa_pma.rda")

test_that("extract_x() extracts values from VPA result", {
  expect_equal(
    extract_x(vpadata, "faa", 2011), vpadata$faa["2011"]
  )
  expect_equal(
    extract_x(vpadata, "faa", 2010:2011), vpadata$faa[as.character(2010:2011)]
  )

  expect_equal(
    extract_x(vpadata, "naa", 2011), vpadata$naa["2011"]
  )
  expect_equal(
    extract_x(vpadata, "naa", 2010:2011), vpadata$naa[as.character(2010:2011)]
  )
})

context("- extract_year_from()")

test_that("vpa object", {
  years <- as.numeric(colnames(vpadata$faa))

  expect_equal(extract_year_from(vpadata = vpadata), years)
  expect_equal(extract_year_from(vpadata = vpadata), years)
  expect_equal(extract_year_from(vpadata = vpadata, var = "naa"), years)
  expect_equal(extract_year_from(vpadata = vpadata, var = "baa"), years)
  expect_equal(extract_year_from(vpadata = vpadata, var = "faa"), years)
  expect_equal(extract_year_from(vpadata = vpadata, var = "saa"), years)

  expect_error(extract_year_from(),
               "Give me either 'vpadata' or 'msydata'")
})

test_that("msy object", {
  expect_error(extract_year_from(msydata = "input_MSY_object_here"),
               "Not implemented")
})

context("- extract_value()")

test_that("future_new", {
  expect_df(extract_value(from = generate_dummy_future_new_object(),
                          what = "ssb"))
  quick_test <- function(name) {
    expect_df(extract_value(from = generate_dummy_future_new_object(),
                            what = name))
  }

  quick_test("biomass")
  quick_test("catch")
})


context("Handling vectors")

test_that("select_from_tail() selects elements from tail", {
  vec <- 1:10
  expect_equal(select_from_tail(vec, relative = -1),    10)
  expect_equal(select_from_tail(vec, relative = -3:-1), 8:10)
  expect_equal(select_from_tail(vec, relative = -1:-3), 8:10)

  msg <- "'relative' should be negative"
  expect_error(select_from_tail(vec, relative = 0), msg)
  expect_error(select_from_tail(vec, relative = 1), msg)
})


test_that("wrap_by_paren() works", {

  expect_equal(wrap_by_paren(1), "(1)")
  expect_equal(wrap_by_paren("a"), "(a)")
  expect_equal(wrap_by_paren(NA), "(NA)")
  expect_equal(wrap_by_paren(NULL), "()")

  expect_equal(wrap_by_paren(1:3), "(1, 2, 3)")
  expect_equal(wrap_by_paren(c(1, 2, 3)), "(1, 2, 3)")
  expect_equal(wrap_by_paren(c(1, NA, 3)), "(1, NA, 3)")
})
