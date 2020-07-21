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
context("- extract_fmsy()")

test_that("extract_fmsy()", {
  msyoutput <- load_data("../inst/extdata/res_MSY_pma_pre.rda")

  fmsy_expected  <- c("0歳"     = 0.2796186,
                      "1歳"     = 0.6052315,
                      "2歳"     = 0.6717214,
                      "3歳以上" = 0.6717213)

  expect_equal(extract_fmsy(msyoutput),
               fmsy_expected,
               tolerance = 10^-7)

  expect_equal(extract_fmsy(msyoutput, mean = TRUE),
               mean(fmsy_expected),
               tolerance = 10^-7)
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
  quick_test("fratio")
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

test_that("add_age_suffix() simply adds age suffixes to vectors", {

  expect_equal(add_age_suffix(0:2),  c("0歳", "1歳", "2歳以上"))
  expect_equal(add_age_suffix(-1:1), c("-1歳", "0歳", "1歳以上"))

  # Single vectors
  expect_equal(add_age_suffix(0),    c("0歳以上"))
  expect_equal(add_age_suffix(-1),   c("-1歳以上"))

  # Characters
  expect_equal(add_age_suffix(""),   c("歳以上"))
  expect_equal(add_age_suffix(c("a", "b", "c")), c("a歳", "b歳", "c歳以上"))
  expect_equal(add_age_suffix(c("a", "b", "")),  c("a歳", "b歳", "歳以上"))

  # NAs
  expect_equal(add_age_suffix(NA),   c("NA歳以上"))
  expect_equal(add_age_suffix(c(NA, "b", NA)), c("NA歳", "b歳", "NA歳以上"))

  # NULLs
  expect_equal(add_age_suffix(NULL), "歳以上")
  expect_equal(add_age_suffix(c(1, 2, NULL)), c("1歳", "2歳以上"))
  expect_equal(add_age_suffix(c(NULL, 1, 2)), c("1歳", "2歳以上"))
})

test_that("give_agename() works", {

  vec       <- 1:3
  named_vec <- give_agename(1:3)
  expect_equal(unname(named_vec), vec)
  expect_equal(names(named_vec),  c("0歳", "1歳", "2歳以上"))

  expect_error(
    # Not mere concats of vector and age suffix
    expect_equal(names(give_agename(1:3)), c("1歳", "2歳", "3歳以上"))
  )

  quick_test <- function(v, name_expected) {
    named <- give_agename(v)
    expect_equal(unname(named), v)
    expect_equal(names(named),  name_expected)
  }

  quick_test(1:3, c("0歳", "1歳", "2歳以上"))
  quick_test(0:3, c("0歳", "1歳", "2歳", "3歳以上"))
  quick_test(2:3, c("0歳", "1歳以上"))
  quick_test(2:4, c("0歳", "1歳", "2歳以上"))
  quick_test(-1:1, c("0歳", "1歳", "2歳以上"))
  quick_test(c(1.234, 2.345, 3.456),
             c("0歳", "1歳", "2歳以上"))
})

test_that("convert_unit() works", {
  expect_equal(convert_unit(tons = 1000, to = "千トン"), 1)
  expect_equal(convert_unit(tons = 1500, to = "千トン"), 1.5)

  expect_equal(convert_unit(tons = 1000, to = "千トン", add_unit = TRUE),
               "1 千トン")
  expect_equal(convert_unit(tons = 1500, to = "千トン", add_unit = TRUE),
               "1.5 千トン")

  expect_equal(convert_unit(tons = c(1500, 2500), to = "千トン"),
               c(1.5, 2.5))
  expect_equal(convert_unit(tons = c(1234, 23456), to = "千トン"),
               c(1.234, 23.456))

  expect_equal(convert_unit(tons     = c(1500, 2500),
                            to       = "千トン",
                            add_unit = TRUE),
               c("1.5 千トン", "2.5 千トン"))

  expect_equal(convert_unit(tons     = c(1234, 23456),
                            to       = "千トン",
                            add_unit = TRUE),
               c("1.234 千トン", "23.456 千トン"))

  expect_equal(convert_unit(tons     = c(1234, 23456),
                            to       = "千トン",
                            add_unit = TRUE,
                            round    = 1),
               c("1.2 千トン", "23.5 千トン"))

  expect_equal(convert_unit(tons     = c(1234, 23456),
                            to       = "千トン",
                            round    = 1),
               c(1.2, 23.5))

  expect_equal(convert_unit(tons = 1000, to = "百トン"), 10)
  expect_equal(convert_unit(tons = 1500, to = "百トン"), 15)

  msg <- "'to' should be either '百トン' or '千トン'"
  expect_error(convert_unit(tons = 1000, to = "万トン"), msg)
  expect_error(convert_unit(tons = 1000, to = "foo"),    msg)
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
