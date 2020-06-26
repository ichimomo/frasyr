context("User interfaces")

test_that("load_data() works", {
  expect_is(load_data("../../inst/extdata/res_vpa_pma.rda"), "list")
  expect_is(load_data("../../inst/extdata/vpa.csv"), "list")
})

context("Reuse input from the result of scientist meeting")

dummy_future_result <-
  list(
    data = "data",
    input = list(model_average_option = 123,
                 res_SR = list(pars = list(sd = 0.123)))
  )

test_that("reuse inputs as-is", {

  retrieved <- retrieve_input(dummy_future_result)

  expect_setequal(names(retrieved), c("model_average_option", "res_SR"))

  expect_equal(retrieved$model_average_option, 123)
  expect_equal(retrieved$res_SR$pars$sd,     0.123)

  list_without_input <- list(data = "data")

  expect_error(
    retrieve_input(list_without_input),
  )
})

test_that("reuse inputs with modification", {

  retrieved <- retrieve_input(dummy_future_result, new_sd = 0)
  expect_equal(retrieved$res_SR$pars$sd, 0)

  retrieved <- retrieve_input(dummy_future_result, new_sd = 0.456)
  expect_equal(retrieved$res_SR$pars$sd, 0.456)

  expect_error(retrieve_input(dummy_future_result, new_sd = "a"),
               "'new_sd' should be numeric")
})
