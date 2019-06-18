context("start up")

test_that("read the data",{
    caa <- read.csv(system.file("extdata","caa_pma.csv",package="frasyr"))
    expect_equal(nrow(caa),4)
})




