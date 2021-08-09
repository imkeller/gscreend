context("3 - Functions used for calculating statistics")

test_that("alphaBeta() returns a numeric value", {
    expect_is(alphaBeta(c(0.1, 0.2, 0.3)), "numeric")
})
