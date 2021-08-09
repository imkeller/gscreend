context("2 - Quality control functions")

test_that("plotReplicateCorrelation() throws error when replicate names are not known", {
    pse_an <- readRDS(
        system.file('extdata', 'gscreend_analysed_experiment.RData',
        package = 'gscreend'))
    expect_error(plotReplicateCorrelation(pse_an, rep1 = "wrong", rep2 = "R2"))
})

test_that("ResultsTable() throws error when direction is wrong", {
    pse_an <- readRDS(
        system.file('extdata', 'gscreend_analysed_experiment.RData',
                    package = 'gscreend'))
    expect_error(ResultsTable(pse_an, direction = "wrong"))
})

