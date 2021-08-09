context("1 - Data input and generate SummarizedExperiment")

test_that("createPoolScreenExp() throws error without valid input", {
    expect_error(createPoolScreenExp(1:10))
})


test_that("createPoolScreenExpFromSE() throws error when sample names are not as correct", {
    raw_counts <- read.table(
        system.file('extdata', 'simulated_counts.txt',
        package = 'gscreend'),
        header=TRUE)

    # 1
    # more than one reference sample T0
    counts_matrix <- cbind(raw_counts$library0, raw_counts$R0_0, raw_counts$R1_0)
    rowData <- data.frame(sgRNA_id = raw_counts$sgrna_id,
        gene = raw_counts$Gene)
    colData_wrong1 <- data.frame(samplename = c('library', 'R1', 'R2'),
        timepoint = c('T0', 'T0', 'T1'))
    se_wrong1 <- SummarizedExperiment(assays=list(counts=counts_matrix),
        rowData=rowData, colData=colData_wrong1)
    # create a PoolScreenExp experiment
    expect_error(createPoolScreenExpFromSE(se_wrong1))

    # 2
    # wrong name
    colData_wrong2 <- data.frame(samplename = c('library', 'R1', 'R2'),
                                 timepoint = c('T0', 'T1', 'Wrong'))
    se_wrong2 <- SummarizedExperiment(assays=list(counts=counts_matrix),
                                      rowData=rowData, colData=colData_wrong2)
    # create a PoolScreenExp experiment
    expect_error(createPoolScreenExpFromSE(se_wrong2))

    # 3
    # no sample T1
    counts_matrix_wrong3 <- cbind(raw_counts$library0)
    colData_wrong3 <- data.frame(samplename = c('library'),
                                 timepoint = c('T0'))
    se_wrong3 <- SummarizedExperiment(assays=list(counts=counts_matrix_wrong3),
                                      rowData=rowData, colData=colData_wrong3)
    # create a PoolScreenExp experiment
    expect_error(createPoolScreenExpFromSE(se_wrong3))

})
