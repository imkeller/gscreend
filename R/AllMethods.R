# Access the sgRNA slot
setMethod("sgRNAData", "PoolScreenExp", function(x) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    se <- x@sgRNAData
    se
})

# Write into the sgRNA slot
setMethod("sgRNAData<-", "PoolScreenExp", function(x, value) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    x@sgRNAData <- value
    x
})

setMethod("normcounts", "PoolScreenExp", function(x)
    assays(sgRNAData(x))$normcounts)

setMethod("normcounts<-", "PoolScreenExp", function(x, value) {
    assays(sgRNAData(x))$normcounts <- value
    x
})

# get the reference counts for LFC calculation
setMethod("refcounts", "PoolScreenExp", function(x) {
    se <- sgRNAData(x)
    # filter to get T0 samples
    assays(se[, se$timepoint == "T0"])$normcounts
})

# get the reference counts for LFC calculation
setMethod("samplecounts", "PoolScreenExp", function(x) {
    se <- sgRNAData(x)
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$normcounts
})

setMethod("samplelfc", "PoolScreenExp", function(x) {
    se <- sgRNAData(x)
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$lfc
})

setMethod("samplelfc", "PoolScreenExp", function(x) {
    se <- sgRNAData(x)
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$lfc
})

setMethod("samplepval", "PoolScreenExp", function(x) {
    se <- sgRNAData(x)
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$pval
})
