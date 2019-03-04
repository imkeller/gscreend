setMethod("normcounts", "PoolScreenExp", function(x)
    assays(x@sgRNAData)$normcounts)
setMethod("normcounts<-", "PoolScreenExp", function(x, value) {
    assays(x@sgRNAData)$normcounts <- value
    x
})

# get the reference counts for LFC calculation
setMethod("refcounts", "PoolScreenExp", function(x) {
    se <- x@sgRNAData
    # filter to get T0 samples
    assays(se[, se$timepoint == "T0"])$normcounts
})

# get the reference counts for LFC calculation
setMethod("samplecounts", "PoolScreenExp", function(x) {
    se <- x@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$normcounts
})

setMethod("samplelfc", "PoolScreenExp", function(x) {
    se <- x@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$lfc
})

setMethod("samplelfc", "PoolScreenExp", function(x) {
    se <- x@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$lfc
})

setMethod("samplepval", "PoolScreenExp", function(x) {
    se <- x@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$pval
})