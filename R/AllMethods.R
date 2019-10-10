# Access the sgRNA slot
#' @export
setMethod("sgRNAData", "PoolScreenExp", function(x) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    se <- x@sgRNAData
    se
})

# Write into the sgRNA slot
#' @export
setMethod("sgRNAData<-", "PoolScreenExp", function(x, value) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    x@sgRNAData <- value
    x
})

# Access the Gene slot
setMethod("GeneData", "PoolScreenExp", function(x) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    se <- x@GeneData
    se
})

# Write into the Gene slot
setMethod("GeneData<-", "PoolScreenExp", function(x, value) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    x@GeneData <- value
    x
})

# Fitting options slot
setMethod("FittingOptions", "PoolScreenExp", function(x) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    se <- x@FittingOptions
    se
})
setMethod("FittingOptions<-", "PoolScreenExp", function(x, value) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    x@FittingOptions <- value
    x
})

# Fitting intervals slot
setMethod("FittingIntervals", "PoolScreenExp", function(x) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    se <- x@FittingIntervals
    se
})
setMethod("FittingIntervals<-", "PoolScreenExp", function(x, value) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    x@FittingIntervals <- value
    x
})

# slot containing fitted parameters
setMethod("LFCModelParameters", "PoolScreenExp", function(x) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    se <- x@LFCModelParameters
    se
})
setMethod("LFCModelParameters<-", "PoolScreenExp", function(x, value) {
    # I assume I am allowed to use the
    # direct slot access here, because there is no other way?
    x@LFCModelParameters <- value
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
