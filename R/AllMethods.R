setMethod("normcounts", "PoolScreenExp", function(x)
    assays(x@sgRNAData)$normcounts)
setMethod("normcounts<-", "PoolScreenExp", function(x, value) {
    assays(x@sgRNAData)$normcounts <- value
    x
})

# get the reference counts for LFC calculation
setMethod("refcounts", "PoolScreenExp", function(x) {
    se <- object@sgRNAData
    # filter to get T0 samples
    assays(se[, se$timepoint == "T0"])$normcounts
})

# get the reference counts for LFC calculation
setMethod("samplecounts", "PoolScreenExp", function(x) {
    se <- object@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$normcounts
})

setMethod("samplelfc", "PoolScreenExp", function(x) {
    se <- object@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$lfc
})

setMethod("samplelfc", "PoolScreenExp", function(x) {
    se <- object@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$lfc
})

setMethod("samplepval", "PoolScreenExp", function(x) {
    se <- object@sgRNAData
    # filter to get T0 samples
    # names have to be checked upon input
    assays(se[, se$timepoint == "T1"])$pval
})


# Has a Summarized Experiment (takes as input)

# Can be normalized <- normalized counts go back to object as assay

# LFC can be calculated <- lfc go back to object as assay

# model can be associated per cellline (technical replicates)

# p values can be calculated per sgRNA

# another table is there to show gene p values and gene LFCs
