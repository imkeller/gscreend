# sgRNA slot
#' @title sgRNAData: set and retrieve sgRNAData of PoolScreenExp
#' @name sgRNAData
#' @param x PoolScreenExp object
#' @return sgRNA slot of the object
#' @export
#' @examples # import a PoolScreenExp object that has been generated using
#' # RunGscreend()
#' pse_an <- readRDS(
#' system.file('extdata', 'gscreend_analysed_experiment.RData',
#' package = 'gscreend'))
#'
#' sgRNAData(pse_an)
#'
setGeneric("sgRNAData", function(x) standardGeneric("sgRNAData"))

setGeneric("sgRNAData<-", function(x, value) standardGeneric("sgRNAData<-"))

# Gene slot
#' @title GeneData: set and retrieve GeneData of PoolScreenExp
#' @name GeneData
#' @param x PoolScreenExp object
#' @return Gene slot of the object
#' @export
#' @examples # import a PoolScreenExp object that has been generated using
#' # RunGscreend()
#' pse_an <- readRDS(
#' system.file('extdata', 'gscreend_analysed_experiment.RData',
#' package = 'gscreend'))
#'
#' GeneData(pse_an)
#'
setGeneric("GeneData", function(x) standardGeneric("GeneData"))

setGeneric("GeneData<-", function(x, value) standardGeneric("GeneData<-"))

# FittingOptions slot
setGeneric("FittingOptions",
        function(x) standardGeneric("FittingOptions"))
setGeneric("FittingOptions<-",
        function(x, value) standardGeneric("FittingOptions<-"))

# FittingIntervals slot
setGeneric("FittingIntervals",
        function(x) standardGeneric("FittingIntervals"))
setGeneric("FittingIntervals<-",
        function(x, value) standardGeneric("FittingIntervals<-"))

# LFCModelParameters slot
setGeneric("LFCModelParameters",
        function(x) standardGeneric("LFCModelParameters"))
setGeneric("LFCModelParameters<-",
        function(x, value) standardGeneric("LFCModelParameters<-"))

# set and get normalized counts
setGeneric("normcounts", function(x) standardGeneric("normcounts"))
setGeneric("normcounts<-", function(x, value) standardGeneric("normcounts<-"))

# get reference(library) counts
setGeneric("refcounts", function(x) standardGeneric("refcounts"))

# get sample counts
setGeneric("samplecounts", function(x) standardGeneric("samplecounts"))

# get sample lfc
setGeneric("samplelfc", function(x) standardGeneric("samplelfc"))

# get sample pvalues
setGeneric("samplepval", function(x) standardGeneric("samplepval"))
