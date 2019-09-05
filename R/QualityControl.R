#' Plot replicate correlation
#'
#' @param object PoolScreenExp object
#' @param rep1 Name of replicate 1
#' @param rep2 Name of replicate 2
#'
#' @return replicate_plot
#' @export
#'
#' @examples
#' # import a PoolScreenExp object that has been generated using RunGscreend()
#' pse_an <- readRDS(
#' system.file("extdata", "gscreend_analysed_experiment.RData",
#' package = "gscreend"))
#' plotReplicateCorrelation(pse_an, rep1 = "R1", rep2 = "R2")
#'
plotReplicateCorrelation <- function(object, rep1 = "R1", rep2 = "R2") {
    se <- object@sgRNAData
    counts1 <- assays(se[, se$samplename == rep1])$normcounts
    counts2 <- assays(se[, se$samplename == rep2])$normcounts
    plot(log2(counts1+1), log2(counts2+1))
}


#' Plot model parameters from the fitting
#'
#' @param object  PoolScreenExp object
#'
#' @return plot
#' @export
#'
#' @examples # import a PoolScreenExp object that has been generated using
#' # RunGscreend()
#' pse_an <- readRDS(
#' system.file("extdata", "gscreend_analysed_experiment.RData",
#' package = "gscreend"))
#' plotModelParameters(pse_an)
plotModelParameters <- function(object) {
    limits <- object@FittingIntervals
    limits <- limits[1:length(limits) - 1]
    parameters <- object@LFCModelParameters
    graphics::par(mfrow=c(2,2))
    plot(limits, parameters[,1], main = "Mean")
    plot(limits, parameters[,2], main = "Sd")
    plot(limits, parameters[,3], main = "Xi")
}
