#' Plot replicate correlation
#'
#' @param object 
#' @param rep1 
#' @param rep2 
#'
#' @return plot
#' @export
#'
#' @examples load(system.file("data", "poolscreen_experiment.RData", package = "poolscreen"))
#' # pse is a poolscreen Experiment with two replicates R1 and R2
#' plotReplicateCorrelation(pse, rep1 = "R1", rep2 = "rep2")
plotReplicateCorrelation <- function(object, rep1 = "R1", rep2 = "R2") {
    se <- object@sgRNAData
    counts1 <- assays(se[, se$samplename == rep1])$normcounts
    counts2 <- assays(se[, se$samplename == rep2])$normcounts
    plot(log2(counts1+1), log2(counts2+1))
}


#' Plot model parameters from the fitting
#'
#' @param object 
#'
#' @return plot
#' @export
#'
#' @examples load(system.file("data", "poolscreen_experiment.RData", package = "poolscreen"))
#' # pse is a poolscreen Experiment
#' plotModelParameters(pse)
plotModelParameters <- function(object) {
    limits <- object@FittingIntervals
    limits <- limits[1:length(limits) - 1]
    parameters <- object@LFCModelParameters
    par(mfrow=c(2,2))
    plot(limits, parameters[,1], main = "Mean")
    plot(limits, parameters[,2], main = "Sd")
    plot(limits, parameters[,3], main = "Xi")
}