library(SummarizedExperiment)
library(nloptr)
library(fGarch)

# wrapper funtion to run complete analysis
#' run PoolScreen
#'
#' @param object
#'
#' @return
#' @export
#' @import magrittr
#'
#' @examples
runPoolScreen <- function(object) {
    # normalize
    object <- normalizePoolScreenExp(object)
    # improve this unsing setMethod
    object <- calculateLFC(object)

    #sgRNA fitting
    object <- defineFittingIntervals(object)
    # rank sgRNAs
    object <- calculateIntervalFits(object)
    # limits determination must be somewhere outside of function!
    object <- calculatePValues(object)
    # rank genes
    object <- calculateGeneRank(object)
    object
}
