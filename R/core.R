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
runPoolScreen <- function(object, quant1, quant2, alphacutoff) {
    # normalize
    object <- normalizePoolScreenExp(object)
    # improve this unsing setMethod
    object <- calculateLFC(object)

    #sgRNA fitting
    object <- defineFittingIntervals(object)
    object <- calculateIntervalFits(object, quant1, quant2)
    # pvalues and rank sgRNAs
    object <- calculatePValues(object)
    # rank genes
    object <- assignGeneData(object, alphacutoff)
    object
}
