library(SummarizedExperiment)
library(nloptr)
library(fGarch)

# wrapper funtion to run complete analysis
#' run PoolScreen
#'
#' @param object
#'
#' @return object
#' @export
#' @import magrittr
#'
#' @examples
runPoolScreen <- function(object, 
                          quant1 = 0.1, 
                          quant2 = 0.9, 
                          alphacutoff = 0.05) {
    # normalize
    normalizePoolScreenExp(object) %>%
    # improve this unsing setMethod
    calculateLFC() %>%
    #sgRNA fitting
    defineFittingIntervals() %>%
    calculateIntervalFits(quant1, quant2) %>%
    # pvalues and rank sgRNAs
    calculatePValues() %>%
    # rank genes
    assignGeneData(alphacutoff)
}
