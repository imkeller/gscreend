#' Normalize raw count
#'
#' @param object PoolScreenExp object
#'
#' @return object
#' @export
#'
#' @examples pse <- readRDS(system.file("extdata", "gscreend_experiment.RData",
#' package = "gscreend"))
#' # pse is a PoolScreenExp
#' normalizePoolScreenExp(pse)
normalizePoolScreenExp <- function(object) {
    sgRNAse <- object@sgRNAData
    sizefactors <- colSums(assays(sgRNAse)$counts)
    maxfact <- max(sizefactors)
    # want to divide each row of matrix by vector elements.
    # only works with t()
    normcounts(object) <- t(t(assays(sgRNAse)$counts)/sizefactors) * maxfact
    object
}

#' Calculate log fold changes
#'
#' @param object PoolScreenExp object
#'
#' @return object
#' @export
#'
#' @examples pse <- readRDS(system.file("extdata", "gscreend_experiment.RData",
#' package = "gscreend"))
#' # pse is a of PoolScreenExp class
#' norm_pse <- normalizePoolScreenExp(pse)
#' calculateLFC(norm_pse)
calculateLFC <- function(object) {
    sgRNAse <- object@sgRNAData
    assays(object@sgRNAData)$lfc <- log2((normcounts(object) + 1)/
                                             (refcounts(object)[,1] + 1))
    object
}
