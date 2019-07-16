#' Normalize raw count
#'
#' @param object
#'
#' @return object
#' @export
#'
#' @examples load(system.file("data", "gscreend_experiment.RData", package = "gscreend"))
#' # pse is a of PoolScreenExp class
#' normalizePoolScreenExp(pse)
normalizePoolScreenExp <- function(object) {
    sgRNAse <- object@sgRNAData
    sizefactors <- colSums(assays(sgRNAse)$counts)
    maxfact <- max(sizefactors)
    # want to divide each row of matrix by vector elements.
    # only works with t() trick
    normcounts(object) <- t(t(assays(sgRNAse)$counts)/sizefactors) * maxfact
    object
}

#' Calculate log fold changes
#'
#' @param object
#'
#' @return object
#' @export
#'
#' @examples load(system.file("data" "gscreend_experiment.RData", package = "gscreend"))
#' # pse is a of PoolScreenExp class
#' calculateLFC(pse)
calculateLFC <- function(object) {
    sgRNAse <- object@sgRNAData
    assays(object@sgRNAData)$lfc <- log2((normcounts(object) + 1)/
                               (refcounts(object)[,1] + 1))
    object
}
