#' Normalize raw count
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples load(system.file("data", "poolscreen_experiment.RData", package = "poolscreen"))
#' # pse is a poolscreen Experiment
#' normalizePoolScreenExp(pse)
normalizePoolScreenExp <- function(object) {
    sgRNAse <- object@sgRNAData
    sizefactors <- colSums(assays(sgRNAse)$counts)
    normcounts(object) <- assays(sgRNAse)$counts/sizefactors * max(sizefactors)
    object
}

#' Calculate log fold changes
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples load(system.file("data", "poolscreen_experiment.RData", package = "poolscreen"))
#' # pse is a poolscreen Experiment
#' calculateLFC(pse)
calculateLFC <- function(object) {
    sgRNAse <- object@sgRNAData
    assays(object@sgRNAData)$lfc <- log2((normcounts(object) + 1)/
                               (refcounts(object)[,1] + 1))
    object
}
