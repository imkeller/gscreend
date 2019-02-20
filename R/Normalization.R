#' Normalize raw count
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
normalizePoolScreenExp <- function(object) {
    sgRNAse <- object@sgRNAData
    sizefactors <- colSums(assays(sgRNAse)$counts)
    normcounts(object) <- assays(sgRNAse)$counts/sizefactors * max(sizefactors)
    # how can you modify the object withour returning and copying?
    object
}

#' Calculate log fold changes
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
calculateLFC <- function(object) {
    sgRNAse <- object@sgRNAData
    assays(object@sgRNAData)$lfc <- log2((normcounts(object) + 1)/
                               (refcounts(object)[,1] + 1))
    object
}
