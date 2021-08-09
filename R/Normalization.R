#' Normalize raw count
#'
#' @param object PoolScreenExp object
#'
#' @return object
#' @keywords internal
normalizePoolScreenExp <- function(object) {
    sgRNAse <- sgRNAData(object)
    sizefactors <- colSums(assays(sgRNAse)$counts)
    maxfact <- max(sizefactors)
    # want to divide each row of matrix by vector elements.  only works with t()
    normcounts(object) <- t(t(assays(sgRNAse)$counts)/sizefactors) * maxfact
    message("Size normalized count data.")
    object
}

#' Calculate log fold changes
#'
#' @param object PoolScreenExp object
#'
#' @return object
#' @keywords internal
calculateLFC <- function(object) {
    sgRNAse <- sgRNAData(object)
    assays(sgRNAData(object))$lfc <-log2(
        (normcounts(object) + 1)/(refcounts(object)[, 1] + 1)
        )
    message("Calculated LFC.")
    object
}
