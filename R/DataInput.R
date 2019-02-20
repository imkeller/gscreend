
#' Create PoolScreen Experiment
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
createPoolScreenExp <- function(data) {
    if(is(data, "SummarizedExperiment")) {
        message("Creating PoolScreenExp object from a SummarizedExperiment object.")
        object <- createPoolScreenExpFromSE(data)
    } else {
        message("Error: input data is not of the suported input format: SummarizedExperiment.")
    }
    object
}

#' Create PoolScreen Experiment from summarized experiment
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
createPoolScreenExpFromSE <- function(data) {
    object <- new("PoolScreenExp")
    object@sgRNAData <- data
    object@FittingOptions <- list(IntervalFraction = 0.1,
                                  alphaCutoff = 0.05)
    object
}
